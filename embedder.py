import os
import time
import scanpy as sc
import numpy as np
import anndata as ad
import sys

from argument import printConfig, config2string
from misc.utils import drop_data
from math import sqrt
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from misc.utils import imputation_error, cluster_acc
from sklearn.metrics import confusion_matrix, mean_absolute_error, mean_squared_error, jaccard_score
from sklearn.metrics.cluster import adjusted_rand_score, normalized_mutual_info_score
from sklearn.preprocessing import LabelEncoder

class embedder:
    def __init__(self, args):
        self.args = args
        printConfig(args)
        # self.config_str = config2string(args)
        self.device = f'cuda:{args.device}' if self.args.gpu else "cpu"
        self._init_dataset()
        # self.start_time = time.time()

    def _init_dataset(self):

        self.adata = sc.read(f'./dataset/{self.args.name}.h5ad')
        # if self.args.eval_clustering and self.adata.obs[self.args.cell_type_label].dtype != int:
        #     self.label_encoding()

        # self.preprocess(filter=self.args.filter, cell_min=self.args.cell_min, gene_min=self.args.gene_min, hvg=self.args.hvg, n_hvg=self.args.n_hvg, size_factors=self.args.sf, logtrans_input=self.args.log)
        # self.adata = drop_data(self.adata, rate=self.args.drop_rate)


    def label_encoding(self):
        label_encoder = LabelEncoder()
        celltype = self.adata.obs[self.args.cell_type_label]
        celltype = label_encoder.fit_transform(celltype)
        self.adata.obs[self.args.cell_type_label] = celltype

    # def preprocess(self, filter, cell_min=1, gene_min=1, hvg=True, n_hvg=2000, size_factors=True, logtrans_input=True):

    #     if filter:
    #         sc.pp.filter_cells(self.adata, min_counts=cell_min)
    #         sc.pp.filter_genes(self.adata, min_counts=gene_min)

    #     if hvg:
    #         variance = np.array(self.adata.X.todense().var(axis=0))[0]
    #         hvg_gene_idx = np.argsort(variance)[-int(n_hvg):]
    #         self.adata = self.adata[:,hvg_gene_idx]

    #     self.adata.raw = self.adata.copy()        
    #     if size_factors:
    #         sc.pp.normalize_per_cell(self.adata)
    #         self.adata.obs['size_factors'] = self.adata.obs.n_counts / np.median(self.adata.obs.n_counts)
    #     else:
    #         self.adata.obs['size_factors'] = 1.0

    #     if logtrans_input:
    #         sc.pp.log1p(self.adata)
    def preprocess_citeseq(self, citeseq):
        adata_GEX = citeseq[:, citeseq.var["feature_types"] == "GEX"].copy()
        adata_ADT = citeseq[:, citeseq.var["feature_types"] == "ADT"].copy()
        sc.pp.normalize_total(adata_GEX, target_sum=1e4)
        sc.pp.normalize_total(adata_ADT, target_sum=1e4)
        sc.pp.log1p(adata_GEX)
        sc.pp.log1p(adata_ADT)
        sc.pp.highly_variable_genes(
            adata_GEX,
            n_top_genes=2000,
            subset=True
        )
        citeseq = ad.concat([adata_GEX, adata_ADT], axis=1, merge="first")   # X(:,1): GEX, X(:,2): ADT
        print(f"Finish preprocessing\n")
        return citeseq


    def evaluate(self):
        SITE1_CELL = 16311
        SITE2_CELL = 25171
        SITE3_CELL = 32029
        SITE4_CELL = 16750
        GEX = 2000

        adata_imputed = self.adata.copy()
        adata_imputed.X = self.adata.obsm['imputation']

        gt = self.preprocess_citeseq(self.adata)
        adata_imputed = self.preprocess_citeseq(adata_imputed)
        # similarity
        X_true = gt.X.toarray()
        X_imputed_np = adata_imputed.X
        nonzero_mask31 = (X_true[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX] != 0)
        MAE = mean_absolute_error(X_imputed_np[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31], X_true[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31])
        RMSE = sqrt(mean_squared_error(X_imputed_np[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31], X_true[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31]))
        PCC = pearsonr(X_imputed_np[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31], X_true[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31])[0]
        
        # clustering
        if self.args.eval_clustering:
            sc.pp.pca(adata_imputed)
            sc.pp.neighbors(adata_imputed, use_rep="X_pca")
            true_labels = self.adata.obs[self.args.cell_type_label].values

            sc.tl.leiden(adata_imputed, resolution=0.2, flavor="igraph", n_iterations=2)
            predicted_labels = adata_imputed.obs["leiden"]
    
            ari = adjusted_rand_score(true_labels, predicted_labels)
            nmi = normalized_mutual_info_score(true_labels, predicted_labels)
            true_labels = [str(label) for label in true_labels]
            predicted_labels = [str(label) for label in predicted_labels]
            contingency_matrix = confusion_matrix(true_labels, predicted_labels)
            purity = np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix)

        print(f" ==================== Dataset: {self.args.name} ==================== ")
        if self.args.drop_rate != 0.0:
            print("Drop Rate {} -> RMSE : {:.4f} / Median L1 Dist : {:.4f}\n".format(self.args.drop_rate, rmse, median_l1_distance))

        if self.args.eval_clustering:
            print("Imputed --> MAE : {} / RMSE : {} / PCC : {} / ARI : {} / NMI : {} / PURITY : {}\n".format(MAE, RMSE, PCC, ari, nmi, purity))


            

