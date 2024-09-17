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
import matplotlib.pyplot as plt
import seaborn as sns


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
        print(self.adata)


    def evaluate(self):
        if self.args.name == "citeseq":
            SITE1_CELL = 16311
            SITE2_CELL = 25171
            SITE3_CELL = 32029
            SITE4_CELL = 16750
            GEX = 2000
        elif self.args.name == "multiome":
            SITE1_CELL = 17243
            SITE2_CELL = 15226
            SITE3_CELL = 14556
            SITE4_CELL = 22224
            GEX = 4000 # TODO: wrong code

        adata_imputed = self.adata.copy()
        adata_imputed.X = self.adata.obsm['imputation']

        gt = self.adata
        # similarity
        X_true = gt.X.toarray()
        X_imputed_np = adata_imputed.X
        nonzero_mask31 = (X_true[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX] != 0)
        MAE = mean_absolute_error(X_imputed_np[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31], X_true[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31])
        RMSE = sqrt(mean_squared_error(X_imputed_np[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31], X_true[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31]))
        PCC = pearsonr(X_imputed_np[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31], X_true[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX][nonzero_mask31])[0]
        
        # clustering
        sc.pp.pca(adata_imputed)
        sc.pp.neighbors(adata_imputed, use_rep="X_pca")
        true_labels = self.adata.obs[self.args.cell_type_label].values

        resolution = 0.2 if self.args.name == "citeseq" else 0.65
        sc.tl.leiden(adata_imputed, resolution=resolution, n_iterations=2)
        predicted_labels = adata_imputed.obs["leiden"]
    
        ari = adjusted_rand_score(true_labels, predicted_labels)
        nmi = normalized_mutual_info_score(true_labels, predicted_labels)
        true_labels = [str(label) for label in true_labels]
        predicted_labels = [str(label) for label in predicted_labels]
        contingency_matrix = confusion_matrix(true_labels, predicted_labels)
        purity = np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix)

        # visualization
        sc.tl.umap(adata_imputed)
        plt.figure(figsize=(8, 6))
        sc.pl.umap(adata_imputed, color='leiden', show=False)
        plt.savefig(f"/workspace/scBFP/visualization/umap_{self.args.name}.png", format="png")
        plt.close()
        print(f"UMAP plot saved as 'umap_{self.args.name}.png'.")

        print(f" ==================== Dataset: {self.args.name} ==================== ")
        if self.args.drop_rate != 0.0:
            print("Drop Rate {} -> RMSE : {:.4f} / Median L1 Dist : {:.4f}\n".format(self.args.drop_rate, rmse, median_l1_distance))

        if self.args.eval_clustering:
            print("Imputed --> MAE : {} / RMSE : {} / PCC : {} / ARI : {} / NMI : {} / PURITY : {}\n".format(MAE, RMSE, PCC, ari, nmi, purity))


