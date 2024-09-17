import anndata as ad
import pandas as pd
import scanpy as sc
import sys
from scipy.sparse import csr_matrix


SITE1_CELL = 16311
SITE2_CELL = 25171
SITE3_CELL = 32029
SITE4_CELL = 16750
GEX = 2000


adata = ad.read_h5ad("/workspace/scBFP/dataset/citeseq_processed.h5ad")
adata.var_names_make_unique()

adata_GEX = adata[:, adata.var["feature_types"] == "GEX"].copy()
adata_ADT = adata[:, adata.var["feature_types"] == "ADT"].copy()
sc.pp.normalize_total(adata_GEX, target_sum=1e4)
sc.pp.normalize_total(adata_ADT, target_sum=1e4)
sc.pp.log1p(adata_GEX)
sc.pp.log1p(adata_ADT)
sc.pp.highly_variable_genes(
    adata_GEX,
    n_top_genes=2000,
    subset=True
)
adata = ad.concat([adata_GEX, adata_ADT], axis=1, merge="first")   # X(:,1): GEX, X(:,2): ADT
print(f"Finish preprocessing\n")

X = adata.X.copy()
X = X.toarray()
X[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX] = 0
adata.obsm["missing"] = csr_matrix(X)
adata.write_h5ad("/workspace/scBFP/dataset/citeseq.h5ad")
