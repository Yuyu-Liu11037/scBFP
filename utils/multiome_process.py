import anndata as ad
import pandas as pd
import scanpy as sc
import sys
from scipy.sparse import csr_matrix


SITE1_CELL = 17243
SITE2_CELL = 15226
SITE3_CELL = 14556
SITE4_CELL = 22224

adata = ad.read_h5ad("/workspace/scBFP/dataset/multiome_processed.h5ad")
adata.var_names_make_unique()

### preprocess
adata_GEX = adata[:, adata.var['feature_types'] == 'GEX'].copy()
adata_ATAC = adata[:, adata.var['feature_types'] == 'ATAC'].copy()
### step 1: normalize
sc.pp.normalize_total(adata_GEX, target_sum=1e4)
sc.pp.normalize_total(adata_ATAC, target_sum=1e4)
### step 2: log transform
sc.pp.log1p(adata_GEX)
sc.pp.log1p(adata_ATAC)
### step 3: select highly variable features
sc.pp.highly_variable_genes(adata_GEX, subset=True)
sc.pp.highly_variable_genes(
    adata_ATAC,
    n_top_genes=4000,
    subset=True
)

num_atac = adata_ATAC.X.shape[1]
adata = ad.concat([adata_ATAC, adata_GEX], axis=1, merge="first")   # left num_atac: ATAC, right 2832: GEX
print(f"Finish preprocessing\n")

X = adata.X.copy()
X = X.toarray()
X[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :num_atac] = 0
adata.obsm["missing"] = csr_matrix(X)
adata.write_h5ad("/workspace/scBFP/dataset/multiome.h5ad")
