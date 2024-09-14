import anndata as ad
import pandas as pd
import scanpy as sc
import sys
from scipy.sparse import csr_matrix


SITE1_CELL = 16311
SITE2_CELL = 25171
SITE3_CELL = 32029
SITE4_CELL = 16750
GEX = 13953


adata = ad.read_h5ad("/workspace/scBFP/dataset/citeseq_processed.h5ad")
adata.var_names_make_unique()

X = adata.X.toarray()
X[SITE1_CELL + SITE2_CELL: SITE1_CELL + SITE2_CELL + SITE3_CELL, :GEX] = 0

adata.obsm["missing"] = csr_matrix(X)
adata.write_h5ad("/workspace/scBFP/dataset/citeseq_processed.h5ad")
