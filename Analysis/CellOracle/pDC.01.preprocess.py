#!/usr/bin/env python


import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.figsize"] = [6, 4.5]

# cell
fname = "/work/home/project/scESCA/200227_6samples/output/04.rm_cells/seurat_pDC.h5ad"
cell_group = "level_3"
cell_name = "pDC"

# file
adata = sc.read_h5ad(fname)
adata

# Only consider genes with more than 1 count
sc.pp.filter_genes(adata, min_counts=1)
# Normalize gene expression matrix with total UMI count per cell
sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

# Select top 2000 highly-variable genes
filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                              flavor='cell_ranger',
                                              n_top_genes=2000,
                                              log=False)

# Subset the genes
adata = adata[:, filter_result.gene_subset]

# keep raw cont data before log transformation
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()

# Renormalize after filtering
sc.pp.normalize_per_cell(adata)

# Log transformation and scaling
sc.pp.log1p(adata)
sc.pp.scale(adata)

# PCA
sc.tl.pca(adata, svd_solver='arpack')

# Diffusion map
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)

sc.tl.diffmap(adata)
# Calculate neihbors again based on diffusionmap
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

# PAGA graph construction
sc.tl.paga(adata, groups=cell_group)
plt.rcParams["figure.figsize"] = [10, 8]
sc.pl.paga(adata)
sc.tl.draw_graph(adata, init_pos='paga', random_state=123)
sc.pl.draw_graph(adata, color=cell_group, legend_loc='on data')

adata.write_h5ad(f"data/{cell_name}.h5ad")
os.getcwd()




