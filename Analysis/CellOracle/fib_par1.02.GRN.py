#!/usr/bin/env python

# 0. Import

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import celloracle as co
co.__version__

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

save_folder = "figures"
os.makedirs(save_folder, exist_ok=True)

cell_name = "fib_par1"
cell_group = "level_3"
adata = sc.read_h5ad(f"data/{cell_name}.h5ad")
adata

print(f"Cell number is :{adata.shape[0]}")
print(f"Gene number is :{adata.shape[1]}")

# Load TF info which was made from mouse cell atlas dataset.
base_GRN = co.data.load_human_promoter_base_GRN(version="hg38_gimmemotifsv5_fpr2")

# Check data
base_GRN.head()
#base_GRN.to_csv('test.csv')

# Instantiate Oracle object
oracle = co.Oracle()

# Check data in anndata
print("Metadata columns :", list(adata.obs.columns))
print("Dimensional reduction: ", list(adata.obsm.keys()))

# In this notebook, we use the unscaled mRNA count for the nput of Oracle object.
adata.X = adata.layers["raw_count"].copy()

# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name=cell_group,
                                   embedding_name="X_umap")

# You can load TF info dataframe with the following code.
oracle.import_TF_data(TF_info_matrix=base_GRN)

# Perform PCA
oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)

# Save oracle object.
oracle.to_hdf5(f"data/{cell_name}.celloracle.oracle")

# Load file.
oracle = co.load_hdf5(f"data/{cell_name}.celloracle.oracle")

# Check clustering data
sc.pl.scatter(oracle.adata, color="level_3", basis="umap")

# Calculate GRN for each population in "louvain_annot" clustering unit.
# This step may take some time.(~30 minutes)
links = oracle.get_links(cluster_name_for_GRN_unit="level_3", alpha=10,
                         verbose_level=10)


links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
plt.rcParams["figure.figsize"] = [9, 4.5]
links.plot_degree_distributions(plot_model=True,
                                               #save=f"{save_folder}/degree_distribution/",
                                               )

plt.rcParams["figure.figsize"] = [10, 8]
# Calculate network scores.
links.get_network_score()
links.merged_score.head()

# Save Links object.
links.to_hdf5(file_path=f"data/{cell_name}.links.celloracle.links")

# You can load files with the following command.
links = co.load_hdf5(file_path=f"data/{cell_name}.links.celloracle.links")

# Check cluster name
print(links.cluster)

# Visualize Gata2 network score dynamics
links.plot_score_per_cluster(goi="SOX4", save=f"{save_folder}/{cell_name}.network_score_per_gene/")





