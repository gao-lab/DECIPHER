# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Run other spatial methods in downsampled mouse brain 3D atlas
#
# We compare with following methods:
#
# - Banksy
# - SLAT
# - STAGTE

# %%
import scanpy as sc
import pandas as pd
import numpy as np
import torch
import os
from torch_geometric.nn import SimpleConv

from decipher.graphic.build import build_graph
from decipher.utils import scanpy_viz, gex_embedding

# %% [markdown]
# ## Downsample data

# %%
adata = sc.read_h5ad('./data/zhuang_dataset/abca1.h5ad')

# %%
subsample_idx = np.random.choice(adata.shape[0], 200_000, replace=False)
adata_sub = adata[subsample_idx, :].copy()

# %%
np.save('./results/other_methods_compare/subsample_idx.npy', subsample_idx)

# %%
adata_sub.write_h5ad('./results/other_methods_compare/adata_sub_200k.h5ad')

# %%
adata_sub.layers['counts'] = adata_sub.X.copy()

# %%
adata_sub.obsm['X_spatial_3d'] = adata_sub.obs[['x_ccf', 'y_ccf', 'z_ccf']].values

# %% [markdown]
# ## SLAT

# %%
adata_sub.X = adata_sub.layers["counts"].copy()
adata_sub = gex_embedding(adata_sub, rapids_after_scale=True, viz=False)

# %%
edge_index = build_graph(adata_sub.obsm['X_spatial_3d'], k=24, mode='knn')
torch.save(edge_index, './results/other_methods_compare/edge_index_k24.pt')

# %%
x = torch.from_numpy(adata_sub.obsm['X_pca'])
gcn = SimpleConv(aggr="mean")
slat_emb = gcn(x, edge_index)
adata_sub.obsm['X_slat'] = slat_emb.detach().numpy()

os.makedirs('results/other_methods_compare/slat', exist_ok=True)
np.save('./results/other_methods_compare/slat/slat_emb.npy', slat_emb)

# %%
adata_sub = scanpy_viz(adata_sub, keys = ['slat'])

# %%
adata_sub

# %%
adata_sub.obsm['X_umap'] = adata_sub.obsm['X_umap_slat'].copy()
sc.pl.umap(adata_sub, color=['class', 'parcellation_division'], ncols=1)

# %% [markdown]
# ## STAGATE

# %%
import STAGATE_pyG


# %%
# Spatial graph 
# STAGATE_pyG.Cal_Spatial_Net(adata_sub, k_cutoff=24, model='KNN')

# %% [markdown]
# STAGATE do not suport 3D data, we need build 3d graph by ourself.

# %%
Spatial_Net = pd.DataFrame(
    {'Cell1' : adata_sub.obs.index[edge_index[0]],
     'Cell2' : adata_sub.obs.index[edge_index[1]],
     'Distance' : edge_index.shape[1]*[0.1]
     }
    )
adata_sub.uns['Spatial_Net'] = Spatial_Net

# %%
adata_sub.X = adata_sub.layers["counts"].copy()

# %%
# Preprocessing
sc.pp.normalize_total(adata_sub, target_sum=1e4)
sc.pp.log1p(adata_sub)
sc.pp.highly_variable_genes(adata_sub)

# Train model
adata_sub = STAGATE_pyG.train_STAGATE(adata_sub)

# %%
os.makedirs('./results/other_methods_compare/stagate', exist_ok=True)
np.save('./results/other_methods_compare/stagate/stagate_emb.npy', adata_sub.obsm['STAGATE'])

# %%
adata_sub.obsm['X_stagate'] = np.load('./results/other_methods_compare/stagate/stagate_emb.npy')
adata_sub = scanpy_viz(adata_sub, keys = ['stagate'])

# %%
adata_sub.obsm['X_umap'] = adata_sub.obsm['X_umap_stagate'].copy()
sc.pl.umap(adata_sub, color=['class', 'parcellation_division'], ncols=1)

# %% [markdown]
# ## Banksy

# %%
import time
import random
from banksy.main import median_dist_to_nearest_neighbour
from banksy.banksy_utils.filter_utils import normalize_total, filter_hvg
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.run_banksy import run_banksy_multiparam

# %%
os.makedirs('results/other_methods_compare/banksy', exist_ok=True)

adata_sub.X = adata_sub.layers["counts"].copy()

coord_keys = ['x', 'y', 'spatial']
adata_sub.obs['x'] = adata_sub.obsm['spatial'][:, 0]
adata_sub.obs['y'] = adata_sub.obsm['spatial'][:, 1]

# Preprocessing
adata_sub = normalize_total(adata_sub)
adata_sub, adata_allgenes = filter_hvg(adata_sub,
                                   n_top_genes=2000,
                                   flavor="seurat")

# set
plot_graph_weights = True
k_geom = 15  # only for fixed type
max_m = 1  # azumithal transform up to kth order
nbr_weight_decay = "scaled_gaussian"  # can also be "reciprocal", "uniform" or "ranked"

# Find median distance to closest neighbours, the median distance will be `sigma`
nbrs = median_dist_to_nearest_neighbour(adata_sub, key=coord_keys[2])



banksy_dict = initialize_banksy(adata_sub,
                                coord_keys,
                                k_geom,
                                nbr_weight_decay=nbr_weight_decay,
                                max_m=max_m,
                                plt_edge_hist=False,
                                plt_nbr_weights=False,
                                plt_agf_angles=False
                                )


# The following are the main hyperparameters for BANKSY
resolutions = [0.7]  # clustering resolution for UMAP
pca_dims = [20]  # Dimensionality in which PCA reduces to
lambda_list = [0.2]  # list of lambda parameters

banksy_dict, banksy_matrix = generate_banksy_matrix(adata_sub,
                                                    banksy_dict,
                                                    lambda_list,
                                                    max_m)
number_of_colors = 100
colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

results_df = run_banksy_multiparam(
    adata_sub,
    banksy_dict,
    lambda_list,
    resolutions,
    color_list=colors,
    max_m=max_m,
    filepath='/tmp',
    key=coord_keys,
    pca_dims=pca_dims,
    annotation_key=None,
    max_labels=None,
    cluster_algorithm='leiden',
    match_labels=False,
    savefig=False,
    add_nonspatial=True,
    variance_balance=False
)


banksy_emb = results_df['adata'][1].obsm['reduced_pc_20'].copy()
np.save('./results/other_methods_compare/banksy/banksy_emb.npy', banksy_emb)


# %%
adata_sub.obsm['X_banksy'] = np.load('./results/other_methods_compare/banksy/banksy_emb.npy')
adata_sub = scanpy_viz(adata_sub, keys = ['banksy'])

# %%
adata_sub.obsm['X_umap'] = adata_sub.obsm['X_umap_banksy'].copy()
sc.pl.umap(adata_sub, color=['class', 'parcellation_division'], ncols=1)
