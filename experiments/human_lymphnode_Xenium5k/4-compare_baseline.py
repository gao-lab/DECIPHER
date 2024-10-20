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
# # Compare with simple baseline in find LR
#
# We run with Banksy/UTAG to define spatial cluster and then find the different ligand-receptors across the clusters.

# %%
import os

import scanpy as sc
import pandas as pd
import torch
import numpy as np
from torch_geometric.nn import SimpleConv

from decipher.explain.lr import get_lr_expr
from decipher.graphic.build import build_graph
from decipher.utils import scanpy_viz, gex_embedding


# %%
adata = sc.read_h5ad('./data/lymph_node.h5ad')
lr = pd.read_csv('../../resource/lr_data/cellchat.csv')

# %% [markdown]
# ## Ligand-receptor pairs

# %%
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# %%
expr, lr = get_lr_expr(adata, lr, 20)

# %%
np.save('./results/lr_baseline/lr_expr.npy', expr)
lr.to_csv('./results/lr_baseline/lr_filter.csv')

# %% [markdown]
# ## SLAT embedding

# %%
adata.X = adata.layers["counts"].copy()
adata = gex_embedding(adata, rapids_after_scale=True, viz=False)

# %%
edge_index = build_graph(adata.obsm['spatial'], mode='knn', k=20)

# %%
x = torch.from_numpy(adata.obsm['X_pca'])
gcn = SimpleConv(aggr="mean")
slat_emb = gcn(x, edge_index)
adata.obsm['X_slat'] = slat_emb.detach().numpy()

os.makedirs('results/slat', exist_ok=True)
np.save('./results/slat/slat_emb.npy', slat_emb)

# %%
adata = scanpy_viz(adata, keys = ['slat'])

# %% [markdown]
# ## Banksy embedding

# %%
import time
import random
from banksy.main import median_dist_to_nearest_neighbour
from banksy.banksy_utils.filter_utils import normalize_total, filter_hvg
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.run_banksy import run_banksy_multiparam

# %%
os.makedirs('results/banksy', exist_ok=True)

# %%
start = time.time()
adata.X = adata.layers["counts"].copy()

coord_keys = ['x', 'y', 'spatial']
adata.obs['x'] = adata.obsm['spatial'][:, 0]
adata.obs['y'] = adata.obsm['spatial'][:, 1]

# Preprocessing
adata = normalize_total(adata)
adata, adata_allgenes = filter_hvg(adata,
                                   n_top_genes=2000,
                                   flavor="seurat")

# set
plot_graph_weights = True
k_geom = 15  # only for fixed type
max_m = 1  # azumithal transform up to kth order
nbr_weight_decay = "scaled_gaussian"  # can also be "reciprocal", "uniform" or "ranked"

# Find median distance to closest neighbours, the median distance will be `sigma`
nbrs = median_dist_to_nearest_neighbour(adata, key=coord_keys[2])



banksy_dict = initialize_banksy(adata,
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

banksy_dict, banksy_matrix = generate_banksy_matrix(adata,
                                                    banksy_dict,
                                                    lambda_list,
                                                    max_m)
number_of_colors = 100
colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

results_df = run_banksy_multiparam(
    adata,
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


run_time = str(time.time() - start)
print('Runtime: ' + run_time)

banksy_emb = results_df['adata'][1].obsm['reduced_pc_20'].copy()
np.save('./results/banksy/banksy_emb.npy', banksy_emb)

# %% [markdown]
# ## Define spatial domain based on embeddings

# %%
adata.obsm['X_slat'] = np.load('./results/slat/slat_emb.npy')
adata.obsm['X_banksy'] = np.load('./results/banksy/banksy_emb.npy')

# %%
adata = scanpy_viz(adata, keys = ['banksy', 'slat'], resolution=0.3)

# %%
adata.write_h5ad('./results/lr_baseline/slat_banksy.h5ad')

# %%
sc.set_figure_params(dpi=150)
# sc.pl.spatial(adata, color=['leiden_banksy'], spot_size=6, ncols=1)
# sc.pl.spatial(adata, color=['leiden_slat'], spot_size=6, ncols=1)

# %% [markdown]
# ## Find differential expressed LR across spatial domains

# %%
lr = pd.read_csv('./results/lr_baseline/lr_filter.csv')
lr = pd.read_csv('./results/lr_baseline/lr_filter.csv')

# %%
lr.index= lr['interaction_name']
adata_lr = sc.AnnData(
    X = expr,
    obs = adata.obs,
    obsm = adata.obsm,
    var = lr,
)
bcell = adata_lr[adata_lr.obs['cell_type'] == 'Bcell'].copy()

# %%
sc.pp.normalize_total(bcell, target_sum=1e4)
sc.pp.log1p(bcell)

# %%
# sc.pl.violin(bcell, keys=['CXCL12_CXCR4', 'CXCL13_CXCR5', 'FCER2A_CR2'], size = 0, rotation=45, ylabel='Activity')
sc.pl.violin(bcell, keys=['CXCL12_CXCR4', 'CXCL13_CXCR5'], size = 0, rotation=45, ylabel='Activity')

# %%
sc.tl.rank_genes_groups(bcell, "leiden_banksy", method="wilcoxon")
sc.pl.rank_genes_groups(bcell, n_genes=10, sharey=False)

# %%
# filter leiden_slat > 100
cluster = bcell.obs['leiden_slat'].value_counts().index[bcell.obs['leiden_slat'].value_counts() > 100]
bcell = bcell[bcell.obs['leiden_slat'].isin(cluster.to_list())]

# %%
bcell = bcell[bcell.obs['leiden_slat'].isin(cluster.to_list())]

# %%
sc.tl.rank_genes_groups(bcell, "leiden_slat", method="wilcoxon")
sc.pl.rank_genes_groups(bcell, n_genes=10, sharey=False)
