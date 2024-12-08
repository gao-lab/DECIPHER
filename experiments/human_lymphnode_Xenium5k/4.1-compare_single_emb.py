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
#     display_name: conda
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Compare with simple embedding in find LR
#
# We use Spider spatial embedding to define spatial cluster and then find the different ligand-receptors across the clusters.

# %%
import os

import scanpy as sc
import pandas as pd
import torch
import numpy as np
from torch_geometric.nn import SimpleConv

from spider.explain.gene.lr import get_lr_expr
from spider.graphic.build import build_graph
from spider.utils import scanpy_viz, gex_embedding


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
# ## Spider embbeding

# %%
adata.obsm['X_nbr'] =  np.load('./results/spider/nbr_emb.npy')
adata = scanpy_viz(adata, keys = ['nbr'])

# %% [markdown]
# ## Find differential expressed LR across spatial domains

# %%
lr.index= lr['interaction_name']
adata_lr = sc.AnnData(
    X = expr.numpy(),
    obs = adata.obs,
    obsm = adata.obsm,
    var = lr,
)
bcell = adata_lr[adata_lr.obs['cell_type'] == 'Bcell'].copy()

# %%
sc.pp.normalize_total(bcell, target_sum=1e4)
sc.pp.log1p(bcell)

# %%
sc.tl.rank_genes_groups(bcell, "leiden_nbr", method="wilcoxon")
sc.pl.rank_genes_groups(bcell, n_genes=10, sharey=False)

# %%
