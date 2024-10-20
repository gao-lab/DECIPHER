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
# # Analysis mouse brain MEFISH dataset results
# > We use rapids single cell as backend

# %%
from pathlib import Path

import scanpy as sc
import numpy as np
import rapids_singlecell as rsc
import pandas as pd

from decipher.utils import scanpy_viz, clip_umap, manage_gpu

# %%
sc.set_figure_params(dpi=120)

# %%
adata = sc.read_h5ad("./data/zhuang_dataset/abca_processed.h5ad")
adata = adata[adata.obs["feature_matrix_label"].isin(["Zhuang-ABCA-1"])].copy()
adata

# %% [markdown]
# ## Embedding visualization

# %%
center_emb_3d = np.load('./results/decipher_abca-1_3d_0705/center_emb.npy')
nbr_emb_3d = np.load('./results/decipher_abca-1_3d_0705/nbr_emb.npy')
center_emb_2d = np.load('./results/decipher_abca-1_2d_0705/center_emb.npy')
nbr_emb_2d = np.load('./results/decipher_abca-1_2d_0705/nbr_emb.npy')

# %%
adata.obsm['X_center_3d'] = center_emb_3d
adata.obsm['X_nbr_3d'] = nbr_emb_3d
adata.obsm['X_center_2d'] = center_emb_2d
adata.obsm['X_nbr_2d'] = nbr_emb_2d

# %%
adata = scanpy_viz(adata, keys=['center_3d', 'nbr_3d', 'center_2d', 'nbr_2d'])

# %%
color_by = ['class', 'parcellation_division']

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_center_2d'].copy()
sc.pl.umap(adata, color=color_by, ncols=1, wspace=0.3)

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_center_3d'].copy()
sc.pl.umap(adata, color=color_by, ncols=1, wspace=0.3)

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_nbr_2d'].copy()
adata.obsm['X_umap'] = clip_umap(adata.obsm['X_umap'], 0.05)
sc.pl.umap(adata, color=color_by, ncols=1, wspace=0.3)

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_nbr_3d'].copy()
adata.obsm['X_umap'] = clip_umap(adata.obsm['X_umap'])
sc.pl.umap(adata, color=color_by, ncols=1, wspace=0.3)

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_center_3d'].copy()
sc.pl.umap(adata, color=color_by, ncols=1, wspace=0.3)

# %%
vascular = adata[adata.obs['class'] == '33 Vascular']

# %%
vascular.obsm['X_umap'] = vascular.obsm['X_umap_nbr_1'].copy()
sc.pl.umap(vascular, color=['class', 'parcellation_division', 'feature_matrix_label'], ncols=1, wspace=0.3)

# %% [markdown]
# ## Cell type density

# %%
N_SAMPLE = 500_000
adata_dense = sc.AnnData(X=np.arange(adata.shape[0],dtype=np.float32).reshape(-1, 1) , obs=adata.obs, obsm=adata.obsm)
adata_proc_sampled = adata_dense[np.random.choice(adata_dense.shape[0], N_SAMPLE, replace=False), :].copy()


# %%
manage_gpu(7, 'large')
# 4 mins for 0.5 M cells, 8 mins for  0.8 M cells
rsc.get.anndata_to_GPU(adata_proc_sampled)
rsc.tl.embedding_density(adata_proc_sampled, basis='umap', groupby='class')
rsc.get.anndata_to_CPU(adata_proc_sampled)

# %%
adata_proc_sampled.obsm['X_umap'] = adata_proc_sampled.obsm['X_umap_nbr_3d'].copy()
adata_proc_sampled.obsm['X_umap'] = clip_umap(adata_proc_sampled.obsm['X_umap'])
sc.pl.embedding_density(adata_proc_sampled, basis='umap', groupby='class', bg_dotsize=1, fg_dotsize=3, ncols=7)

# %%
adata_proc_sampled.obsm['X_umap'] = adata_proc_sampled.obsm['X_umap_nbr_2d'].copy()
adata_proc_sampled.obsm['X_umap'] = clip_umap(adata_proc_sampled.obsm['X_umap'])
sc.pl.embedding_density(adata_proc_sampled, basis='umap', groupby='class', bg_dotsize=1, fg_dotsize=3, ncols=7)
