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
# # Xeneium breast dataset explain analysis

# %%
import pandas as pd
import scanpy as sc
import numpy as np
import torch
import seaborn as sns
import matplotlib.pyplot as plt

from decipher.utils import scanpy_viz, gex_embedding


# %%
sc.set_figure_params(dpi=120)

# %% [markdown]
# ## Load data

# %%
adata = sc.read_h5ad('./data/adata.h5ad')
print(adata)
print((adata.X >0).sum(1).mean())


center_emb = np.load('./results/decipher_explain/center_emb.npy')
nbr_emb = np.load('./results/decipher_explain/nbr_emb.npy')

adata.obsm['X_center'] = center_emb
adata.obsm['X_nbr'] = nbr_emb


# %%
# adata = scanpy_viz(adata, keys=['pred_nbr', 'center', 'nbr'])
# adata

# %% [markdown]
# ## Select genes

# %%
gene_mask = np.load('./results/decipher_explain/explain/select_celltype:Bcell/gene_mask.npy')
print(gene_mask.shape)

gene_df_bcell = adata.var
gene_df_bcell['cell_type'] = 'B cell'
gene_df_bcell['counts'] = gene_mask.sum(0)
# sort by counts
gene_df_bcell = gene_df_bcell.sort_values('counts', ascending=False)
gene_df_bcell['gene'] = gene_df_bcell.index
gene_df_bcell['score'] = gene_df_bcell['counts'] / gene_mask.shape[0]
gene_df_bcell[['counts', 'score']].head(10)

# %%
gene_mask = np.load('./results/decipher_explain/explain/select_celltype:Tcell/gene_mask.npy')

gene_df_tcell = adata.var
gene_df_tcell['counts'] = gene_mask.sum(0)
gene_df_bcell['cell_type'] = 'T cell'
# sort by counts
gene_df_tcell = gene_df_tcell.sort_values('counts', ascending=False)
gene_df_tcell['gene'] = gene_df_tcell.index
gene_df_tcell['score'] = gene_df_tcell['counts'] / gene_mask.shape[0]
gene_df_tcell[['counts', 'score']].head(10)


# %%
# set figure size as (4, 6)
plt.figure(figsize=(4, 2.5))
sns.barplot(data=gene_df_tcell.head(5), x='score', y='gene')
plt.xlabel('Score of T cell')
plt.ylabel('Gene')
plt.show()

plt.figure(figsize=(4, 2.5))
sns.barplot(data=gene_df_bcell.head(5), x='score', y='gene')
plt.xlabel('Score of B cell')
plt.ylabel('Gene')

# %%
gene_df.head(20).index.to_frame().to_csv('./gene_list_top.csv', index=False)
gene_df.index.to_frame().to_csv('./gene_list.csv', index=False)

# %%
gene_df['counts'].hist(bins=100)

# %% [markdown]
# ## T cell analysis

# %%
subset = (adata.obs["cell_type"] == "Tcell").values
tcell = adata[subset].copy()

# %%
sc.set_figure_params(dpi=200, transparent=True, vector_friendly=True)
sc.pl.spatial(tcell, color='PTGDS', spot_size=40, legend_loc='on data', vmax=5)

# %%
sc.set_figure_params(dpi=200, transparent=True, vector_friendly=True)
sc.pl.spatial(tcell, color='CXCL12', spot_size=40, legend_loc='on data', vmax=10)

# %%
sc.pl.spatial(tcell, color='LYZ', spot_size=40, legend_loc='on data', vmax=10)

# %%
sc.pl.spatial(tcell, color=['CXCR4'], spot_size=50, legend_loc='on data', ncols =2)


# %%
sc.pl.spatial(tcell, color=['SFRP4'], vmax=10, spot_size=50, legend_loc='on data', ncols =2)


# %%
tcell = scanpy_viz(tcell, keys=['center', 'nbr'])

# %%
tcell.obsm['X_umap'] = tcell.obsm['X_umap_nbr'].copy()
# sc.pl.umap(tcell, color = gene_df.head(10).index.to_list(), vmax=2)
sc.pl.umap(tcell, color = ['PTGDS', 'CXCL12'], vmax=2)

# %%
tcell.obsm['X_umap'] = tcell.obsm['X_umap_center'].copy()
# sc.pl.umap(tcell, color = gene_df.head(10).index.to_list(), vmax=2)
sc.pl.umap(tcell, color = ['PTGDS', 'CXCL12'], vmax=2)

# %%
sc.pl.umap(tcell, color = ['PDCD1', 'HAVCR2', 'TIGIT', 'APOBEC3B'], vmax=2)

# %% [markdown]
# ## B cell analysis

# %%
subset = (adata.obs["cell_type"] == "Bcell").values
bcell = adata[subset].copy()

# %%
# sc.pl.spatial(bcell, color=['PTGDS', 'CXCL12', 'LYZ'], spot_size=40, legend_loc='on data', vmax=5, ncols = 1)
sc.pl.spatial(bcell, color=['LYZ'], spot_size=40, legend_loc='on data', vmax=10, ncols = 1)

# %%
sc.pl.spatial(bcell, color=['CD86'], spot_size=40, legend_loc='on data', ncols = 1, vmax = 2)

# %%
bcell = scanpy_viz(bcell, keys=['center', 'nbr'])

# %%
bcell.obsm['X_umap'] = bcell.obsm['X_umap_nbr'].copy()
# sc.pl.umap(bcell, color = gene_df.head(10).index.to_list(), vmax=2)
sc.pl.umap(bcell, color = ['PTGDS', 'LYZ', 'LUM'], vmax=2)

# %%
bcell.obsm['X_umap'] = bcell.obsm['X_umap_center'].copy()
# sc.pl.umap(bcell, color = gene_df.head(10).index.to_list(), vmax=2)
sc.pl.umap(bcell, color = ['PTGDS', 'LYZ', 'LUM'], vmax=2)

# %% [markdown]
# ## Check in whole spatial dataset

# %%
sc.pl.spatial(adata, spot_size=12, color=['PTGDS', 'LUM', 'LYZ'], ncols=1, vmax=5)

# %%
sc.pl.spatial(bcell, spot_size=40, color=['LYZ', 'LUM'], ncols=1, vmax=5)

# %% [markdown]
# ## Check in paired 3' single cell data

# %%
adata_sc = sc.read_10x_h5('./data/Chromium_FFPE_Human_Breast_Cancer_Chromium_FFPE_Human_Breast_Cancer_count_sample_filtered_feature_bc_matrix.h5')
adata_sc

# %%
adata_sc = gex_embedding(adata_sc, n_top_genes=3000)

# %%
sc.pl.umap(adata_sc, color = ['PTGDS', 'LYZ', 'LUM', 'CD79A','CD8A'], vmax=3, ncols=2)

# %%
adata_sc.var[adata_sc.var.index.str.startswith('TR')]

# %%
