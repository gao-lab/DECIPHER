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
# # Xenium 5k results analysis

# %%
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import rapids_singlecell as rsc

from spider.utils import scanpy_viz, manage_gpu, clip_umap

# %%
adata = sc.read_h5ad('./data/lymph_node.h5ad')

center_emb = np.load('./results/spider/center_emb.npy')
nbr_emb = np.load('./results/spider/nbr_emb.npy')

adata.obsm['X_center'] = center_emb
adata.obsm['X_nbr'] = nbr_emb
adata

# %% [markdown]
# ## analysis embedding

# %%
adata = scanpy_viz(adata, keys=['center', 'nbr'])
adata

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_center'].copy()
adata.obsm['X_umap'] = clip_umap(adata.obsm['X_umap'], 0.15)
sc.pl.umap(adata, color='cell_type')

# %%
sc.pl.umap(adata, color = ['CD3E', 'CD79A', 'PLVAP', ])

# %%
sc.pl.umap(adata, color = ['CXCR4', 'CXCR5', 'CXCL12', 'CXCL13', 'CR2', 'FCER2'])

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_nbr'].copy()
adata.obsm['X_umap'] = clip_umap(adata.obsm['X_umap'], 0.01)
sc.pl.umap(adata, color=['cell_type', 'leiden_nbr'], wspace=0.5)

# %%
sc.pl.umap(adata, color = ['CXCR4', 'CXCR5', 'CXCL12', 'CXCL13', 'CR2', 'FCER2'])

# %% [markdown]
# ## analysis gene selection

# %%
gene_mask = np.load('./results/spider/explain/select_batch:None_celltype:Bcell/gene_mask.npy')
print(gene_mask.shape)

lr_df = pd.read_csv('./results/spider/explain/lr.csv')

lr_df['counts'] = gene_mask.sum(0)
# sort by counts
lr_df = lr_df.sort_values('counts', ascending=False)
lr_df[['counts', 'interaction_name']].head(20)

lr_df['score'] = lr_df['counts'] / gene_mask.shape[0]

# %%
# bar plot of top 20 genes x is 'interaction_name' y is 'counts'
sns.barplot(x='score', y='interaction_name', data=lr_df.head(15))
# remove y label
plt.ylabel('')
# change y ticks font size to 12
_ = plt.yticks(fontsize=9)
_ = plt.xticks(fontsize=9)
# change x label font size to 10
_ = plt.xlabel('Score', fontsize=10)

# %%
gene_mask = np.load('./results/spider/explain/select_batch:None_celltype:Plasma/gene_mask.npy')
print(gene_mask.shape)

lr_df = pd.read_csv('./results/spider/explain/lr.csv')

lr_df['counts'] = gene_mask.sum(0)
# sort by counts
lr_df = lr_df.sort_values('counts', ascending=False)
lr_df[['counts', 'interaction_name']].head(20)

# %%
# filter the rows whose interaction_name colum contain 'CD40'
lr_df[lr_df['interaction_name'].str.contains('CD40')]


# %% [markdown]
# ## Visualization

# %%
sc.set_figure_params(dpi=600)
sc.pl.spatial(adata, color='cell_type', spot_size=6)
sc.set_figure_params(dpi=100)

# %%
# sc.set_figure_params(dpi=500, dpi_save=500, vector_friendly=True, transparent=True, figsize=(20,20))
# sc.pl.spatial(adata, color='leiden_nbr', spot_size=4, save='lymph_node_leiden_nbr.pdf')
sc.pl.spatial(adata, color='leiden_nbr', spot_size=5)
sc.set_figure_params(dpi=120)

# %%
adata.obsm['X_spatial'] = adata.obsm['spatial'].copy()
manage_gpu(7, 'large')
# 4 mins for 0.5 M cells, 8 mins for  0.8 M cells
rsc.get.anndata_to_GPU(adata)
rsc.tl.embedding_density(adata, basis='spatial', groupby='cell_type')
rsc.get.anndata_to_CPU(adata)

# %%
bcell = adata[adata.obs['cell_type'] == 'Bcell'].copy()

# %%
sc.set_figure_params(dpi=200)
lr_genes = ['CXCR4', 'CXCL12', 'CXCR5', 'CXCL13', 'CR2', 'FCER2']
sc.pl.spatial(bcell, color=lr_genes, spot_size=8, cmap='viridis', ncols=2)
sc.set_figure_params(dpi=120)

# %% [markdown]
# ## Save results

# %%
adata.write_h5ad('./results/adata_analysis.h5ad')

# %%
