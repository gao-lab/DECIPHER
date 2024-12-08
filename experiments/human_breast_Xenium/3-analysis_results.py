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

# %%
from spider import Spider
import scanpy as sc

# %%
model = Spider(work_dir='./results/spider')

# %%
model.visualize()

# %% [markdown]
# ## Analysis results
# ### T cells

# %%
embedding = sc.read_h5ad('./results/spider/embedding.h5ad')

# %%
adata.obsm.update(embedding.obsm)

# %%
adata

# %%
tcell = adata[adata.obs['cell_type'] == 'Tcell'].copy()

# %%
tcell = scanpy_viz(tcell, n_neighbors=15, resolution=0.3)

# %% [markdown]
# #### Neighbor

# %%
tcell.obsm['X_umap'] = tcell.obsm['X_umap_nbr'].copy()

# %%
sc.pl.umap(tcell, color = ['CD4', 'CD8A', 'CCR7' ,'TCF7'], vmax='3')

# %%
sc.pl.umap(tcell, color = ['HAVCR2', 'FOXP3', 'LAG3', 'TIGIT', 'CTLA4'], vmax=1)

# %%
sc.pl.dotplot(tcell, ['HAVCR2', 'FOXP3'], groupby='leiden_nbr')
sc.pl.dotplot(tcell, ['LAG3', 'TIGIT', 'CTLA4'], groupby='leiden_nbr')

# %%
# sc.pl.violin(tcell, ['HAVCR2', 'PDCD1'], groupby='leiden_nbr', size=0)

# %%
sc.set_figure_params(dpi=300)
sc.pl.spatial(tcell, color='leiden_nbr', spot_size=embedding.uns['spot_size']*1.5)

# %%
sc.pl.umap(tcell, color='leiden_nbr')

# %%
tcell.obs['Environment'] = 'unknown'
tcell.obs.loc[tcell.obs['leiden_nbr'].isin(['0']), 'Environment'] = 'Tumor-invasive'
tcell.obs.loc[tcell.obs['leiden_nbr'].isin(['3']), 'Environment'] = 'Tumor-DICS-infiltated'
tcell.obs.loc[tcell.obs['leiden_nbr'].isin(['1', '2', '4']), 'Environment'] = 'Tumor-DICS'
tcell.obs.loc[tcell.obs['leiden_nbr'].isin(['5']), 'Environment'] = 'Non-infiltrated'

# %%
sc.pl.umap(tcell, color='Environment')

# %%
sc.set_figure_params(dpi=300)
sc.pl.spatial(tcell, color='Environment', spot_size=embedding.uns['spot_size']*1.5)

# %%
sc.pl.dotplot(tcell, ['HAVCR2', 'FOXP3', 'LAG3', 'TIGIT', 'CTLA4'], groupby='Environment')
# sc.pl.dotplot(tcell, ['TIGIT', 'CTLA4'], groupby='Environment')

# %% [markdown]
# #### Gex

# %%
tcell.obsm['X_umap'] = tcell.obsm['X_umap_center'].copy()
# sc.pl.umap(tcell, color = ['CD4', 'CD8A', 'CCR7' ,'TCF7', 'GZMA', 'CD83', 'HAVCR2', 'PDCD1'], vmax='3')

# %%
sc.set_figure_params(dpi=150)
sc.pl.umap(tcell, color = ['CD4', 'CD8A', 'GZMA' ,'TCF7', 'FOXP3', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4'], vmax='3', ncols=3)

# %%
sc.pl.dotplot(tcell, ['FOXP3', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4'], groupby='Subtype')

# %%
sc.tl.embedding_density(tcell, basis='umap', groupby='Subtype')
sc.pl.embedding_density(
    tcell, basis='umap', key='umap_density_Subtype', group='TCD4+'
)

# %%
tcell.obs['Subtype'] = 'unknown'
tcell.obs.loc[tcell.obs['leiden_center'].isin(['0', '1', '2']), 'Subtype'] = 'TCD4+'
tcell.obs.loc[tcell.obs['leiden_center'].isin(['3']), 'Subtype'] = 'TCD8+'

# %%
sc.pl.umap(tcell, color='Subtype')

# %% [markdown]
# #### Merge
# Will dominate by the cell embedding

# %%
tcell.obsm['X_umap'] = tcell.obsm['X_umap_merge'].copy()
sc.pl.umap(tcell, color = ['CD4', 'CD8A', 'GZMA' ,'TCF7', 'HAVCR2', 'FOXP3', 'LAG3', 'TIGIT'], vmax='3')

# %%
sc.pl.umap(tcell, color='leiden_merge')

# %% [markdown]
# #### Union Gex and Neighbor

# %%
tcell.obs.groupby(['Environment', 'Subtype']).size().unstack()

# %% [markdown]
# #### Save

# %%
tcell.write_h5ad('results/tcell.h5ad')

# %%
tcell = sc.read_h5ad('./results/tcell.h5ad')

# %% [markdown]
# ### B cell

# %%
bcell = adata[adata.obs['cell_type'] == 'Bcell'].copy()

# %%
bcell = scanpy_viz(bcell, n_neighbors=15, resolution=0.3)

# %%
bcell.obsm['X_umap'] = bcell.obsm['X_umap_nbr'].copy()

# %%
sc.pl.umap(bcell, color = ['CD27', 'CD79A', 'MZB1'], vmax=3)

# %%
sc.pl.umap(bcell, color = ['leiden_nbr'])

# %%
sc.pl.spatial(bcell, color='leiden_nbr', spot_size=embedding.uns['spot_size'])
