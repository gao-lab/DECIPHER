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
# # Analysis pan-cancer MERSCOPE dataset results
# > Run with hpc A100 80G
#
# > We use rapids single cell as backend
#

# %%
from pathlib import Path
# import python reload 
from importlib import reload

import dask
dask.config.set({'dataframe.query-planning': False})

import scanpy as sc
import rapids_singlecell as rsc
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt

import decipher
from decipher.utils import scanpy_viz, clip_umap, manage_gpu, gex_embedding

# reload scanpy_viz function
reload(decipher.utils)

# %%
def split_umap(adata, split_by, ncol=1, nrow=None, **kwargs):
    r"""
    Split umap by split_by like Seurat

    Parameters
    ----------
    adata
        AnnData object
    split_by
        split by variable
    ncol
        number of columns
    nrow
        number of rows
    **kwargs
        other parameters for `sc.pl.umap`
    """
    categories = adata.obs[split_by].cat.categories
    # print(categories)
    if nrow is None:
        nrow = int(np.ceil(len(categories) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4 * nrow))
    axs = axs.flatten()
    for i, cat in enumerate(categories):
        ax = axs[i]
        sc.pl.umap(adata[adata.obs[split_by] == cat], ax=ax, show=False, title=cat, **kwargs)
    plt.tight_layout()


# %%
adata_path = "./data/pancancer_filter_anno.h5ad"
adata = sc.read_h5ad(adata_path)

# %%
adata.obs.head()

# %% [markdown]
# # RAPIDS analysis

# %%
embs = Path('./results/decipher_6_10/model/lightning_logs/version_3').glob('*.npy')
for file in embs:
    if 'umap' in str(file):
        continue
    emb = np.load(file)
    print(file)
    epoch = file.stem.split('-')[1].split('_')[0]
    mode = file.stem.split('-')[1].split('_')[1]
    print(epoch)
    print(mode)
    adata.obsm[f'X_{mode}_{epoch}'] = emb

# %%
adata.obs['cell_type'].value_counts()

# %%
adata = scanpy_viz(adata, keys=['nbr_0'], approx = True)

# %%
np.save('./results/decipher_6_10/model/lightning_logs/version_3/umap_gex_0.npy', adata.obsm['X_umap_gex_0'])
np.save('./results/decipher_6_10/model/lightning_logs/version_3/umap_nbr_0.npy', adata.obsm['X_umap_nbr_0'])

# %%
adata.obsm['X_umap_gex_0'] = np.load('./results/decipher_6_10/model/lightning_logs/version_3/umap_gex_0.npy')
adata.obsm['X_umap_nbr_0'] = np.load('./results/decipher_6_10/model/lightning_logs/version_3/umap_nbr_0.npy')

# %%
adata_plot = adata[~adata.obs['cell_type'].isin(['unknown', 'Hepatocyte', 'pneumocyte']), :]

# %%
adata_plot.obsm['X_umap'] = adata_plot.obsm['X_umap_gex_0'].copy()
adata_plot.obsm['X_umap'] = clip_umap(adata_plot.obsm['X_umap'],  0.001)
sc.pl.umap(adata_plot, color=['cell_type'], ncols=1, legend_loc='on data', legend_fontsize=6)

# %%
adata_plot.obsm['X_umap'] = adata_plot.obsm['X_umap_gex_0'].copy()
adata_plot.obsm['X_umap'] = clip_umap(adata_plot.obsm['X_umap'],  0.001)
sc.pl.umap(adata_plot, color=['Tissue_1'], ncols=1)

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_nbr_0'].copy()
sc.pl.umap(adata, color=['leiden_nbr_0'], ncols=1, wspace=0.3, legend_loc='on data')

# %%
adata_plot.obsm['X_umap'] = adata_plot.obsm['X_umap_nbr_0'].copy()
sc.pl.umap(adata_plot, color=['cell_type', 'Tissue_1'], ncols=1, wspace=0.3)

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_nbr_0'].copy()
adata_proc = sc.AnnData(X=np.arange(adata.shape[0],dtype=np.float32).reshape(-1, 1) , obs=adata.obs, obsm=adata.obsm)

# random sample 10,000 cells
adata_proc_sampled = adata_proc[np.random.choice(adata_proc.shape[0], 500_000, replace=False), :].copy()
# remove the 'unknown', 'Hepatocyte' and 'pneumocyte' cells
adata_proc_sampled = adata_proc_sampled[~adata_proc_sampled.obs['cell_type'].isin(['unknown', 'Hepatocyte', 'pneumocyte']), :].copy()

# %%
# 4 mins for 0.5 M cells, 8 mins for 0.8 M cells in A100-GPU
manage_gpu(2)
rsc.get.anndata_to_GPU(adata_proc_sampled)
rsc.tl.embedding_density(adata_proc_sampled, basis='umap', groupby='cell_type')
rsc.get.anndata_to_CPU(adata_proc_sampled)

# %%
sc.pl.embedding_density(adata_proc_sampled, basis='umap', groupby='cell_type', bg_dotsize=1, fg_dotsize=3, ncols = 8)

# %%
adata_plot.obs["region"] = "unknown"
adata_plot.obs.loc[adata_plot.obs["leiden_nbr_0"].isin(["1", "5"]), "region"] = "Lymph high infiltrated tumor"
adata_plot.obs.loc[adata_plot.obs["leiden_nbr_0"].isin(["0", "2", "6", "7", "9", "10"]), "region"] = "Lymph low infiltrated tumor"
adata_plot.obs.loc[adata_plot.obs["leiden_nbr_0"].isin(["3", "4", "8"]), "region"] = "Tumor core"
adata_plot.obs["region"].value_counts()

# %%
adata_plot.obsm['X_umap'] = adata_plot.obsm['X_umap_nbr_0'].copy()
sc.pl.umap(adata_plot, color=['region'], ncols=1, legend_loc='on data', legend_fontsize=7)

# %%
sc.pl.dotplot(adata_plot, var_names=['CD274'], groupby = 'region')

# %%
ax = adata_plot.obs.groupby('leiden_nbr_0')['cell_type'].value_counts(normalize=True).unstack('cell_type').plot.bar(stacked=True)
# set the legend to the right
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# %%
ax = adata_plot.obs.groupby('region')['cell_type'].value_counts(normalize=True).unstack('cell_type').plot.bar(stacked=True)
# set the legend to the right
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# remove x label
ax.set_xlabel('')
# set x ticks angle to 45
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, horizontalalignment='right', ha = 'center')

# %%
adata.obs['Tissue_0'].value_counts()

# %%
adata.write_h5ad('./results/adata_analysis.h5ad')

# %%
adata = sc.read_h5ad('./results/adata_analysis.h5ad')

# %% [markdown]
# # Sub region analysis

# %%
adata.obs['region'].value_counts()

# %%
high_inf = adata[adata.obs['region'] == 'Lymphatic infiltrated tumor', :].copy()

# %%
high_inf = scanpy_viz(high_inf, keys=['nbr_0'], approx = True)

# %%
high_inf.obsm['X_umap'] = high_inf.obsm['X_umap_nbr_0'].copy()
sc.pl.umap(high_inf, color=['cell_type'], ncols=1)

# %%
adata_plot.obs.groupby(['Tissue_0', 'Index']).size().unstack().fillna(0)


# %%
# high_inf.obs.columns
# high_inf.obs['Tissue_0'].value_counts()
sub = adata_plot[adata_plot.obs['Index'] == '2', :].copy()

# %%
sc.set_figure_params(dpi=150, fontsize=15)
sc.pl.spatial(sub, color=['region'], ncols=1, spot_size=8)

# %%
sq.gr.spatial_neighbors(sub, n_neighs=30, coord_type='generic')
sq.gr.nhood_enrichment(sub, cluster_key="cell_type",)

# %%
from matplotlib.colors import LinearSegmentedColormap

cmap_colors = [
    (0, 'blue'),  # Color for the low end (-10)
    (0.35, 'white'),   # 0 will be white
    (1, 'red')     # Color for the high end (20)
]

# Create a custom colormap
cmap = LinearSegmentedColormap.from_list('cmap', [(pos, color) for pos, color in cmap_colors])

sq.pl.nhood_enrichment(sub, cluster_key="cell_type", cmap=cmap)

# %% [markdown]
# # Subtype analysis

# %% [markdown]
# ## B cell analysis

# %%
bcell = adata_plot[adata_plot.obs['cell_type'] == 'B/Plasma cell', :].copy()
bcell

# %%
bcell_df = bcell.obs.groupby(['region', 'Tissue_0']).size().unstack().fillna(0)
# print(bcell_df)
tmp = bcell_df.iloc[0].values / bcell_df.sum(axis=0).values
bcell_df.loc['Lymphatic infiltrated tumor'] = tmp
bcell_df

# %%
sc.pl.dotplot(bcell, var_names=['HLA-DMA', 'HLA-DPB1', 'HLA-DQA1'], groupby='Tissue_0', layer = None)


# %%
sc.pl.dotplot(bcell, var_names=['HLA-DMA', 'HLA-DPB1', 'HLA-DQA1'], groupby='region')


# %% [markdown]
# ## T cell analysis

# %%
tcell = adata_plot[adata_plot.obs['cell_type'] == 'T cell', :].copy()

# %%
tcell.X = tcell.layers['counts'].copy()
tcell = gex_embedding(tcell, method='harmony', batch_key='Index', resolution=0.5, rapids_after_scale=True)
tcell.obs['leiden'].value_counts()

# %%
tcell.obs.groupby(['region', 'Tissue_0']).size().unstack().fillna(0)

# %%
tcell_marker = ['CD4', 'CD8A', 'FOXP3', 'PDCD1', 'CCR7', 'GZMK']
sc.pl.dotplot(tcell, var_names=tcell_marker, groupby='leiden')

# %%
sc.pl.umap(tcell, color=['leiden', 'PDCD1'], ncols=4, wspace=0.3)

# %%
sc.pl.dotplot(tcell, var_names=['CXCR4'], groupby='region', layer = 'counts')
sc.pl.dotplot(tcell, var_names=['CCR7'], groupby='region', layer = 'counts')


# %%
sc.pl.dotplot(tcell, var_names=['CXCR4'], groupby='Tissue_0')
sc.pl.dotplot(tcell, var_names=['CCR7'], groupby='Tissue_0')


# %%
# sc.pl.umap(tcell, color=['CCL5', 'ZAP70'], ncols=4, wspace=0.3)
# dp = sc.pl.dotplot(tcell, ['CCL5'], groupby='region', layer = 'counts')
# dp.style(dot_edge_color='black', dot_edge_lw=0.5).show()
# sc.pl.dotplot(tcell, var_names=['FOXP3'], groupby='region', layer = 'counts' )
sc.pl.dotplot(tcell, var_names=['CCL5', 'ZAP70'], groupby='region', layer = 'counts', dot_max=0.65, vmin=1.5, vmax=2.4)
# sc.pl.dotplot(tcell, var_names=['CCL5'], groupby='region', layer = 'counts', figsize=[2, 2], dot_max=0.60)
# sc.pl.dotplot(tcell, var_names=['ZAP70'], groupby='region', layer = 'counts')

# %%
markers = ['PDCD1', 'CTLA4', 'LAG3', 'GZMA', 'GNLY', 'GZMB', 'IFNG', 'CCL5', 'ZAP70']
# markers = [ 'CCL5', 'ZAP70']
for idx in tcell.obs.Index.unique():
    tcell_sub = tcell[tcell.obs['Index'] == idx, :]
    sc.pl.dotplot(tcell_sub, var_names=markers, groupby='region', layer = 'counts', dot_max=0.65)

# %%
tcell.obsm['X_umap'] = tcell.obsm['X_umap_nbr_0'].copy()
sc.pl.umap(tcell, color=['leiden', 'Tissue_1'], ncols=1, wspace=0.3)
sc.pl.umap(tcell, color=['PDCD1', 'CTLA4', 'LAG3'], ncols=1, wspace=0.3)

# %%
sc.pl.dotplot(tcell, var_names=['PDCD1', 'CTLA4', 'LAG3', 'GZMA', 'GNLY', 'GZMB', 'IFNG'], groupby='region', )

# %%
# 4 mins for 0.5 M cells, 8 mins for 0.8 M cells in A100-GPU
manage_gpu(2)
rsc.get.anndata_to_GPU(tcell)
rsc.tl.embedding_density(tcell, basis='umap', groupby='leiden')
rsc.get.anndata_to_CPU(tcell)

# %%
sc.pl.embedding_density(tcell, basis='umap', groupby='leiden', bg_dotsize=1, fg_dotsize=3)

# %%
adata.var.iloc[:, :2].to_csv('./results/genes.csv')
