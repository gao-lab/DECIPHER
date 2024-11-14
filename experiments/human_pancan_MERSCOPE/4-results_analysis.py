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
# # Analysis pan-cancer MERSCOPE dataset results
# > Run with hpc A100 80G
# > We use rapids single cell as backend
#

# %%
from pathlib import Path
# import python reload 
from importlib import reload

import scanpy as sc
import rapids_singlecell as rsc
import numpy as np
import decipher
from decipher.utils import scanpy_viz, clip_umap, manage_gpu, gex_embedding
from decipher.plot import split_umap

# reload scanpy_viz function
reload(decipher.utils)

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

# %% [markdown]
# ## Subtype analysis

# %% [markdown]
# ### T cell analysis

# %%
tcell = adata_plot[adata_plot.obs['cell_type'] == 'T cell', :].copy()

# %%
tcell.X = tcell.layers['counts'].copy()
tcell = gex_embedding(tcell, method='harmony', batch_key='Index', resolution=0.5, rapids_after_scale=True)
tcell.obs['leiden'].value_counts()

# %%
tcell_marker = ['CD4', 'CD8A', 'FOXP3', 'PDCD1', 'CCR7', 'GZMK']
sc.pl.dotplot(tcell, var_names=tcell_marker, groupby='leiden')

# %%
sc.pl.umap(tcell, color=['leiden', 'PDCD1'], ncols=4, wspace=0.3)

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
adata = sc.read_h5ad('./results/adata_analysis.h5ad')
adata.obs

# %%
adata.obs['cell_type'].value_counts()
