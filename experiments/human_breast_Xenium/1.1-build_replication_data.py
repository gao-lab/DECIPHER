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
# # technology replication in Xenium and Visium
#
#

# %%
import scanpy as sc
import pandas as pd
from decipher.utils import gex_embedding, IMMUNE_MARKER, CANCER_MARKER, OTHER_MARKER, scanpy_viz


# %% [markdown]
# ## Xenium rep2

# %%
adata = sc.read_10x_h5('./data/rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_cell_feature_matrix.h5')
df = pd.read_csv('./data/rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_cells.csv.gz')

# %%
df.set_index(adata.obs_names, inplace=True)
adata.obs = df.copy()
adata.layers["counts"] = adata.X.copy()
adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()

# %%
adata = gex_embedding(adata)
adata.X = adata.layers["counts"].copy()
adata.write_h5ad("data/rep2/adata.h5ad")

# %%
del adata
adata1 = sc.read_h5ad('./data/adata.h5ad')
adata2 = sc.read_h5ad('./data/rep2/adata.h5ad')
adata = adata1.concatenate(adata2)

# %%
adata = gex_embedding(adata)

# %%
XENIUM_BREAST_MARKER = [
    'BANK1',
    'KRT15', # myoepithelial
    'ACAT2', # myoepithelial
    'POSTN', # stromal
    'LRRC15', # stromalz
    'FABP4', # adipocytes
    'CEACAM6', # DCIS
    'FASN', # invasive
    'VWF', # endothelial,
    'CD69', # Mast cells
    'FGL2', # macrophages
    'CCPG1',
    ]

markers = IMMUNE_MARKER + CANCER_MARKER + OTHER_MARKER + XENIUM_BREAST_MARKER + ['total_counts']
markers = list(set(markers))
markers = [m for m in markers if m in adata.var_names]
len(markers)
sc.pl.dotplot(adata, markers, groupby='leiden')


# %%
import rapids_singlecell as rsc

rsc.get.anndata_to_GPU(adata)
rsc.tl.rank_genes_groups_logreg(adata, groupby="leiden", use_raw=False)
rsc.get.anndata_to_CPU(adata)
sc.pl.rank_genes_groups(adata, n_genes=20)

# %%
sc.pl.umap(adata, color=['cell_type', 'batch', 'leiden'], vmax=6, ncols=1, legend_loc='on data', legend_fontsize=5)


# %%
adata.obs['cell_type'] = 'unknown'
adata.obs.loc[adata.obs['leiden'].isin(['13']), 'cell_type'] = 'Bcell'
adata.obs.loc[adata.obs['leiden'].isin(['4']), 'cell_type'] = 'Plasma'
adata.obs.loc[adata.obs['leiden'].isin(['0']), 'cell_type'] = 'DC'
adata.obs.loc[adata.obs['leiden'].isin(['19']), 'cell_type'] = 'Breast glandular'

adata.obs.loc[adata.obs['leiden'].isin([ '7', '12', '14',  '15','28',]), 'cell_type'] = 'Stromal'
adata.obs.loc[adata.obs['leiden'].isin(['5', '10']), 'cell_type'] = 'Tcell'
adata.obs.loc[adata.obs['leiden'].isin(['25']), 'cell_type'] = 'Endothelial'
adata.obs.loc[adata.obs['leiden'].isin(['16']), 'cell_type'] = 'Myoepithelial'
adata.obs.loc[adata.obs['leiden'].isin(['1', '2', '3', '6', '8', '9', '21','22', '23' ,'24', '27']), 'cell_type'] = 'Invasive'
adata.obs.loc[adata.obs['leiden'].isin(['11', '17', '18' ,'26']), 'cell_type'] = 'Macrophage'
adata.obs.loc[adata.obs['leiden'].isin(['20']), 'cell_type'] = 'IRF7+ Dendritic' 


# %%
adata.write_h5ad('./data/merged_adata.h5ad')

# %% [markdown]
# ## Visium (depricated)

# %%
visium = sc.read_10x_h5('./data/rep_visium/CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5')
visium.shape

# %%
spatial = pd.read_csv('./data/rep_visium/spatial/tissue_positions.csv')
spatial.shape
visium.obs = spatial.copy()
visium.obsm['spatial'] = visium.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].copy().to_numpy()
visium.write_h5ad('data/rep_visium/adata.h5ad')
visium

# %%
sc.pl.spatial(visium, color=['CD3D'], spot_size=100, vmax=5)

# %%
