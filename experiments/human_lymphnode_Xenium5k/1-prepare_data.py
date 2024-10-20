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
# # Xenium 5K data

# %%
import scanpy as sc
import pandas as pd
import rapids_singlecell as rsc
from decipher.utils import gex_embedding, clip_umap, IMMUNE_MARKER, CANCER_MARKER, OTHER_MARKER

# %%
adata = sc.read_10x_h5('./data/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs/cell_feature_matrix.h5')
adata.layers['counts'] = adata.X.copy()

# %%
meta = pd.read_parquet('./data/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs/cells.parquet' )
adata.obs = meta
adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values

# %%
adata.obs['total_counts'] = adata.layers['counts'].sum(axis=1)

# %%
adata = gex_embedding(adata, hvg_only=False, resolution = 0.5)

# %%
sc.set_figure_params(dpi=300)
sc.pl.spatial(adata, color='leiden', spot_size=5)
sc.set_figure_params(dpi=100)

# %%
adata.obsm['X_umap_raw'] = adata.obsm['X_umap'].copy()
adata.obsm['X_umap'] = clip_umap(adata.obsm['X_umap'])

sc.pl.umap(adata, color='leiden', legend_loc='on data')
sc.pl.umap(adata, color=['CD3E' ,'CD79A' ,'CD14'], ncols=2)

# %%
sc.pl.umap(adata, color=['CCL14' ,'TSPAN7', 'PLVAP'], ncols=2)

# %%
rsc.get.anndata_to_GPU(adata)
rsc.tl.rank_genes_groups_logreg(adata, groupby="leiden", use_raw=False)
rsc.get.anndata_to_CPU(adata)
sc.pl.rank_genes_groups(adata, n_genes=20)

# %%
LYMPH = ['PDCD1', 'ACTA2', 'RGS5', 'CD14', 'CD79A', 'TCF7', 'TPSAB1', 'MMP9','C1QC',
         'CR2', 'CXCL13', 'FDCSP', 'MKI67' ,'IRF7', 'BCL11A', 'GMDS', 'LPP', '"BCL6', 'S100A3']

markers = IMMUNE_MARKER + CANCER_MARKER + OTHER_MARKER  + LYMPH
markers = list(set(markers))
markers = [m for m in markers if m in adata.var_names]
len(markers)
sc.pl.dotplot(adata, markers, groupby='leiden')


# %%
sc.pl.dotplot(adata, ['TOX', 'TOX2' ,'CD200', 'BCL6'], groupby='leiden')


# %%

adata.obs['cell_type'] = 'unknown'
adata.obs.loc[adata.obs['leiden'].isin(['0', '1' ,'2' ,'3', '6']), 'cell_type'] = 'Tcell'
adata.obs.loc[adata.obs['leiden'].isin(['10']), 'cell_type'] = 'Plasma'
adata.obs.loc[adata.obs['leiden'].isin(['4']), 'cell_type'] = 'Macrophage'
adata.obs.loc[adata.obs['leiden'].isin(['7']), 'cell_type'] = 'pDC'
adata.obs.loc[adata.obs['leiden'].isin(['9']), 'cell_type'] = 'DC'
adata.obs.loc[adata.obs['leiden'].isin(['5', '8']), 'cell_type'] = 'Bcell'
# adata.obs.loc[adata.obs['leiden'].isin(['8']), 'cell_type'] = 'Bcell_prolifing'
adata.obs.loc[adata.obs['leiden'].isin(['12']), 'cell_type'] = 'VSMCs'
adata.obs.loc[adata.obs['leiden'].isin(['11']), 'cell_type'] = 'low_quality'
adata.obs.loc[adata.obs['leiden'].isin(['14']), 'cell_type'] = 'Endothelial'
adata.obs.loc[adata.obs['leiden'].isin(['13']), 'cell_type'] = 'fDC'
adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

# %%
adata.obsm['X_umap'] = clip_umap(adata.obsm['X_umap'])
sc.pl.umap(adata, color='cell_type')

# %%
sc.pl.dotplot(adata,
              ['CD79A',
               'ITGAX',
               'PLVAP',
               'MMP9',
               'MZB1',
               'CD3E',
               'NOTCH3',
               'CR2',
               'IRF7'
               ],
            groupby='cell_type')


# %%
sc.pl.violin(adata, ['total_counts'], groupby='cell_type', rotation=90, ylabel='Total counts')

# %%
adata.write_h5ad('data/lymph_node.h5ad')

# %%
adata = sc.read_h5ad('./data/lymph_node.h5ad')
