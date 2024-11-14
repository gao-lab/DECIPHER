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
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import rapids_singlecell as rsc

from decipher.utils import gex_embedding, IMMUNE_MARKER, CANCER_MARKER, OTHER_MARKER, scanpy_viz

# %% [markdown]
# ## Build data

# %%
adata = sc.read_10x_h5('./data/Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5')
df = pd.read_csv('./data/Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz')

# %%
df.set_index(adata.obs_names, inplace=True)
adata.obs = df.copy()
adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()


# %%
adata.obs

# %%
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)


# %%
cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")

# %%
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=5)

# %%
adata.layers["counts"] = adata.X.copy()

# %%
adata = gex_embedding(adata)

# %%
# sc.pp.normalize_total(adata, inplace=True)
# sc.pp.log1p(adata)
# sc.pp.pca(adata)
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)
# sc.tl.leiden(adata)

# %%
adata.X = adata.layers["counts"].copy()
adata.write_h5ad("data/adata.h5ad")

# %% [markdown]
# ## Cell type anno

# %%
adata = sc.read_h5ad("data/adata.h5ad")

# %%
sc.pl.umap(adata, color=["leiden"], legend_loc="on data")

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
markers

# %%
sc.pl.dotplot(adata, markers, groupby='leiden')

# %%
sc.pl.umap(adata, color=markers, vmax=6, ncols=6)

# %%
# from https://www.nature.com/articles/s41467-023-43458-x/figures/3
# paper_marker = [
#     'GZMB', # IRF7 + Dendritic
#     'SDC4', # Stromal
#     'LRRC15', # Stromal
#     'MMP1', # Myoepithelial
#     'EDN1', # Myoepithelial
#     'MYBPC1', # Myoepithelial
# ]
paper_tumor_marker = [
    'CCPG1',
    'DAPK3',
    'CTTN',
    'SQLE',
]
sc.pl.umap(adata, color='DNAAF1', vmax=3, ncols=2)

# %%
rsc.get.anndata_to_GPU(adata)
rsc.tl.rank_genes_groups_logreg(adata, groupby="leiden", use_raw=False)
rsc.get.anndata_to_CPU(adata)


# %%
sc.pl.rank_genes_groups(adata, n_genes=20)

# %%
sc.pl.spatial(adata, color='leiden', spot_size=10)

# %%
adata.obs['cell_type'] = 'unknown'
adata.obs.loc[adata.obs['leiden'].isin(['23']), 'cell_type'] = 'Bcell'
adata.obs.loc[adata.obs['leiden'].isin(['4']), 'cell_type'] = 'Plasma'

adata.obs.loc[adata.obs['leiden'].isin(['5', '7', '14', '21', '22']), 'cell_type'] = 'Stromal'
adata.obs.loc[adata.obs['leiden'].isin(['15', '18']), 'cell_type'] = 'Tcell'
adata.obs.loc[adata.obs['leiden'].isin(['19']), 'cell_type'] = 'Endothelial'
adata.obs.loc[adata.obs['leiden'].isin(['20']), 'cell_type'] = 'Myoepithelial'
adata.obs.loc[adata.obs['leiden'].isin(['0', '1', '2', '3', '6', '8', '9', '10','12', '16', '24']), 'cell_type'] = 'Invasive'
adata.obs.loc[adata.obs['leiden'].isin(['17']), 'cell_type'] = 'Macrophage'
adata.obs.loc[adata.obs['leiden'].isin(['11']), 'cell_type'] = 'IRF7+ Dendritic' 
adata.obs.loc[adata.obs['leiden'].isin(['13']), 'cell_type'] = 'DNAAF1+ cell' 

# %%
sc.pl.umap(adata, color='cell_type')

# %%
adata.X = adata.layers["counts"].copy()
adata.write_h5ad("data/adata.h5ad")

# %%
adata = sc.read_h5ad("data/adata.h5ad")
sc.pl.spatial(adata, color='cell_type', spot_size=10)
