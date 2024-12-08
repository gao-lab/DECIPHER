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
# # Compare with COMMOT

# %%
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np

# %%
adata = sc.read_h5ad('./data/lymph_node.h5ad')
adata.X = adata.layers["counts"].copy()
adata.var_names_make_unique()

# %%
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

# %%
adata.X = adata.X.astype(np.float32)

# %%
# random downsample to 10k cells
adata = adata[np.random.choice(adata.obs.index, 10000, replace=False), :]

# %%
df_ligrec=ct.pp.ligand_receptor_database(database='CellChat', species='human')

# %%
ct.tl.spatial_communication(adata,
    database_name='cellchat', df_ligrec=df_ligrec, dis_thr=200, heteromeric=True, pathway_sum=True)

# %%
adata

# %%
ct.tl.communication_direction(adata, database_name='cellchat', pathway_name='CXCL12-CXCR4', k=5)
ct.pl.plot_cell_communication(adata, database_name='cellchat', pathway_name='CXCL12-CXCR4', plot_method='grid', background_legend=True,
    scale=0.00003, ndsize=3, grid_density=0.4, summary='sender', clustering='leiden',
    normalize_v = True, normalize_v_quantile=0.995)

# %%
ct.tl.communication_direction(adata, database_name='cellchat', pathway_name='CXCL13-CXCR5', k=5)
ct.pl.plot_cell_communication(adata, database_name='cellchat', pathway_name='CXCL13-CXCR5', plot_method='grid', background_legend=True,
    scale=0.00003, ndsize=3, grid_density=0.4, summary='sender', clustering='leiden',
    normalize_v = True, normalize_v_quantile=0.995)

# %%
