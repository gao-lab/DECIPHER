# ---
# jupyter:
#   jupytext:
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
# # Plot UMAP
#
# Plot UMAPs of STADIA, Graph

import numpy as np
import pandas as pd

# %%
import scanpy as sc

from decipher.utils import scanpy_viz

# %%
adata = sc.read_h5ad("../../input/mimic/merged.h5ad")

# %% [markdown]
# ## STADIA

# %%
stadia_emb = pd.read_csv("./results/mimic_STADIA_emb.csv")
cell_id = stadia_emb.columns
stadia_emb = stadia_emb.T
print(stadia_emb.shape)

# %%
adata.obs.index = adata.obs_names.str.rsplit("-", n=1).str[0]
adata_sub = adata[adata.obs_names.isin(cell_id)].copy()
print(adata_sub.shape)

# %%
# reorder the emb by the adata_sub.obs.index
stadia_emb_reorder = stadia_emb.reindex(adata_sub.obs.index)
adata_sub.obsm["X_stadia"] = stadia_emb_reorder.values

# %%
scanpy_viz(adata_sub, ["stadia"], leiden=False, rapids=True)

# %%
dataset = "10x_mimic"
method = "stadia"

adata_sub.obsm["X_umap"] = adata_sub.obsm[f"X_umap_{method}"].copy()
sc.pl.umap(adata_sub, color="cell_type", title=method, save=f"_{dataset}_{method}_omics.pdf")
sc.pl.umap(adata_sub, color="region", title=method, save=f"_{dataset}_{method}_spatial.pdf")
sc.pl.umap(adata_sub, color="batch", title=method, save=f"_{dataset}_{method}_batch.pdf")

# %% [markdown]
# ## GraphST

# %%
graphst_emb = np.load("../../results/human_pbmc_10x_mimic/seed:0/graphst/center_emb.npy")

# %%
adata.obsm["X_graphst"] = graphst_emb

# %%
scanpy_viz(adata, ["graphst"], leiden=False, rapids=True)

# %%
dataset = "10x_mimic"
for method in ["graphst"]:
    adata.obsm["X_umap"] = adata.obsm[f"X_umap_{method}"].copy()
    sc.pl.umap(adata, color="cell_type", title=method, save=f"_{dataset}_{method}_omics.pdf")
    sc.pl.umap(adata, color="region", title=method, save=f"_{dataset}_{method}_spatial.pdf")
    sc.pl.umap(adata, color="batch", title=method, save=f"_{dataset}_{method}_batch.pdf")
