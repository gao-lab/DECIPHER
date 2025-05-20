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
# # None spatial methods on pan-cancer dataset
#
# We compare with following methods:
#  - Scanpy
#  - Harmony
#  - scVI

# %%
import time

import scanpy as sc
import numpy as np

from spider.utils import scanpy_viz

# %%
adata2 = sc.read_h5ad("./results/adata_analysis.h5ad")
adata2.obs.head()

# %%
adata = sc.read_h5ad("./data/pancancer_filter_anno.h5ad")

# %%
adata.obs['region'] = adata2.obs['region'].copy()

# %%
adata.layers['counts'] = adata.layers['counts'].astype(np.int32)

# %%
adata.X = adata.layers['counts'].copy()

# %% [markdown]
# ## Run scVI

# %%
import scvi

# %%
start = time.time()
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["dataset"],
)
model = scvi.model.SCVI(adata)
model.train(max_epochs = 100, batch_size = 512)

run_time = str(time.time() - start)
print('Runtime: ' + run_time)


# %%
emb = model.get_latent_representation()
np.save('./results/scvi_emb.npy', emb)

# %%
emb = np.load('./results/scvi_emb.npy')

# %%
adata.obsm['X_scVI'] = emb
adata = scanpy_viz(adata, keys=['scVI'], approx = True)

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_scVI'].copy()
sc.pl.umap(adata, color=['cell_type'], ncols=1, title=['Cell type'])

# %%
sc.pl.umap(adata, color=['Tissue_1'], ncols=1, title=['Tissue'])

# %%
sc.pl.umap(adata, color=['region'], ncols=1, title=['Region'])

# %%
adata.obs['Tissue_1'].value_counts()

# %% [markdown]
# ## Run harmony

# %%
from importlib import reload
import spider
reload(spider.utils)

# %%
from spider.utils import gex_embedding

# %%
# del adata.layers

# %%
harm_emb = gex_embedding(adata, method ='harmony', batch_key='dataset', emb_only=True, memory_strategy='large')

# %%
np.save('./results/harmony_emb.npy', harm_emb)

# %%
adata.obsm['X_harmony'] = harm_emb
adata = scanpy_viz(adata, keys=['harmony'], approx = True)

# %%
adata.obsm['X_umap'] = adata.obsm['X_umap_harmony'].copy()
sc.pl.umap(adata, color=['cell_type'], ncols=1, title=['Cell type'])

# %%
sc.pl.umap(adata, color=['Tissue_1'], ncols=1, title=['Tissue'])

# %%
sc.pl.umap(adata, color=['region'], ncols=1, title=['Region'])

# %% [markdown]
# ## Run Scanpy
#
# We have run scanpy in previous notebook

# %%
sc.pl.umap(adata, color=['region'], ncols=1, title=['Region'])

# %%
sc.pl.umap(adata, color=['cell_type'], ncols=1, title=['Cell type'])

# %%
sc.pl.umap(adata, color=['Tissue_1'], ncols=1, title=['Tissue'])
