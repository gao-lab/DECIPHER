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
# ## T cell analysis in pan-cancer dataset
#
# > Run with hpc A100 80G
#
# > We use rapids single cell as backend

# %%
from pathlib import Path
# import python reload 
from importlib import reload

import scanpy as sc
import rapids_singlecell as rsc
import numpy as np
import spider
from spider.utils import scanpy_viz, clip_umap, manage_gpu, gex_embedding
from spider.plot import split_umap

# reload scanpy_viz function
reload(spider.utils)

# %%
adata_path = "./data/pancancer_filter_anno.h5ad"
adata = sc.read_h5ad(adata_path)
