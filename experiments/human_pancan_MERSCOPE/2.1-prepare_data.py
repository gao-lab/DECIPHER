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
from decipher import DECIPHER, CFG

# %%
CFG.omics.model.augment.dropout_gex = 0.4
CFG.omics.model.batch_size = 512
# disable validation plot
CFG.omics.model.plot = False
CFG.device_num = 4
CFG.omics.model.max_steps = 20_000
CFG.omics.spatial_graph.k = 15

# %%
model = DECIPHER(work_dir="./results/decipher_6_10", user_cfg=CFG)

# %%
adata_path = "./data/pancancer_filter_anno.h5ad"
adata = sc.read_h5ad(adata_path)

# %%
adata.obs['_batch'] = adata.obs['dataset'].astype(int)
model.register_data(adata, preprocess=False)
