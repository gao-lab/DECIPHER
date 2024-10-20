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
import time

import yaml
import scanpy as sc
from spider import Spider, CFG

# %% tags=["parameters"]
# parameter tag
input_file = ""
work_dir = ""
out_file = ""
time_file = ""

# %%
CFG.omics.pretrain.force = True
if '10x' in input_file:
    CFG.omics.model.augment.dropout_gex = 0.6
    CFG.omics.spatial_graph.k = 25
    CFG.omics.model.emb_dim = 32
    CFG.omics.model.prj_dims[-1] = 32
    CFG.omics.model.epochs = 8
    CFG.omics.model.dropout = 0.0
    CFG.omics.model.transformer_layers = 3
    CFG.omics.model.temperature_center = 0.1
    CFG.omics.model.temperature_nbr = 0.01
    CFG.omics.loader.batch_size = 128
if 'xenium' in input_file:
    CFG.omics.model.augment.dropout_gex = 0.3
    CFG.omics.spatial_graph.k = 25
    CFG.omics.model.emb_dim = 32
    CFG.omics.model.prj_dims[-1] = 32
    CFG.omics.model.epochs = 6
    CFG.omics.model.dropout = 0.2
    CFG.omics.model.transformer_layers = 2
    CFG.omics.model.temperature_center = 0.05
    CFG.omics.model.temperature_nbr = 0.07
    CFG.omics.loader.batch_size = 64
if 'merfish' in input_file:
    CFG.omics.model.augment.dropout_gex = 0.3
    CFG.omics.spatial_graph.k = 40
    CFG.omics.model.emb_dim = 32
    CFG.omics.model.prj_dims[-1] = 32
    CFG.omics.model.epochs = 2
    CFG.omics.model.dropout = 0.0
    CFG.omics.model.transformer_layers = 1
    CFG.omics.model.temperature_center = 0.05
    CFG.omics.model.temperature_nbr = 0.1
    CFG.omics.loader.batch_size = 256

# %%
adata = sc.read_h5ad(input_file)
if 'batch' in adata.obs.columns:
    split_by = "batch"
elif 'graph_batch' in adata.obs.columns:
    CFG.omics.ignore_batch = True
    split_by = "graph_batch"
else:
    CFG.omics.ignore_batch = True
    split_by = None

# %%
start = time.time()
model = Spider(work_dir=work_dir, user_cfg=CFG)
model.register_data(adata, split_by=split_by)
model.fit_omics()
run_time = str(time.time() - start)
print("Runtime: " + run_time)

# %%
time_dic = {"run_time": run_time}
with open(time_file, "w") as f:
    yaml.dump(time_dic, f)
