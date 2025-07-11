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
from decipher import DECIPHER, CFG

# %% tags=["parameters"]
# parameter tag
input_file = ""
work_dir = ""
out_file = ""
time_file = ""

# %%
CFG.omics.pretrain.force = True

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
model = DECIPHER(work_dir=work_dir, user_cfg=CFG)
model.register_data(adata, split_by=split_by)
model.fit_omics()
run_time = str(time.time() - start)
print("Runtime: " + run_time)

# %%
time_dic = {"run_time": run_time}
with open(time_file, "w") as f:
    yaml.dump(time_dic, f)
