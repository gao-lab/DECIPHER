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
import torch
from GraphST import GraphST
import numpy as np

# %% tags=["parameters"]
# parameter tag
input_file = ""
center_emb_file = ""
nbr_emb_file = ""
time_file = ""

# %%
adata = sc.read_h5ad(input_file)

# %%
start = time.time()

device = 'cuda' if torch.cuda.is_available() else 'cpu'
model = GraphST.GraphST(adata, device=device)

# train model
adata = model.train()

run_time = str(time.time() - start)
print('Runtime: ' + run_time)

# %%
emb = adata.obsm['emb']
np.save(nbr_emb_file, emb)
np.save(center_emb_file, emb)

# %%
time_dic = {"run_time": run_time}
with open(time_file, "w") as f:
    yaml.dump(time_dic, f)


