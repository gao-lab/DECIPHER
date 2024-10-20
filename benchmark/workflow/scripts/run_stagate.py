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
import STAGATE_pyG
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

# Spatial graph
STAGATE_pyG.Cal_Spatial_Net(adata, k_cutoff=20, model='KNN')

# Preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)

# Train model
adata = STAGATE_pyG.train_STAGATE(adata)

run_time = str(time.time() - start)
print('Runtime: ' + run_time)

# %%
emb = adata.obsm['STAGATE']
np.save(nbr_emb_file, emb)
np.save(center_emb_file, emb)

# %%
time_dic = {"run_time": run_time}
with open(time_file, "w") as f:
    yaml.dump(time_dic, f)
