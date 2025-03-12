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

# %%
import time
import yaml

import numpy as np
import scniche as sn
import scanpy as sc
import torch
import warnings
warnings.filterwarnings('ignore')

print("Last run with scNiche version:", sn.__version__)

# set seed
sn.pp.set_seed()

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

cutoff = 20
lr = 0.01
epochs = 100
batch_num = 100

# prepare
adata = sn.pp.cal_spatial_neighbors(adata=adata, celltype_key='cell_type', mode='KNN', k_cutoff=cutoff, verbose=False)
adata = sn.pp.cal_spatial_exp(adata=adata, mode='KNN', k_cutoff=cutoff, is_pca=True, n_comps=50, verbose=False)
adata = sn.pp.prepare_data_batch(adata=adata, verbose=False, batch_num=batch_num)
                
# training
if torch.cuda.is_available():
    model = sn.tr.Runner_batch(adata=adata, device='cuda:0', verbose=False)
else:
    model = sn.tr.Runner_batch(adata=adata, device='cpu', verbose=False)
adata = model.fit(lr=lr, epochs=epochs)
                
emb = adata.obsm['X_scniche']

run_time = str(time.time() - start)
print('Runtime: ' + run_time)

# %%
np.save(nbr_emb_file, emb)
np.save(center_emb_file, emb)

# %%
time_dic = {"run_time": run_time}
with open(time_file, "w") as f:
    yaml.dump(time_dic, f)

