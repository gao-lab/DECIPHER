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
import torch
import scanpy as sc
import numpy as np
from torch_geometric.nn import SimpleConv
from decipher.graphic.build import build_graph
from sklearn.preprocessing import LabelEncoder

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

# preocessing
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",
)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata)

# build spatial edges
if 'batch' in adata.obs.columns:
    batch = LabelEncoder().fit_transform(adata.obs['batch'])
elif 'graph_batch' in adata.obs.columns:
    batch = LabelEncoder().fit_transform(adata.obs['graph_batch'])
else:
    batch = None
edge_index = build_graph(
    adata.obsm['spatial'],
    batch=batch,
    mode='knn',
    k = 15,
)

# run
x = torch.Tensor(adata.obsm['X_pca'])
gcn = SimpleConv(aggr="mean", combine_root='cat')
emb = gcn(x, edge_index)

run_time = str(time.time() - start)
print('Runtime: ' + run_time)

# %%
np.save(nbr_emb_file, emb)
np.save(center_emb_file, emb)

# %%
time_dic = {"run_time": run_time}
with open(time_file, "w") as f:
    yaml.dump(time_dic, f)