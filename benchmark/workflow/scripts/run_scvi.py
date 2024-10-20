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
import scvi
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

# preocessing
adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
)

# set data
if 'batch' in adata.obs.columns:
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        categorical_covariate_keys=["batch"],
    )
else:
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
    )

# init model
model = scvi.model.SCVI(adata)
model

# train model
model.train()

run_time = str(time.time() - start)
print('Runtime: ' + run_time)

# %%
emb = model.get_latent_representation()
np.save(nbr_emb_file, emb)
np.save(center_emb_file, emb)

# %%
time_dic = {"run_time": run_time}
with open(time_file, "w") as f:
    yaml.dump(time_dic, f)
