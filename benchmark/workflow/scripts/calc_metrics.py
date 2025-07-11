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
import os
from pathlib import Path

import yaml
import numpy as np
import scanpy as sc
import scib
import numpy as np
from sklearn.neighbors import KNeighborsTransformer
from addict import Dict

from decipher.utils import scanpy_viz, estimate_spot_size

import ray

sc._settings.settings.n_jobs = 1
# %% tags=["parameters"]
# parameter tag
center_emb = "./center_emb.npy"
nbr_emb = "./nbr_emb.npy"
metrics_file = "./metrics.yaml"
time = "./time.yaml"

# %%
time = yaml.load(open(time), Loader=yaml.FullLoader)
time = float(time["run_time"])

# %% [markdown]
# ## Build anndata object

# %%
if 'search_decipher' not in str(Path(center_emb).absolute()):
    h5ad_file = Path(center_emb).absolute().parent.parent / "data.h5ad"
else:
    h5ad_file = Path(center_emb).absolute().parent.parent.parent / "data.h5ad"
adata = sc.read_h5ad(h5ad_file)

# %%
if 'decipher' not in str(Path(center_emb).absolute()):
    adata.obsm["X_center"] = np.load(center_emb)
    CLUSTER_KEYS = ["center"]
else:
    adata.obsm["X_center"] = np.load(center_emb)
    adata.obsm["X_nbr"] = np.load(nbr_emb)
    CLUSTER_KEYS = ["center", "nbr"]
adata

# %%
if "batch" in adata.obs.columns:
    print("The 'batch' column is found in the adata.obs.")
    BATCH_FALG = True
else:
    BATCH_FALG = False

# %% [markdown]
# ## Clustering

# ## Calculate metrics
#
# we calculate the following metrics for each embedding:
# - ARI
# - NMI
# - ASW
#
# For batch removal, we calculate the following metrics:
# - Batch ASE
# - Graph connectivity
# - iLSI
# - kBET

# %%

def calc_metric_once(adata):
    metric_dict = Dict()
    transformer = KNeighborsTransformer(n_neighbors=15, metric = 'euclidean', n_jobs=8)
    nbr_slot = "X_nbr" if len(CLUSTER_KEYS) > 1 else "X_center"

    metric_dict.gex_asw = scib.me.silhouette(adata, label_key="cell_type", embed="X_center")
    metric_dict.nbr_asw = scib.me.silhouette(adata, label_key="region", embed=nbr_slot)

    # batch metric
    if "batch" in adata.obs.columns:
        print("The 'batch' column is found in the adata.obs.")
        metric_dict.gex_batch_asw = scib.me.silhouette_batch(adata, batch_key="batch", label_key="cell_type", embed="X_center")
        metric_dict.gex_batch_ilisi = scib.me.ilisi_graph(adata, batch_key="batch", type_="embed", use_rep="X_center", n_cores=8)
        sc.pp.neighbors(adata, use_rep="X_center", transformer=transformer)
        metric_dict.gex_batch_gc = scib.me.graph_connectivity(adata, label_key="cell_type")
        metric_dict.gex_batch_kbet = scib.me.kBET(
            adata, batch_key="batch", label_key="cell_type", type_="embed", embed="X_center"
            )
        # nbr part

        metric_dict.nbr_batch_asw = scib.me.silhouette_batch(adata, batch_key="batch", label_key="region", embed=nbr_slot)
        metric_dict.nbr_batch_ilisi = scib.me.ilisi_graph(adata, batch_key="batch", type_="embed", use_rep=nbr_slot, n_cores=8)
        sc.pp.neighbors(adata, use_rep=nbr_slot, transformer=transformer)
        metric_dict.nbr_batch_gc= scib.me.graph_connectivity(adata, label_key="region")
        metric_dict.nbr_batch_kbet = scib.me.kBET(
            adata, batch_key="batch", label_key="region", type_="embed", embed=nbr_slot
            )

    else:
        print('Do not have batch info')
        metric_dict.gex_batch_asw = 0.0
        metric_dict.gex_batch_ilisi = 0.0
        metric_dict.gex_batch_gc = 0.0
        metric_dict.gex_batch_kbet = 0.0
        metric_dict.nbr_batch_asw = 0.0
        metric_dict.nbr_batch_ilisi = 0.0
        metric_dict.nbr_batch_gc = 0.0
        metric_dict.nbr_batch_kbet = 0.0
    return metric_dict.to_dict()


def calc_metric(adata, resolution:float):
    adata = scanpy_viz(adata, resolution=resolution, keys=CLUSTER_KEYS)
    print('Finish scanpy_viz')
    metric_dict = Dict()
    # omics metric
    metric_dict.gex_ari = scib.me.ari(adata, 'leiden_center', 'cell_type')
    metric_dict.gex_nmi = scib.me.nmi(adata, 'leiden_center', 'cell_type')

    # spatial metric
    if len(CLUSTER_KEYS) > 1:
        metric_dict.nbr_ari = scib.me.ari(adata, 'leiden_nbr', 'region')
        metric_dict.nbr_nmi = scib.me.nmi(adata, 'leiden_nbr', 'region')
    else:
        metric_dict.nbr_ari = scib.me.ari(adata, 'leiden_center', 'region')
        metric_dict.nbr_nmi = scib.me.nmi(adata, 'leiden_center', 'region')

    # batch_metric_dict = calc_batch_metric(adata)
    # metric_dict.update(batch_metric_dict)
    # convert all values to float
    return metric_dict.to_dict()

@ray.remote(num_cpus=5)
def remote_calc_metric(adata, resolution:float):
    return calc_metric(adata, resolution)

# %%
# metric without parallel
metric_dict_once = calc_metric_once(adata)

# %%
# for resolution in np.linspace(0.1, 1, 10):
#     adata = scanpy_viz(adata, resolution=0.3)
ray.init()
metrics_dicts = ray.get([remote_calc_metric.remote(adata, resolution) for resolution in np.linspace(0.1, 1, 10)])
result_dict = {}
for resolution, metrics_dict in zip(np.linspace(0.1, 1, 10), metrics_dicts):
    metrics_dict.update(metric_dict_once)
    for k, v in metrics_dict.items():
        metrics_dict[k] = float(v)
    result_dict[f'resolution:{resolution:.1f}'] = metrics_dict
ray.shutdown()

# %% [markdown]
# ## Save results

# %%
with open(metrics_file, "w") as f:
    yaml.dump(result_dict, f)


# %% [markdown]
# ## Visualization

# %%
# vars_plot = ["leiden_center", "cell_type", "region"]
# if BATCH_FALG:
#     vars_plot.append("batch")

# if len(CLUSTER_KEYS) > 1:
#     vars_plot.append("leiden_nbr")

# adata.obsm["X_umap"] = adata.obsm["X_umap_center"]
# sc.pl.umap(adata, color = vars_plot, ncols=1)

# if len(CLUSTER_KEYS) > 1:
#     adata.obsm["X_umap"] = adata.obsm["X_umap_nbr"]
#     sc.pl.umap(adata, color = vars_plot, ncols=1)

# spot_size = estimate_spot_size(adata.obsm["spatial"]) * 0.5
# sc.pl.spatial(adata, color = vars_plot, spot_size=spot_size)