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

# %% [markdown]
# # Compute the metrics based on the labels
#
# The STADIA and BASS only provide the labels

import numpy as np
import pandas as pd
import ray
import scanpy as sc
import scib

# %%
import yaml
from addict import Dict

from decipher.utils import scanpy_viz

# %%
adata = sc.read_h5ad("../../input/mimic/merged.h5ad")


# %% [markdown]
# ## BASS


# %%
def calc_metric(adata):
    metric_dict = Dict()
    # omics metric
    metric_dict.gex_ari = scib.me.ari(adata, "leiden_center", "cell_type")
    metric_dict.gex_nmi = scib.me.nmi(adata, "leiden_center", "cell_type")

    # spatial metric
    metric_dict.nbr_ari = scib.me.ari(adata, "leiden_nbr", "region")
    metric_dict.nbr_nmi = scib.me.nmi(adata, "leiden_nbr", "region")
    return metric_dict.to_dict()


# %%
df = pd.read_csv("./results/mimic_BASS_res.csv")
df.head(3)

# %%
adata.obs["leiden_nbr"] = df["zlabels"].values
adata.obs["leiden_center"] = df["clabels"].values


# %%
bass_metrics = calc_metric(adata)

# %% [markdown]
# ## STADIA

# %%
# from sklearn.neighbors import KNeighborsTransformer

CLUSTER_KEYS = ["center"]


def calc_metric_once(adata):
    metric_dict = Dict()
    # transformer = KNeighborsTransformer(n_neighbors=15, metric = 'euclidean', n_jobs=8)
    nbr_slot = "X_nbr" if len(CLUSTER_KEYS) > 1 else "X_center"

    metric_dict.gex_asw = scib.me.silhouette(adata, label_key="cell_type", embed="X_center")
    metric_dict.nbr_asw = scib.me.silhouette(adata, label_key="region", embed=nbr_slot)

    # batch metric
    if "batch" in adata.obs.columns:
        print("The 'batch' column is found in the adata.obs.")
        metric_dict.gex_batch_asw = scib.me.silhouette_batch(
            adata, batch_key="batch", label_key="cell_type", embed="X_center"
        )
        metric_dict.gex_batch_ilisi = scib.me.ilisi_graph(
            adata, batch_key="batch", type_="embed", use_rep="X_center", n_cores=8
        )
        sc.pp.neighbors(adata, use_rep="X_center")
        metric_dict.gex_batch_gc = scib.me.graph_connectivity(adata, label_key="cell_type")
        metric_dict.gex_batch_kbet = scib.me.kBET(
            adata, batch_key="batch", label_key="cell_type", type_="embed", embed="X_center"
        )

        # nbr part
        metric_dict.nbr_batch_asw = scib.me.silhouette_batch(
            adata, batch_key="batch", label_key="region", embed=nbr_slot
        )
        metric_dict.nbr_batch_ilisi = scib.me.ilisi_graph(
            adata, batch_key="batch", type_="embed", use_rep=nbr_slot, n_cores=8
        )
        if nbr_slot == "X_nbr":
            sc.pp.neighbors(adata, use_rep=nbr_slot)
        metric_dict.nbr_batch_gc = scib.me.graph_connectivity(adata, label_key="region")
        metric_dict.nbr_batch_kbet = scib.me.kBET(
            adata, batch_key="batch", label_key="region", type_="embed", embed=nbr_slot
        )

    else:
        print("Do not have batch info")
        metric_dict.gex_batch_asw = 0.0
        metric_dict.gex_batch_ilisi = 0.0
        metric_dict.gex_batch_gc = 0.0
        metric_dict.gex_batch_kbet = 0.0
        metric_dict.nbr_batch_asw = 0.0
        metric_dict.nbr_batch_ilisi = 0.0
        metric_dict.nbr_batch_gc = 0.0
        metric_dict.nbr_batch_kbet = 0.0
    return metric_dict.to_dict()


def calc_metric(adata, resolution: float):
    adata = scanpy_viz(adata, resolution=resolution, keys=CLUSTER_KEYS)
    print("Finish scanpy_viz")
    metric_dict = Dict()
    # omics metric
    metric_dict.gex_ari = scib.me.ari(adata, "leiden_center", "cell_type")
    metric_dict.gex_nmi = scib.me.nmi(adata, "leiden_center", "cell_type")

    # spatial metric
    if len(CLUSTER_KEYS) > 1:
        metric_dict.nbr_ari = scib.me.ari(adata, "leiden_nbr", "region")
        metric_dict.nbr_nmi = scib.me.nmi(adata, "leiden_nbr", "region")
    else:
        metric_dict.nbr_ari = scib.me.ari(adata, "leiden_center", "region")
        metric_dict.nbr_nmi = scib.me.nmi(adata, "leiden_center", "region")

    # batch_metric_dict = calc_batch_metric(adata)
    # metric_dict.update(batch_metric_dict)
    # convert all values to float
    return metric_dict.to_dict()


@ray.remote(num_cpus=5)
def remote_calc_metric(adata, resolution: float):
    return calc_metric(adata, resolution)


# %%
emb = pd.read_csv("./results/mimic_STADIA_emb.csv")

# %%
adata.obs.index = adata.obs_names.str.rsplit("-", n=1).str[0]

# %%
cell_id = emb.columns
adata_sub = adata[adata.obs_names.isin(cell_id)].copy()
adata_sub

# %%
emb = emb.T
emb.head()

# %%
# reorder the emb by the adata_sub.obs.index
emb_reorder = emb.reindex(adata_sub.obs.index)
emb_reorder.head()


# %%
adata_sub.obsm["X_center"] = emb_reorder.values

# %%
metric_dict_once = calc_metric_once(adata_sub)

# %%
metric_dict_once

# %%
ray.init()
metrics_dicts = ray.get(
    [remote_calc_metric.remote(adata_sub, resolution) for resolution in np.linspace(0.1, 1, 10)]
)
result_dict = {}
for resolution, metrics_dict in zip(np.linspace(0.1, 1, 10), metrics_dicts):
    metrics_dict.update(metric_dict_once)
    for k, v in metrics_dict.items():
        metrics_dict[k] = float(v)
    result_dict[f"resolution:{resolution:.1f}"] = metrics_dict
ray.shutdown()

# %%
metrics_file = "./results/mimic_STADIA_metrics.yaml"
with open(metrics_file, "w") as f:
    yaml.dump(result_dict, f)
