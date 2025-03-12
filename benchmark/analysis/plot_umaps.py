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
# # UMAP plot of benchmark datasets

# %%
from pathlib import Path

import numpy as np
import ray
import scanpy as sc

from decipher.utils import scanpy_viz

# %% [markdown]
# ## Simulation dataset


# %%
def plot_umap(adata: sc.AnnData, dataset: str, method: str, force: bool = False):
    root_dir = Path("../results") / dataset / "seed:3" / method
    print(root_dir)
    if not root_dir.is_dir():
        return
    if "gnn" in method or "search" in method:
        return
    if not (root_dir / "nbr_emb.npy").exists():
        print(f"{method} does not have embeddings")
        return
    if (
        Path(f"./figures/umap_{dataset}_{method}_omics.pdf").exists()
        and Path(f"./figures/umap_{dataset}_{method}_spatial.pdf").exists()
        and not force
    ):
        print(f"{method} of dataset :{dataset} already plotted")
        return

    center_emb = np.load(root_dir / "center_emb.npy")
    adata.obsm[f"X_{method}_center_emb"] = center_emb
    scanpy_viz(adata, f"{method}_center_emb", leiden=False, gpu_id=0, rapids=True)
    adata.obsm["X_umap"] = adata.obsm[f"X_umap_{method}_center_emb"].copy()
    sc.pl.umap(adata, color="cell_type", title=method, save=f"_{dataset}_{method}_omics.pdf")
    if "mimic" in dataset:
        sc.pl.umap(adata, color="batch", title=method, save=f"_{dataset}_{method}_batch.pdf")

    if "decipher" in method:
        nbr_emb = np.load(root_dir / "nbr_emb.npy")
        adata.obsm[f"X_{method}_nbr_emb"] = nbr_emb
        scanpy_viz(adata, f"{method}_nbr_emb", leiden=False, gpu_id=0)
        adata.obsm["X_umap"] = adata.obsm[f"X_umap_{method}_nbr_emb"].copy()
    sc.pl.umap(adata, color="region", title=method, save=f"_{dataset}_{method}_spatial.pdf")


@ray.remote(num_cpus=4, num_gpus=0.5)
def ray_plot_umap(*args):
    plot_umap(*args)


# %%
tasks = []
ray.shutdown()
ray.init(ignore_reinit_error=True)
for dataset in Path("../results").glob("human*"):
    adata = sc.read_h5ad(dataset / "seed:0" / "data.h5ad")
    for method in (dataset / "seed:0").iterdir():
        tasks.append(ray_plot_umap.remote(adata, dataset.name, method.name, False))
ray.get(tasks)
ray.shutdown()


# %%
# for dataset in Path('../results').glob('human*'):
#     adata = sc.read_h5ad(dataset / 'seed:0' / 'data.h5ad')
#     for method in (dataset / 'seed:0' ).iterdir():
#         plot_umap(adata, dataset.name, method.name)
