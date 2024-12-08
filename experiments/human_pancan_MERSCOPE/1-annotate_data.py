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

# %% [markdown]
# # Annotate the MERSCOPE pan-cancer data
#
# > WARNING: we use `rapids-singlecell` which is based on GPU but have less reproducibility than CPU-based scanpy.
#
# We found different tissue seems have batch effect while same tissue do not (it just technological replication)

# %%
import glob
from pathlib import Path

import cupy as cp
import joblib
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc
from anndata import AnnData
from cupyx.scipy import sparse
from rapids_singlecell.preprocessing._scale import (
    _check_gpu_X,
    _check_mask,
    _get_obs_rep,
    _set_obs_rep,
    view_to_actual,
)

from spider.utils import gex_embedding

# %% [markdown]
# ## Define some useful functions


# %%
def scale_by_batch(
    adata: AnnData,
    batch_key: str,
    *,
    zero_center: bool = True,
    max_value: float | None = None,
    copy: bool = False,
    layer: str | None = None,
    obsm: str | None = None,
    mask_obs: np.ndarray | str | None = None,
    inplace: bool = True,
) -> None | cp.ndarray:
    """
    Scales matrix to unit variance and clips values

    Parameters
    ----------
        adata
            AnnData object

        zero_center
            If `False`, omit zero-centering variables, which allows to handle sparse
            input efficiently.

        max_value
            Clip (truncate) to this value after scaling. If `None`, do not clip.

        copy
            Whether this function should be performed inplace. If an AnnData object
            is passed, this also determines if a copy is returned.

        layer
            If provided, which element of layers to scale.

        obsm
            If provided, which element of obsm to scale.

        mask_obs
            Restrict both the derivation of scaling parameters and the scaling itself
            to a certain set of observations. The mask is specified as a boolean array
            or a string referring to an array in :attr:`~anndata.AnnData.obs`. If the matrix is in csc format and a mask is provided, the matrix will be transformed to csr format.

        inplace
            If True, update AnnData with results. Otherwise, return results. See below for details of what is returned.

    Returns
    -------
    Returns a sacled copy or updates `adata` with a scaled version of the original `adata.X` and `adata.layers['layer']`, \
    depending on `inplace`.

    """
    if copy:
        if not inplace:
            raise ValueError("`copy=True` cannot be used with `inplace=False`.")
        adata = adata.copy()

    if isinstance(adata, AnnData):
        view_to_actual(adata)

    X = _get_obs_rep(adata, layer=layer, obsm=obsm)
    _check_gpu_X(X, require_cf=True)

    if mask_obs is not None:
        mask_obs = _check_mask(adata, mask_obs, "obs")

    batch_idx = adata.obs[batch_key].unique()
    X_list = []
    for batch in batch_idx:
        X_batch = X[adata.obs[batch_key] == batch]
        if isinstance(X_batch, cp.ndarray):
            X_batch, means, std = scale_array(
                X_batch, mask_obs=mask_obs, zero_center=zero_center, inplace=inplace
            )
        else:
            X_batch, means, std = scale_sparse(
                X_batch, mask_obs=mask_obs, zero_center=zero_center, inplace=inplace
            )
        X_list.append(X_batch)

    X = cp.vstack(X_list)

    if max_value:
        if zero_center:
            X = cp.clip(X, a_min=-max_value, a_max=max_value)
        else:
            if isinstance(X, sparse.spmatrix):
                X.data[X.data > max_value] = max_value
            else:
                X[X > max_value] = max_value

    if inplace:
        _set_obs_rep(adata, X, layer=layer, obsm=obsm)

    if copy:
        return adata
    elif not inplace:
        return X


# %%
def clip_to_percentile(array, percent: float = 0.1):
    assert 0 < percent < 50
    half_percent = percent / 2
    percentile_down = np.percentile(array, half_percent, axis=0)
    percentile_up = np.percentile(array, 100 - half_percent, axis=0)
    return np.clip(array, a_min=percentile_down, a_max=percentile_up)


# %%
def read_adata(
    path, min_genes: int = 50, max_scale_value: float = None, rapids: bool = False
) -> AnnData:
    adata = sc.read_h5ad(path)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata.var.index = adata.var["gene"]
    # sc.pp.scale(adata, max_value=max_scale_value)
    return adata.copy()


# %%
root_dir = "../../../SpatialFormer-dev/data/processed/RNA/Human/"
meta = pd.read_csv("../../data/spatial_meta.csv")
# select 'Disease' is 'Cancer' and Technology is 'MERSCOPE'
meta = meta[(meta["Disease"] == "Cancer") & (meta["Technology"] == "MERSCOPE")]
out_dir = Path("../../data/processed/human/pancancer/MERSCOPE/nine_cancers/")

# %%
marker_list = [
    "FOXM1",
    "BIRC5",
    "ERBB3",
    "MKI67",
    "VCAM1",
    "LAMB3",  # cancer
    "CD207",
    "ITGAX",  # dc
    "COL4A1",
    "PLVAP",  # endothelial
    "COL1A1",
    "COL5A1",
    "COL11A1",  # fibroblast
    "CD79A",
    "XBP1",
    "MZB1",
    "POU2AF1",  # b cell
    "CD3D",
    "CD3E",
    "CD3G",  # t cell
    "CD14",
    "FCGR3A",
    "CYBB",
    "C1QC",
    "SPP1",  # macrophage
    "KIT",
    "CD22",  # mast cell,
    "MYH11",
    "ACTA2",  # smooth muscle
    "PLA2G2A",  # epithelial
    "S100A9",
    "PTPRC",  # other
]

# %%
# all_gene = pd.DataFrame(adata.var.index)

# %% [markdown]
# ## Breast

# %%
file_list = []
for i, row in meta.iterrows():
    path = (
        Path(root_dir)
        / row["Tissue_0"]
        / f"{row['Tissue_1']}_{row['Disease']}_{row['Platform']}_{row['Technology']}_{str(row['Index'])}"
        / "processed.h5ad"
    )
    if path.exists():
        if "breast" in str(path).lower():
            print(path)
            file_list.append(path)
print(len(file_list))
if len(file_list) > 1:
    adata_list = joblib.Parallel(n_jobs=20)(
        joblib.delayed(read_adata)(file) for file in file_list
    )
    adata = adata_list[0].concatenate(adata_list[1:], batch_key="batch", join="inner")
else:
    adata = read_adata(file_list[0], min_genes=50)

# %%
df = pd.DataFrame({"sum": np.array(adata.X.sum(1)).flatten()})
# set x-axis range 0 to 1000
df["sum"].hist(bins=100, range=(0, 1800))

# %%
adata_raw = adata.copy()

# %%
adata = adata_raw.copy()

# %%
adata = gex_embedding(adata, call_hvg=False, viz=True, resolution=0.8)
adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()

# %%
adata.obsm["X_umap_clip"] = clip_to_percentile(adata.obsm["X_umap_raw"])

# %%
rsc.tl.leiden(adata, resolution=0.8)

# %%
adata.obsm["X_umap"] = adata.obsm["X_umap_clip"]
sc.pl.umap(adata, color=["leiden"], legend_loc="on data")

# %%
rsc.get.anndata_to_CPU(adata)

# %%
marker_list = [
    "FOXM1",
    "BIRC5",
    "ERBB3",
    "MKI67",
    "VCAM1",
    "LAMB3",  # cancer
    "CD207",  # dc
    "COL4A1",
    "PLVAP",  # endothelial
    "COL1A1",
    "COL11A1",  # fibroblast
    "CD79A",
    "XBP1",
    "MZB1",
    "POU2AF1",  # b cell
    "CD3D",  # t cell
    "CD14",
    "FCGR3A",
    "CYBB",
    "C1QC",
    "SPP1",  # macrophage
    "KIT",
    "CD22",  # mast cell
]

# %%
sc.pl.dotplot(adata, marker_list, groupby="leiden", dendrogram=False)

# %%
# annotate
adata.obs["cell_type"] = "unknown"
adata.obs.loc[
    adata.obs["leiden"].isin(["2", "5", "6", "9", "11"]), "cell_type"
] = "Cancer"
adata.obs.loc[adata.obs["leiden"].isin(["4", "7"]), "cell_type"] = "TAM"
adata.obs.loc[adata.obs["leiden"].isin(["12"]), "cell_type"] = "DC"
adata.obs.loc[adata.obs["leiden"].isin(["0"]), "cell_type"] = "B/Plasma cell"
adata.obs.loc[adata.obs["leiden"].isin(["8"]), "cell_type"] = "Mast cell"
adata.obs.loc[adata.obs["leiden"].isin(["1"]), "cell_type"] = "T cell"
adata.obs.loc[adata.obs["leiden"].isin(["3"]), "cell_type"] = "Endothelial"
adata.obs.loc[adata.obs["leiden"].isin(["10"]), "cell_type"] = "Fibroblast"
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color=["cell_type", "n_genes"], legend_loc="on data", ncols=1)

# %%
sc.pl.dotplot(adata, marker_list, groupby="cell_type", dendrogram=False, vmax=2)

# %%
adata.layers["counts"] = adata_raw.X.copy()
adata.write(out_dir / "breast_filter_anno.h5ad")

# %% [markdown]
# ## Colon

# %%
file_list = []
for i, row in meta.iterrows():
    path = (
        Path(root_dir)
        / row["Tissue_0"]
        / f"{row['Tissue_1']}_{row['Disease']}_{row['Platform']}_{row['Technology']}_{str(row['Index'])}"
        / "processed.h5ad"
    )
    if path.exists():
        if "colon" in str(path).lower():
            print(path)
            file_list.append(path)
print(len(file_list))
if len(file_list) > 1:
    adata_list = joblib.Parallel(n_jobs=len(file_list))(
        joblib.delayed(read_adata)(file) for file in file_list
    )
    adata = adata_list[0].concatenate(adata_list[1:], batch_key="batch", join="inner")
else:
    adata = read_adata(file_list[0], min_genes=50)

# %%
adata_raw = adata.copy()

# %%
adata = adata_raw.copy()

# %%
df = pd.DataFrame({"sum": np.array(adata.X.sum(1)).flatten()})
# set x-axis range 0 to 1000
df["sum"].hist(bins=100, range=(0, 1800))

# %%
adata = gex_embedding(adata, batch_key="batch", call_hvg=False, viz=True, resolution=0.8)
adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()

# %%
adata.obsm["X_umap_clip"] = clip_to_percentile(adata.obsm["X_umap_raw"])

# %%
rsc.tl.leiden(adata, resolution=0.8)

# %%
adata.obsm["X_umap"] = adata.obsm["X_umap_clip"]
sc.pl.umap(adata, color=["leiden"], legend_loc="on data")

# %%
sc.pl.umap(
    adata,
    color=["Index"],
    legend_loc="on data",
)

# %%
rsc.get.anndata_to_CPU(adata)

# %%
sc.pl.dotplot(adata, marker_list, groupby="leiden", dendrogram=False)

# %%
# annotate
adata.obs["cell_type"] = "unknown"
adata.obs.loc[
    adata.obs["leiden"].isin(["0", "3", "6", "7", "9", "10", "13"]), "cell_type"
] = "Cancer"
adata.obs.loc[adata.obs["leiden"].isin(["4"]), "cell_type"] = "TAM"
adata.obs.loc[adata.obs["leiden"].isin(["12"]), "cell_type"] = "DC"
adata.obs.loc[adata.obs["leiden"].isin(["5", "11"]), "cell_type"] = "B/Plasma cell"
adata.obs.loc[adata.obs["leiden"].isin(["8"]), "cell_type"] = "Mast cell"
adata.obs.loc[adata.obs["leiden"].isin(["12"]), "cell_type"] = "T cell"
adata.obs.loc[adata.obs["leiden"].isin(["1"]), "cell_type"] = "Endothelial"
adata.obs.loc[adata.obs["leiden"].isin(["2", "8"]), "cell_type"] = "Fibroblast"
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color=["cell_type", "n_genes"], legend_loc="on data", ncols=1)

# %%
adata.layers["counts"] = adata_raw.X.copy()
adata.write_h5ad(out_dir / "colon_filter_anno.h5ad")

# %% [markdown]
# ## Liver

# %%
file_list = []
for i, row in meta.iterrows():
    path = (
        Path(root_dir)
        / row["Tissue_0"]
        / f"{row['Tissue_1']}_{row['Disease']}_{row['Platform']}_{row['Technology']}_{str(row['Index'])}"
        / "processed.h5ad"
    )
    if path.exists():
        if "liver" in str(path).lower():
            print(path)
            file_list.append(path)
print(len(file_list))
if len(file_list) > 1:
    adata_list = joblib.Parallel(n_jobs=len(file_list))(
        joblib.delayed(read_adata)(file) for file in file_list
    )
    adata = adata_list[0].concatenate(adata_list[1:], batch_key="batch", join="inner")
else:
    adata = read_adata(file_list[0], min_genes=50)

# %%
adata_raw = adata.copy()

# %%
adata = adata_raw.copy()

# %%
df = pd.DataFrame({"sum": np.array(adata.X.sum(1)).flatten()})
# set x-axis range 0 to 1000
df["sum"].hist(bins=100, range=(0, 1800))

# %%
adata = gex_embedding(adata, batch_key="batch", call_hvg=False, viz=True, resolution=0.8)
adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()

# %%
adata.obsm["X_umap_clip"] = clip_to_percentile(adata.obsm["X_umap_raw"])

# %%
rsc.tl.leiden(adata, resolution=0.8)

# %%
adata.obsm["X_umap"] = adata.obsm["X_umap_clip"]
sc.pl.umap(adata, color=["leiden", "batch"], legend_loc="on data")

# %%
rsc.tl.rank_genes_groups_logreg(adata, groupby="leiden", use_raw=False)

# %%
sc.pl.rank_genes_groups(adata)

# %%
rsc.get.anndata_to_CPU(adata)

# %%
marker_liver = [
    # 'CYP3A4',  # hepatocyte
    "CD163",  # Kupffer cell
    "ACTA2",
    "SPP1",  # stellate
    "VWF",  # endothelial
    "LRP6",
    "SERPINA1",  # hepatocyte
    "EPCAM",  # tumor
    # 'KRT7', 'KRT19', # cholagiocyte
]
marker_list_liver = marker_list + marker_liver

# %%
sc.pl.dotplot(adata, marker_list_liver, groupby="leiden", dendrogram=False)

# %%
# # annotate
adata.obs["cell_type"] = "unknown"
adata.obs.loc[
    adata.obs["leiden"].isin(["3", "4", "6", "7", "9", "18", "19"]), "cell_type"
] = "Cancer"
adata.obs.loc[adata.obs["leiden"].isin(["1", "10"]), "cell_type"] = "TAM"
adata.obs.loc[adata.obs["leiden"].isin(["17"]), "cell_type"] = "Hepatocyte"
adata.obs.loc[adata.obs["leiden"].isin(["0"]), "cell_type"] = "DC"
adata.obs.loc[adata.obs["leiden"].isin(["8", "16"]), "cell_type"] = "B/Plasma cell"
adata.obs.loc[adata.obs["leiden"].isin(["15"]), "cell_type"] = "Mast cell"
adata.obs.loc[adata.obs["leiden"].isin(["5", "13"]), "cell_type"] = "T cell"
adata.obs.loc[adata.obs["leiden"].isin(["2", "11"]), "cell_type"] = "Endothelial"
adata.obs.loc[adata.obs["leiden"].isin(["12"]), "cell_type"] = "Fibroblast"
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color=["cell_type", "n_genes"], legend_loc="on data", ncols=1)

# %%
adata.layers["counts"] = adata_raw.X.copy()
adata.write_h5ad(out_dir / "liver_filter_anno.h5ad")

# %% [markdown]
# ## Lung

# %%
file_list = []
for i, row in meta.iterrows():
    path = (
        Path(root_dir)
        / row["Tissue_0"]
        / f"{row['Tissue_1']}_{row['Disease']}_{row['Platform']}_{row['Technology']}_{str(row['Index'])}"
        / "processed.h5ad"
    )
    if path.exists():
        if "lung" in str(path).lower():
            print(path)
            file_list.append(path)
print(len(file_list))
if len(file_list) > 1:
    adata_list = joblib.Parallel(n_jobs=len(file_list))(
        joblib.delayed(read_adata)(file) for file in file_list
    )
    adata = adata_list[0].concatenate(adata_list[1:], batch_key="batch", join="inner")
else:
    adata = read_adata(file_list[0], min_genes=50)

# %%
adata_raw = adata.copy()

# %%
adata = adata_raw.copy()

# %%
df = pd.DataFrame({"sum": np.array(adata.X.sum(1)).flatten()})
# set x-axis range 0 to 1000
df["sum"].hist(bins=100, range=(0, 1800))

# %%
adata = gex_embedding(adata, batch_key="batch", call_hvg=False, viz=True, resolution=0.8)
adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()

# %%
adata.obsm["X_umap_clip"] = clip_to_percentile(adata.obsm["X_umap_raw"])

# %%
rsc.tl.leiden(adata, resolution=0.8)

# %%
adata.obsm["X_umap"] = adata.obsm["X_umap_clip"]
sc.pl.umap(adata, color=["leiden", "batch"], legend_loc="on data")

# %%
rsc.tl.rank_genes_groups_logreg(adata, groupby="leiden", use_raw=False)

# %%
sc.pl.rank_genes_groups(adata)

# %%
rsc.get.anndata_to_CPU(adata)

# %%
marker_lung = [
    # 'CYP3A4',  # hepatocyte
    "CD163",  # Kupffer cell
    "ACTA2",
    "SPP1",  # stellate
    "VWF",  # endothelial
    "LRP6",
    "SERPINA1",  # pneumocyte
    "EPCAM",  # tumor,
    "CXCL1"
    # 'KRT7', 'KRT19', # cholagiocyte
]
marker_list_lung = marker_list + marker_lung

# %%
sc.pl.dotplot(adata, marker_list_lung, groupby="leiden", dendrogram=False)

# %%
# # annotate
adata.obs["cell_type"] = "unknown"
adata.obs.loc[
    adata.obs["leiden"].isin(["0", "3", "4", "7", "10", "11", "17"]), "cell_type"
] = "Cancer"
adata.obs.loc[adata.obs["leiden"].isin(["1", "15"]), "cell_type"] = "TAM"
adata.obs.loc[adata.obs["leiden"].isin(["9"]), "cell_type"] = "pneumocyte"
adata.obs.loc[adata.obs["leiden"].isin(["6", "18"]), "cell_type"] = "B/Plasma cell"
adata.obs.loc[adata.obs["leiden"].isin(["13"]), "cell_type"] = "Mast cell"
adata.obs.loc[adata.obs["leiden"].isin(["5", "8", "16"]), "cell_type"] = "T cell"
adata.obs.loc[adata.obs["leiden"].isin(["2", "14"]), "cell_type"] = "Endothelial"
adata.obs.loc[adata.obs["leiden"].isin(["12"]), "cell_type"] = "Fibroblast"
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color=["cell_type", "n_genes"], legend_loc="on data", ncols=1)

# %%
# import matplotlib.pyplot as plt
# plt.figure(figsize=(15, 15))
# sc.pl.spatial(adata, color = ['cell_type'], spot_size=3)

# %%
adata.layers["counts"] = adata_raw.X.copy()
adata.write_h5ad(out_dir / "lung_filter_anno.h5ad")

# %% [markdown]
# ## Ovarian

# %%
file_list = []
for i, row in meta.iterrows():
    path = (
        Path(root_dir)
        / row["Tissue_0"]
        / f"{row['Tissue_1']}_{row['Disease']}_{row['Platform']}_{row['Technology']}_{str(row['Index'])}"
        / "processed.h5ad"
    )
    if path.exists():
        if "ovarian" in str(path).lower():
            print(path)
            file_list.append(path)
print(len(file_list))
if len(file_list) > 1:
    adata_list = joblib.Parallel(n_jobs=len(file_list))(
        joblib.delayed(read_adata)(file) for file in file_list
    )
    adata = adata_list[0].concatenate(adata_list[1:], batch_key="batch", join="inner")
else:
    adata = read_adata(file_list[0], min_genes=50)

# %%
adata_raw = adata.copy()

# %%
adata = adata_raw.copy()

# %%
df = pd.DataFrame({"sum": np.array(adata.X.sum(1)).flatten()})
# set x-axis range 0 to 1000
df["sum"].hist(bins=100, range=(0, 1800))

# %%
adata = gex_embedding(adata, batch_key="batch", call_hvg=False, viz=True, resolution=0.8)
adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()

# %%
adata.obsm["X_umap_clip"] = clip_to_percentile(adata.obsm["X_umap_raw"])

# %%
rsc.tl.leiden(adata, resolution=0.8)

# %%
adata.obsm["X_umap"] = adata.obsm["X_umap_clip"]
sc.pl.umap(adata, color=["leiden", "batch"], legend_loc="on data")

# %%
rsc.tl.rank_genes_groups_logreg(adata, groupby="leiden", use_raw=False)

# %%
sc.pl.rank_genes_groups(adata)

# %%
rsc.get.anndata_to_CPU(adata)

# %%
marker_ov = [
    # 'CYP3A4',  # hepatocyte
    "CD163",  # Kupffer cell
    "ACTA2",
    "SPP1",  # stellate
    "VWF",  # endothelial
    "LRP6",
    "SERPINA1",  # pneumocyte
    "EPCAM",  # tumor,
    "CXCL1"
    # 'KRT7', 'KRT19', # cholagiocyte
]
marker_list_ov = marker_list + marker_ov

# %%
sc.pl.dotplot(adata, marker_list_ov, groupby="leiden", dendrogram=False)

# %%
# # annotate
adata.obs["cell_type"] = "unknown"
adata.obs.loc[
    adata.obs["leiden"].isin(["4", "5", "6", "9", "12", "15", "16"]), "cell_type"
] = "Cancer"
adata.obs.loc[adata.obs["leiden"].isin(["2", "10", "13"]), "cell_type"] = "TAM"
adata.obs.loc[adata.obs["leiden"].isin(["0"]), "cell_type"] = "DC"
adata.obs.loc[adata.obs["leiden"].isin(["1"]), "cell_type"] = "Mast cell"
adata.obs.loc[adata.obs["leiden"].isin(["17"]), "cell_type"] = "T cell"
adata.obs.loc[adata.obs["leiden"].isin(["7", "14"]), "cell_type"] = "Endothelial"
adata.obs.loc[adata.obs["leiden"].isin(["3", "8", "11"]), "cell_type"] = "Fibroblast"
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color=["cell_type", "n_genes"], legend_loc="on data", ncols=1)

# %%
adata.layers["counts"] = adata_raw.X.copy()
adata.write_h5ad(out_dir / "ovarian_filter_anno.h5ad")

# %% [markdown]
# ## Prostate

# %%
file_list = []
for i, row in meta.iterrows():
    path = (
        Path(root_dir)
        / row["Tissue_0"]
        / f"{row['Tissue_1']}_{row['Disease']}_{row['Platform']}_{row['Technology']}_{str(row['Index'])}"
        / "processed.h5ad"
    )
    if path.exists():
        if "prostate" in str(path).lower():
            print(path)
            file_list.append(path)
print(len(file_list))
if len(file_list) > 1:
    adata_list = joblib.Parallel(n_jobs=len(file_list))(
        joblib.delayed(read_adata)(file) for file in file_list
    )
    adata = adata_list[0].concatenate(adata_list[1:], batch_key="batch", join="inner")
else:
    adata = read_adata(file_list[0], min_genes=50)

# %%
adata_raw = adata.copy()

# %%
adata = adata_raw.copy()

# %%
df = pd.DataFrame({"sum": np.array(adata.X.sum(1)).flatten()})
# set x-axis range 0 to 1000
df["sum"].hist(bins=100, range=(0, 1800))

# %%
adata = gex_embedding(adata, batch_key="batch", call_hvg=False, viz=True, resolution=0.8)
adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()

# %%
adata.obsm["X_umap_clip"] = clip_to_percentile(adata.obsm["X_umap_raw"])

# %%
rsc.tl.leiden(adata, resolution=1)

# %%
adata.obsm["X_umap"] = adata.obsm["X_umap_clip"]
sc.pl.umap(adata, color=["leiden", "batch"], legend_loc="on data")

# %%
rsc.tl.rank_genes_groups_logreg(adata, groupby="leiden", use_raw=False)

# %%
sc.pl.rank_genes_groups(adata)

# %%
rsc.get.anndata_to_CPU(adata)

# %%
marker_prostate = [
    # 'CYP3A4',  # hepatocyte
    "CD163",  # Kupffer cell
    "ACTA2",
    "SPP1",  # stellate
    "VWF",  # endothelial
    "LRP6",
    "SERPINA1",  # pneumocyte
    "EPCAM",  # tumor,
    "CXCL1"
    # 'KRT7', 'KRT19', # cholagiocyte
]
marker_list_prostate = marker_list + marker_prostate

# %%
sc.pl.dotplot(adata, marker_list_prostate, groupby="leiden", dendrogram=False)

# %%
# # annotate
adata.obs["cell_type"] = "unknown"
adata.obs.loc[adata.obs["leiden"].isin(["12", "15"]), "cell_type"] = "B/Plasma cell"
adata.obs.loc[
    adata.obs["leiden"].isin(["1", "2", "6", "7", "8", "10", "14", "16"]), "cell_type"
] = "Cancer"
adata.obs.loc[adata.obs["leiden"].isin(["0"]), "cell_type"] = "TAM"
adata.obs.loc[adata.obs["leiden"].isin(["13"]), "cell_type"] = "Mast cell"
adata.obs.loc[adata.obs["leiden"].isin(["5"]), "cell_type"] = "T cell"
adata.obs.loc[adata.obs["leiden"].isin(["9"]), "cell_type"] = "Endothelial"
adata.obs.loc[adata.obs["leiden"].isin(["3", "4", "11"]), "cell_type"] = "Fibroblast"
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color=["cell_type", "n_genes"], legend_loc="on data", ncols=1)

# %%
adata.layers["counts"] = adata_raw.X.copy()
adata.write_h5ad(out_dir / "prosate_filter_anno.h5ad")

# %% [markdown]
# ## Skin

# %%
file_list = []
for i, row in meta.iterrows():
    path = (
        Path(root_dir)
        / row["Tissue_0"]
        / f"{row['Tissue_1']}_{row['Disease']}_{row['Platform']}_{row['Technology']}_{str(row['Index'])}"
        / "processed.h5ad"
    )
    if path.exists():
        if "skin" in str(path).lower():
            print(path)
            file_list.append(path)
print(len(file_list))
if len(file_list) > 1:
    adata_list = joblib.Parallel(n_jobs=len(file_list))(
        joblib.delayed(read_adata)(file) for file in file_list
    )
    adata = adata_list[0].concatenate(adata_list[1:], batch_key="batch", join="inner")
else:
    adata = read_adata(file_list[0], min_genes=50)

# %%
adata_raw = adata.copy()

# %%
adata = adata_raw.copy()

# %%
df = pd.DataFrame({"sum": np.array(adata.X.sum(1)).flatten()})
# set x-axis range 0 to 1000
df["sum"].hist(bins=100, range=(0, 1800))

# %%
adata = gex_embedding(adata, batch_key="batch", call_hvg=False, viz=True, resolution=0.8)
adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()

# %%
adata.obsm["X_umap_clip"] = clip_to_percentile(adata.obsm["X_umap_raw"])

# %%
rsc.tl.leiden(adata, resolution=0.8)

# %%
adata.obsm["X_umap"] = adata.obsm["X_umap_clip"]
sc.pl.umap(adata, color=["leiden", "batch"], legend_loc="on data")

# %%
rsc.tl.rank_genes_groups_logreg(adata, groupby="leiden", use_raw=False)

# %%
sc.pl.rank_genes_groups(adata)

# %%
rsc.get.anndata_to_CPU(adata)

# %%
marker_skin = [
    # 'CYP3A4',  # hepatocyte
    "CD163",  # Kupffer cell
    "ACTA2",
    "SPP1",  # stellate
    "VWF",  # endothelial
    "LRP6",
    "SERPINA1",  # pneumocyte
    "EPCAM",  # tumor,
    "CXCL1"
    # 'KRT7', 'KRT19', # cholagiocyte
]
marker_list_skin = marker_list + marker_skin

# %%
sc.pl.dotplot(adata, marker_list_skin, groupby="leiden", dendrogram=False)

# %%
# # annotate
adata.obs["cell_type"] = "unknown"
adata.obs.loc[adata.obs["leiden"].isin(["7", "11"]), "cell_type"] = "B/Plasma cell"
adata.obs.loc[
    adata.obs["leiden"].isin(["2", "4", "3", "5", "9", "14", "15", "16"]), "cell_type"
] = "Cancer"
adata.obs.loc[adata.obs["leiden"].isin(["1", "13"]), "cell_type"] = "TAM"
adata.obs.loc[adata.obs["leiden"].isin(["8"]), "cell_type"] = "Mast cell"
adata.obs.loc[adata.obs["leiden"].isin(["12"]), "cell_type"] = "T cell"
adata.obs.loc[adata.obs["leiden"].isin(["0", "10"]), "cell_type"] = "Endothelial"
adata.obs.loc[adata.obs["leiden"].isin(["6"]), "cell_type"] = "Fibroblast"
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color=["cell_type", "n_genes"], legend_loc="on data", ncols=1)

# %%
adata.layers["counts"] = adata_raw.X.copy()
adata.write_h5ad(out_dir / "skin_filter_anno.h5ad")

# %% [markdown]
# ## Uterine

# %%
file_list = []
for i, row in meta.iterrows():
    path = (
        Path(root_dir)
        / row["Tissue_0"]
        / f"{row['Tissue_1']}_{row['Disease']}_{row['Platform']}_{row['Technology']}_{str(row['Index'])}"
        / "processed.h5ad"
    )
    if path.exists():
        if "uterine" in str(path).lower():
            print(path)
            file_list.append(path)
print(len(file_list))
if len(file_list) > 1:
    adata_list = joblib.Parallel(n_jobs=len(file_list))(
        joblib.delayed(read_adata)(file) for file in file_list
    )
    adata = adata_list[0].concatenate(adata_list[1:], batch_key="batch", join="inner")
else:
    adata = read_adata(file_list[0], min_genes=50)

# %%
adata_raw = adata.copy()

# %%
adata = adata_raw.copy()

# %%
df = pd.DataFrame({"sum": np.array(adata.X.sum(1)).flatten()})
# set x-axis range 0 to 1000
df["sum"].hist(bins=100, range=(0, 1800))

# %%
adata = gex_embedding(adata, batch_key="batch", call_hvg=False, viz=True, resolution=0.8)
adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()

# %%
adata.obsm["X_umap_clip"] = clip_to_percentile(adata.obsm["X_umap_raw"])

# %%
rsc.tl.leiden(adata, resolution=1)

# %%
adata.obsm["X_umap"] = adata.obsm["X_umap_clip"]
sc.pl.umap(adata, color=["leiden", "batch"], legend_loc="on data")

# %%
rsc.tl.rank_genes_groups_logreg(adata, groupby="leiden", use_raw=False)

# %%
sc.pl.rank_genes_groups(adata)

# %%
rsc.get.anndata_to_CPU(adata)

# %%
marker_ut = [
    # 'CYP3A4',  # hepatocyte
    "CD163",  # Kupffer cell
    "ACTA2",
    "SPP1",  # stellate
    "VWF",  # endothelial
    "LRP6",
    "SERPINA1",  # pneumocyte
    "EPCAM",  # tumor,
    "CXCL1"
    # 'KRT7', 'KRT19', # cholagiocyte
]
marker_list_ut = marker_list + marker_ut

# %%
sc.pl.dotplot(adata, marker_list_skin, groupby="leiden", dendrogram=False)

# %%
# # annotate
adata.obs["cell_type"] = "unknown"
adata.obs.loc[adata.obs["leiden"].isin(["1", "5"]), "cell_type"] = "B/Plasma cell"
adata.obs.loc[
    adata.obs["leiden"].isin(["3", "4", "6", "7", "10", "11", "12", "13", "16"]),
    "cell_type",
] = "Cancer"
adata.obs.loc[adata.obs["leiden"].isin(["9"]), "cell_type"] = "TAM"
adata.obs.loc[adata.obs["leiden"].isin(["8"]), "cell_type"] = "T cell"
adata.obs.loc[adata.obs["leiden"].isin(["0"]), "cell_type"] = "Endothelial"
adata.obs.loc[adata.obs["leiden"].isin(["14"]), "cell_type"] = "Fibroblast"
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color=["cell_type", "n_genes"], legend_loc="on data", ncols=1)

# %%
adata.layers["counts"] = adata_raw.X.copy()
adata.write_h5ad(out_dir / "uterine_filter_anno.h5ad")

# %% [markdown]
# # Merge the annotated data

# %%
adata_name = glob.glob(str(out_dir / "*filter*"))
# remove str contrain 'pancancer_fiter'
adata_name = [name for name in adata_name if "pancancer_filter" not in name]
adata_name

# %%
adata_list = joblib.Parallel(n_jobs=8)(
    joblib.delayed(read_adata)(file) for file in adata_name
)
adata = adata_list[0].concatenate(adata_list[1:], batch_key="dataset", join="inner")

# %%
adata

# %%
sc.pp.pca(adata, n_comps=50)

# %%
# rsc.tl.pca(adata, n_comps=50) # can not run this on GPU because memory limitation
rsc.get.anndata_to_GPU(adata)
rsc.pp.neighbors(adata, algorithm = 'cagra')
rsc.tl.umap(adata)
rsc.get.anndata_to_CPU(adata)
# rsc.tl.leiden(adata, resolution=1)

# %%
adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()

# %%
adata_plot = adata[~adata.obs['cell_type'].isin(['unknown', 'Hepatocyte', 'pneumocyte']), :].copy()
adata_plot.obsm["X_umap"] = clip_to_percentile(adata_plot.obsm["X_umap_raw"], percent=0.001)
sc.pl.umap(adata_plot, color=["cell_type", "Tissue_0", "Index"], ncols=1)
del adata_plot

# %%
adata.write_h5ad( "./data/pancancer_filter_anno.h5ad")

# %%
adata = sc.read_h5ad('./data/pancancer_filter_anno.h5ad')

# %%
adata.obs['dataset'].astype(int).value_counts()

# %% [markdown]
# # Build edges

# %%
from spider.graphic.build import build_graph
batch = adata.obs["Index"].astype(int).to_numpy()
edge_index = build_graph(adata.obsm["spatial"], batch, mode="knn", k=20)

# %%
import torch
torch.save(edge_index, "./edge_index.pt")
