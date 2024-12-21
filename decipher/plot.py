r"""
Plot functions
"""
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from anndata import AnnData
from loguru import logger


def plot_sc(adata: AnnData, color_vars: list[str], suffix: str = "") -> None:
    r"""
    Plot single cell data

    Parameters
    ----------
    adata
        AnnData object
    color_vars
        variables to plot
    suffix
        suffix of the plot name
    """
    color_vars = [color_vars] if isinstance(color_vars, str) else color_vars
    suffix = suffix if suffix == "" else f"__{suffix}"
    for color in color_vars:
        name = f"coords:spatial__color:{color}" + suffix
        kwargs = dict(adata=adata, color=color, save=f"__{name}.pdf", title=name)
        sc.pl.spatial(**kwargs, spot_size=adata.uns["spot_size"])
        for key in ["center", "nbr", "merge"]:
            umap_key = f"X_umap_{key}"
            if umap_key not in adata.obsm.keys():
                logger.warning(f"{umap_key} not in adata.obsm")
                continue
            adata.obsm["X_umap"] = adata.obsm[umap_key].copy()
            name = f"coords:{key}umap_color:{color}" + suffix
            kwargs["save"] = f"__{name}.pdf"
            kwargs["title"] = name
            sc.pl.umap(**kwargs)


def split_umap(adata, split_by, ncol=1, nrow=None, **kwargs):
    r"""
    Split umap by split_by like Seurat

    Parameters
    ----------
    adata
        AnnData object
    split_by
        split by variable
    ncol
        number of columns
    nrow
        number of rows
    **kwargs
        other parameters for `sc.pl.umap`
    """
    categories = adata.obs[split_by].cat.categories
    # print(categories)
    if nrow is None:
        nrow = int(np.ceil(len(categories) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4 * nrow))
    axs = axs.flatten()
    for i, cat in enumerate(categories):
        ax = axs[i]
        sc.pl.umap(adata[adata.obs[split_by] == cat], ax=ax, show=False, title=cat, **kwargs)
    plt.tight_layout()
