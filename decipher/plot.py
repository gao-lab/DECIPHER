r"""
Plot functions
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
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


def Sankey(
    matching_table: pd.DataFrame,
    filter_num: int = 50,
    color: list[str] = "red",
    title: str = "",
    layout: list[int] = [1300, 900],
    font_size: float = 15,
    font_color: str = "Black",
    save_name: str = None,
    format: str = "png",
    width: int = 1200,
    height: int = 1000,
    return_fig: bool = False,
) -> None:
    r"""
    Sankey plot of celltype

    Parameters
    ----------
    matching_table
        matching table
    filter_num
        filter number of matches
    color
        color of node
    title
        plot title
    layout
        layout size of picture
    font_size
        font size in plot
    font_color
        font color in plot
    save_name
        save file name (None for not save)
    format
        save picture format (see https://plotly.com/python/static-image-export/ for more details)
    width
        save picture width
    height
        save picture height
    return_fig
        if return plotly figure
    """
    source, target, value = [], [], []
    label_ref = matching_table.columns.to_list()
    label_query = matching_table.index.to_list()
    label_all = label_query + label_ref
    label2index = dict(zip(label_all, list(range(len(label_all)))))

    for i, query in enumerate(label_query):
        for j, ref in enumerate(label_ref):
            if int(matching_table.iloc[i, j]) > filter_num:
                target.append(label2index[query])
                source.append(label2index[ref])
                value.append(int(matching_table.iloc[i, j]))

    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=50,
                    thickness=50,
                    line=dict(color="green", width=0.5),
                    label=label_all,
                    color=color,
                ),
                link=dict(
                    source=source,  # indices correspond to labels, eg A1, A2, A1, B1, ...
                    target=target,
                    value=value,
                ),
            )
        ],
        layout=go.Layout(autosize=False, width=layout[0], height=layout[1]),
    )

    fig.update_layout(title_text=title, font_size=font_size, font_color=font_color)
    fig.show()
    if save_name is not None:
        fig.write_image(save_name + f".{format}", width=width, height=height)
    if return_fig:
        return fig
