r"""
Ligand and receptor prior
"""
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from loguru import logger
from scipy.sparse import isspmatrix
from torch_geometric.nn import SimpleConv

from ...graphic.build import build_graph


def get_lr_expr(
    adata: sc.AnnData,
    lr_df: pd.DataFrame,
    radius: float,
    aggr: str = "mean",
    binary: bool = False,
    threshold: float = 0.01,
) -> tuple[torch.Tensor, pd.DataFrame]:
    r"""
    Get ligand-receptor activity from spatial data.

    We multiply mean ligand expression in given radius and receptor expression in the cell as the activity of a ligand-receptor pair.

    Args:
        adata: AnnData object, adata.X is log-norm data
        lr_df: Ligand receptor pairs, must contain columns `ligand.symbol` and `receptor.symbol`
        radius: float radius for valid ligand
        aggr: str aggregation method of ligand expression, default is `mean`
        binary: bool if convert the expression to binary, default is `False`
        threshold: float threshold to filter low-expressed LR pairs, default is `0.01`
    Returns:
        LR activity: Ligand-receptor pair activity
        lr_df: filtered ligand receptor pairs
    Note:
        Not support batched data yet.
    """
    # check
    assert "ligand.symbol" in lr_df.columns
    assert "receptor.symbol" in lr_df.columns

    # filter non-LR genes (reduce memory usage)
    lr_genes = []
    for _, row in lr_df.iterrows():
        ligand = row["ligand.symbol"]
        if "," in ligand:
            lr_genes += ligand.split(",")
        else:
            lr_genes.append(ligand)

        receptor = row["receptor.symbol"]
        if "," in receptor:
            lr_genes += receptor.split(",")
        else:
            lr_genes.append(receptor)
    lr_genes = list(set(lr_genes))
    lr_idx = [x in lr_genes for x in adata.var.index]
    adata = adata[:, lr_idx]

    logger.info(f"Find {adata.n_vars} LR-related genes.")

    # merge neighbor expr
    edge_index = build_graph(adata.obsm["spatial"], radius=radius, mode="radius")
    aggr_net = SimpleConv(aggr=aggr, combine_root=None)
    expr = torch.from_numpy(adata.X.toarray()) if isspmatrix(adata.X) else torch.from_numpy(adata.X)
    # expr = expr.to(torch.float32)
    expr_neighbors = aggr_net(expr, edge_index)

    # get the ligand and receptor gene index in adata.var
    lr_activity_list = []
    gene_set = set(adata.var.index)
    lr_filter = np.zeros(len(lr_df), dtype=bool)
    for i, lr_pair in lr_df.iterrows():
        ligand = lr_pair["ligand.symbol"]
        receptor = lr_pair["receptor.symbol"]
        ligand = [ligand] if "," not in ligand else ligand.split(",")
        receptor = [receptor] if "," not in receptor else receptor.split(",")
        if len(set(ligand) & gene_set) == 0 or len(set(receptor) & gene_set) == 0:
            continue
        # get gene id
        ligand_idx = [adata.var.index.get_loc(x) for x in ligand if x in gene_set]
        receptor_idx = [adata.var.index.get_loc(x) for x in receptor if x in gene_set]
        ligand_expr = expr_neighbors[:, ligand_idx].mean(dim=1)  # from neighbor
        receptor_expr = expr[:, receptor_idx].mean(dim=1)  # from cell
        lr_activity = ligand_expr * receptor_expr
        if (lr_activity > 0).float().mean().item() < threshold:
            continue
        lr_activity_list.append(ligand_expr * receptor_expr)
        lr_filter[i] = True
    if len(lr_activity_list) == 0:
        logger.warning(f"No ligand-receptor pairs (threshold: {threshold}).")
        return None, None
    lr_activity = torch.stack(lr_activity_list, dim=1)
    if binary:
        lr_activity = (lr_activity > 0).float()
    lr_df = lr_df[lr_filter]
    return lr_activity, lr_df
