r"""
Graph building
"""

import numpy as np
import torch
from loguru import logger
from torch import Tensor
from torch_geometric.nn import knn_graph, radius_graph
from torch_geometric.utils import to_undirected

from ..utils import estimate_spot_distance
from .knn import knn


def build_graph(
    coords: np.ndarray,
    batch: np.ndarray = None,
    mode: str = "radius",
    radius: float = None,
    max_num_neighbors: int = 64,
    k: int = 10,
    pyg_backend: bool = True,
    self_loop: bool = False,
    undirected: bool = False,
    num_workers: bool = 1,
) -> Tensor:
    r"""
    Build cell graph of spatial slices

    Parameters
    ----------
    coords
        coordinates of cells
    batch
        batch of cells
    mode
        mode of building graph, 'knn' or 'radius', default is 'radius'
    radius
        radius of radius graph
    max_num_neighbors
        max number of neighbors for radius graph
    k
        number of nearest neighbors
    pyg_backend
        whether to use PyG backend to build graph
    self_loop
        whether to add self loop
    undirected
        whether to build undirected graph
    num_workers
        number of workers (only for PyG backend)

    Returns
    ----------
    edge_index of PyG graph
    """
    assert mode.lower() in ["radius", "knn"], f"Unsupported mode: {mode}"

    if mode == "radius" or batch is not None:
        if not pyg_backend:
            logger.warning("Other backend is not available, switch to PyG backend.")
        pyg_backend = True

    if not pyg_backend and mode == "knn":
        knn_result, _ = knn(coords, k=k, metric="euclidean")
        edge_index = knn_to_edge_index(knn_result)
        edge_index = to_undirected(edge_index) if undirected else edge_index
        assert edge_index.shape[1] > 0, "No edges contained in spatial graph."
        return edge_index

    coords = torch.from_numpy(coords)
    if batch is not None:
        logger.debug("Build graph with batch.")
        batch = torch.tensor(batch, dtype=torch.long)
        assert batch.shape[0] == coords.shape[0], "Batch length != coords length."
    if mode == "knn":
        edge_index = knn_graph(
            x=coords,
            k=k,
            batch=batch,
            loop=self_loop,
            num_workers=num_workers,
        )
    elif mode == "radius":
        assert isinstance(radius, (int, float)), "Radius must be a number."
        min_distance = estimate_spot_distance(coords)
        if radius < min_distance:
            logger.warning(f"Radius: {radius} < min_distance: {min_distance}, reset min_distance.")
            radius = float(min_distance) * 1.01
        edge_index = radius_graph(
            x=coords,
            r=radius,
            batch=batch,
            loop=self_loop,
            num_workers=num_workers,
            max_num_neighbors=max_num_neighbors,
        )
    edge_index = to_undirected(edge_index) if undirected else edge_index
    assert edge_index.shape[1] > 0, "No edges contained in spatial graph."
    per_edge = edge_index.shape[1] / coords.shape[0]
    logger.success(
        f"Built {mode} graph with {coords.shape[0]} nodes and {edge_index.shape[1]} edges, {per_edge:.2f} edges/node."
    )
    if mode == "radius":
        num_neighbors = edge_index[0].bincount().float()
        logger.info(f"Mean number of neighbors: {num_neighbors.mean().item():.2f}")
        logger.info(f"Max number of neighbors: {num_neighbors.max().item()}")
        logger.info(f"Min number of neighbors: {num_neighbors.min().item()}")
        logger.info(f"Median number of neighbors: {num_neighbors.median().item()}")
        for percentile in [0.05, 0.25, 0.75, 0.95]:
            tile = num_neighbors.kthvalue(int(num_neighbors.numel() * percentile)).values.item()
            logger.info(f"{percentile * 100}th percentile of number of neighbors: {tile}")
    return edge_index


def knn_to_edge_index(knn_result: np.ndarray) -> Tensor:
    r"""
    Convert KNN result to edge index

    Parameters
    ----------
    knn_result
        KNN result in shape (num_nodes, k)
    """
    num_nodes = knn_result.shape[0]
    k = knn_result.shape[1]

    src_nodes = torch.repeat_interleave(torch.arange(num_nodes), k)
    dst_nodes = torch.tensor(knn_result.flatten())

    edge_index = torch.stack([src_nodes, dst_nodes], dim=0)

    return edge_index
