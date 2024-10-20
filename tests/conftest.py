r"""
Shared data
"""
import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from pytest import fixture
from torch_geometric.data import Data
from torch_geometric.utils import remove_self_loops


@fixture
def adata() -> AnnData:
    r"""
    Generate omics data
    """
    N_CELL = 1000
    N_GENES = 200
    return AnnData(
        X=np.random.randn(N_CELL, N_GENES),
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(100)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(100)]),
        obsm={"spatial": np.random.rand(N_CELL, 2)},
    )


@fixture
def img_graph() -> Data:
    r"""
    Generate image graph
    """
    N_CELL = 500
    IMAGE_DIM = 16
    N_NEIGHBOR_PER_CELL = 10
    edge_index = torch.randint(0, N_CELL, (2, N_CELL * N_NEIGHBOR_PER_CELL))
    return Data(
        x=torch.randint(0, 255, (N_CELL, 3, IMAGE_DIM, IMAGE_DIM)),
        edge_index=remove_self_loops(edge_index),
    )


@fixture
def wsi() -> np.ndarray:
    r"""
    Generate WSI image
    """
    WEIGHT = 1000
    HEIGHT = 800
    return np.random.randint(0, 255, (HEIGHT, WEIGHT, 3))
