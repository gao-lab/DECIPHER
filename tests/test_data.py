import numpy as np
import torch
from addict import Dict
from torch_geometric.data import Data

from decipher.data.augment import OmicsSpatialAugment, ScAugment
from decipher.data.mnn_dataset import findMNN


def test_OmicsSpatialAugment():
    # build a graph
    edge_index = torch.tensor([[0, 1, 2, 3, 4, 5], [1, 2, 0, 3, 5, 13]], dtype=torch.long)
    batch = torch.tensor([0, 1, 0, 0, 1, 1], dtype=torch.long)
    n_id = torch.arange(0, 5, dtype=torch.long)
    graph = Data(x=torch.rand(6, 10), batch=batch, edge_index=edge_index, batch_size=2, n_id=n_id)

    config = Dict(dropout_gex=0.5, dropout_nbr_prob=-1, mask_hop=-1, max_neighbor=-1)
    aug = OmicsSpatialAugment(config)
    x1, x2, batch_mask = aug(graph)
    assert x1.shape == x2.shape
    assert (x1 - x2).sum() != 0
    if batch_mask is not None:  # TODO: add sc_batch in the test
        assert x1.shape[0] == batch_mask.shape[0]


def test_ScAugment():
    config = Dict(dropout_gex=0.1)
    aug = ScAugment(config)

    x = torch.rand(6, 10)
    x1, x2 = aug(x)

    assert x1.shape == x2.shape
    assert (x1 - x2).sum() != 0


def test_findMNN():
    x = np.random.randn(100, 50)
    y = np.random.randn(100, 50)
    name_x = np.array([f"cell_{i}" for i in range(100)])
    name_y = np.array([f"cell_{i}" for i in range(100)])
    df = findMNN(x, y, name_x, name_y, 5, k_components=10)
    print(df.head())
    print(df.shape)


if __name__ == "__main__":
    test_OmicsSpatialAugment()
    test_ScAugment()
    test_findMNN()
