import torch
from torch_geometric.data import Data

from decipher.data.augment import OmicsSpatialAugment, ScAugment


def test_OmicsSpatialAugment():
    # build a graph
    edge_index = torch.tensor([[0, 1, 2, 3, 4, 5], [1, 2, 0, 3, 5, 13]], dtype=torch.long)
    batch = torch.tensor([0, 1, 0, 0, 1, 1], dtype=torch.long)
    n_id = torch.arange(0, 5, dtype=torch.long)
    graph = Data(x=torch.rand(6, 10), batch=batch, edge_index=edge_index, batch_size=2, n_id=n_id)

    aug = OmicsSpatialAugment()
    x1, x2, batch_mask = aug(graph)
    assert x1.shape == x2.shape
    assert (x1 - x2).sum() != 0
    if batch_mask is not None:  # TODO: add sc_batch in the test
        assert x1.shape[0] == batch_mask.shape[0]


def test_ScAugment():
    aug = ScAugment()

    x = torch.rand(6, 10)
    x1, x2 = aug(x)

    assert x1.shape == x2.shape
    assert (x1 - x2).sum() != 0


if __name__ == "__main__":
    test_OmicsSpatialAugment()
    test_ScAugment()
