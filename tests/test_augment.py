import torch
from addict import Dict
from loguru import logger
from pytest import mark
from torch_geometric.data import Data

from spider.data.augment import OmicsSpatialAugment

config = Dict(
    {
        "dropout_gex": 0.1,
        "mask_hop": -1,
        "dropout_nbr_prob": -1,
    }
)
edge_index = torch.tensor([[0, 1, 2, 3, 4, 5], [1, 2, 0, 3, 5, 13]], dtype=torch.long)
batch = torch.tensor([0, 1, 0, 0, 1, 1], dtype=torch.long)
n_id = torch.range(0, 5, dtype=torch.long)
g_omics = Data(
    x=torch.rand(6, 10), batch=batch, edge_index=edge_index, batch_size=2, n_id=n_id
)


@mark.parametrize("config,edge_index,x", [(config, edge_index, g_omics)])
def test_OmicsSpatialAugment(config, edge_index, g):
    aug = OmicsSpatialAugment(config, edge_index)
    x1, x2, batch_mask = aug(g)


if __name__ == "__main__":
    test_OmicsSpatialAugment(config, edge_index, g_omics)
    logger.success("test passed")
