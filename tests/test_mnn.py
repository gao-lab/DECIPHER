import numpy as np

from spider.data.mnn_dataset import MNNDataset
from spider.utils import CFG


def test_MNNDataset(*args, **kwargs):
    MNNDataset(*args, **kwargs)


if __name__ == "__main__":
    N_CELL = 1000_000
    N_DIM = 400
    x = np.random.randn(3 * N_CELL, N_DIM)
    batch = np.array([0] * N_CELL + [1] * N_CELL + [2] * N_CELL)
    test_MNNDataset(x, batch, CFG.omics.mnn)
