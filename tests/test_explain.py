import torch
from torch_geometric.data import Batch, Data

from decipher.explain.gene.gene_selection import train_GAE
from decipher.explain.regress.mixin import train_regress
from decipher.utils import GENESELECT_CFG, REGRESS_CFG


def test_GAE_mimic():
    GENESELECT_CFG.center_dim = 32
    GENESELECT_CFG.expr_dim = 300
    GENESELECT_CFG.work_dir = "./results/explain"
    GENESELECT_CFG.gae_epochs = 10
    GENESELECT_CFG.fit.epochs = 10

    N_CELL1 = 200
    graph1 = Data(
        x=torch.randn(N_CELL1, 32),
        edge_index=torch.randint(0, N_CELL1, (2, 1000)),
        expr=torch.randn(N_CELL1, 300),
    )

    N_CELL2 = 300
    graph2 = Data(
        x=torch.randn(N_CELL2, 32),
        edge_index=torch.randint(0, N_CELL2, (2, 1000)),
        expr=torch.randn(N_CELL2, 300),
    )
    graph_all = Batch.from_data_list([graph1, graph2])
    train_GAE(graph1, GENESELECT_CFG, save_dir="test_GAE")
    train_GAE(graph_all, GENESELECT_CFG, save_dir="test_GAE_batched")


def test_regress_mimic():
    REGRESS_CFG.center_dim = 32
    REGRESS_CFG.nbr_dim = 32
    REGRESS_CFG.work_dir = "./results/explain"
    REGRESS_CFG.fit.epochs = 5

    x = torch.randn(100, 32)
    y = torch.randn(100, 32)
    train_regress(x, y, REGRESS_CFG, save_dir="test_regress")


# NOTE: the test of LR is included in test_cls.py


if __name__ == "__main__":
    test_GAE_mimic()
    test_regress_mimic()
