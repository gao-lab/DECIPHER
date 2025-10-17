import torch
from torch_geometric.data import Batch, Data

from decipher import CFG
from decipher.explain.gene.gene_selection import train_GAE
from decipher.explain.regress.mixin import train_regress


def test_GAE_mimic():
    CFG.gene_select.center_dim = 8
    CFG.gene_select.expr_dim = 100
    CFG.gene_select.work_dir = "./results/explain"
    CFG.gene_select.gae_epochs = 2

    N_CELL1 = 100
    graph1 = Data(
        x=torch.randn(N_CELL1, 8),
        edge_index=torch.randint(0, N_CELL1, (2, 500)),
        expr=torch.randn(N_CELL1, 100),
    )

    N_CELL2 = 120
    graph2 = Data(
        x=torch.randn(N_CELL2, 8),
        edge_index=torch.randint(0, N_CELL2, (2, 600)),
        expr=torch.randn(N_CELL2, 100),
    )
    graph_all = Batch.from_data_list([graph1, graph2])
    train_GAE(graph1, CFG.gene_select, save_dir="test_GAE")
    train_GAE(graph_all, CFG.gene_select, save_dir="test_GAE_batched")


def test_regress_mimic():
    CFG.regress.center_dim = 8
    CFG.regress.nbr_dim = 8
    CFG.regress.work_dir = "./results/explain"
    CFG.regress.trainer.epochs = 2

    x = torch.randn(100, 8)
    y = torch.randn(100, 8)
    train_regress(x, y, CFG.regress, save_dir="test_regress")


if __name__ == "__main__":
    test_GAE_mimic()
    test_regress_mimic()
