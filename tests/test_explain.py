import numpy as np
import pandas as pd
import scanpy as sc
import torch
from torch_geometric.data import Batch, Data

from decipher.explain.gene.gene_selection import train_GAE
from decipher.explain.gene.lr import get_lr_expr
from decipher.explain.regress.mixin import train_regress
from decipher.utils import GENESELECT_CFG, REGRESS_CFG


def test_GAE_mimic():
    GENESELECT_CFG.center_dim = 64
    GENESELECT_CFG.expr_dim = 300
    GENESELECT_CFG.work_dir = "./results"
    GENESELECT_CFG.loader.batch_size = 32

    N_CELL1 = 100
    graph1 = Data(
        x=torch.randn(N_CELL1, 64),
        edge_index=torch.randint(0, N_CELL1, (2, 1000)),
        expr=torch.randn(N_CELL1, 300),
    )

    N_CELL2 = 250
    graph2 = Data(
        x=torch.randn(N_CELL2, 64),
        edge_index=torch.randint(0, N_CELL2, (2, 1000)),
        expr=torch.randn(N_CELL2, 300),
    )
    graph = Batch.from_data_list([graph1, graph2])
    train_GAE(graph, GENESELECT_CFG, save_dir="test_GAE_batch")
    train_GAE(graph1, GENESELECT_CFG, save_dir="test_GAE")
    train_GAE(graph1, GENESELECT_CFG, batched=True, save_dir="test_GAE_mini_batch")


def test_regress_mimic():
    REGRESS_CFG.center_dim = 128
    REGRESS_CFG.nbr_dim = 32
    REGRESS_CFG.work_dir = "./results"

    x = torch.randn(100, 128)
    y = torch.randn(100, 32)
    train_regress(x, y, REGRESS_CFG, save_dir="test_regress")


def test_get_lr_expr_mimic():
    adata = sc.AnnData(
        X=np.random.randint(0, 10, (100, 50)),
        obsm={"spatial": np.random.randn(100, 2)},
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    lr = pd.read_csv("../resource/lr_data/cellchat.csv")
    get_lr_expr(adata, lr, radius=0.5)


if __name__ == "__main__":
    # test_GAE_mimic()
    test_regress_mimic()
    # test_get_lr_expr_mimic()
    # test_mine_mimic()
