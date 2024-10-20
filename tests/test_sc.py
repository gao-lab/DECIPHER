import scanpy as sc

from spider.emb import sc_emb
from spider.utils import CFG, sync_config


def test_sc_plain(adata):
    CFG.work_dir = "./results/sc_pbmc_mimic_plain"
    sync_config(CFG)
    sc_emb(adata, CFG.omics, preprocess=False)


def test_sc_mnn(adata):
    CFG.work_dir = "./results/sc_pbmc_mimic_mnn"
    sync_config(CFG)
    sc_emb(adata, CFG.omics, "batch", preprocess=False)


if __name__ == "__main__":
    adata = sc.read_h5ad("../data/processed/human/pbmc/mimic/pbmc3k/adata.h5ad")
    test_sc_plain(adata)
    test_sc_mnn(adata)
