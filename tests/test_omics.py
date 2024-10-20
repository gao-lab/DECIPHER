import scanpy as sc

from spider.emb import spatial_emb
from spider.graphic.build import build_graph
from spider.utils import CFG, sync_config


def test_omics_contrastive(adata):
    CFG.work_dir = "./results/omics_pbmc_mimic"
    CFG.omics.fit.device_num = 1
    CFG.omics.model.epochs = 2
    CFG.omics.loader.batch_size = 128
    sync_config(CFG)
    spatial_edge = build_graph(adata.obsm["spatial"], **CFG.omics.spatial_graph)
    adata = spatial_emb(adata, spatial_edge, CFG.omics)
    print(adata)


if __name__ == "__main__":
    test_omics_contrastive(
        sc.read_h5ad("../data/processed/human/pbmc/mimic/pbmc3k/adata.h5ad")
    )
