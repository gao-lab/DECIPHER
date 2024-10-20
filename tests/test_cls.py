import numpy as np
import scanpy as sc

from spider import Spider
from spider.utils import CFG, GENESELECT_CFG, REGRESS_CFG, GetRunTime, global_seed


@GetRunTime
def test_spider_omics_single(h5ad_path: str = None):
    global_seed(0)
    CFG.omics.model.augment.dropout_gex = 0.6
    CFG.omics.model.epochs = 3
    CFG.omics.model.plot = True
    CFG.omics.loader.batch_size = 128
    CFG.omics.pretrain.force = True
    CFG.omics.pretrain.epochs = 1
    GENESELECT_CFG.gae_epochs = 3
    REGRESS_CFG.fit.epochs = 3

    adata = (
        sc.datasets.visium_sge("V1_Human_Lymph_Node") if h5ad_path is None else sc.read(h5ad_path)
    )
    adata.layers["counts"] = adata.X.copy()
    # random choose cell type from ['a', 'b', 'c']
    adata.obs["cell_type"] = np.random.choice(["a", "b", "c"], adata.n_obs).tolist()

    model = Spider(work_dir="./results/spider_omics_single", user_cfg=CFG)
    model.register_data(adata)
    model.fit_omics()
    model.visualize()
    adata.X = adata.layers["counts"].copy()
    model.train_gene_select(adata, cell_type="cell_type")
    model.train_regress_explain(adata, cell_type="cell_type", reverse_regress=True)


@GetRunTime
def test_spider_omics_multi():
    CFG.omics.model.epochs = 3
    work_dir = "./results/spider_omics_multi"
    # adata1 = sc.datasets.visium_sge("V1_Breast_Cancer_Block_A_Section_1")
    # adata2 = sc.datasets.visium_sge("V1_Breast_Cancer_Block_A_Section_2")
    model = Spider(work_dir=work_dir, recover=True)
    # model.register_data([adata1, adata2], group_list=["Section1", "Section2"])
    model.fit_omics()
    model.visualize()


if __name__ == "__main__":
    test_spider_omics_single()
    # test_spider_omics_multi()
