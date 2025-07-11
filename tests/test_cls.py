import numpy as np
import scanpy as sc

from decipher import DECIPHER
from decipher.utils import CFG, GENESELECT_CFG, REGRESS_CFG, GetRunTime, global_seed


def get_adata() -> sc.AnnData:
    adata = sc.datasets.visium_sge("V1_Breast_Cancer_Block_A_Section_1")
    # sample first 200 cells and 500 genes
    adata = adata[:200, :500]
    adata.layers["counts"] = adata.X.copy()
    adata.obs["cell_type"] = np.random.choice(["a", "b"], adata.n_obs).tolist()
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata.copy()


def get_adata_2() -> sc.AnnData:
    adata = sc.datasets.visium_sge("V1_Breast_Cancer_Block_A_Section_2")
    # sample first 210 cells and 500 genes
    adata = adata[:210, :500]
    adata.layers["counts"] = adata.X.copy()
    adata.obs["cell_type"] = np.random.choice(["a", "b"], adata.n_obs).tolist()
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata.copy()


@GetRunTime
def test_decipher_single_slice():
    adata = get_adata()

    global_seed(0)
    CFG.omics.model.augment.dropout_gex = 0.6
    CFG.omics.model.epochs = 2
    CFG.omics.loader.batch_size = 64
    CFG.omics.pretrain.force = True
    CFG.omics.pretrain.epochs = 2
    GENESELECT_CFG.gae_epochs = 2
    REGRESS_CFG.fit.epochs = 2
    work_dir = "./results/decipher_single_slice"

    model = DECIPHER(work_dir=work_dir, user_cfg=CFG)
    model.register_data(adata)

    model.fit_omics()

    adata.X = adata.layers["counts"].copy()
    model.train_gene_select(adata, cell_type="cell_type", min_cells=1, n_jobs=1)
    model.train_regress_explain(adata, cell_type="cell_type", reverse_regress=True, n_jobs=1)

    # test the recovery of the model
    model_recover = DECIPHER(work_dir=work_dir, recover=True)
    CFG.omics.pretrain.force = False
    model_recover.fit_omics()


@GetRunTime
def test_decipher_multi_slices():
    adata_1, adata_2 = get_adata(), get_adata_2()
    adata = adata_1.concatenate(adata_2)

    global_seed(0)
    CFG.omics.model.augment.dropout_gex = 0.6
    CFG.omics.model.epochs = 2
    CFG.omics.loader.batch_size = 64
    CFG.omics.pretrain.force = True
    CFG.omics.pretrain.epochs = 2
    GENESELECT_CFG.gae_epochs = 2
    REGRESS_CFG.fit.epochs = 2
    work_dir = "./results/decipher_multi_slices"

    model = DECIPHER(work_dir=work_dir, user_cfg=CFG)
    model.register_data([adata_1, adata_2], group_list=["batch_1", "batch_2"])

    model.fit_omics()

    adata.X = adata.layers["counts"].copy()
    model.train_gene_select(adata, cell_type="cell_type", min_cells=1, n_jobs=1)
    model.train_regress_explain(adata, cell_type="cell_type", reverse_regress=True, n_jobs=1)

    # test the recovery of the model
    model_recover = DECIPHER(work_dir=work_dir, recover=True)
    CFG.omics.pretrain.force = False
    model_recover.fit_omics()


if __name__ == "__main__":
    test_decipher_single_slice()
    test_decipher_multi_slices()
