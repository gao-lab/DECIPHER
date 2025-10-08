import scanpy as sc
import torch
from loguru import logger

from decipher import CFG, DECIPHER


def test_ddp():
    r"""
    DDP can not be run in pytest pipeline, only for manual test
    """
    if torch.cuda.device_count() < 2:
        logger.error("Only have < 2 gpus, Skip DDP test")
        return

    work_dir = "./results/decipher_ddp"
    adata = sc.datasets.visium_sge("V1_Breast_Cancer_Block_A_Section_1")
    adata = adata[:400, :500].copy()

    CFG.sp_model.model.epochs = 2
    CFG.sc_model.model.epochs = 2

    model = DECIPHER(work_dir, overwrite=True)
    model.register_data(adata)
    model.fit_ddp(gpus=2)

    # test model recovery
    model_recover = DECIPHER(work_dir, recover=True)
    model_recover.fit_ddp(gpus=2)


if __name__ == "__main__":
    test_ddp()
