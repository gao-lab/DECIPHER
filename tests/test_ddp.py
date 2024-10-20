import scanpy as sc
from addict import Dict

from decipher import DECIPHER


def test_run_sc_emb(work_dir, config):
    model = DECIPHER(work_dir, config)
    model = model.fit_sc()


def test_run_spatial_emb(work_dir, config):
    model = DECIPHER(work_dir, config)
    model = model.fit_spatial()


def test_run_spatial_emb_infer(work_dir, config):
    model = DECIPHER(work_dir, config)
    model.inference_spaital()


def test_ddp(adata, work_dir, config):
    model = DECIPHER(work_dir, config)
    model.register_data(adata)
    model.fit_ddp(gpus=2)


if __name__ == "__main__":
    args = Dict()
    work_dir = "./results/ddp"
    adata = sc.datasets.visium_sge("V1_Human_Lymph_Node")
    test_ddp(adata, work_dir, args)
