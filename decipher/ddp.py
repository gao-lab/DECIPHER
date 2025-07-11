r"""
DDP mixin for multi-GPU training
"""
import subprocess

import numpy as np
import torch
from click import command, option
from loguru import logger
from torch_geometric.data import Data

from .data.mnn_dataset import MNNDataset
from .emb import (
    _update_config,
    get_graph_datamodule,
    load_sc_model,
    load_spatial_model,
    sc_emb,
    spatial_emb,
)
from .nn.trainer import inference


class DDPMixin:
    r"""
    DDP mixin
    """

    def fit_ddp(self, gpus: int = -1, ddp_pretrain: bool = False) -> None:
        r"""
        DDP training wrapper

        Parameters
        ----------
        gpus
            gpu number, -1 means using all available gpus
        ddp_pretain
            if using DDP strategy for pretrain, default is False
        """
        if self.x.shape[0] < 500_000:
            logger.warning("Using DDP with < 500k cells is not recommended.")

        max_gpus = torch.cuda.device_count()
        assert max_gpus > 1, "DDP requires at least 2 GPUs."
        gpus = min(gpus, max_gpus) if gpus > 0 else max_gpus

        if ddp_pretrain:
            self.update_config({"device_num": gpus})
            # DDP fit omics
            subprocess.run(
                [f"decipher_ddp_sc --work_dir {str(self.work_dir)}"], shell=True, check=True
            )
        else:
            self.fit_sc()
            self.update_config({"device_num": gpus})
        # DDP fit spatial
        subprocess.run(
            [f"decipher_ddp_spatial --work_dir {str(self.work_dir)}"], shell=True, check=True
        )
        self.inference_spaital()
        logger.success("DDP training finished.")

    def fit_sc(self) -> None:
        r"""
        Only fit on single cell data (for DDP)
        """
        # mnn dataset
        mnn_dataset = None
        if self.batch is not None and not self.cfg.omics.ignore_batch:
            mnn_dataset = MNNDataset(self.x, self.valid_cellidx, self.mnn_dict)
            logger.info(f"Using MNN with {len(np.unique(self.batch))} batches.")
        # train model
        _, self.center_emb = sc_emb(self.x, self.cfg.omics, mnn_dataset, self.meta, self.batch)

    def fit_spatial(self) -> None:
        r"""
        Only fit on spatial cell data (for DDP)
        """
        # mnn dataset
        mnn_dataset = None
        if self.batch is not None and not self.cfg.omics.ignore_batch:
            mnn_dataset = MNNDataset(self.x, self.valid_cellidx, self.mnn_dict)
            logger.info(f"Using MNN with {len(np.unique(self.batch))} batches.")
        # train model
        mnn_flag = True if mnn_dataset is not None else False
        sc_model = load_sc_model(self.cfg.omics, mnn_flag, self.meta)
        spatial_emb(
            self.x,
            self.edge_index,
            self.cfg.omics,
            mnn_dataset,
            self.meta,
            sc_model,
            self.batch,
            DDP=True,
        )

    def inference_spaital(self) -> None:
        r"""
        Inference spatial embeddings
        """
        config = self.cfg.omics
        config = _update_config(self.x.shape[0], self.x.shape[1], config)
        # mnn dataset
        mnn_dataset = None
        mnn_flag = False
        if self.batch is not None and self.cfg.mnn:
            mnn_dataset = MNNDataset(self.x, self.valid_cellidx, self.mnn_dict)
            logger.info(f"Using MNN with {len(np.unique(self.batch))} batches.")
            mnn_flag = True
        # data
        graph = Data(x=torch.Tensor(self.x), edge_index=self.edge_index)
        datamodule = get_graph_datamodule(graph, config, mnn_dataset)
        # trained model
        model = load_spatial_model(self.cfg.omics, mnn_flag)
        # inference
        inference(model, datamodule, config=config.model)
        center_emb, self.nbr_emb = model.gather_output()
        if self.center_emb is None:
            self.center_emb = center_emb
        np.save(self.work_dir / "center_emb.npy", self.center_emb)
        np.save(self.work_dir / "nbr_emb.npy", self.nbr_emb)


@command()
@option("--work_dir", help="work directory", required=True, type=str)
def decipher_ddp_spatial(work_dir):
    from decipher import DECIPHER

    model = DECIPHER(work_dir, recover=True)
    model.fit_spatial()


@command()
@option("--work_dir", help="work directory", required=True, type=str)
def decipher_ddp_sc(work_dir):
    from decipher import DECIPHER

    model = DECIPHER(work_dir, recover=True)
    model.fit_sc()
