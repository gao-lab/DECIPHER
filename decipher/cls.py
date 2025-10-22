r"""
Main class of `decipher` package
"""
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import torch
from anndata import AnnData
from loguru import logger
from omegaconf import OmegaConf
from pytorch_lightning import seed_everything

from . import CFG
from .data.mnn_dataset import MNNDataset
from .data.process import omics_data_process
from .emb import ScModel, SpatialModel, load_model
from .explain.gene.mixin import GeneSelectMixin
from .explain.regress.mixin import RegressMixin
from .graphic.build import build_graph
from .graphic.knn import cal_mnn


class DECIPHER(RegressMixin, GeneSelectMixin):
    r"""
    Base class of model definition and training

    Parameters
    ----------
    work_dir
        working directory
    recover
        if recover from a previous run
    """

    def __init__(
        self, work_dir: str = "DECIPHER", recover: bool = False, overwrite: bool = False
    ) -> None:
        self.work_dir = Path(work_dir)
        now = datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
        if self.work_dir.exists():
            if recover:
                self.load(work_dir)
            else:
                if not overwrite:
                    raise FileExistsError
        else:
            self.work_dir.mkdir(parents=True, exist_ok=True)
        CFG.sc_model.trainer.model_dir = str(self.work_dir / "sc_model")
        CFG.sp_model.trainer.model_dir = str(self.work_dir / "sp_model")
        logger.add(self.work_dir / "logs" / f"run_{now}.log")
        seed_everything(CFG.seed)

    def register_data(
        self,
        adata: list[AnnData] | AnnData,
        split_by: str = None,
        preprocess: bool = True,
        edge_index: np.ndarray = None,
        save_data: bool = True,
    ) -> None:
        r"""
        Register spatial omics data

        Parameters
        ----------
        adata:
            `AnnData` or a list of `AnnData`, each element is a batch
        split_by:
            Only works when `adata` is an `AnnData` object, indicates batch column in `.obs`
        preprocess:
            if `False`, will skip the preprocess steps, only for advanced users
        edge_index
            use self-defined spatial neighbor edges (PyG format), only for advanced users
        """
        # process spatial omics data
        self.x, coords, self.batch_idx = omics_data_process(adata, split_by, preprocess)
        del adata
        if save_data:
            np.save(self.work_dir / "x.npy", self.x)

        # mnn correction
        if self.batch_idx is not None:
            self.mnn_dict, self.mnn_cell_idx = cal_mnn(self.x, self.batch_idx, CFG.mnn)
            np.save(self.work_dir / "mnn_cell_idx.npy", self.mnn_cell_idx)
            torch.save(self.mnn_dict, self.work_dir / "mnn_dict.pt", pickle_protocol=4)
        torch.save(self.batch_idx, self.work_dir / "batch.pt")

        # build spatial graph
        if edge_index is None:
            self.edge_index = build_graph(coords, self.batch_idx, **CFG.graph)
        else:
            logger.info("Use self-defined edge index.")
            assert edge_index.max() < self.x.shape[0], "Edge index out of range."
            self.edge_index = edge_index
        np.save(self.work_dir / "edge_index.npy", self.edge_index.numpy())

        # Save hyperparameters
        OmegaConf.save(CFG, self.work_dir / "hyperparams.yaml")

    def fit_sc(self) -> None:
        r"""
        Fit on single cell data
        """
        CFG.sc_model.model.gex_dims[0] = self.x.shape[1]
        mnn_dataset = (
            MNNDataset(self.x, self.mnn_cell_idx, self.mnn_dict)
            if self.batch_idx is not None
            else None
        )
        self.sc_model = ScModel(self.x, mnn_dataset)
        self.sc_model.train()

    def fit_spatial(self) -> None:
        r"""
        Fit on spatial cell data
        """
        CFG.sp_model.model.gex_dims[0] = self.x.shape[1]
        mnn_dataset = (
            MNNDataset(self.x, self.mnn_cell_idx, self.mnn_dict)
            if self.batch_idx is not None
            else None
        )
        # train model
        sc_encoder = load_model(self.work_dir / "sc_model", CFG.sc_model.model, spatial=False)
        self.sp_model = SpatialModel(
            self.x, self.edge_index, mnn_dataset, self.batch_idx, sc_encoder
        )
        self.sp_model.train()

    def fit_omics(self) -> None:
        r"""
        Fit model on spatial omics data
        """
        self.fit_sc()
        self.fit_spatial()
        self.center_emb, self.nbr_emb = self.sp_model.infer()
        # save embeddings
        np.save(self.work_dir / "center_emb.npy", self.center_emb)
        np.save(self.work_dir / "nbr_emb.npy", self.nbr_emb)
        logger.info(f"Results saved to {self.work_dir}")

    def load(self, from_dir: str) -> None:
        r"""
        Load saved results, should run after `register_data`

        Args:
            from_dir: directory to load from, default is `self.work_dir`
        """
        logger.info(f"Loading from {from_dir}")
        assert Path(from_dir).is_dir(), f"Directory {from_dir} not exists."
        self.work_dir = Path(from_dir)
        CFG = OmegaConf.load(self.work_dir / "hyperparams.yaml")
        # load data
        self.x = np.load(self.work_dir / "x.npy").astype(np.float32)
        self.edge_index = torch.from_numpy(np.load(self.work_dir / "edge_index.npy"))
        self.batch_idx = torch.load(self.work_dir / "batch.pt", weights_only=False)
        if self.batch_idx is not None and not CFG.pp.ignore_batch:
            self.mnn_cell_idx = np.load(self.work_dir / "mnn_cell_idx.npy")
            self.mnn_dict = torch.load(self.work_dir / "mnn_dict.pt", weights_only=False)
        for var in ["center_emb", "nbr_emb"]:
            if (self.work_dir / f"{var}.npy").exists():
                setattr(self, var, np.load(self.work_dir / f"{var}.npy").astype(np.float32))

    def fit_ddp(self, gpus: int = -1) -> None:
        r"""
        DDP training

        Args:
            gpus: number of gpus to use, -1 for all available gpus
        """
        if self.x.shape[0] < 500_000:
            logger.warning("Using DDP with < 500k cells is not recommended.")
        max_gpus = torch.cuda.device_count()
        assert max_gpus > 1, "DDP requires at least 2 GPUs."
        gpus = min(gpus, max_gpus) if gpus > 0 else max_gpus
        logger.info(f"Using {gpus} GPUs for DDP training.")

        # DDP fit omics
        subprocess.run(
            [f"decipher_ddp_sc --work_dir {str(self.work_dir)} --gpus {gpus}"],
            shell=True,
            check=True,
        )
        # DDP fit spatial
        subprocess.run(
            [f"decipher_ddp_spatial --work_dir {str(self.work_dir)} --gpus {gpus}"],
            shell=True,
            check=True,
        )
        logger.success("DDP training finished.")
