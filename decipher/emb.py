r"""
Core embedding models
"""
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import torch
from loguru import logger
from rui_utils.torch.trainer import init_trainer
from rui_utils.utils import l2norm
from torch import Tensor
from torch.utils.data import Dataset, TensorDataset
from torch_geometric.data import Data

from . import CFG
from .data.mnn_dataset import LightningScMNNData, get_graph_datamodule
from .nn import ScSimCLR, SpatialSimCLR


class ScModel:
    r"""
    Single cell embedding model

    Args:
        x: gene expression matrix
        mnn_dataset: MNN dataset
    """

    def __init__(self, x, mnn_dataset):
        train_dataset = TensorDataset(torch.from_numpy(x))
        self.datamodule = LightningScMNNData(CFG.sc_model.loader, train_dataset, mnn_dataset)
        self.model = ScSimCLR(CFG.sc_model.model)

        self.trainer = init_trainer(CFG.sc_model.trainer)

    def train(self) -> None:
        self.trainer.fit(self.model, self.datamodule)

    def infer(self, norm: bool = True) -> np.ndarray:
        self.trainer.test(self.model, self.datamodule)
        emb = self.model.gather_output()
        if norm:
            emb = l2norm(emb.astype(np.float32))
        return emb


class SpatialModel:
    r"""
    Spatial embedding model

    Args:
        x: gene expression matrix
        spatial_edge: spatial graph edge index
        mnn_dataset: MNN dataset
        batch: batch labels
        pretrained_model: pre-trained single cell model
    """

    def __init__(
        self,
        x: np.ndarray,
        spatial_edge: Tensor,
        mnn_dataset: Dataset = None,
        batch: np.ndarray = None,
        pretrained_model: ScSimCLR = None,
    ):
        graph = Data(x=torch.from_numpy(x), edge_index=spatial_edge)
        if batch is not None:
            graph.sc_batch = torch.tensor(batch, dtype=int)
        self.datamodule = get_graph_datamodule(graph, mnn_dataset)
        self.model = SpatialSimCLR(CFG.sp_model.model)

        if pretrained_model is not None:
            self.model.center_encoder = deepcopy(pretrained_model.center_encoder)
        self.trainer = init_trainer(CFG.sp_model.trainer)

    def train(self) -> None:
        self.trainer.fit(self.model, self.datamodule)

    def infer(self, norm: bool = True) -> tuple[np.ndarray, np.ndarray]:
        self.trainer.test(self.model, self.datamodule)
        center_emb, nbr_emb = self.model.gather_output()
        if norm:
            center_emb = l2norm(center_emb.astype(np.float32))
            nbr_emb = l2norm(nbr_emb.astype(np.float32))
        return center_emb, nbr_emb


def load_model(dir: str, model_cfg, spatial: bool = False) -> ScSimCLR | SpatialSimCLR:
    r"""
    Load omics encoder model

    Args:
        dir: model directory
        model_cfg: model config
        spatial: whether load spatial model
    """
    model_cls = SpatialSimCLR if spatial else ScSimCLR
    # sort by modification time
    model_path = sorted(Path(dir).glob("*.ckpt"), key=os.path.getmtime)[-1]
    logger.info(f"Loading model from {model_path}")
    model = model_cls.load_from_checkpoint(model_path, config=model_cfg, weights_only=False)
    logger.success(f"Pre-trained model loaded from {model_path}.")
    return model
