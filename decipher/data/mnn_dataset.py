r"""
Mutual nearest neighbor (MNN) torch dataset
"""
from copy import deepcopy

import numpy as np
import torch
from omegaconf import OmegaConf
from pytorch_lightning import LightningDataModule
from pytorch_lightning.utilities import CombinedLoader
from torch import Tensor
from torch.utils.data import DataLoader, Dataset
from torch_geometric.data import Data
from torch_geometric.data.lightning import LightningNodeData

from .. import CFG


class LightningSpatialMNNData(LightningNodeData):
    r"""
    Wrapper of `LightningNodeData` to support extra MNN pairs.

    Parameters
    ----------
    graph:
        spatial graph
    loader_config:
        loader configuration
    mnn_dataset:
        MNN dataset
    **kwargs:
        extra keyword arguments for loader
    """

    def __init__(
        self,
        graph: Data,
        loader_config: OmegaConf,
        mnn_dataset: DataLoader = None,
        **kwargs,
    ):
        super().__init__(graph, **loader_config, **kwargs)
        self.loader_config = loader_config
        self.mnn_dataset = mnn_dataset

    def train_dataloader(self):
        graph_train_loader = super().train_dataloader()
        if self.mnn_dataset is None:
            return graph_train_loader
        else:
            mnn_loader = DataLoader(self.mnn_dataset, **self.loader_config)
            loaders = {"graph": graph_train_loader, "mnn": mnn_loader}
            return CombinedLoader(loaders, mode="max_size_cycle")


class LightningScMNNData(LightningDataModule):
    r"""
    Single cell lighting datamodule support MNN

    Parameters
    ----------
    loader_config
        DataLoader configuration
    train_dataset
        training dataset
    mnn_dataset
        MNN dataset for batch correction
    """

    def __init__(
        self,
        loader_config: OmegaConf,
        train_dataset: Dataset,
        mnn_dataset: Dataset = None,
    ) -> None:
        super().__init__()
        self.loader_config = loader_config
        self.train_dataset = train_dataset
        self.mnn_dataset = mnn_dataset

    def train_dataloader(self):
        train_dataset = DataLoader(self.train_dataset, **self.loader_config)
        if self.mnn_dataset is None:
            return train_dataset
        else:
            # with batch
            mnn_loader = DataLoader(self.mnn_dataset, **self.loader_config)
            loaders = {"x": train_dataset, "mnn": mnn_loader}
            return CombinedLoader(loaders, mode="max_size_cycle")

    def test_dataloader(self):
        val_cfg = deepcopy(self.loader_config)
        val_cfg.update({"batch_size": 1024, "shuffle": False, "drop_last": False})
        return DataLoader(self.train_dataset, **val_cfg)


class MNNDataset(Dataset):
    r"""
    Mutual nearest neighbors (MNNs) dataset to provide positive samples

    Parameters
    ----------
    X:
        scaled single cell expression matrix (cell x gene)
    valid_cellidx:
        valid cell index (cells included in MNN pairs)
    mnn_dict:
        MNN pairs dict
    """

    def __init__(self, x: np.ndarray, valid_cellidx: np.ndarray, mnn_dict: dict) -> None:
        super().__init__()
        self.x = torch.from_numpy(x)
        self.valid_cellidx = valid_cellidx
        self.mnn_dict = mnn_dict

    def __len__(self) -> int:
        return len(self.valid_cellidx)

    def __getitem__(self, idx: int) -> tuple[Tensor]:
        idx = self.valid_cellidx[idx]
        x1 = self.x[idx]
        x2 = self.x[np.random.choice(self.mnn_dict[idx])]
        return x1, x2


def get_graph_datamodule(graph: Data, mnn_dataset=None) -> LightningSpatialMNNData:
    r"""
    Build `LightningNodeMNNData` datamodule

    Parameters
    ----------
    graph
        spatial graph
    config
        model config
    mnn_dataset
        mnn dataset
    """
    # del config.loader["shuffle"]
    mask = torch.ones(graph.num_nodes, dtype=torch.bool)
    datamodule = LightningSpatialMNNData(
        graph,
        CFG.sp_model.loader,
        mnn_dataset,
        num_neighbors=CFG.sp_model.num_neighbors,
        input_train_nodes=mask,
        input_val_nodes=mask,
        input_test_nodes=mask,
        subgraph_type="bidirectional",
        disjoint=True,
        eval_loader_kwargs={"drop_last": False},
    )
    return datamodule
