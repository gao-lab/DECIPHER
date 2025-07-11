r"""
Mutual nearest neighbor (MNN) for batch correction.
"""

from collections import defaultdict
from copy import deepcopy

import numpy as np
import pandas as pd
import torch
from addict import Dict
from loguru import logger
from numba import njit, prange
from pytorch_lightning import LightningDataModule
from pytorch_lightning.utilities import CombinedLoader
from sklearn.utils.extmath import randomized_svd
from torch import Tensor
from torch.utils.data import DataLoader, Dataset
from torch_geometric.data import Data
from torch_geometric.data.lightning import LightningNodeData

from ..graphic.knn import knn
from ..utils import l2norm


class MNNMixin:
    r"""
    MNN mixin
    """

    def initMNN(self):
        r"""
        Init nearest neighbors pairs
        """
        self.mnn_flag = True
        mnn_config = self.cfg.omics.mnn
        n_cell = self.x.shape[0]
        cell_names = np.array([f"batch_{b}_{idx}" for idx, b in enumerate(self.batch)])
        name2idx = dict(zip(cell_names, range(n_cell)))

        # get MNN pairs
        n_pcs, k_mnn = mnn_config.k_components, mnn_config.k_anchor
        logger.info(f"Computing MNN pairs with k = {k_mnn}")
        batchs, cells_per_batch = np.unique(self.batch, return_counts=True)
        n_batch = len(batchs)

        all_mnn_pairs = []
        # TODO: parallel with multi gpus
        if mnn_config.ref_based:
            # ref batch is batch with most cells
            logger.debug("Use reference-based MNN")
            ref_batch = batchs[np.argmax(cells_per_batch)]
            ref_batch_idx = self.batch == ref_batch
            X_i, name_i = self.x[ref_batch_idx], cell_names[ref_batch_idx]
            for batch_j in batchs:
                if batch_j == ref_batch:
                    continue
                batch_j_idx = self.batch == batch_j
                X_j, name_j = self.x[batch_j_idx], cell_names[batch_j_idx]
                mnn_pairs = findMNN(X_i, X_j, name_i, name_j, k_mnn, n_pcs)
                all_mnn_pairs.append(mnn_pairs)
        else:
            for i, batch_i in enumerate(batchs):
                batch_i_idx = self.batch == batch_i
                X_i, name_i = self.x[batch_i_idx], cell_names[batch_i_idx]
                for j in range(i + 1, n_batch):
                    batch_j_idx = self.batch == batchs[j]
                    X_j, name_j = self.x[batch_j_idx], cell_names[batch_j_idx]
                    mnn_pairs = findMNN(X_i, X_j, name_i, name_j, k_mnn, n_pcs)
                    all_mnn_pairs.append(mnn_pairs)
        pairs = pd.concat(all_mnn_pairs)
        # convert to global cell index
        pairs["cell1"] = pairs["cellname1"].apply(lambda x: name2idx[x])
        pairs["cell2"] = pairs["cellname2"].apply(lambda x: name2idx[x])
        pairs = pairs[["cell1", "cell2"]].values
        self.pairs = np.unique(pairs, axis=0)

        # get MNN dict
        self.mnn_dict = defaultdict(list)
        for r, c in self.pairs:
            self.mnn_dict[r].append(c)
            self.mnn_dict[c].append(r)

        # exclude cells without MNN
        exclude_fn = True
        if exclude_fn:
            self.valid_cellidx = np.unique(self.pairs.ravel())
        else:
            self.valid_cellidx = np.arange(n_cell)


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
        loader_config: Dict,
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
            combined_loader = CombinedLoader(loaders, mode="max_size_cycle")
            return combined_loader


class LightningScMNNData(LightningDataModule):
    r"""
    Single cell lighting datamodule support MNN

    Parameters
    ----------
    loader_config
        DataLoader configuration
    train_dataset
        training dataset
    val_dataset
        validation dataset
    mnn_dataset
        MNN dataset for batch correction
    """

    def __init__(
        self,
        loader_config: Dict,
        train_dataset: Dataset,
        val_dataset: Dataset,
        mnn_dataset: Dataset = None,
    ) -> None:
        super().__init__()
        self.loader_config = loader_config
        self.train_dataset = train_dataset
        self.val_dataset = val_dataset
        self.mnn_dataset = mnn_dataset

    def train_dataloader(self):
        train_dataset = DataLoader(self.train_dataset, **self.loader_config)
        # without batch
        if self.mnn_dataset is None:
            return train_dataset
        else:
            # with batch
            mnn_loader = DataLoader(self.mnn_dataset, **self.loader_config)
            loaders = {"x": train_dataset, "mnn": mnn_loader}
            combined_loader = CombinedLoader(loaders, mode="max_size_cycle")
            return combined_loader

    def test_dataloader(self):
        val_cfg = deepcopy(self.loader_config)
        val_cfg.update({"batch_size": 1024, "shuffle": False, "drop_last": False})
        return DataLoader(self.val_dataset, **val_cfg)


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


def findMNN(
    x: np.ndarray,
    y: np.ndarray,
    name_x: list[str],
    name_y: list[str],
    mnn_k: int = 5,
    k_components: int = 20,
) -> pd.DataFrame:
    r"""
    Find mutual nearest neighbors (MNN) pairs

    Parameters
    ----------
    x:
        scaled expression of batch A (gene x cell)
    y:
        scaled expression of batch b (gene x cell)
    name_x:
        cell names of batch A
    name_y:
        cell names of batch B
    mnn_k:
        number of nearest neighbors in computing MNN
    k_components:
        SVD components
    """
    # x, y = x.copy(), y.copy()
    z_norm = svd(x, y, k_components)
    z_x, z_y = z_norm[: x.shape[0]], z_norm[x.shape[0] :]

    knn_a2b, _ = knn(z_x, z_y, k=mnn_k)
    knn_b2a, _ = knn(z_y, z_x, k=mnn_k)

    pairs_a2b = create_pairs(knn_a2b)
    pairs_b2a = create_pairs(knn_b2a, reverse=True)
    pairs = pairs_a2b & pairs_b2a
    pairs = np.array(list(pairs))
    mnns = pd.DataFrame(pairs, columns=["cell1", "cell2"])
    mnns["cellname1"] = mnns.cell1.apply(lambda x: name_x[x])
    mnns["cellname2"] = mnns.cell2.apply(lambda x: name_y[x])
    logger.info(f"Found {mnns.shape[0]} MNN pairs")
    return mnns


def svd(x: np.ndarray, y: np.ndarray, k_components: int = 20) -> tuple[np.ndarray]:
    r"""
    Fast SVD

    Parameters
    ----------
    x:
        scaled expression of batch A (gene x cell)
    y:
        scaled expression of batch b (gene x cell)
    k_components:
        number of components to keep

    Returns
    ----------
    z_norm:
        normalized embedding
    """
    logger.debug(f"x shape: {x.shape}, y shape: {y.shape}")
    if x.shape[0] > 1_000_000 or y.shape[0] > 1_000_000:
        logger.debug("Use harmony-based SVD for large dataset.")
        from harmony import harmonize

        # batch
        batch = [0] * x.shape[0] + [1] * y.shape[0]
        batch = pd.DataFrame(batch, columns=["batch"])
        # pca
        z = np.vstack([x, y])
        z, _, _ = randomized_svd(z, n_components=k_components, random_state=0)
        # harmonize
        z_norm = harmonize(z, batch, "batch", use_gpu=True)
        return z_norm

    try:
        dot = torch.from_numpy(x).cuda().half() @ torch.from_numpy(y).T.cuda().half()
        dot = dot.cpu().float().numpy()
        logger.info("Use CUDA for small dataset")
    except:  # noqa
        logger.error(f"CUDA failed: {x.shape}, {y.shape}, use CPU instead.")
        dot = torch.from_numpy(x) @ torch.from_numpy(y).T
        dot = dot.numpy()
    torch.cuda.empty_cache()
    u, s, vh = randomized_svd(dot, n_components=k_components, random_state=0)
    z = np.vstack([u, vh.T])  # gene x k_components
    z = z @ np.sqrt(np.diag(s))  # will reduce the MNN pairs number greatly
    z_norm = l2norm(z)  # follow Seurat
    return z_norm


@njit
def create_pairs(knn_result: np.ndarray, reverse=False) -> tuple:
    r"""
    Create MNN pairs

    Parameters
    ----------
    knn_result:
        knn results
    reverse:
        reverse the pairs
    """
    num_rows = knn_result.shape[0]
    pairs = set()
    for i in prange(num_rows):
        for j in knn_result[i]:
            if reverse:
                pairs.add((j, i))
            else:
                pairs.add((i, j))
    return pairs
