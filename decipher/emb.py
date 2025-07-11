r"""
Get omics embedding and spatial embedding
"""
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import torch
from addict import Dict
from loguru import logger
from torch import Tensor
from torch.utils.data import Dataset, TensorDataset
from torch_geometric.data import Data

from .data.mnn_dataset import LightningScMNNData, LightningSpatialMNNData
from .nn.models.sc import ScSimCLR, ScSimCLRMNN
from .nn.models.spatial import OmicsSpatialSimCLR, OmicsSpatialSimCLRMNN
from .nn.trainer import fit, fit_and_inference


def _update_config(n_obs: int, gene_dim: int, config: Dict) -> Dict:
    r"""
    Update the config before training

    Parameters
    ----------
    n_obs
        number of train elements
    gene_dim
        number of genes
    config
        model config
    """
    config.model.gex_dims[0] = gene_dim
    logger.info(f"Using {gene_dim} genes as model input.")

    # set max steps
    batch_size = config.loader.batch_size * config.model.device_num
    step_per_batch = n_obs // batch_size
    max_steps = config.model.epochs * step_per_batch
    max_steps = min(max_steps, config.model.max_steps)
    max_steps = max(max_steps, step_per_batch)
    config.model.max_steps = max_steps

    # set scheduler config
    if max_steps < 500:
        logger.warning(f"Too few steps {max_steps}, try train more epochs.")
    config.model.warmup_steps = min(int(0.1 * max_steps), config.model.warmup_steps)
    config.model.first_cycle_steps = min(config.model.first_cycle_steps, max_steps)

    return config.copy()


def spatial_emb(
    x: np.ndarray,
    spatial_edge: Tensor,
    config: Dict,
    mnn_dataset: Dataset = None,
    pretrained_model: ScSimCLR = None,
    batch: np.ndarray = None,
    DDP: bool = False,
) -> tuple[np.ndarray, np.ndarray] | None:
    r"""
    Spatial omics embedding

    Parameters
    ----------
    x
        omics expression matrix
    spatial_edge
        spatial edge index
    config
        model config
    mnn_dataset
        mnn dataset
    pretrained_model
        pre-trained single cell model
    batch
        batch index
    DDP
        whether use DDP

    Returns
    ----------
    center_emb
        omics embedding
    nbr_emb
        spatial embedding
    """
    config = _update_config(x.shape[0], x.shape[1], config)
    graph = Data(x=torch.from_numpy(x), edge_index=spatial_edge)
    if batch is not None:
        graph.sc_batch = torch.tensor(batch, dtype=int)
    datamodule = get_graph_datamodule(graph, config, mnn_dataset)

    if mnn_dataset is None:
        model = OmicsSpatialSimCLR(config.model)
    else:
        model = OmicsSpatialSimCLRMNN(config.model)

    if pretrained_model is not None:
        model.center_encoder = deepcopy(pretrained_model.center_encoder)
        # logger.error("NOT support pretrained single cell model.")
    else:
        logger.warning("NOT use pretrained single cell model.")

    if DDP:
        fit(model, datamodule, config.model, show_name="spatial omics DDP")
    else:
        fit_and_inference(model, datamodule, config.model, show_name="spatial omics")
        center_emb, nbr_emb = model.gather_output()
        return center_emb, nbr_emb


def sc_emb(
    x: np.ndarray,
    config: Dict,
    mnn_dataset: Dataset = None,
    batch: np.ndarray = None,
) -> tuple[ScSimCLR, np.ndarray | None]:
    r"""
    Pre-train omics encoder

    Parameters
    ----------
    x:
        omics expression matrix
    config:
        model config
    mnn_dataset:
        mnn dataset
    batch:
        batch index

    Returns
    -----------
    model
        Pre-trained omics encoder
    center_emb
        omics embedding
    """
    config = deepcopy(config)
    config.model.update(config.pretrain)
    mnn_flag = True if mnn_dataset is not None else False
    if not config.pretrain.force:
        try:
            return load_sc_model(config, mnn_flag), None
        except Exception as e:  # noqa
            logger.info(f"Not found pre-trained model: {e}")

    config = _update_config(x.shape[0], x.shape[1], config)

    if batch is not None:
        train_dataset = TensorDataset(torch.from_numpy(x))  # FIXME: add batch
    else:
        train_dataset = TensorDataset(torch.from_numpy(x))
    val_dataset = TensorDataset(torch.from_numpy(x), torch.arange(x.shape[0], dtype=torch.int32))
    datamodule = LightningScMNNData(config.loader, train_dataset, val_dataset, mnn_dataset)

    if mnn_flag:
        model = ScSimCLRMNN(config.model)
    else:
        model = ScSimCLR(config.model)

    if config.model.fix_sc:
        fit_and_inference(model, datamodule, config.model, show_name="single cell")
        center_emb, _ = model.gather_output()
    else:
        fit(model, datamodule, config=config.model, show_name="single cell")
        center_emb = None
    return model, center_emb


def load_sc_model(config, mnn_flag: bool):
    r"""
    Load omics encoder model

    Parameters
    ----------
    config
        model config
    mnn_flag
        whether use mnn
    """
    model_path = Path(config.model.work_dir) / "pretrain"
    # sort by modification time
    model_path = sorted(model_path.glob("*.ckpt"), key=os.path.getmtime)[-1]
    logger.info(f"Loading model from {model_path}")
    kwargs = {"config": config.model}
    if mnn_flag:
        sc_model = ScSimCLRMNN.load_from_checkpoint(model_path, **kwargs)
    else:
        sc_model = ScSimCLR.load_from_checkpoint(model_path, **kwargs)
    logger.success(f"Pre-trained sc model loaded from {model_path}.")
    return sc_model


def load_spatial_model(config, mnn_flag: bool):
    r"""
    Load decipher spatial model

    Parameters
    ----------
    config
        model config
    mnn_flag
        whether use mnn
    """
    model_path = Path(config.model.work_dir) / "model"
    model_path = sorted(model_path.glob("*.ckpt"), key=os.path.getmtime)[-1]
    logger.info(f"Loading model from {model_path}")
    config.model.device_num = 1
    kwargs = {"config": config.model}
    if mnn_flag:
        model = OmicsSpatialSimCLRMNN.load_from_checkpoint(model_path, **kwargs)
    else:
        model = OmicsSpatialSimCLR.load_from_checkpoint(model_path, **kwargs)
    logger.success(f"Pre-trained spatial model loaded from {model_path}.")
    return model


def get_graph_datamodule(graph: Data, config: Dict, mnn_dataset=None) -> LightningSpatialMNNData:
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
        config.loader,
        mnn_dataset,
        num_neighbors=config.num_neighbors,
        input_train_nodes=mask,
        input_val_nodes=mask,
        input_test_nodes=mask,
        subgraph_type="bidirectional",
        disjoint=True,
        eval_loader_kwargs={"drop_last": False},
    )
    return datamodule
