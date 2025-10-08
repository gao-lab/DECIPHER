r"""
Utility functions
"""
from copy import deepcopy
from subprocess import run

import numpy as np
import requests
import scanpy as sc
import torch
from addict import Dict
from anndata import AnnData
from bs4 import BeautifulSoup
from loguru import logger
from pytorch_lightning import seed_everything as global_seed  # noqa
from rui_utils.gpu import manage_gpu, select_free_gpu  # noqa
from rui_utils.sc import clip_umap, gex_embedding, scanpy_viz  # noqa
from rui_utils.utils import l2norm  # noqa
from scipy.spatial.distance import cdist
from torch import Tensor
from torch_geometric.nn import SimpleConv

sc.set_figure_params(dpi=120, dpi_save=300, format="png", transparent=True)
IMMUNE_MARKER = [
    "CD3G",
    "CD4",
    "CD8A",
    "NKG7",  # T cell and NKs
    "CD14",
    "FCGR3A",
    "SPP1",
    "ITGAX",  # Meoncytes and macrophages
    "MS4A1",
    "POU2AF1",
    "MZB1",  # B cell
    "CD22",  # Mast cells
]
CANCER_MARKER = ["MKI67", "LAMB3"]
OTHER_MARKER = [
    "PLVAP",  # endothelial
    "COL1A1",  # fibroblast
    "PLA2G2A",  # epithelial
    "ACTA2",  # smooth muscle
]

CFG = Dict()
CFG.seed = 0
CFG.work_dir = "DECIPHER"
CFG.device = "gpu" if torch.cuda.is_available() else "cpu"
CFG.device_num = 1

# base config
model_cfg = Dict(
    model_dir="model",
    fix_sc=False,
    # for model
    spatial_emb="attn",
    transformer_layers=3,
    num_heads=1,
    dropout=0.1,
    prj_dims=None,
    temperature_center=0.07,
    temperature_nbr=0.07,
    lr_base=1e-4,
    lr_min=1e-5,
    weight_decay=1e-5,
    first_cycle_steps=99999,
    warmup_steps=200,
    epochs=6,
    nbr_loss_weight=0.5,
    plot=False,
    plot_hist=False,
    # for fit
    device="auto",
    select_gpu=True,
    device_num=1,
    fp16=True,
    patient=10,
    log_every_n_steps=1,
    gradient_clip_val=5.0,
    check_val_every_n_epoch=1,
    max_steps=10_000,
)
loader_cfg = Dict(batch_size=256, shuffle=True, num_workers=8, drop_last=True, pin_memory=True)

contrast_cfg = Dict(model=model_cfg, loader=loader_cfg)

# Config for omics
CFG.omics = deepcopy(contrast_cfg)
CFG.omics.ignore_batch = False
CFG.omics.spatial_graph = Dict(k=20, mode="knn", max_num_neighbors=30)
CFG.omics.mnn = Dict(k_anchor=5, k_components=50, ref_based=True)
CFG.omics.pp = Dict(
    hvg=2000,
    normalize=True,
    log=True,
    scale=True,
    min_genes=0,
    min_cells=0,
    per_batch_scale=True,
)
CFG.omics.num_neighbors = [-1]
CFG.omics.model.update(
    augment=Dict(dropout_gex=0.5, dropout_nbr_prob=-1, mask_hop=-1, max_neighbor=-1),
    emb_dim=128,
    gex_dims=[-1, 256, 32],
    prj_dims=[32, 32, 32],
)
CFG.omics.pretrain = Dict(lr_base=1e-2, lr_min=1e-3, epochs=3, model_dir="pretrain", force=False)

fit_cfg = Dict(
    device="auto",
    select_gpu=True,
    device_num=1,
    fp16=True,
    patient=100,
    log_every_n_steps=1,
    gradient_clip_val=5.0,
    check_val_every_n_epoch=1,
    max_steps=10_000,
)

REGRESS_CFG = Dict(
    select_gpu=False,
    lr_base=1e-3,
    shuffle=False,
    hidden_dim=64,
    val_ratio=0.1,
    test_ratio=0.3,
    fit=deepcopy(fit_cfg),
)
REGRESS_CFG.fit.update(epochs=50)

GENESELECT_CFG = Dict(
    k=30,
    lr_base=1e-3,
    l1_weight=1.0,
    gae_epochs=300,
    test_ratio=0.3,
    gumbel_threshold=0.5,
    num_neighbors=[-1],
    select_gpu=False,
)


def sync_config(cfg: Dict) -> None:
    r"""
    Sync the config
    """
    cfg.omics.model.work_dir = cfg.work_dir
    cfg.omics.model.device_num = cfg.device_num

    cfg.omics.model.gex_dims[-1] = CFG.omics.model.emb_dim
    cfg.omics.model.prj_dims[0] = CFG.omics.model.emb_dim

    if cfg.omics.spatial_graph.mode == "knn":
        cfg.omics.model.augment.max_neighbor = cfg.omics.spatial_graph.k + 1
    elif cfg.omics.spatial_graph.mode == "radius":
        cfg.omics.model.augment.max_neighbor = cfg.omics.spatial_graph.max_num_neighbors + 1


sync_config(CFG)


def install_pyg_dep(torch_version: str | None = None, cuda_version: str | None = None) -> None:
    r"""
    Automatically install PyG dependencies

    Parameters
    ----------
    torch_version
        torch version, e.g. 2.2.1
    cuda_version
        cuda version, e.g. 12.1
    """
    if torch_version is None:
        torch_version = torch.__version__
        torch_version = torch_version.split("+")[0]

    if cuda_version is None:
        cuda_version = torch.version.cuda
        cuda_version = cuda_version.replace(".", "")

    gpu_version = f"torch-{torch_version}+cu{cuda_version}"
    cpu_version = f"torch-{torch_version}+cpu"

    url = "https://data.pyg.org/whl/"
    response = requests.get(url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")
    whl_links = []
    for link in soup.find_all("a"):
        href = link.get("href")
        if href and "torch" in href:
            # replace the + with %2B for URL encoding
            href = href.replace("%2B", "+").replace(".html", "")
            whl_links.append(href)

    if gpu_version in whl_links:
        logger.info(f"Install PyG deps for {gpu_version}")
        version = gpu_version
    elif cpu_version in whl_links:
        logger.warning(f"PyG deps for {gpu_version} not found, use {cpu_version} instead")
        version = cpu_version
    else:
        help_url = "https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html"
        raise ValueError(
            f"PyG deps for {torch_version} not found, please install manually, see {help_url}"
        )

    cmd = f"pip --no-cache-dir install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/{version}.html"
    run(cmd, shell=True)


def estimate_spot_distance(coords: np.ndarray, n_sample: int = 50) -> float:
    r"""
    Estimate the minimum distance between spots

    Parameters
    ----------
    coords
        2D coordinates of spots
    n_sample
        Number of samples to estimate the distance
    """
    n_sample = min(n_sample, coords.shape[0])
    sample_idx = np.random.choice(coords.shape[0], n_sample)
    sample_coords = coords[sample_idx]
    distance = cdist(sample_coords, coords)
    # sort the distance by each row
    distance.sort(axis=1)
    est_distance = np.mean(distance[:, 1])
    return est_distance


def estimate_spot_size(coords: np.ndarray) -> float:
    r"""
    Estimate proper spot size for visualization

    Parameters
    ----------
    coords
        2D coordinates of spots
    """
    x_min, x_max = np.min(coords[:, 0]), np.max(coords[:, 0])
    y_min, y_max = np.min(coords[:, 1]), np.max(coords[:, 1])
    region_area = (x_max - x_min) * (y_max - y_min)
    region_per_cell = region_area / coords.shape[0]
    spot_size = np.sqrt(region_per_cell)
    return spot_size


def nbr_embedding(
    adata: AnnData,
    edge_index: Tensor,
    X_gex: str,
    viz: bool = True,
    n_neighbors: int = 15,
    resolution: float = 0.3,
) -> AnnData:
    r"""
    Get neighbor embedding by aggregating the spatial neighbors
    """
    logger.info("Spatial neighbor embedding...")
    x = torch.tensor(adata.obsm[X_gex], dtype=torch.float32)
    gcn = SimpleConv(aggr="mean")
    embd = gcn(x, edge_index)
    adata.obsm["X_nbr"] = embd.cpu().detach().numpy()

    if viz:
        adata = scanpy_viz(adata, keys=["nbr"], resolution=resolution, n_neighbors=n_neighbors)
    return adata.copy()
