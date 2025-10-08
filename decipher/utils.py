r"""
Utility functions
"""
from subprocess import run

import numpy as np
import requests
import scanpy as sc
import torch
from anndata import AnnData
from bs4 import BeautifulSoup
from click import command, option
from loguru import logger
from pytorch_lightning import seed_everything as global_seed  # noqa
from rui_utils.gpu import manage_gpu, select_free_gpu  # noqa
from rui_utils.sc import clip_umap, gex_embedding, scanpy_viz  # noqa
from rui_utils.utils import l2norm  # noqa
from scipy.spatial.distance import cdist
from torch import Tensor

sc.set_figure_params(dpi=120, dpi_save=300, format="png", transparent=True)


@command()
@option("--work_dir", help="work directory", required=True, type=str)
@option("--gpus", help="number of gpus to use, -1 for all available gpus", default=-1, type=int)
def decipher_ddp_spatial(work_dir, gpus):
    from decipher import CFG, DECIPHER

    CFG.sp_model.trainer.device_num = gpus
    model = DECIPHER(work_dir, recover=True)
    model.fit_spatial()


@command()
@option("--work_dir", help="work directory", required=True, type=str)
@option("--gpus", help="number of gpus to use, -1 for all available gpus", default=-1, type=int)
def decipher_ddp_sc(work_dir, gpus):
    from decipher import CFG, DECIPHER

    CFG.sc_model.trainer.device_num = gpus
    model = DECIPHER(work_dir, recover=True)
    model.fit_sc()


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
    try:
        from torch_geometric.nn import SimpleConv
    except ImportError:
        raise ImportError("Please install PyG to use this function.")

    logger.info("Spatial neighbor embedding...")
    x = torch.tensor(adata.obsm[X_gex], dtype=torch.float32)
    gcn = SimpleConv(aggr="mean")
    embd = gcn(x, edge_index)
    adata.obsm["X_nbr"] = embd.cpu().detach().numpy()

    if viz:
        adata = scanpy_viz(adata, keys=["nbr"], resolution=resolution, n_neighbors=n_neighbors)
    return adata.copy()
