r"""
Utility functions
"""
import os
import random
import time
from copy import deepcopy
from heapq import nlargest
from subprocess import run

import numpy as np
import pandas as pd
import pynvml
import scanpy as sc
import torch
import yaml
from addict import Dict
from anndata import AnnData
from joblib import Parallel, delayed
from loguru import logger
from torch import Tensor
from torch_geometric.nn import SimpleConv

try:
    import cupy
    import rapids_singlecell as rsc
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator

    assert torch.cuda.is_available()
    RSC_FLAG = True
    sc_ = rsc
except:  # noqa
    logger.warning("Rapids not avaliable, use Scanpy")
    RSC_FLAG = False
    sc_ = sc


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


def save_dict(dic: dict | Dict, path: str) -> None:
    r"""
    Save the dict to a yaml file
    """
    if isinstance(dic, Dict):
        dic = dic.to_dict()
    with open(path, "w") as f:
        yaml.dump(dic, f, default_flow_style=False)


def GetRunTime(func):
    r"""
    Decorator to get the run time of a function
    """

    def call_func(*args, **kwargs):
        begin_time = time.time()
        ret = func(*args, **kwargs)
        end_time = time.time()
        Run_time = end_time - begin_time
        logger.debug(f"{func.__name__} run time: {Run_time:.2f}s")
        return ret

    return call_func


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

    if torch_version < "2.0":
        raise ValueError(f"PyG only support torch>=2.0, but get {torch_version}")
    elif "2.0" <= torch_version < "2.1":
        torch_version = "2.0.0"
    elif "2.1" <= torch_version < "2.2":
        torch_version = "2.1.0"
    elif "2.2" <= torch_version < "2.3":
        torch_version = "2.2.0"
    elif "2.3" <= torch_version < "2.4":
        torch_version = "2.3.0"
    else:
        raise ValueError(f"Automatic install only support torch<=2.3, but get {torch_version}")

    if "cu" in cuda_version and not torch.cuda.is_available():
        logger.warning("CUDA is not available, try install CPU version, but may raise error.")
        cuda_version = "cpu"
    elif cuda_version >= "12.1":
        cuda_version = "cu121"
    elif "11.8" <= cuda_version < "12.1":
        cuda_version = "cu118"
    elif "11.7" <= cuda_version < "11.8":
        cuda_version = "cu117"
    else:
        raise ValueError(f"PyG only support cuda>=11.7, but get {cuda_version}")

    url = "https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html"
    if torch_version in ["2.2.0", "2.1.0"] and cuda_version == "cu117":
        raise ValueError(
            f"PyG not support torch-{torch_version} with cuda-11.7, please check {url}"
        )
    if torch_version == "2.0.0" and cuda_version == "cu121":
        raise ValueError(f"PyG not support torch-2.0.* with cuda-12.1, please check {url}")

    logger.info(f"Installing PyG dependencies for torch-{torch_version} and cuda-{cuda_version}")
    cmd = f"pip --no-cache-dir install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-{torch_version}+{cuda_version}.html"
    run(cmd, shell=True)


def select_free_gpu(n: int = 1) -> list[int | None]:
    r"""
    Get torch computation device automatically

    Parameters
    ----------
    n
        Number of GPUs to request

    Returns
    -------
    n_devices
        list of devices index
    """
    assert n > 0
    try:
        pynvml.nvmlInit()
        devices = os.environ.get("CUDA_VISIBLE_DEVICES", None)
        devices = (
            range(pynvml.nvmlDeviceGetCount())
            if devices is None
            else [int(d.strip()) for d in devices.split(",") if d != ""]
        )
        free_mems = {
            i: pynvml.nvmlDeviceGetMemoryInfo(pynvml.nvmlDeviceGetHandleByIndex(i)).free
            for i in devices
        }
        n_devices = nlargest(n, free_mems, free_mems.get)
        if len(n_devices) == 0:
            raise pynvml.NVMLError("GPU disabled.")
        logger.info(f"Using GPU {n_devices} as computation device.")
        # os.environ["CUDA_VISIBLE_DEVICES"] = ",".join(map(str, n_devices))
        return n_devices
    except pynvml.NVMLError:  # pragma: no cover
        logger.warning("No GPU available.")
        return [None] * n


def global_seed(seed: int, cuda_deterministic: bool = False) -> None:
    r"""
    Set seed in global scope

    Parameters
    ----------
    seed
        int
    cuda_deterministic
        set `True` can make cuda more deterministic, but more slower
    """
    seed = seed if seed != -1 else torch.seed()
    if seed > 2**32 - 1:
        seed = seed >> 32

    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    # Speed-reproducibility tradeoff https://pytorch.org/docs/stable/notes/randomness.html
    if cuda_deterministic:  # slower, more reproducible
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = True
    else:  # faster, less reproducible
        torch.backends.cudnn.deterministic = False
        torch.backends.cudnn.benchmark = True
    logger.success(f"Global seed set to {seed}.")


def euclidean_distance(
    x: np.ndarray | Tensor,
    y: np.ndarray | Tensor | None = None,
    cuda: bool = False,
) -> np.ndarray:
    r"""
    Calculate the euclidean distance between two arrays

    Parameters
    ----------
    x
        array
    y
        array (if `None`, `y` will be set to `x`)
    cuda
        whether to use GPU
    """
    if y is None:
        y = x
    assert x.shape[1] == y.shape[1], "The dimension of two arrays should be the same"
    x, y = torch.tensor(x, dtype=torch.float), torch.tensor(y, dtype=torch.float)
    if cuda and torch.cuda.is_available():
        x, y = x.cuda(), y.cuda()

    m, n = x.shape[0], y.shape[0]
    x = x.view(m, 1, -1)
    y = y.view(1, n, -1)
    return torch.norm(x - y, dim=-1).numpy()


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
    distance = euclidean_distance(sample_coords, coords)
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


def manage_gpu(gpu_id: int, memory_strategy: str | None = None):
    r"""
    Manage Rapids GPU index and memory strategy
    """
    assert memory_strategy in ["large", "fast", "auto", None]
    if memory_strategy is not None:
        if memory_strategy == "large":
            managed_memory, pool_allocator = True, False
        elif memory_strategy == "fast":
            managed_memory, pool_allocator = False, True
        rmm.reinitialize(
            managed_memory=managed_memory,
            pool_allocator=pool_allocator,
            devices=gpu_id,
        )
    else:
        rmm.reinitialize(devices=gpu_id)
    cupy.cuda.set_allocator(rmm_cupy_allocator)
    cupy.cuda.Device(gpu_id).use()
    logger.info(f"Using GPU {gpu_id} and {memory_strategy} memory strategy.")


def harmony_pytorch(emb: np.ndarray, obs: pd.DataFrame, batch_key: str, seed: int) -> np.ndarray:
    r"""
    Run `harmony` in pytorch version

    Parameters
    ----------
    emb
        Embedding matrix (such as PCA)
    obs
        `obs` of AnnData
    batch_key
        Batch information column in `obs`
    seed
        Random seed
    """
    try:
        from harmony import harmonize
    except ImportError:
        raise ImportError("Please install `harmony-pytorch` package first.")
    GPU_FLAG = True if torch.cuda.is_available() and emb.shape[0] > 2e3 else False
    if GPU_FLAG:
        logger.warning("Use GPU for harmony")
    return harmonize(
        emb,
        obs,
        batch_key=batch_key,
        random_state=seed,
        max_iter_harmony=30,
        use_gpu=GPU_FLAG,
    )


def gex_embedding(
    adata: AnnData,
    method: str = "pca",
    filter: bool = False,
    min_gene: int = 100,
    min_cell: int = 30,
    batch_key: str | None = None,
    call_hvg: bool = True,
    n_top_genes: int = 2000,
    hvg_by_batch: bool = False,
    hvg_only: bool = True,
    n_comps: int = 50,
    seed: int = 0,
    viz: bool = True,
    forbid_rapids: bool = False,
    rapids_after_scale: bool = True,
    resolution: float = 0.8,
    harmony_version: str = "torch",
    emb_only: bool = False,
    gpu_id: str | int = "auto",
    memory_strategy: str | None = None,
) -> AnnData:
    r"""
    Gene expression embedding via Scanpy pipeline

    Parameters
    ----------
    adata
        AnnData object
    method
        method for embedding, `pca` or `harmony`
    filter
        If filter cells and genes, default is False
    min_gene
        Minimum number of genes for each cell
    min_cell
        Minimum number of cells for each gene
    batch_key
        Key for batch information
    call_hvg
        If call highly variable genes, default is True
    n_top_genes
        Number of top genes, default is 2000
    hvg_by_batch
        If call HVG by batch, default is False
    hvg_only
        If subet HVG of adata only, default is True
    n_comps
        Number of components for PCA, default is 50
    seed
        Random seed
    viz
        If visualize the embedding, default is True
    forbid_rapids
        If forbid rapids, default is False
    rapids_after_scale
        If use rapids after scale, default is False, should set to True on large dataset
    resolution
        Resolution for leiden clustering, default is 0.8
    harmony_version
        Version of harmony, `torch` or `rapids`, default is `torch`
    emb_only
        If return embedding only, default is False
    gpu_id
        GPU index, default is `auto`
    memory_strategy
        Memory strategy for Rapids, `large` or `fast`, default is None

    Warning
    ----------
    Input can not be the View of anndata
    """
    logger.info("Gene expression embedding...")
    assert method.lower() in ["pca", "harmony"], f"Method {method} not supported"
    if batch_key is not None:
        method = "harmony"
        logger.info(f"Use {method} for batch correction.")

    forbid_rapids = forbid_rapids if RSC_FLAG else True
    if not forbid_rapids:
        if gpu_id == "auto":
            gpu_id = select_free_gpu(1)[0]
        else:
            assert isinstance(gpu_id, int), "Invalid gpu_id"
        manage_gpu(gpu_id, memory_strategy)

    if filter:
        raw_cell, raw_genes = adata.n_obs, adata.n_vars
        sc.pp.filter_cells(adata, min_genes=min_gene)
        sc.pp.filter_genes(adata, min_cells=min_cell)
        logger.info(f"Filter {raw_cell-adata.n_obs} cells and {raw_genes-adata.n_vars} genes")

    sc_ = sc  # default backend is scanpy
    if not rapids_after_scale:
        rsc.get.anndata_to_GPU(adata)
        sc_ = rsc

    call_hvg = call_hvg and adata.n_vars > n_top_genes
    if call_hvg:
        if isinstance(n_top_genes, int):
            sc_.pp.highly_variable_genes(
                adata,
                n_top_genes=n_top_genes,
                flavor="seurat_v3",
                batch_key=batch_key if hvg_by_batch else None,
            )
        elif isinstance(n_top_genes, list):
            adata.var["highly_variable"] = False
            n_top_genes = list(set(adata.var.index).intersection(set(n_top_genes)))
            adata.var.loc[n_top_genes, "highly_variable"] = True
            logger.warning(f"{len(n_top_genes)} genes detected given lists.")
    else:
        logger.warning("All genes are seen as highly variable genes.")
        adata.var["highly_variable"] = True

    if hvg_only:
        adata = adata[:, adata.var["highly_variable"]]

    sc_.pp.normalize_total(adata, target_sum=1e4)
    sc_.pp.log1p(adata)
    sc_.pp.scale(adata, max_value=10)  # Use a lot GPU memory
    if rapids_after_scale and not forbid_rapids:
        rsc.get.anndata_to_GPU(adata)
        sc_ = rsc

    sc_.pp.pca(adata, n_comps=n_comps)
    if method.lower() == "pca":
        adata.obsm["X_gex"] = adata.obsm["X_pca"]
    elif method.lower() == "harmony":
        if harmony_version == "rapids":
            rsc.pp.harmony_integrate(adata, key=batch_key)
            adata.obsm["X_gex"] = adata.obsm["X_pca_harmony"]
        elif harmony_version in ["pytorch", "torch"]:
            pca_emb = adata.obsm["X_pca"]
            if not isinstance(pca_emb, np.ndarray):
                pca_emb = pca_emb.to_numpy()
            adata.obsm["X_gex"] = harmony_pytorch(adata.obsm["X_pca"], adata.obs, batch_key, seed)
    if emb_only:
        return adata.obsm["X_gex"]

    if viz:
        sc_.pp.neighbors(adata, use_rep="X_gex")  # TODO: support approx neighbor, add n_neighbors
        sc_.tl.leiden(adata, resolution=resolution)
        sc_.tl.umap(adata)

    if sc_ is not sc:
        sc_.get.anndata_to_CPU(adata)
        rmm.reinitialize()
    return adata.copy()


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


def _scanpy_viz(
    adata: AnnData,
    gpu_id: int,
    rapids: bool = True,
    memory_strategy: str = None,
    approx: bool = False,
    leiden: bool = True,
    resolution: float = 0.5,
    n_neighbors: int = 15,
) -> tuple[np.ndarray, np.ndarray | None]:
    neighbor_kwargs = dict(n_neighbors=n_neighbors, use_rep="X")
    if rapids and RSC_FLAG:
        manage_gpu(gpu_id, memory_strategy)
        rsc.get.anndata_to_GPU(adata)
        logger.debug("Use rapids for visualization.")
        neighbor_kwargs["algorithm"] = "cagra" if approx else "brute"

    sc_.pp.neighbors(adata, **neighbor_kwargs)
    sc_.tl.umap(adata)
    umap = adata.obsm["X_umap"]
    if leiden:
        sc_.tl.leiden(adata, resolution=resolution)
    leiden = adata.obs["leiden"].values if leiden else None
    return umap, leiden


def scanpy_viz(
    adata: AnnData,
    keys: list[str] = ["center", "nbr"],
    gpu_id: int | list = None,
    resolution: float = 0.5,
    rapids: bool = True,
    approx: bool | str = False,
    memory_strategy: str | None = None,
    leiden: bool = True,
    n_neighbors: int = 15,
    n_jobs=-1,
) -> sc.AnnData:
    r"""
    Fast clustering and visualization via scanpy/rapids

    Parameters
    -----------
    adata
        AnnData object
    keys
        Keys for visualization, should match the keys in `adata.obsm[X_{key}]`
    gpu_id
        GPU index, default is None
    resolution
        Resolution for leiden clustering, default is 0.5
    rapids
        If use rapids for visualization, default is True
    approx
        If use approximate nearest neighbors, default is False
    memory_strategy
        Memory strategy for Rapids, `large` or `fast`, default is None
    leiden
        If use leiden clustering, default is True
    n_neighbors
        Number of neighbors, default is 15
    n_jobs
        jobs for parallel computation, default is -1
    """
    kwargs = dict(
        resolution=resolution,
        rapids=rapids,
        memory_strategy=memory_strategy,
        approx=approx,
        leiden=leiden,
        n_neighbors=n_neighbors,
    )

    # make agency adata
    viz_adatas = []
    keys = [keys] if isinstance(keys, str) else keys
    for key in keys:
        if key.lower() == "x":
            viz_adatas.append(adata.copy())
        elif f"X_{key}" in adata.obsm.keys():
            viz_adatas.append(sc.AnnData(X=adata.obsm[f"X_{key}"].copy()))
        else:
            logger.warning(f"Key {key} not found in adata.obsm, skip.")

    # select gpu
    if RSC_FLAG:
        n_jobs = len(viz_adatas) if n_jobs == -1 else n_jobs
        if gpu_id is None:
            gpu_ids = select_free_gpu(len(viz_adatas))
            if len(gpu_ids) == 1:
                gpu_ids = gpu_ids * len(viz_adatas)
        elif isinstance(gpu_id, int):
            gpu_ids = [gpu_id] * len(viz_adatas)
        elif isinstance(gpu_id, list):
            assert len(gpu_id) == len(viz_adatas)
            gpu_ids = gpu_id
        else:
            raise ValueError(f"Invalid gpu_id: {gpu_id}")
    else:
        gpu_ids = [None] * len(viz_adatas)
        n_jobs = 1

    # parallel computation
    if n_jobs == 1:
        results = [
            _scanpy_viz(_adata, _gpu_id, **kwargs) for _adata, _gpu_id in zip(viz_adatas, gpu_ids)
        ]
    else:
        results = Parallel(n_jobs=n_jobs)(
            delayed(_scanpy_viz)(_adata, gpu_id, **kwargs)
            for _adata, gpu_id in zip(viz_adatas, gpu_ids)
        )
    umaps, leidens = zip(*results)

    for key, umap, leiden in zip(keys, umaps, leidens):
        adata.obsm[f"X_umap_{key}"] = umap
        if leiden is not None:
            adata.obs[f"leiden_{key}"] = leiden
    if RSC_FLAG:
        rmm.reinitialize()
    return adata


def clip_umap(array, percent: float = 0.1):
    r"""
    Clip the outlier to the percentile range of UMAP

    Parameters
    ----------
    array
        UMAP array
    percent
        Percentile for clipping, default is 0.1
    """
    assert 0 < percent < 50
    half_percent = percent / 2
    percentile_down = np.percentile(array, half_percent, axis=0)
    percentile_up = np.percentile(array, 100 - half_percent, axis=0)
    return np.clip(array, a_min=percentile_down, a_max=percentile_up)


def l2norm(mat: np.ndarray) -> np.ndarray:
    r"""
    L2 norm of numpy array
    """
    stats = np.sqrt(np.sum(mat**2, axis=1, keepdims=True)) + 1e-9
    mat = mat / stats
    return mat
