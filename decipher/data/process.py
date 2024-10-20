r"""
Single cell data process
"""
import time

import numpy as np
import scanpy as sc
from addict import Dict
from anndata import AnnData
from loguru import logger
from sklearn.preprocessing import LabelEncoder


def _filter_gene(x: np.ndarray, idx: np.ndarray) -> np.ndarray:
    r"""
    Filter the genes with zero variance
    """
    gene_var = x[:, idx].std(axis=0)
    return idx[gene_var != 0]


def omics_data_process(
    adatas: list[AnnData] | AnnData,
    cell_type: str = None,
    batch_list: list[str] = None,
    group_list: list[str] = None,
    split_by: str = None,
    config: Dict = None,
) -> AnnData:
    r"""
    Process single cell data

    Parameters
    ----------
    adatas
        AnnData or list of slice, each is a spatial slice
    cell_type
        Cell type column name in `obs` of AnnData object
    batch_list
        Batch list for each AnnData
    group_list
        Group list for each AnnData
    split_by
        Split by column name in `obs` of each AnnData object
    config
        Single cell data preprocess config

    Returns
    ----------
    Preprocessed AnnData
    """
    if isinstance(adatas, AnnData):
        adatas = [adatas]

    if isinstance(group_list, list):
        assert len(adatas) == len(group_list), "Length of group_list != adatas"
    else:
        group_list = [None] * len(adatas)

    for i, (adata, group) in enumerate(zip(adatas, group_list)):
        assert "spatial" in adata.obsm.keys(), f"Missing spatial info in {i}th adata"
        if cell_type is not None:
            assert cell_type in adata.obs.keys(), f"Missing cell type in {i}th adata.obs"
            adata.obs["_celltype"] = adata.obs[cell_type].astype(str)
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        # adata.layers["counts"] = adata.X.copy()
        logger.debug(f"adata {i} has {adata.n_obs} cells, {adata.n_vars} genes.")
        adata.obs["_group"] = group

    if len(adatas) == 1:
        adata = adatas[0]
        if split_by in adata.obs.columns:
            batch_int = LabelEncoder().fit_transform(adata.obs[split_by])
            adata.obs["_batch"] = batch_int
    else:
        adata = adatas[0].concatenate(adatas[1:], batch_key="_batch", uns_merge="same")
        if batch_list is not None:
            assert len(batch_list) == len(adatas), "Length of batch_list != adatas"
            le = LabelEncoder().fit(batch_list)
            batch_id = np.hstack([np.repeat(i, adatas[i].n_obs) for i in batch_list])
            adata.obs["_batch"] = le.transform(batch_id)

    batch_col = "_batch" if "_batch" in adata.obs.columns and not config.ignore_batch else None
    n_batch = len(np.unique(adata.obs[batch_col].values)) if batch_col else 1
    logger.info(f"Input: {len(adatas)} slice(s) with {adata.n_obs} cells and {n_batch} batches.")
    return _preprocess_adata(adata, config.pp, batch=batch_col)


def _preprocess_adata(
    adata: AnnData,
    config: Dict,
    batch: str | None = None,
) -> AnnData:
    r"""
    Preprocess single cell data
    """
    logger.info(f"Preprocessing {adata.n_obs} cells.")
    start_time = time.time()
    adata.X = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X

    cell_num, gene_num = adata.X.shape
    if config.min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=config.min_genes)
        logger.warning(f"Filte {cell_num - adata.n_obs} cells with < {config.min_genes} genes.")
    if config.min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=config.min_cells)
        logger.warning(f"Filte {gene_num - adata.n_vars} genes with < {config.min_cells} cells.")

    if config.hvg > 0:
        if config.hvg >= adata.n_vars:
            adata.var["highly_variable"] = True
            logger.warning(f"hvg:{config.hvg} >= n_vars:{adata.n_vars}, set all genes as HVGs.")
        else:
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=config.hvg,
                # layer="counts",
                subset=True,  # NOTE: Must be True for gene selection
                batch_key=batch,
                flavor="seurat_v3",
            )
    # data normalization and scaling
    if config.normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
    if config.log:
        sc.pp.log1p(adata)

    gene_idx = np.arange(adata.n_vars)
    if config.scale and batch is not None and config.per_batch_scale:
        logger.debug("Per batch scaling.")
        batch_name = np.unique(adata.obs[batch].values)
        adatas = [adata[adata.obs[batch] == i] for i in batch_name]
        for ad in adatas:
            gene_idx = _filter_gene(ad.X, gene_idx)
        for i, ad in enumerate(adatas):
            ad = ad[:, gene_idx]
            sc.pp.scale(ad)
            adatas[i] = ad
        adata = adatas[0].concatenate(adatas[1:], batch_key=batch)
        adata.obs[batch] = adata.obs[batch].astype(int)
    elif config.scale:
        gene_idx = _filter_gene(adata.X, gene_idx)
        adata = adata[:, gene_idx]
        sc.pp.scale(adata, max_value=10)
    else:
        logger.warning("No scaling is not recommended !")
    logger.success(f"Preprocessing finished in {time.time() - start_time:.2f} seconds.")
    return adata.copy()
