r"""
Single cell data process
"""
import time

import numpy as np
import scanpy as sc
from anndata import AnnData
from loguru import logger
from sklearn.preprocessing import LabelEncoder

from .. import CFG


def omics_data_process(
    adata_raw: list[AnnData] | AnnData,
    split_by: str = None,
    preprocess: bool = True,
) -> tuple[AnnData, np.ndarray | None]:
    r"""
    Process single cell data

    Parameters
    ----------
    adatas
        AnnData or list of slice, each is a spatial slice
    split_by
        Split by column name in `obs` of each AnnData object
    config
        Single cell data preprocess config

    Returns
    ----------
    Preprocessed AnnData
    """
    batch_idx = None
    if isinstance(adata_raw, AnnData):
        if split_by is not None and split_by in adata_raw.obs.columns:
            batch_idx = LabelEncoder().fit_transform(adata_raw.obs[split_by])
        adata = adata_raw.copy()
    else:
        adata = adata_raw[0].concatenate(adata_raw[1:], batch_key="batch", uns_merge="same")
        batch_idx = LabelEncoder().fit_transform(adata.obs["batch"])
    del adata_raw

    if batch_idx is not None and not CFG.pp.ignore_batch:
        logger.info(f"Detected {np.unique(batch_idx).shape[0]} batches.")
        adata.obs["_batch"] = batch_idx

    if not preprocess:
        logger.warning("Skip preprocessing steps, only for advanced users.")
        return adata.X.astype(np.float32), adata.obsm["spatial"], batch_idx

    logger.info(f"Preprocessing {adata.n_obs} cells.")
    start_time = time.time()

    cell_num, gene_num = adata.X.shape
    if CFG.pp.min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=CFG.pp.min_genes)
        logger.info(f"Filter {cell_num - adata.n_obs} cells with < {CFG.pp.min_genes} genes.")
    if CFG.pp.min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=CFG.pp.min_cells)
        logger.info(f"Filter {gene_num - adata.n_vars} genes with < {CFG.pp.min_cells} cells.")

    if CFG.pp.hvg > 0:
        if CFG.pp.hvg >= adata.n_vars:
            adata.var["highly_variable"] = True
            logger.info(f"hvg:{CFG.pp.hvg} >= n_vars:{adata.n_vars}, set all genes as HVGs.")
        else:
            logger.info(f"Selecting {CFG.pp.hvg} highly variable genes.")
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=CFG.pp.hvg,
                subset=True,  # NOTE: Must be True for gene selection
                batch_key="_batch" if "_batch" in adata.obs.columns else None,
                flavor="seurat_v3",
            )

    # data normalization and scaling
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # NOTE: we found per-batch scaling is necessary to find more knn pairs
    if batch_idx is not None and not CFG.pp.ignore_batch:
        logger.info("Per batch scaling.")
        gene_idx = np.arange(adata.n_vars)
        batch_name = np.unique(adata.obs["_batch"].values)
        adatas = [adata[adata.obs["_batch"] == i] for i in batch_name]
        for i, ad in enumerate(adatas):
            ad = ad[:, gene_idx]
            sc.pp.scale(ad)
            adatas[i] = ad
        adata = adatas[0].concatenate(adatas[1:], batch_key="_batch")
        del adatas
        adata.obs["_batch"] = adata.obs["_batch"].astype(int)
    else:
        sc.pp.scale(adata, max_value=10)
    logger.success(f"Preprocessing finished in {time.time() - start_time:.2f} seconds.")
    return adata.X.astype(np.float32), adata.obsm["spatial"], batch_idx
