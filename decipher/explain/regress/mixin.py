r"""
Mixin for explain embeddings
"""
import glob
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import ray
import scanpy as sc
import torch
import yaml
from addict import Dict
from loguru import logger

from ...utils import REGRESS_CFG
from .regression import ray_train_regress, train_regress


class RegressMixin:
    r"""
    Mixin for explain train results
    """

    def train_regress_explain(
        self,
        adata: sc.AnnData,
        cell_type: str,
        explain_dir: str = "explain",
        subsets: list | str = None,
        batch: str = None,
        regress_per_celltype: bool = True,
        reverse_regress: bool = False,
        min_cells: int = 100,
        user_cfg: Dict = {},
        n_jobs: int = -1,
    ):
        r"""
        Train explain model

        Specify `cell_type` and `subsets` to train explain model on specific cell types.

        Parameters
        ----------
        adata
            AnnData object
        cell_type
            Cell type column name in `obs` of AnnData
        explain_dir
            Explain work directory, will be a subdirectory of `work_dir`
        subsets
            Cell type subsets to train explain model, default use all cell types
        batch
            Batch column name in `obs` of AnnData
        regress_per_celltype
            If train regress model per cell type, otherwise train on all cell types (not recommended)
        reverse_regress
            If train reverse regress model
        min_cells
            Minimum cells for each cell type
        user_cfg
            User config for explain model
        n_jobs
            Number of jobs to run in parallel
        """
        start_time = time.time()
        assert hasattr(self, "nbr_emb"), "Lack self.nbr_emb, Please run `fit_omics()` first."

        # set config
        cfg_explain = deepcopy(REGRESS_CFG)
        explain_work_dir = self.work_dir / explain_dir
        explain_work_dir.mkdir(exist_ok=True, parents=True)
        cfg_explain.work_dir = str(explain_work_dir)

        cfg_explain.update(user_cfg)
        cfg_explain.center_dim = self.center_emb.shape[1]
        cfg_explain.nbr_dim = self.nbr_emb.shape[1]
        cfg_explain.work_dir = str(explain_work_dir)
        with open(explain_work_dir / "regress_config.yaml", "w") as f:
            yaml.dump(cfg_explain.to_dict(), f)
        logger.debug(f"Regress config: {cfg_explain}")

        # prepare data
        nbr_emb = torch.from_numpy(self.nbr_emb)
        center_emb = torch.from_numpy(self.center_emb)
        if cfg_explain.shuffle:
            logger.warning("Shuffle center_emb")
            nbr_emb = nbr_emb[torch.randperm(nbr_emb.size(0))]

        # get cell type data
        cell_types = adata.obs[cell_type].value_counts().index.tolist()
        if isinstance(subsets, str):
            assert subsets in cell_types, f"Cell type {subsets} not found."
            cell_types = [subsets]
        elif isinstance(subsets, list):
            assert all([x in cell_types for x in subsets]), f"Cell type {subsets} not found."
            cell_types = subsets
        logger.info(f"Run on {len(cell_types)} cell types: {cell_types}")

        # split to batch
        batches = adata.obs[batch].value_counts().index.tolist() if batch else [None]

        # prepare data
        xy_list, regress_dirs = [], []
        for _batch in batches:
            if _batch is None:
                subset_batch = np.ones(adata.shape[0], dtype=bool)
            else:
                subset_batch = (adata.obs[batch] == _batch).values

            for i, _cell_type in enumerate(cell_types):
                subset_celltype = (adata.obs[cell_type] == _cell_type).values

                _name = f"batch:{_batch}_celltype:{_cell_type}"
                _name = _name.replace(" ", "_").replace("/", "_")
                subset = subset_batch & subset_celltype

                if subset.sum() < min_cells:
                    logger.warning(f"Skip {_name} with {subset.sum()} cells (<{min_cells}).")
                    continue
                if regress_per_celltype:
                    xy_list.append((center_emb[subset], nbr_emb[subset]))
                    regress_dirs.append(f"regress_{_name}")
                else:
                    if i == 0:
                        xy_list.append((center_emb, nbr_emb))
                        regress_dirs.append(f"regress_batch:{_batch}_all")

        # train models
        train_regress_parallel(
            xy_list, cfg_explain, regress_dirs, n_jobs, reverse_regress, adata.obs[cell_type]
        )
        self.explain_df = merge_regress_results(explain_work_dir, start_time)
        logger.success(f"Explain model training finished in {time.time() - start_time:.2f}s.")


def merge_regress_results(work_dir: Path, start_time) -> pd.DataFrame | None:
    r"""
    Merge explain results to a DataFrame
    """
    metrics = glob.glob(str(work_dir / "*regress*/*.json"))
    if len(metrics) == 0:
        logger.warning("No explain results found.")
    else:
        results = {}
        for file in metrics:
            if "all" in str(file):
                continue
            # if time of file is earlier than start_time, skip
            if start_time > Path(file).stat().st_mtime:
                continue
            with open(file) as f:
                metric = yaml.load(f, Loader=yaml.FullLoader)
            name = file.split("regress_")[-1].split("/")[0]
            results[name] = metric
        df = pd.DataFrame(results).T
        df.to_csv(work_dir / "explain_results.csv")
        return df


def train_regress_parallel(
    xy_list, cfg_explain, regress_dirs, n_jobs, reverse_regress, cell_type
) -> list:
    r"""
    train regress models in parallel
    """
    n_datasets = len(xy_list)
    if reverse_regress:
        logger.info("Also train reverse regression models...")
        cfg_explain_rev = deepcopy(cfg_explain)
        cfg_explain_rev.nbr_dim = cfg_explain.center_dim
        cfg_explain_rev.center_dim = cfg_explain.nbr_dim
        # reverse xy_list
        xy_list_rev = [(y, x) for x, y in xy_list]
        regress_dirs_rev = [x.replace("regress", "regress_reverse") for x in regress_dirs]
        xy_list += xy_list_rev
        regress_dirs += regress_dirs_rev
        cfg_explain_list = [cfg_explain] * n_datasets + [cfg_explain_rev] * n_datasets
    else:
        cfg_explain_list = [cfg_explain] * n_datasets

    if n_jobs == 1:
        for xy, cfg, save_dir in zip(xy_list, cfg_explain_list, regress_dirs):
            train_regress(*xy, cfg, save_dir, cell_type)
    else:
        ray.init(ignore_reinit_error=True)
        queue = [
            ray_train_regress.remote(*xy, cfg, save_dir, cell_type)
            for xy, cfg, save_dir in zip(xy_list, cfg_explain_list, regress_dirs)
        ]
        ray.get(queue)
        ray.shutdown()
