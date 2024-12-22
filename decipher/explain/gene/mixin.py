r"""
Mixin for select genes
"""
import time
from copy import deepcopy

import pandas as pd
import ray
import scanpy as sc
import torch
import yaml
from addict import Dict
from loguru import logger
from scanpy._utils import check_nonnegative_integers
from scipy.sparse import issparse
from torch_geometric.data import Batch, Data

from ...graphic.build import build_graph
from ...utils import GENESELECT_CFG
from .gene_selection import ray_train_GAE, train_GAE
from .lr import get_lr_expr


class GeneSelectMixin:
    r"""
    Mixin for gene selection
    """

    def train_gene_select(
        self,
        adata: sc.AnnData,
        cell_type: str,
        sub_dir: str = "explain",
        subsets: list | str = None,
        batch: str = None,
        per_batch: bool = False,
        lr_mode: bool = False,
        lr_data: pd.DataFrame = None,
        lr_radius: float = 0.5,
        check_data: bool = True,
        normalize: bool = True,
        min_cells: int = 1000,
        user_cfg: Dict = {},
        disable_gpu: bool = False,
        n_jobs: int = -1,
    ):
        r"""
        Train gene selection model by cell type

        Parameters
        ----------
        adata
            AnnData object
        cell_type
            Cell type column name in `obs` of AnnData
        sub_dir
            sub directory to save explain model
        subsets
            Cell type subsets to train explain model, default use all cell types
        batch
            Batch column name in `obs` of AnnData
        lr_mode
            If run in Ligand-Receptor mode
        lr_data
            Ligand-Receptor pairs
        lr_radius
            Radius for ligand-receptor pairs
        check_data
            If check adata.X is non-negative integers
        normalize
            If normalize adata.X
        min_cells
            Minimum cells for each cell type
        user_cfg
            User config for explain model
        disable_gpu
            If disable GPU (only when GPU memory is not enough)
        n_jobs
            Number of jobs to run in parallel
        """
        start_time = time.time()
        assert hasattr(self, "nbr_emb"), "Lack self.nbr_emb, Please run `fit_omics()` first."
        # set config
        cfg = deepcopy(GENESELECT_CFG)
        work_dir = self.work_dir / sub_dir
        work_dir.mkdir(exist_ok=True, parents=True)
        cfg.update(user_cfg)
        cfg.center_dim = self.center_emb.shape[1]
        cfg.nbr_dim = self.nbr_emb.shape[1]
        cfg.work_dir = str(work_dir)
        with open(work_dir / "gene_select_config.yaml", "w") as f:
            yaml.dump(cfg.to_dict(), f)
        logger.debug(f"Gene select config: {cfg}")

        # prepare data
        if check_data:
            assert check_nonnegative_integers(adata.X), "Please ensure adata.X is raw count."

        if not lr_mode and adata.n_vars > 3000:
            logger.warning("Using all genes is not recommended, subset top 3000 HVG only.")
            sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", subset=True)

        if normalize:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        if lr_mode:
            logger.info("Run in Ligand-Receptor mode")
            expr, lr = get_lr_expr(adata, lr_data, lr_radius)
            cfg.expr_dim = expr.shape[1]
            lr.to_csv(work_dir / "lr.csv")
        else:
            cfg.expr_dim = adata.n_vars
            expr = torch.from_numpy(adata.X.A) if issparse(adata.X) else torch.from_numpy(adata.X)
        nbr_emb = torch.from_numpy(self.nbr_emb)
        center_emb = torch.from_numpy(self.center_emb)
        assert expr.shape[0] == nbr_emb.shape[0] == center_emb.shape[0], "Data shape not match."

        # get cell type data
        cell_types = adata.obs[cell_type].value_counts().index.tolist()
        if isinstance(subsets, str):
            assert subsets in cell_types, f"Cell type {subsets} not found."
            cell_types = [subsets]
        elif isinstance(subsets, list):
            assert all([x in cell_types for x in subsets]), f"Cell type {subsets} not found."
            cell_types = subsets
        logger.info(f"Run on {len(cell_types)} cell types: {cell_types}")

        batches = adata.obs[batch].value_counts().index.tolist() if batch else [None]
        graph_list, save_dir_list = [], []
        for _cell_type in cell_types:
            subset_cell_type = (adata.obs[cell_type] == _cell_type).values
            if subset_cell_type.sum() < min_cells:
                logger.warning(
                    f"Skip {_cell_type} with {subset_cell_type.sum()} cells (<{min_cells})."
                )
                continue
            name = f"select_celltype_{_cell_type}"
            name = name.replace(" ", "_").replace("/", "_")
            if not per_batch:
                save_dir_list.append(name)

            if len(batches) == 1:
                edge_index = build_graph(
                    nbr_emb[subset_cell_type].numpy(), mode="knn", k=cfg.k, pyg_backend=False
                )
                graph = Data(
                    x=center_emb[subset_cell_type],
                    expr=expr[subset_cell_type],
                    edge_index=edge_index,
                )
                graph_list.append(graph)
                if per_batch:
                    save_dir = f"{name}_batch_{batches[0]}"
                    save_dir_list.append(save_dir)
            else:
                batch_graph_list = []
                for _batch in batches:
                    subset_batch = (adata.obs[batch] == _batch).values
                    subset = subset_batch & subset_cell_type

                    if subset.sum() < min_cells:
                        logger.warning(
                            f"Skip cell type {_cell_type} in batch {_batch} with {subset.sum()} cells (<{min_cells})."
                        )
                        continue

                    edge_index = build_graph(
                        nbr_emb[subset].numpy(), mode="knn", k=cfg.k, pyg_backend=False
                    )
                    batch_graph = Data(
                        x=center_emb[subset],
                        expr=expr[subset],
                        edge_index=edge_index,
                    )
                    batch_graph_list.append(batch_graph)
                    if per_batch:
                        save_dir = f"{name}_batch_{_batch}"
                        save_dir_list.append(save_dir)
                if not per_batch:
                    graph = Batch.from_data_list(batch_graph_list)
                    graph_list.append(graph)
                else:
                    graph_list.extend(batch_graph_list)

        # train models
        len_graph, len_save_dir = len(graph_list), len(save_dir_list)
        assert len_graph == len_save_dir, f"Graphs: {graph_list} != save_dirs: {len_save_dir}"

        if n_jobs == 1:
            for graph, save_dir in zip(graph_list, save_dir_list):
                train_GAE(graph, cfg, save_dir, disable_gpu)
        else:
            ray.init(ignore_reinit_error=True)
            q = [
                ray_train_GAE.remote(graph, cfg, save_dir, disable_gpu)
                for graph, save_dir in zip(graph_list, save_dir_list)
            ]
            ray.get(q)
            ray.shutdown()
        logger.info(f"Gene selection finished in {time.time() - start_time:.2f}s")
