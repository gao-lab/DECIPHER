r"""
Mixin for omics data
"""

import numpy as np
import pandas as pd
import scanpy as sc
from addict import Dict
from anndata import AnnData
from loguru import logger
from torch_geometric.nn import MLP

from ...utils import estimate_spot_size, sc_


class OmicsMixin:
    r"""
    Mixin class for spatial omics data
    """

    def register_omics(self, meta: pd.DataFrame, config: Dict) -> None:
        self.batched = False
        self.center_encoder = MLP(config.gex_dims, dropout=config.dropout)
        if meta is not None:
            self.adata = AnnData(
                X=np.arange(meta.shape[0], dtype=np.float32).reshape(-1, 1),
                obs=meta.astype(str),
                obsm={"spatial": meta[["x", "y"]].values},
            )
        else:
            self.adata = None

    def save_embedding(
        self,
        x: np.ndarray,
        name: str = "emb",
        resolution: float = 0.3,
        color_by: list[str] = ["leiden", "_celltype", "_batch"],
    ) -> None:
        # Save to numpy
        logger.debug("Save embedding to disk")
        epoch = self.current_epoch
        if not hasattr(self, "version_dir"):
            self.version_dir = (
                self.work_dir
                / self.config.model_dir
                / "lightning_logs"
                / f"version_{self.logger.version}"
            )
        file_name = self.version_dir / f"epoch-{epoch}_{name}_emb.npy"
        np.save(file_name, x)

        # Save to adata
        if not isinstance(self.adata, sc.AnnData):
            return
        logger.debug("Save embedding to adata")
        curr_name = f"X_{name}_{epoch}"
        self.adata.obsm[curr_name] = x

        # Ploting
        if not self.config.plot or self.adata is None:
            return
        logger.debug("Plot umap...")
        sc_.pp.neighbors(self.adata, use_rep=curr_name)
        sc_.tl.leiden(self.adata, resolution=resolution)
        sc_.tl.umap(self.adata)

        self.adata.obsm[f"X_umap_{name}_{epoch}"] = self.adata.obsm["X_umap"]
        self.adata.obs[f"leiden_{name}_{epoch}"] = self.adata.obs["leiden"]
        for var in color_by:
            if var in self.adata.obs.keys():
                umap = sc.pl.umap(self.adata, color=var, return_fig=True)
                self.logger.experiment.add_figure(f"{name}/{var}", umap, epoch)
        if isinstance(self.config.marker_list, list):
            markers = list(set(self.config.marker_list) & set(self.adata.var_names))
            marker_umap = sc.pl.umap(self.adata, color=markers, return_fig=True, vmax=5)
            self.logger.experiment.add_figure(f"{name}/marker", marker_umap, epoch)
        # split dataset by batch to plot spatial
        if "spatial" in self.adata.obsm.keys():
            if "_batch" in self.adata.obs.keys() and self.batched:
                for batch in self.adata.obs["_batch"].unique():
                    batch_idx = self.adata.obs["_batch"] == batch
                    adata = self.adata[batch_idx]
                    spot_size = estimate_spot_size(adata.obsm["spatial"])
                    spatial_fig = sc.pl.spatial(
                        adata,
                        color="leiden",
                        return_fig=True,
                        spot_size=spot_size,
                    )
                    self.logger.experiment.add_figure(
                        f"{name}/sp_{batch}", spatial_fig, epoch
                    )
            else:
                spot_size = estimate_spot_size(self.adata.obsm["spatial"])
                spatial_fig = sc.pl.spatial(
                    self.adata,
                    color="leiden",
                    return_fig=True,
                    spot_size=spot_size,
                )
                self.logger.experiment.add_figure(f"{name}/sp", spatial_fig, epoch)
