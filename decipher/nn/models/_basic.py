r"""
Basic classes for embedding
"""
from abc import ABC

import numpy as np
import torch
import torch.nn.functional as F
from loguru import logger
from omegaconf import OmegaConf
from pytorch_lightning import LightningModule
from torch import Tensor, nn
from torch_geometric.nn import MLP
from transformers.optimization import get_cosine_schedule_with_warmup

from ..layers import AttentionPooling, ViT1D


class EmbeddingModel(LightningModule, ABC):
    r"""
    Basic class embedding
    """

    def __init__(self, config: OmegaConf) -> None:
        super().__init__()
        self.config = config
        logger.debug(config)
        self.save_hyperparameters(config)
        self.z_center_list = []
        self.z_nbr_list = []

    def _reset_prams(self) -> None:
        r"""
        xavier init (unavailable for lazy init), skip frozen parameters
        """
        for p in self.parameters():
            if not isinstance(p, nn.UninitializedParameter):
                if p.dim() > 1 and p.requires_grad:
                    nn.init.xavier_uniform_(p)

    def configure_optimizers(self):
        r"""
        pl optimizer and scheduler
        """
        optimizer = torch.optim.AdamW(
            params=filter(lambda p: p.requires_grad, self.parameters()),
            lr=self.config.lr,
            weight_decay=self.config.weight_decay,
        )
        if not self.config.lr_scheduler:
            return optimizer
        scheduler = get_cosine_schedule_with_warmup(
            optimizer=optimizer,
            num_warmup_steps=self.config.warmup_steps,
            num_training_steps=self.config.num_training_steps,
        )
        return [optimizer], [
            {
                "scheduler": scheduler,
                "interval": "step",
                "frequency": 1,
                "monitor": "train/total_loss",
            }
        ]

    def on_exception(self, exception: Exception) -> None:
        logger.error(f"Training failed: {exception}")

    def gather_output(self) -> tuple[np.ndarray, np.ndarray] | np.ndarray:
        z_center = np.vstack(self.z_center_list)
        if len(self.z_nbr_list) > 0:
            z_nbr = np.vstack(self.z_nbr_list)
        self.z_center_list = []
        self.z_nbr_list = []
        if len(z_nbr) > 0:
            return z_center, z_nbr
        else:
            return z_center


class NeighborEmbeddingModel(EmbeddingModel, ABC):
    r"""
    Transformer-based neighbor embedding model
    """

    def __init__(self, config: OmegaConf) -> None:
        super().__init__(config)
        self.center_encoder = MLP(list(config.gex_dims), dropout=config.dropout)
        self.nbr_encoder = ViT1D(
            config.emb_dim, config.transformer_layers, config.num_heads, config.dropout
        )
        self.projection_head = MLP(list(config.prj_dims), dropout=config.dropout)
        if config.spatial_emb == "attn":
            self.attn_pool = AttentionPooling(config.emb_dim, config.emb_dim // 2, config.dropout)

        self.augment = None  # need to be set in subclass
        self.center_encoder = None  # need to be set in subclass

    def create_attn_mask(self, x: Tensor) -> Tensor:
        r"""
        Create attention mask for padding cells
        """
        # x is Tensor in shape (batch_size, max_len, emb_dim)
        sum_along_feature = torch.sum(x, dim=-1)
        return sum_along_feature == 0  # (batch_size, max_len)

    def center_forward(self, x: Tensor) -> tuple[Tensor, Tensor]:
        r"""
        Center node forward

        Parameters
        ----------
        x
           Tensor in (batch_size, max_neighbor, *feature_dim)
        coords
            spatial coordinates of cells, in shape (batch_size, max_len, 2)
        """
        # flatten 0 and 1 dim only
        batch, max_neighbor = x.size(0), x.size(1)
        x = torch.flatten(x, start_dim=0, end_dim=1)
        z: Tensor = self.center_encoder(x)
        z = z.reshape(batch, max_neighbor, -1)
        z_center = z[:, 0, :].clone()  # center node
        return z_center, z

    def nbr_forward(self, z: Tensor, mask: Tensor) -> Tensor:
        r"""
        Neighbor forward
        """
        cls_emb, cell_emb = self.nbr_encoder(z[:, 1:, :], key_padding_mask=mask[:, 1:])
        if self.config.spatial_emb == "cls":
            return cls_emb
        elif self.config.spatial_emb == "attn":
            attn = self.attn_pool(cell_emb)
            attn = F.softmax(attn, dim=1)
            attn = torch.transpose(attn, 1, 2)
            cell_emb = torch.matmul(attn, cell_emb)
            cell_emb = cell_emb.squeeze(1)
            return cell_emb
        elif self.config.spatial_emb == "mean":
            return cell_emb.mean(dim=1)
        else:
            raise ValueError(f"Unknown spatial embedding type: {self.config.spatial_emb}")
