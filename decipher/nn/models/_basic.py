r"""
Basic classes for embedding
"""
from abc import ABC

import numpy as np
import torch
import torch.nn.functional as F
from addict import Dict
from loguru import logger
from pytorch_lightning import LightningModule
from pytorch_lightning.strategies import DDPStrategy
from torch import Tensor, nn
from torch_geometric.nn import MLP

from ..layers.attn_pooling_layer import AttentionPooling
from ..layers.transformer_layer import ViT1D
from ..scheduler import CosineAnnealingWarmupRestarts


class EmbeddingModel(LightningModule, ABC):
    r"""
    Basic class embedding
    """

    def __init__(self, config: Dict) -> None:
        super().__init__()
        self.config = config
        self.batched = False
        self.center_encoder = MLP(config.gex_dims, dropout=config.dropout)
        logger.debug(config.to_dict())
        self.save_hyperparameters(config.to_dict())

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
            lr=self.config.lr_base,
            weight_decay=self.config.weight_decay,
        )
        scheduler = CosineAnnealingWarmupRestarts(
            optimizer=optimizer,
            first_cycle_steps=self.config.first_cycle_steps,
            warmup_steps=self.config.warmup_steps,
            max_lr=self.config.lr_base,
            min_lr=self.config.lr_min,
            gamma=0.25,
        )
        return [optimizer], [
            {
                "scheduler": scheduler,
                "interval": "step",
                "frequency": 1,
                "reduce_on_plateau": False,
                "monitor": "train/total_loss",
            }
        ]

    def on_exception(self, exception: Exception) -> None:
        logger.error(f"Training failed: {exception}")

    # def on_after_backward(self):  # only for check unused parameters in DDP
    #     for name, param in self.named_parameters():
    #         if param.grad is None:
    #             print(name)

    def gather_output(self) -> tuple[np.ndarray, np.ndarray]:
        z_center = torch.vstack(self.val_z_center_list)
        z_nbr = torch.vstack(self.val_z_nbr_list) if len(self.val_z_nbr_list) else None
        z_order = torch.hstack(self.val_z_order_list)
        # if using DDP
        if isinstance(self.trainer.strategy, DDPStrategy):
            center_dim = z_center.size(-1)
            z_center = self.all_gather(z_center).reshape(-1, center_dim).detach().cpu().numpy()
            if z_nbr is not None:
                nbr_dim = z_nbr.size(-1)
                z_nbr = self.all_gather(z_nbr).reshape(-1, nbr_dim).detach().cpu().numpy()
            z_order = self.all_gather(z_order).flatten().detach().cpu().numpy()
            # NOTE: pytorch lightning automaticly padding latest epoch in DDP
            _, unique_idx = np.unique(z_order, return_index=True)
            if not self.trainer.sanity_checking:
                z_center = z_center[unique_idx]
                z_nbr = z_nbr[unique_idx] if z_nbr is not None else None
        else:
            z_center = z_center.detach().cpu().numpy()
            z_nbr = z_nbr.detach().cpu().numpy() if z_nbr is not None else None
            # z_order = z_order.detach().cpu().numpy()
        self.val_z_center_list = []
        self.val_z_nbr_list = []
        self.val_z_order_list = []
        return z_center, z_nbr

    def on_test_start(self):
        self.val_z_center_list = []
        self.val_z_nbr_list = []
        self.val_z_order_list = []  # order for multi-gpu inference


class NeighborEmbeddingModel(EmbeddingModel, ABC):
    r"""
    Transformer-based neighbor embedding model
    """

    def __init__(self, config: Dict) -> None:
        super().__init__(config)

        self.nbr_encoder = ViT1D(
            config.emb_dim, config.transformer_layers, config.num_heads, config.dropout
        )
        self.projection_head = MLP(config.prj_dims, dropout=config.dropout)
        if config.spatial_emb == "attn":
            self.attn_pool = AttentionPooling(config.emb_dim, config.emb_dim // 2, config.dropout)

        self.augment = None  # need to be set in subclass
        self.center_encoder = None  # need to be set in subclas

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
