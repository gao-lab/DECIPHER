r"""
Spatial omics contrastive learning model
"""
import torch
from omegaconf import OmegaConf
from torch import Tensor
from torch_geometric.data import Data

from ...data.augment import OmicsSpatialAugment
from ..loss import NTXentLoss
from ._basic import NeighborEmbeddingModel


class SpatialSimCLR(NeighborEmbeddingModel):
    r"""
    SimCLR framework for spatial omics data
    """

    def __init__(self, config: OmegaConf) -> None:
        super().__init__(config)
        self.center_criterion = NTXentLoss(config.temperature_center)
        self.nbr_criterion = NTXentLoss(config.temperature_nbr)
        self.augment = OmicsSpatialAugment()
        self._reset_prams()

    def forward(self, x: Tensor) -> tuple[Tensor, Tensor]:
        mask = self.create_attn_mask(x)
        z_center, z = self.center_forward(x)
        z_nbr = self.nbr_forward(z, mask)
        z_nbr = self.projection_head(z_nbr)
        return z_center, z_nbr

    def training_step(self, data: Data | dict, batch_idx: int) -> Tensor:
        if isinstance(data, dict):
            mnn = data["mnn"]
            data = data["graph"]
        else:
            mnn = None

        # x1 and x2: (bz, nbr, *feature_dims)
        # For image data, augment has done in torch dataset, so here is a trick
        x1, x2, batch_mask = self.augment(data)
        z_center1, z_nbr1 = self.forward(x1)
        z_center2, z_nbr2 = self.forward(x2)

        center_loss = self.center_criterion(z_center1, z_center2, batch_mask=batch_mask)
        nbr_loss_weight = self.config.nbr_loss_weight

        nbr_loss = self.nbr_criterion(z_nbr1, z_nbr2, batch_mask=batch_mask)
        loss = nbr_loss * nbr_loss_weight + center_loss * (1 - nbr_loss_weight)

        loss_dict = {"train/center_loss": center_loss, "train/nbr_loss": nbr_loss}
        if mnn is not None:
            mnn_loss = self.get_mnn_loss(mnn)
            loss = loss + mnn_loss
            loss_dict["train/mnn_loss"] = mnn_loss
        loss_dict["train/total_loss"] = loss
        self.log_dict(loss_dict, prog_bar=True)
        return loss

    def test_step(self, data: Data, batch_idx: int) -> None:
        xs_raw = self.augment(data, train=False)
        z_center, z = self.center_forward(xs_raw)
        mask = self.create_attn_mask(xs_raw)
        z_nbr = self.nbr_forward(z, mask)
        z_nbr = self.projection_head(z_nbr)
        self.z_center_list.append(z_center.detach().to(torch.float32).cpu().numpy())
        self.z_nbr_list.append(z_nbr.detach().to(torch.float32).cpu().numpy())

    def get_mnn_loss(self, mnn: tuple[Tensor]) -> Tensor:
        r"""
        MNN loss for batch correction

        Args:
            mnn: MNN pairs (x1, x2)
        """
        x1, x2 = mnn
        z1 = self.center_encoder(x1)
        z2 = self.center_encoder(x2)
        loss = self.center_criterion(z1, z2)
        return loss
