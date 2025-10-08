r"""
Single cell contrastive learning model
"""
import torch
from omegaconf import OmegaConf
from torch import Tensor
from torch_geometric.nn import MLP

from ...data.augment import ScAugment
from ..loss import NTXentLoss
from ._basic import EmbeddingModel


class ScSimCLR(EmbeddingModel):
    r"""
    SimCLR framework for single cell data

    Args:
        config: model configuration
    """

    def __init__(self, config: OmegaConf) -> None:
        super().__init__(config)
        self.center_encoder = MLP(list(config.gex_dims), dropout=config.dropout)
        self.augment = ScAugment()
        self.criterion = NTXentLoss(config.temperature_center)
        self._reset_prams()

    def forward(self, batch: list[Tensor], batch_idx: int) -> Tensor:
        x1, x2 = self.augment(batch[0])
        z1 = self.center_encoder(x1)
        z2 = self.center_encoder(x2)
        loss = self.criterion(z1, z2)
        self.log("train/total_loss", loss, prog_bar=True)
        return loss

    def test_step(self, data: list[Tensor], batch_idx: int) -> None:
        x = data
        z = self.center_encoder(x)
        self.z_center_list.append(z.detach().to(torch.float32).cpu().numpy())

    def training_step(self, batch: dict | list[Tensor], batch_idx: int) -> Tensor:
        if isinstance(batch, list):
            return self.forward(batch, batch_idx)

        # with batch
        assert isinstance(batch, dict)
        x, mnn = batch["x"], batch["mnn"]
        contrast_loss = self.forward(x, batch_idx)

        mnn_loss = self.get_mnn_loss(mnn)
        loss = contrast_loss + mnn_loss
        self.log_dict(
            {
                "train/total_loss": loss,
                "train/contrast_loss": contrast_loss,
                "train/mnn_loss": mnn_loss,
            },
            prog_bar=True,
        )
        return loss

    def get_mnn_loss(self, mnn: tuple[Tensor]) -> Tensor:
        r"""
        MNN loss for batch correction

        Args:
            mnn: MNN pairs (x1, x2)
        """
        x1, x2 = mnn
        z1 = self.center_encoder(x1)
        z2 = self.center_encoder(x2)
        loss = self.criterion(z1, z2)
        return loss
