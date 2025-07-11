r"""
Single cell contrastive learning model
"""
from addict import Dict
from torch import Tensor

from ...data.augment import ScAugment
from ..loss import NTXentLoss
from ._basic import EmbeddingModel


class ScSimCLR(EmbeddingModel):
    r"""
    Single cell embedding

    Parameters
    ----------
    config:
        model configuration
    """

    def __init__(self, config: Dict) -> None:
        super().__init__(config)
        self.augment = ScAugment(config.augment)
        self.criterion = NTXentLoss(config.temperature_center)
        self._reset_prams()

    def training_step(self, batch: list[Tensor], batch_idx: int) -> Tensor:
        x1, x2 = self.augment(batch[0])
        z1 = self.center_encoder(x1)
        z2 = self.center_encoder(x2)
        loss = self.criterion(z1, z2)
        self.log("train/total_loss", loss, prog_bar=True)
        return loss

    def test_step(self, data: list[Tensor], batch_idx: int) -> None:
        x, order = data
        z = self.center_encoder(x)
        self.val_z_center_list.append(z)
        self.val_z_order_list.append(order)


class ScSimCLRMNN(ScSimCLR):
    r"""
    Single cell MNN-based embedding for batch correction

    Parameters
    ----------
    config:
        model configuration
    """

    def __init__(self, config: Dict) -> None:
        super().__init__(config)

    def training_step(self, batch: list[Tensor] | dict, batch_idx: int) -> Tensor:
        if isinstance(batch, dict):  # with extra data
            x, mnn = batch["x"], batch["mnn"]
            contrast_loss = super().training_step(x, batch_idx)
        elif isinstance(batch, list):  # without extra data
            mnn = batch
            contrast_loss1 = super().training_step(mnn[0], batch_idx)
            contrast_loss2 = super().training_step(mnn[1], batch_idx)
            contrast_loss = (contrast_loss1 + contrast_loss2) * 0.5

        mnn_loss = self.mnn_step(mnn)
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

    def mnn_step(self, mnn: tuple[Tensor]) -> Tensor:
        r"""
        MNN forward step for batch correction

        Parameters
        ----------
        mnn:
            MNN pairs (x1, x2)
        """
        x1, x2 = mnn
        z1 = self.center_encoder(x1)
        z2 = self.center_encoder(x2)
        loss = self.criterion(z1, z2)
        return loss
