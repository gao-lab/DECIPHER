r"""
Variance explained by regression model
"""
from pathlib import Path

import pandas as pd
import ray
import torch
from loguru import logger
from omegaconf import OmegaConf
from pytorch_lightning import LightningDataModule, LightningModule
from rui_utils.torch.trainer import fit_and_inference
from rui_utils.utils import save_dict
from torch import Tensor, nn
from torch.utils.data import DataLoader, TensorDataset
from torch_geometric.nn import MLP
from torchmetrics.functional.regression import mean_absolute_error, r2_score


class SimpleLightningDataModule(LightningDataModule):
    r"""
    Simple LightningDataModule
    """

    def __init__(
        self, train_dataset: TensorDataset, test_dataset: TensorDataset, batch_size: int = 1024
    ):
        super().__init__()
        self.train_dataset = train_dataset
        self.test_dataset = test_dataset
        self.batch_size = batch_size

    def train_dataloader(self):
        return DataLoader(
            self.train_dataset,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=8,
            pin_memory=True,
        )

    def test_dataloader(self):
        return DataLoader(
            self.test_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=8,
            pin_memory=True,
        )


class DeepRegression(LightningModule):
    r"""
    Deep regression model

    Parameters
    ----------
    config
        config dict
    test_cell_type
        cell type annotation of test set (only used when training model in whole dataset)
    """

    def __init__(self, config: OmegaConf, test_cell_type: pd.Series = None):
        super().__init__()
        self.model = MLP([config.center_dim, config.hidden_dim, config.nbr_dim])
        self.criterion = nn.MSELoss()
        self.test_y_list = []
        self.test_pred_y_list = []
        self.test_cell_type = test_cell_type

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=1e-3)

    def training_step(self, batch, batch_idx):
        x, y = batch
        pred_y = self.model(x)
        loss = self.criterion(pred_y, y)
        self.log("train/total_loss", loss)
        return loss

    def test_step(self, batch, batch_idx):
        x, y = batch
        pred_y = self.model(x)
        self.test_y_list.append(y)
        self.test_pred_y_list.append(pred_y)

    def on_test_epoch_end(self):
        y = torch.cat(self.test_y_list)
        pred_y = torch.cat(self.test_pred_y_list)
        if self.test_cell_type is None:
            self.calc_metric(pred_y, y)
        else:
            self.calc_metric_split(pred_y, y)
        self.test_x_list = []
        self.test_y_list = []

    def calc_metric_split(self, pred_y, y):
        r"""
        Calculate test metric for each cell type
        """
        self.test_metric = {}
        for _cell_type in self.test_cell_type.unique():
            idx = self.test_cell_type == _cell_type
            y_tmp = y[idx.values]
            pred_y_tmp = pred_y[idx.values]
            mse = self.criterion(pred_y_tmp, y_tmp)
            r2 = r2_score(pred_y_tmp, y_tmp)
            mae = mean_absolute_error(pred_y_tmp, y_tmp)
            self.test_metric.update(
                {
                    f"test_{_cell_type}:mse_loss": float(mse),
                    f"test_{_cell_type}:r2_score": float(r2),
                    f"test_{_cell_type}:mae_score": float(mae),
                }
            )

    def calc_metric(self, pred_y, y):
        r"""
        Calculate test metric
        """
        mse = self.criterion(pred_y, y)
        r2 = r2_score(pred_y, y)
        mae = mean_absolute_error(pred_y, y)
        self.test_metric = {
            "test_mse_loss": float(mse),
            "test_r2_score": float(r2),
            "test_mae_score": float(mae),
        }
        try:
            self.log_dict(self.test_metric, prog_bar=True)
        except Exception as e:  # noqa
            logger.error(f"Error logging test metrics: {e}")
            print(self.test_metric)


def train_regress(
    x: Tensor,
    y: Tensor,
    config: OmegaConf,
    save_dir: str,
    cell_type: pd.Series = None,
) -> None:
    r"""
    Train regression model

    Parameters
    ----------
    x
        x data
    y
        y data
    config
        config dict
    save_dir
        save directory (sub dir of work_dir)
    cell_type
        cell type anno in pandas.Series
    """
    test_size = int(config.test_ratio * len(x))
    # val_size = int(config.val_ratio * len(x))
    train_size = len(x) - test_size
    # get the index of train, test
    train_idx = torch.randperm(len(x))[:train_size]
    test_idx = torch.randperm(len(x))[train_size:]
    train_x, test_x = x[train_idx], x[test_idx]
    train_y, test_y = y[train_idx], y[test_idx]
    train_data = TensorDataset(train_x, train_y)
    test_data = TensorDataset(test_x, test_y)

    if "_all" not in save_dir:
        test_cell_type = None
    else:
        test_cell_type = cell_type[test_idx.numpy()]

    data = SimpleLightningDataModule(train_data, test_data, batch_size=1024)

    config.trainer.model_dir = str(Path(config.work_dir) / save_dir)
    regress_model = DeepRegression(config, test_cell_type)
    fit_and_inference(regress_model, data, config.trainer)
    if not hasattr(regress_model, "test_metric"):
        regress_model.on_test_epoch_end()
    test_metric = regress_model.test_metric

    metric_path = Path(config.work_dir) / save_dir / "test_metric.json"
    save_dict(test_metric, metric_path)


@ray.remote(num_gpus=0.25, num_cpus=2)
def ray_train_regress(
    x: Tensor, y: Tensor, config: OmegaConf, save_dir: Path, cell_type: pd.Series
) -> None:
    return train_regress(x, y, config, save_dir, cell_type)
