r"""
Trainer
"""
import time
from pathlib import Path

import torch
from addict import Dict
from loguru import logger
from pytorch_lightning import LightningDataModule, LightningModule, Trainer
from pytorch_lightning.callbacks import (
    Callback,
    EarlyStopping,
    LearningRateMonitor,
    ModelCheckpoint,
    TQDMProgressBar,
)
from pytorch_lightning.loggers import TensorBoardLogger
from torch.utils.data import DataLoader

from ..utils import select_free_gpu


def init_trainer(config: Dict, callback: list[Callback] = None, plugins=None) -> Trainer:
    r"""
    Init pytorch lightning trainer
    """
    # save dir
    model_dir = Path(config.work_dir) / config.model_dir
    # device
    if torch.cuda.is_available() and config.device in ["auto", "gpu", "cuda"]:
        accelerator = "gpu"
        if config.device_num == 1 and config.select_gpu:
            devices = select_free_gpu(config.device_num)
        else:  # for DDP and rya, it will raise error if specify DDP gpu ids
            devices = config.device_num
    else:
        accelerator = "cpu"
        devices = "auto"
    # callbacks
    progress_bar = TQDMProgressBar(refresh_rate=5)
    early_stop = EarlyStopping(monitor="train/total_loss", patience=config.patient, verbose=True)
    model_checkpoint = ModelCheckpoint(
        dirpath=model_dir,
        monitor="train/total_loss",
        mode="min",
        save_top_k=2,
    )
    lr_monitor = LearningRateMonitor(logging_interval="step")
    callback_list = [progress_bar, early_stop, model_checkpoint, lr_monitor]
    if callback is not None:
        callback_list.extend(callback)

    # precision
    if config.fp16 and torch.cuda.is_available():
        precision = "16-mixed"
    else:
        precision = None  # auto select
    logger.debug(f"Run model with precision {precision}")
    # strategy
    strategy = "ddp" if accelerator == "gpu" and config.device_num > 1 else "auto"
    # logger
    pl_logger = TensorBoardLogger(save_dir=model_dir, default_hp_metric=False)

    # init trainer
    return Trainer(
        accelerator=accelerator,
        devices=devices,
        callbacks=callback_list,
        precision=precision,
        default_root_dir=model_dir,
        strategy=strategy,
        gradient_clip_val=config.gradient_clip_val,
        logger=pl_logger,
        log_every_n_steps=config.log_every_n_steps,
        max_epochs=config.epochs,
        check_val_every_n_epoch=config.check_val_every_n_epoch,
        plugins=plugins,
        max_steps=config.max_steps,
    )


def fit(
    model: LightningModule,
    data: DataLoader | list[DataLoader] | LightningDataModule,
    config: Dict = None,
    return_model_path: bool = False,
    show_name: str = "",
    callback: list[Callback] | None = None,
    plugins=None,
    trainer: Trainer = None,
) -> str | None:
    r"""
    Model training
    """
    trainer = init_trainer(config, callback, plugins) if trainer is None else trainer
    start_time = time.time()
    logger.info(f"Start training {show_name} model.")

    if isinstance(data, LightningDataModule):
        trainer.fit(model, datamodule=data)
    elif isinstance(data, list):
        if len(data) == 3:
            trainer.fit(model, train_dataloaders=data[0], val_dataloaders=data[1])
        else:
            trainer.fit(model, train_dataloaders=data[0])
    elif isinstance(data, DataLoader):
        trainer.fit(model, train_dataloaders=data)
    else:
        raise TypeError(f"Data type {type(data)} is not supported yet.")

    logger.success(f"Train finished in {time.time() - start_time:.2f}s.")
    if return_model_path:
        return trainer.checkpoint_callback.best_model_path


def inference(
    model: LightningModule,
    data: DataLoader | LightningDataModule,
    config: Dict = None,
    trainer: Trainer = None,
) -> None:
    r"""
    Model inference only
    """
    trainer = init_trainer(config) if trainer is None else trainer
    if isinstance(data, LightningDataModule):
        trainer.test(model, datamodule=data)
    elif isinstance(data, DataLoader):
        trainer.test(model, dataloaders=data)
    elif isinstance(data, list):
        assert len(data) > 1
        trainer.test(model, dataloaders=data[-1])
    else:
        raise TypeError(f"Type {type(data)} is not supported yet.")


def fit_and_inference(
    model: LightningModule,
    data: list[DataLoader] | LightningDataModule,
    config: Dict,
    show_name: str = "",
    callback: list[Callback] | None = None,
    plugins=None,
) -> None:
    r"""
    Model training and inference
    """
    trainer = init_trainer(config, callback, plugins)
    fit(model, data, trainer=trainer, show_name=show_name)
    inference(model, data, trainer=trainer)
