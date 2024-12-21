r"""
Spatial location related gene selection model
"""
from pathlib import Path
from time import time

import numpy as np
import pandas as pd
import ray
import torch
import torch_geometric.transforms as T
from addict import Dict
from loguru import logger
from torch import Tensor, nn
from torch_geometric.data import Batch, Data
from torch_geometric.nn import GAE, MLP

from ...utils import select_free_gpu


class SelectEncoder(nn.Module):
    r"""
    Omics encoder via gene selective
    """

    def __init__(self, config: Dict) -> None:
        super().__init__()
        self.config = config
        self.emb2mask = MLP([config.center_dim, config.expr_dim], bias=False)

    def forward(
        self, emb: Tensor, expr: Tensor, return_mask: bool = False, hard: bool = True
    ) -> tuple:
        gene_mask = self.emb2mask(emb)
        gene_mask = gumbel_sigmoid(
            gene_mask, tau=1, hard=hard, threshold=self.config.gumbel_threshold
        )
        if return_mask:
            return gene_mask
        l1_loss = gene_mask.mean()
        return expr * gene_mask, l1_loss


@ray.remote(num_gpus=1, num_cpus=8)
def ray_train_GAE(graph: Data, config: Dict, save_dir: Path, disable_gpu: bool = False):
    config["select_gpu"] = False
    return train_GAE(graph, config, save_dir, disable_gpu)


def train_GAE(
    graph: Data | Batch, config: Dict, save_dir: Path, disable_gpu: bool = False, fp16: bool = True
) -> None:
    r"""
    Train graph autoencoder with whole graph

    Parameters
    ----------
    graph
        k-nn graph of spatial embedding
    config
        config for train GAE
    save_dir
        save directory
    disable_gpu
        if disable gpu
    fp16
        if use fp16
    """
    # Step1. set up model and optimizer
    encoder = SelectEncoder(config)
    model = GAE(encoder)

    optimizer = torch.optim.Adam(model.parameters(), lr=config.lr_base)

    if torch.cuda.is_available() and not disable_gpu:
        if config.select_gpu:
            cuda_idx = select_free_gpu()[0]
            device = f"cuda:{cuda_idx}"
        else:
            device = "cuda"
        model = model.to(device)
    else:
        device = "cpu"

    # Step2. split edge into train and test (not node)
    transform = T.Compose(
        [
            T.ToDevice(device),
            T.RandomLinkSplit(
                num_val=0.0,
                num_test=config.test_ratio,
                is_undirected=True,
                split_labels=True,
                add_negative_train_samples=False,
            ),
        ]
    )
    if fp16 and device != "cpu":
        logger.debug("Train GAE with fp16")
        model = model.bfloat16()
        graph.x, graph.expr = graph.x.bfloat16(), graph.expr.bfloat16()

    if isinstance(graph, Batch):
        graph_list = Batch.to_data_list(graph)
        train_data, test_edge_pos, test_edge_neg = [], [], []
        for g in graph_list:
            train, _, test = transform(g)
            train_data.append(train)
            test_edge_pos.append(test.pos_edge_label_index)
            test_edge_neg.append(test.neg_edge_label_index)
        del graph_list
    else:
        train_data, _, test_data = transform(graph)
        test_edge_pos = test_data.pos_edge_label_index
        test_edge_neg = test_data.neg_edge_label_index
        del test_data
    del graph
    torch.cuda.empty_cache()

    # Step3. model training
    start_time = time()
    log_list = []

    for i in range(config.gae_epochs):
        # train
        model.train()
        optimizer.zero_grad()
        if isinstance(train_data, list):
            loss = 0.0
            for train in train_data:
                z, l1_loss = model.encode(train.x, train.expr)
                recon_loss = model.recon_loss(z, train.pos_edge_label_index)
                loss = loss + recon_loss + l1_loss * config.l1_weight
            loss = loss / len(train_data)
        else:
            z, l1_loss = model.encode(train_data.x, train_data.expr)
            recon_loss = model.recon_loss(z, train_data.pos_edge_label_index)
            loss = recon_loss + l1_loss * config.l1_weight
        loss.backward()
        optimizer.step()

        logger.info(f"Epoch: {i}: Loss: {loss:.4f}, Recon: {recon_loss:.4f}, L1: {l1_loss:.4f}")

        with torch.inference_mode():
            model.eval()
            if isinstance(train_data, list):
                auc, ap = 0.0, 0.0
                for train, pos_edge, neg_edge in zip(train_data, test_edge_pos, test_edge_neg):
                    z, _ = model.encode(train.x, train.expr)
                    auc_, ap_ = model.test(z.float(), pos_edge, neg_edge)
                    auc += auc_
                    ap += ap_
                auc /= len(train_data)
                ap /= len(train_data)
            else:
                z, _ = model.encode(train_data.x, train_data.expr)
                auc, ap = model.test(z.float(), test_edge_pos, test_edge_neg)
            logger.info(f"Test: AUC: {auc:.4f}, AP: {ap:.4f}\n")
        log_list.append(
            dict(
                loss=loss.item(),
                recon_loss=recon_loss.item(),
                l1_loss=l1_loss.item(),
                auc=auc,
                ap=ap,
            )
        )
    df = pd.DataFrame(log_list)
    logger.success(f"Training finished in {time() - start_time:.2f}s")

    # inference on whole dataset
    with torch.inference_mode():
        if isinstance(train_data, list):
            gene_mask = torch.cat(
                [model.encode(g.x, g.expr, return_mask=True) for g in train_data], dim=0
            )
        else:
            gene_mask = model.encode(train_data.x, train_data.expr, return_mask=True)

    # save
    save_dir = Path(config.work_dir) / save_dir
    logger.info(f"Save explain model to {save_dir}")
    save_dir.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), save_dir / "explain_model.pth")
    np.save(save_dir / "gene_mask.npy", gene_mask.float().detach().cpu().numpy())
    df.to_csv(save_dir / "train_log.csv", index=False)


def gumbel_sigmoid(
    logits: Tensor, tau: float = 1, hard: bool = False, threshold: float = 0.5
) -> Tensor:
    r"""
    Samples from the Gumbel-Sigmoid distribution and optionally discretizes.

    The discretization converts the values greater than `threshold` to 1 and the rest to 0.
    The code is adapted from the official PyTorch implementation of gumbel_softmax:
    https://pytorch.org/docs/stable/_modules/torch/nn/functional.html#gumbel_softmax

    Parameters
    ----------
    logits
        `[..., num_features]` unnormalized log probabilities
    tau
        non-negative scalar temperature
    hard
        if ``True``, the returned samples will be discretized,
        but will be differentiated as if it is the soft sample in autograd
    threshold
        threshold for the discretization, values greater than this will be set to 1 and the rest to 0

    Returns
    ----------
    Sampled tensor of same shape as `logits` from the Gumbel-Sigmoid distribution.
    If ``hard=True``, the returned samples are descretized according to `threshold`, otherwise they will
    be probability distributions.

    References
    ----------
    https://github.com/AngelosNal/PyTorch-Gumbel-Sigmoid/blob/main/gumbel_sigmoid.py
    """
    gumbels = (
        -torch.empty_like(logits, memory_format=torch.legacy_contiguous_format).exponential_().log()
    )  # ~Gumbel(0, 1)
    gumbels = (logits + gumbels) / tau  # ~Gumbel(logits, tau)
    y_soft = gumbels.sigmoid()

    if hard:
        # Straight through.
        indices = (y_soft > threshold).nonzero(as_tuple=True)
        y_hard = torch.zeros_like(logits, memory_format=torch.legacy_contiguous_format)
        y_hard[indices[0], indices[1]] = 1.0
        ret = y_hard - y_soft.detach() + y_soft
    else:
        # Reparametrization trick.
        ret = y_soft
    return ret
