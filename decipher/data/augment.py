r"""
Data augment
"""
import numpy as np
import torch
import torch.nn.functional as F
from torch import Tensor
from torch_geometric.data import Batch, Data
from torch_geometric.utils import unbatch

from .. import CFG


def feature_augment(x: Tensor, dropout: float) -> tuple[Tensor, Tensor]:
    r"""
    Feature augmentation

    Args:
        x: input feature tensor (batch_size, num_features)
        dropout: dropout probability
    """
    x1, x2 = x.clone().detach(), x.clone().detach()
    if dropout > 0:
        x1 = F.dropout(x1, p=dropout, training=True)
        x2 = F.dropout(x2, p=dropout, training=True)
    return x1, x2


def dropout_nodes(x: Tensor, dropout_nbr_prob: float) -> Tensor:
    r"""
    Randomly drop nodes

    Args:
        x: input feature tensor (num_nodes, num_features)
        dropout_nbr_prob: dropout probability for nodes
    """
    if 0 < dropout_nbr_prob < 1:
        drop_prob = np.random.uniform(0, 1, x.shape[0])
        drop_mask = drop_prob > dropout_nbr_prob
        # avoid mask too much nodes
        if drop_mask.sum() < 0.5 * x.shape[0]:
            drop_mask = drop_prob < dropout_nbr_prob
        drop_mask[0] = True  # avoid dropping center node
        x = x[drop_mask]
    return x


def pad_nbr_size(x: Tensor, max_neighbor: int) -> Tensor:
    r"""
    Pad / clip neighbors to a fixed size

    Args:
        x: input feature tensor (num_nodes, num_features)
        max_neighbor: maximum neighbor size
    Returns:
        Padded / clipped tensor in shape (max_neighbor, num_features)
    """
    if x.shape[0] < max_neighbor:  # padding to max neighbor size
        pad = torch.zeros(max_neighbor - x.shape[0], *x.shape[1:], device=x.device)
        x = torch.cat([x, pad])
    elif x.shape[0] > max_neighbor:  # select first max_neighbor
        x = x[:max_neighbor]
    return x  # (max_neighbor, *feature_dims)


class OmicsSpatialAugment:
    r"""
    Spatial omics data augment

    Args:
        config: model configuration
    Note:
        Input is batched sub-graph, augment each sub-graph separately
    """

    def __call__(self, graph: Data | Batch, train: bool = True) -> tuple:
        x, batch = graph.x, graph.batch

        if train and hasattr(graph, "sc_batch"):
            batch_id = graph.sc_batch[: graph.batch_size]
            batch_mask = batch_id.unsqueeze(0) != batch_id.unsqueeze(1)
        else:
            batch_mask = None

        # sort the batch and apply the sort to the x (need stable sort)
        sorted_idx = torch.argsort(batch, stable=True)
        batch = batch[sorted_idx]
        x = x[sorted_idx]
        xs = unbatch(x, batch=batch)

        xs_aug1, xs_aug2 = [], []
        # x means nodes in subgraph
        for x in xs:
            if train:
                # Feature augmentation
                x1, x2 = feature_augment(x, CFG.augment.dropout_gex)
                # Neighbor augmentation
                x1 = dropout_nodes(x1, CFG.augment.dropout_nbr_prob)
                x1 = pad_nbr_size(x1, CFG.augment.max_neighbor)
                x2 = dropout_nodes(x2, CFG.augment.dropout_nbr_prob)
                x2 = pad_nbr_size(x2, CFG.augment.max_neighbor)
                xs_aug1.append(x1)
                xs_aug2.append(x2)
            else:
                x = pad_nbr_size(x, CFG.augment.max_neighbor)
                xs_aug1.append(x)

        if train:
            assert len(xs_aug1) == len(xs_aug2)
            return torch.stack(xs_aug1), torch.stack(xs_aug2), batch_mask
        else:
            return torch.stack(xs_aug1)


class ScAugment:
    r"""
    Single cell data augmentation

    Args:
        config: model configuration
    """

    def __call__(self, x: Tensor, train: bool = True) -> tuple:
        if train:
            x1, x2 = feature_augment(x, CFG.augment.dropout_gex)
            return x1, x2
        else:
            return x, None
