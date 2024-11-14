r"""
Data augment
"""
import numpy as np
import torch
import torch.nn.functional as F
from addict import Dict
from torch import Tensor
from torch_geometric.data import Batch, Data
from torch_geometric.utils import unbatch


class SpatialAugmentMixin:
    r"""
    Mixin class for spatial data augment
    """

    def dropout_nodes(self, x: Tensor) -> Tensor:
        r"""
        Randomly drop nodes

        Parameters
        ----------
        x
            input tensor
        """
        if 0 < self.config.dropout_nbr_prob < 1:
            drop_prob = np.random.uniform(0, 1, x.shape[0])
            drop_mask = drop_prob > self.config.dropout_nbr_prob
            # avoid mask too much nodes
            if drop_mask.sum() < 0.5 * x.shape[0]:
                drop_mask = drop_prob < self.config.dropout_nbr_prob
            drop_mask[0] = True  # avoid dropping center node
            x = x[drop_mask]
        return x

    def pad_nbr_size(self, x: Tensor, max_neighbor: int) -> Tensor:
        r"""
        Pad / clip neighbors to a fixed size

        Parameters
        ----------
        x
            input tensor
        max_neighbor:
            max neighbor length for padding
        """
        if x.shape[0] < max_neighbor:  # padding to max neighbor size
            pad = torch.zeros(max_neighbor - x.shape[0], *x.shape[1:], device=x.device)
            x = torch.cat([x, pad])
        elif x.shape[0] > max_neighbor:  # select first max_neighbor
            x = x[:max_neighbor]
        return x  # (max_neighbor, *feature_dims)


class OmicsSpatialAugment(SpatialAugmentMixin):
    r"""
    Spatial omics data augment

    Parameters
    ----------
    dim
        input dimension
    config
        spatial data augmentation configuration

    Returns
    ----------
    List of augmented sets

    Note
    ----------
    Input is batched graph, we need to augment each graph separately
    """

    def __init__(self, config: Dict) -> None:
        self.config = config

    def __call__(self, graph: Data | Batch, train: bool = True) -> tuple:
        x, batch = graph.x, graph.batch

        if train and hasattr(graph, "sc_batch"):
            batch_mask = self.get_mask_batch(graph.sc_batch[: graph.batch_size])
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
                x1, x2 = self.feature_augment(x)
                # Neighbor augmentation
                x1 = self.dropout_nodes(x1)
                x1 = self.pad_nbr_size(x1, self.config.max_neighbor)
                x2 = self.dropout_nodes(x2)
                x2 = self.pad_nbr_size(x2, self.config.max_neighbor)
                xs_aug1.append(x1)
                xs_aug2.append(x2)
            else:
                x = self.pad_nbr_size(x, self.config.max_neighbor)
                xs_aug1.append(x)

        if train:
            assert len(xs_aug1) == len(xs_aug2)
            return torch.stack(xs_aug1), torch.stack(xs_aug2), batch_mask
        else:
            node_idx = graph.n_id[: graph.batch_size]  # center node index
            return torch.stack(xs_aug1), node_idx

    def feature_augment(self, x: Tensor) -> tuple[Tensor, Tensor]:
        r"""
        Feature augmentation

        Parameters
        ----------
        x
            input tensor
        """
        x1, x2 = x.clone().detach(), x.clone().detach()
        if self.config.dropout_gex > 0:
            # FIXME: test dropout on latent space
            x1 = F.dropout(x1, p=self.config.dropout_gex, training=True)
            x2 = F.dropout(x2, p=self.config.dropout_gex, training=True)
        return x1, x2

    def get_mask_batch(self, batch: Tensor) -> Tensor:
        r"""
        Mask cells from different data batches

        Parameters
        ----------
        batch
            batch index
        """
        mask = batch.unsqueeze(0) != batch.unsqueeze(1)
        return mask


class ScAugment:
    r"""
    Single cell data augmentation

    Parameters
    ----------
    config
        single cell augmentation configuration
    """

    def __init__(self, config: Dict) -> None:
        self.config = config

    def __call__(self, x: Tensor, train: bool = True) -> tuple:
        if train:
            x1, x2 = x.clone().detach(), x.clone().detach()
            x1 = F.dropout(x1, p=self.config.dropout_gex, training=True)
            x2 = F.dropout(x2, p=self.config.dropout_gex, training=True)
            return x1, x2
        else:
            return x, None
