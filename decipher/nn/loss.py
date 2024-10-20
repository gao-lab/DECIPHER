r"""
Loss function
"""
import torch
import torch.nn.functional as F
from torch import Tensor, nn


class NTXentLoss(nn.Module):
    r"""
    NT-Xent loss for self-supervised learning in SimCLR.

    Parameters
    ----------
    temperature
        Temperature for the softmax in InfoNCE loss.
    mask_fill
        The value to fill the mask with.

    References
    ----------
    - https://github.com/sthalles/SimCLR/blob/master/simclr.py
    """

    def __init__(self, temperature: float = 0.07, mask_fill: float = -10.0):
        super().__init__()
        self.temperature = temperature
        self.mask_fill = mask_fill
        self.n_views = 2
        self.criterion = nn.CrossEntropyLoss()

    def forward(self, z1: Tensor, z2: Tensor, batch_mask: Tensor | None = None) -> Tensor:
        r"""
        Compute the contrastive loss.

        Parameters
        ----------
        z1
            The output of the first view. (N, D)
        z2
            The output of the second view. (N, D)
        batch_mask
            The mask for the batch, only for spatial loss. (N, N)

        Returns
        ----------
        The contrastive loss.
        """
        batch_size = z1.shape[0]
        labels = torch.cat([torch.arange(batch_size)] * self.n_views, dim=0)
        labels = (labels.unsqueeze(0) == labels.unsqueeze(1)).float()
        labels = labels.to(z1.device)  # (2N, 2N)

        features = torch.cat([z1, z2], dim=0)
        features = F.normalize(features, dim=1)

        similarity_matrix = torch.matmul(features, features.T)  # (2N, 2N)
        # assert similarity_matrix.shape == (
        #     self.n_views * batch_size, self.n_views * batch_size)
        # assert similarity_matrix.shape == labels.shape

        if batch_mask is not None:
            # batch_mask: (batch_size , batch_size) -> (self.n_views * batch_size, self.n_views * batch_size)
            expanded_mask = batch_mask.repeat(self.n_views, self.n_views)
            # fill the mask positions
            similarity_matrix = similarity_matrix.masked_fill(expanded_mask, self.mask_fill)

        # discard the main diagonal from both: labels and similarities matrix
        mask = torch.eye(labels.shape[0], dtype=torch.bool).to(z1.device)  # (2N, 2N)
        labels = labels[~mask].view(labels.shape[0], -1)  # (2N, 2N - 1), rm diagonal
        similarity_matrix = similarity_matrix[~mask].view(
            similarity_matrix.shape[0], -1
        )  # (2N, 2N - 1), remove the diagonal
        # assert similarity_matrix.shape == labels.shape

        # select and combine multiple positives, (2N, 1)
        positives = similarity_matrix[labels.bool()].view(labels.shape[0], -1)

        # select only the negatives the negatives, (2N, 2N - 2), remove positive pairs
        negatives = similarity_matrix[~labels.bool()].view(similarity_matrix.shape[0], -1)

        logits = torch.cat([positives, negatives], dim=1)
        labels = torch.zeros(logits.shape[0], dtype=torch.long).to(z1.device)

        logits = logits / self.temperature

        return self.criterion(logits, labels)


if __name__ == "__main__":
    loss = NTXentLoss()
    z1 = torch.randn(128, 16)
    z2 = torch.randn(128, 16)
    batch_mask = None
    batch_mask = torch.randn(128, 128) > 0
    # set diagonal to True
    batch_mask = batch_mask.fill_diagonal_(True)
    output = loss(z1, z2, batch_mask)
    print(output)
