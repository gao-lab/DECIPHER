from torch import Tensor, nn


class AttentionPooling(nn.Module):
    r"""
    Attention pooling network with sigmoid gating (3 fc layers)

    Parameters
    ----------
    L (int)
        input feature dimension
    D (int)
        hidden layer dimension
    dropout (bool)
        whether to apply dropout (p = 0.25)
    n_classes (int)
        number of classes
    """

    def __init__(
        self, L: int = 1024, D: int = 256, dropout: float = 0.0, n_classes: int = 1
    ) -> None:
        super().__init__()
        self.attention_a = [nn.Linear(L, D), nn.Tanh()]
        self.attention_b = [nn.Linear(L, D), nn.Sigmoid()]

        if dropout > 0.0:
            self.attention_a.append(nn.Dropout(dropout))
            self.attention_b.append(nn.Dropout(dropout))

        self.attention_a = nn.Sequential(*self.attention_a)
        self.attention_b = nn.Sequential(*self.attention_b)
        self.attention_c = nn.Linear(D, n_classes)

    def forward(self, x: Tensor) -> Tensor:
        a = self.attention_a(x)
        b = self.attention_b(x)
        A = a.mul(b)
        A = self.attention_c(A)  # N x n_classes
        return A
