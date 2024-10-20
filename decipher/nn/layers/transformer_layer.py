import torch
from einops import repeat
from torch import Tensor, nn
from torch_geometric.typing import OptTensor


class TransformerLayer(nn.Module):
    r"""
    Vanilla transformer layer for ViT

    Parameters
    ----------
    embed_dim:
        embedding dimension
    num_heads:
        number of attention heads
    dropout:
        dropout rate
    """

    def __init__(self, embed_dim: int, num_heads: int, dropout: float = 0.0):
        super().__init__()
        self.attn = nn.MultiheadAttention(
            embed_dim=embed_dim, num_heads=num_heads, batch_first=True
        )
        self.feed_forward = nn.Sequential(
            nn.Linear(embed_dim, embed_dim),
            nn.ReLU(),
            nn.Linear(embed_dim, embed_dim),
        )
        self.norm1 = nn.LayerNorm(embed_dim)
        self.norm2 = nn.LayerNorm(embed_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(
        self,
        x: Tensor,
        key_padding_mask: OptTensor,
        attn_mask: OptTensor,
    ) -> Tensor:
        attn_out = self.attn(x, x, x, key_padding_mask=key_padding_mask, attn_mask=attn_mask)[0]
        x = x + self.dropout(attn_out)
        x = self.norm1(x)
        ff_out = self.feed_forward(x)
        x = x + self.dropout(ff_out)
        x = self.norm2(x)
        return x


class Transformer(nn.Module):
    r"""
    Vanilla transformer for ViT

    Parameters
    ----------
    embed_dim:
        embedding dimension
    num_layers:
        number of transformer layers
    num_heads:
        number of attention heads
    dropout:
        dropout rate
    """

    def __init__(self, embed_dim: int, num_layers: int, num_heads: int = 1, dropout: float = 0.0):
        super().__init__()
        self.layers = nn.ModuleList(
            [TransformerLayer(embed_dim, num_heads, dropout) for _ in range(num_layers)]
        )

    def forward(
        self,
        x: Tensor,
        key_padding_mask: OptTensor = None,
        attn_mask: OptTensor = None,
    ) -> Tensor:
        for layer in self.layers:
            x = layer(x, key_padding_mask, attn_mask)
        return x


class ViT1D(nn.Module):
    r"""
    Vision Transformer for 1D data

    Parameters
    ----------
    embed_dim:
        embedding dimension
    num_layers:
        number of transformer layers
    num_heads:
        number of attention heads
    dropout:
        dropout rate

    References
    ----------
    https://github.com/lucidrains/vit-pytorch/blob/main/vit_pytorch/vit_1d.py
    """

    def __init__(
        self,
        embed_dim: int,
        num_layers: int,
        num_heads: int = 1,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.cls_token = nn.Parameter(torch.randn(1, 1, embed_dim))
        self.dropout = nn.Dropout(dropout) if dropout > 0.0 else nn.Identity()
        self.pre_patch_norm = nn.LayerNorm(embed_dim)
        self.transformer = Transformer(embed_dim, num_layers, num_heads, dropout)

    def forward(
        self,
        x: Tensor,
        key_padding_mask: OptTensor = None,
        attn_mask: OptTensor = None,
    ) -> tuple[Tensor, Tensor]:
        b = x.shape[0]  # batch size

        cls_tokens = repeat(self.cls_token, "1 1 d -> b 1 d", b=b)
        # add a False row on key_padding_mask to avoid masking cls token
        key_padding_mask = torch.cat(
            (torch.zeros(b, 1, dtype=bool, device=x.device), key_padding_mask), dim=1
        )
        # TODO: test to add a linear layer before transformer
        x = self.pre_patch_norm(x)
        x = torch.cat((cls_tokens, x), dim=1)

        x = self.dropout(x)
        x = self.transformer(x, key_padding_mask, attn_mask)
        cls_tokens, emb = x[:, 0], x[:, 1:]
        return cls_tokens, emb
