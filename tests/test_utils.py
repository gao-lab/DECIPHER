import scanpy as sc

from decipher.utils import scanpy_viz


def test_scanpy_viz(adata: sc.AnnData):
    adata = scanpy_viz(adata)
    print(adata)


if __name__ == "__main__":
    adata = sc.datasets.pbmc3k()
    test_scanpy_viz(adata)
