import numpy as np
import scanpy as sc

from spider.utils import scanpy_viz


def test_scanpy_viz(adata):
    adata = scanpy_viz(adata)
    print(adata)


if __name__ == "__main__":
    X = np.random.rand(100, 1000)
    X_center = np.random.rand(100, 20)
    X_nbr = np.random.rand(100, 20)
    adata = sc.AnnData(
        X,
        obsm={"X_center": X_center, "X_nbr": X_nbr},
    )
    test_scanpy_viz(adata)
