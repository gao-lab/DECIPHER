import scanpy as sc

from decipher.utils import gex_embedding, scanpy_viz


def test_gex_embedding_scanpy_viz():
    adata = sc.datasets.visium_sge("V1_Breast_Cancer_Block_A_Section_1")
    adata = adata[:200, :2500].copy()
    keys = ["gex"]

    adata = gex_embedding(adata)
    adata = scanpy_viz(adata, keys)
    print(adata)


if __name__ == "__main__":
    test_gex_embedding_scanpy_viz()
