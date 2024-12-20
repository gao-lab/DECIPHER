import numpy as np

from decipher.graphic.knn import knn
from decipher.utils import GetRunTime


@GetRunTime
def test_knn(n_cells=1000, n_pcs=50, n_neighbors=10):
    center_emb = np.random.rand(n_cells, n_pcs)
    indices, distances = knn(center_emb, k=n_neighbors)


if __name__ == "__main__":
    test_knn()
