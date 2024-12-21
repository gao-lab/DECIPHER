import numpy as np
from pytest import mark

from decipher.graphic.build import build_graph
from decipher.graphic.knn import knn


@mark.parametrize("method", ["auto", "faiss", "annoy"])
def test_knn(method, n_cells=1000, n_pcs=50, k=10):
    center_emb = np.random.rand(n_cells, n_pcs)
    if method == "annoy":
        approx = True
    else:
        approx = False
    indices, distances = knn(center_emb, k=k, method=method, approx=approx)
    assert indices.shape == (n_cells, k)
    assert distances.shape == (n_cells, k)


def test_build_graph():
    n_cells = 100
    center_emb = np.random.rand(n_cells, 10)
    edge_index = build_graph(center_emb, mode="knn", k=10)
    assert edge_index.shape[1] == n_cells * 10

    edge_index = build_graph(center_emb, mode="radius", radius=0.5)
    print(edge_index.shape)


if __name__ == "__main__":
    test_knn("auto")
    test_build_graph()
