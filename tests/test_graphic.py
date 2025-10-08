import numpy as np

from decipher.graphic.build import build_graph


def test_build_graph():
    n_cells = 100
    center_emb = np.random.rand(n_cells, 10)
    edge_index = build_graph(center_emb, mode="knn", k=5)
    assert edge_index.shape[1] == n_cells * 5

    edge_index = build_graph(center_emb, mode="radius", radius=0.5)
    print(edge_index.shape)


if __name__ == "__main__":
    test_build_graph()
