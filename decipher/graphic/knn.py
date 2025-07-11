r"""
k-NN related functions
"""
import faiss
import numpy as np
from annoy import AnnoyIndex
from loguru import logger

from ..utils import RSC_FLAG


def knn(
    query: np.ndarray,
    ref: np.ndarray = None,
    k: int = 30,
    metric: str = "cosine",
    approx: bool = False,
    method: list[str] | str = "auto",
    method_params: dict = {},
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Build k-NN graph, support multiple backends

    Parameters
    ----------
    query
        query embeddings
    ref
        reference embeddings
    k
        number of neighbors
    metric
        distance metric, can be one of ['euclidean', 'cosine']
    approx
        whether to use approximate method
    method
        which method to use, can be one of ['annoy', 'faiss', 'cuml', 'auto']
    method_params
        parameters to pass to the KNN algorithm

    Note
    ----------
    if ref is None, use self as ref, will ignore the itself in the result.
    """
    # select method
    if method == "auto":
        method = ["cuml", "faiss", "annoy"]
    method = method if isinstance(method, list) else [method]
    if not RSC_FLAG and "cuml" in method:
        method.remove("cuml")
    if not approx and "annoy" in method:
        method.remove("annoy")
    method = method[0]
    approx_str = " approx " if approx else " "
    logger.debug(f"Use {method} to compute{approx_str}KNN graph.")

    # if ref is None, use self as ref
    if ref is None:
        self_as_ref = True
        ref = query
        k = k + 1
    else:
        self_as_ref = False

    if "cuml" in method:
        distances, indices = knn_cuml(ref, query, k, metric)
    elif "faiss" in method:
        distances, indices = knn_faiss(ref, query, k, metric, approx)
    elif "annoy" in method:
        distances, indices = knn_annoy(ref, query, k, metric, **method_params)

    if self_as_ref:
        distances, indices = distances[:, 1:], indices[:, 1:]

    return indices, distances


def knn_faiss(
    ref: np.ndarray,
    query: np.ndarray,
    k: int = 30,
    metric: str = "euclidean",
    approx: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Build k-NN graph by faiss
    """
    shape = ref.shape[0]
    ref = ref.astype(np.float32)
    query = query.astype(np.float32)
    faiss.normalize_L2(ref)
    faiss.normalize_L2(query)
    if metric == "euclidean":
        metric = faiss.METRIC_L2
    elif metric == "cosine":
        metric = faiss.METRIC_INNER_PRODUCT

    if approx:  # 1M
        logger.warning(f"Use HNSW64 index with {shape} elements.")
        index = faiss.index_factory(ref.shape[1], "HNSW64", metric)
        index.train(ref)
    else:  # 500k
        logger.warning(f"Use Flat index with {shape} elements.")
        index = faiss.index_factory(ref.shape[1], "Flat", metric)
    index.add(ref)

    distances, indices = index.search(query, k)
    return distances, indices


def knn_cuml(
    ref: np.ndarray,
    query: np.ndarray,
    k: int = 30,
    metric: str = "euclidean",
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Build k-NN graph by cuML
    """
    from cuml.neighbors import NearestNeighbors as cuNearestNeighbors

    model = cuNearestNeighbors(n_neighbors=k, metric=metric)
    model.fit(ref)
    distances, indices = model.kneighbors(query)
    return distances, indices


def knn_annoy(
    ref: np.ndarray,
    query: np.ndarray = None,
    k: int = 30,
    metric: str = "euclidean",
    n_trees: int = 10,
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Build k-NN graph by annoy
    """
    metric = "angular" if metric == "cosine" else metric
    if metric == "cosine":
        # normalize the vectors
        ref = ref.astype(np.float32)
        query = query.astype(np.float32)
        faiss.normalize_L2(ref)
        faiss.normalize_L2(query)

    index = AnnoyIndex(ref.shape[1], metric)
    for i in np.arange(ref.shape[0]):
        index.add_item(i, ref[i])
    index.build(n_trees)

    ind_list, dist_list = [], []
    for i in np.arange(query.shape[0]):
        holder = index.get_nns_by_vector(query[i], k, include_distances=True)
        ind_list.append(holder[0])
        dist_list.append(holder[1])
    return np.asarray(dist_list), np.asarray(ind_list)
