r"""
k-nearest neighbor searching
"""
from collections import defaultdict

import numpy as np
import pandas as pd
import torch
from loguru import logger
from numba import njit, prange
from rui_utils.knn import knn
from rui_utils.utils import l2norm
from sklearn.utils.extmath import randomized_svd


def cal_mnn(x, batch, mnn_config):
    r"""
    Init nearest neighbors pairs
    """
    n_cell = x.shape[0]
    cell_names = np.array([f"batch_{b}_{idx}" for idx, b in enumerate(batch)])
    name2idx = dict(zip(cell_names, range(n_cell)))

    # get MNN pairs
    n_pcs, k_mnn = mnn_config.k_components, mnn_config.k_anchor
    logger.info(f"Computing MNN pairs with k = {k_mnn}")
    batches, cells_per_batch = np.unique(batch, return_counts=True)
    n_batch = len(batches)

    all_mnn_pairs = []
    if mnn_config.ref_based:
        # ref batch is batch with most cells
        logger.debug("Use reference-based MNN")
        ref_batch = batches[np.argmax(cells_per_batch)]
        ref_batch_idx = batch == ref_batch
        X_i, name_i = x[ref_batch_idx], cell_names[ref_batch_idx]
        for batch_j in batches:
            if batch_j == ref_batch:
                continue
            batch_j_idx = batch == batch_j
            X_j, name_j = x[batch_j_idx], cell_names[batch_j_idx]
            mnn_pairs = findMNN(X_i, X_j, name_i, name_j, k_mnn, n_pcs)
            all_mnn_pairs.append(mnn_pairs)
    else:
        for i, batch_i in enumerate(batches):
            batch_i_idx = batch == batch_i
            X_i, name_i = x[batch_i_idx], cell_names[batch_i_idx]
            for j in range(i + 1, n_batch):
                batch_j_idx = batch == batches[j]
                X_j, name_j = x[batch_j_idx], cell_names[batch_j_idx]
                mnn_pairs = findMNN(X_i, X_j, name_i, name_j, k_mnn, n_pcs)
                all_mnn_pairs.append(mnn_pairs)
    pairs = pd.concat(all_mnn_pairs)
    # convert to global cell index
    pairs["cell1"] = pairs["cellname1"].apply(lambda x: name2idx[x])
    pairs["cell2"] = pairs["cellname2"].apply(lambda x: name2idx[x])
    pairs = pairs[["cell1", "cell2"]].values
    pairs = np.unique(pairs, axis=0)

    # get MNN dict
    mnn_dict = defaultdict(list)
    for r, c in pairs:
        mnn_dict[r].append(c)
        mnn_dict[c].append(r)

    # exclude cells without MNN
    mnn_cell_idx = np.unique(pairs.ravel())
    return mnn_dict, mnn_cell_idx


def findMNN(
    x: np.ndarray,
    y: np.ndarray,
    name_x: list[str],
    name_y: list[str],
    mnn_k: int = 5,
    k_components: int = 20,
) -> pd.DataFrame:
    r"""
    Find mutual nearest neighbors (MNN) pairs

    Parameters
    ----------
    x:
        scaled expression of batch A (gene x cell)
    y:
        scaled expression of batch b (gene x cell)
    name_x:
        cell names of batch A
    name_y:
        cell names of batch B
    mnn_k:
        number of nearest neighbors in computing MNN
    k_components:
        SVD components
    """
    # x, y = x.copy(), y.copy()
    z_norm = svd(x, y, k_components)
    z_x, z_y = z_norm[: x.shape[0]], z_norm[x.shape[0] :]

    knn_a2b, _ = knn(z_x, z_y, k=mnn_k)
    knn_b2a, _ = knn(z_y, z_x, k=mnn_k)

    pairs_a2b = create_pairs(knn_a2b)
    pairs_b2a = create_pairs(knn_b2a, reverse=True)
    pairs = pairs_a2b & pairs_b2a
    pairs = np.array(list(pairs))
    mnns = pd.DataFrame(pairs, columns=["cell1", "cell2"])
    mnns["cellname1"] = mnns.cell1.apply(lambda x: name_x[x])
    mnns["cellname2"] = mnns.cell2.apply(lambda x: name_y[x])
    logger.info(f"Found {mnns.shape[0]} MNN pairs")
    return mnns


def svd(x: np.ndarray, y: np.ndarray, k_components: int = 20) -> tuple[np.ndarray]:
    r"""
    Fast SVD

    Parameters
    ----------
    x:
        scaled expression of batch A (gene x cell)
    y:
        scaled expression of batch b (gene x cell)
    k_components:
        number of components to keep

    Returns
    ----------
    z_norm:
        normalized embedding
    """
    logger.debug(f"x shape: {x.shape}, y shape: {y.shape}")
    if x.shape[0] > 1_000_000 or y.shape[0] > 1_000_000:
        logger.debug("Use harmony-based SVD for large dataset.")
        from harmony import harmonize

        # batch
        batch = [0] * x.shape[0] + [1] * y.shape[0]
        batch = pd.DataFrame(batch, columns=["batch"])
        # pca
        z = np.vstack([x, y])
        z, _, _ = randomized_svd(z, n_components=k_components, random_state=0)
        # harmonize
        z_norm = harmonize(z, batch, "batch", use_gpu=True)
        return z_norm

    try:
        dot = torch.from_numpy(x).cuda().half() @ torch.from_numpy(y).T.cuda().half()
        dot = dot.cpu().float().numpy()
        logger.info("Use CUDA for small dataset")
    except:  # noqa
        logger.error(f"CUDA failed: {x.shape}, {y.shape}, use CPU instead.")
        dot = torch.from_numpy(x) @ torch.from_numpy(y).T
        dot = dot.numpy()
    torch.cuda.empty_cache()
    u, s, vh = randomized_svd(dot, n_components=k_components, random_state=0)
    z = np.vstack([u, vh.T])  # gene x k_components
    z = z @ np.sqrt(np.diag(s))  # will reduce the MNN pairs number greatly
    z_norm = l2norm(z)  # follow Seurat
    return z_norm


@njit
def create_pairs(knn_result: np.ndarray, reverse=False) -> tuple:
    r"""
    Create MNN pairs

    Parameters
    ----------
    knn_result:
        knn results
    reverse:
        reverse the pairs
    """
    num_rows = knn_result.shape[0]
    pairs = set()
    for i in prange(num_rows):
        for j in knn_result[i]:
            if reverse:
                pairs.add((j, i))
            else:
                pairs.add((i, j))
    return pairs
