# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import time
import random
import yaml
import scanpy as sc
import numpy as np

from banksy.main import median_dist_to_nearest_neighbour
from banksy.banksy_utils.filter_utils import normalize_total, filter_hvg
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.run_banksy import run_banksy_multiparam
# from banksy.banksy_utils.color_lists import spagcn_color


# %% tags=["parameters"]
# parameter tag
input_file = ""
center_emb_file = ""
nbr_emb_file = ""
time_file = ""

# %%
adata = sc.read_h5ad(input_file)

# %%
start = time.time()

coord_keys = ['x', 'y', 'spatial']
adata.obs['x'] = adata.obsm['spatial'][:, 0]
adata.obs['y'] = adata.obsm['spatial'][:, 1]

# Preprocessing
adata = normalize_total(adata)
adata, adata_allgenes = filter_hvg(adata,
                                   n_top_genes=2000,
                                   flavor="seurat")

# set
plot_graph_weights = True
k_geom = 15  # only for fixed type
max_m = 1  # azumithal transform up to kth order
nbr_weight_decay = "scaled_gaussian"  # can also be "reciprocal", "uniform" or "ranked"

# Find median distance to closest neighbours, the median distance will be `sigma`
nbrs = median_dist_to_nearest_neighbour(adata, key=coord_keys[2])



banksy_dict = initialize_banksy(adata,
                                coord_keys,
                                k_geom,
                                nbr_weight_decay=nbr_weight_decay,
                                max_m=max_m,
                                plt_edge_hist=False,
                                plt_nbr_weights=False,
                                plt_agf_angles=False
                                )


# The following are the main hyperparameters for BANKSY
resolutions = [0.7]  # clustering resolution for UMAP
pca_dims = [20]  # Dimensionality in which PCA reduces to
lambda_list = [0.2]  # list of lambda parameters

banksy_dict, banksy_matrix = generate_banksy_matrix(adata,
                                                    banksy_dict,
                                                    lambda_list,
                                                    max_m)
number_of_colors = 50
colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

results_df = run_banksy_multiparam(
    adata,
    banksy_dict,
    lambda_list,
    resolutions,
    color_list=colors,
    max_m=max_m,
    filepath='/tmp',
    key=coord_keys,
    pca_dims=pca_dims,
    annotation_key=None,
    max_labels=None,
    cluster_algorithm='leiden',
    match_labels=False,
    savefig=False,
    add_nonspatial=True,
    variance_balance=False
)


run_time = str(time.time() - start)
print('Runtime: ' + run_time)

# %%
emb = results_df['adata'][1].obsm['reduced_pc_20'].copy()
np.save(nbr_emb_file, emb)
np.save(center_emb_file, emb)

# %%
time_dic = {"run_time": run_time}
with open(time_file, "w") as f:
    yaml.dump(time_dic, f)
