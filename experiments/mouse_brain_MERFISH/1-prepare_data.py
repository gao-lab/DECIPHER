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

# %% [markdown]
# # Mouse brain atlas
# There are two main mouse brain data resources:
# - Zeng dataset: 500 genes x 4M cells, 58 slices
# - Zhuang datasets: 1177 genes x 9M cells, contain 4 sub-parts (zhuang-1, zhuang-2, zhuang-3 and zhuang-4)
#
# This [page](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas) describe two projects

# %%
from pathlib import Path

import torch
import scanpy as sc
import pandas as pd
import numpy as np
import plotly.express as px
from sklearn.preprocessing import LabelEncoder

from decipher.utils import gex_embedding
from decipher.graphic.knn import knn
from decipher.graphic.build import knn_to_edge_index


# %%
def build_h5ad(data_dir:str, cluster_anno:pd.DataFrame, region_anno:pd.DataFrame):
    data_dir = Path(data_dir)
    assert data_dir.exists()
    adata_file = list(data_dir.glob('*.h5ad'))[0]
    print(f'Loading {adata_file}')
    adata = sc.read_h5ad(adata_file)
    
    # load meta and filter cells
    meta = pd.read_csv(data_dir / 'cell_metadata.csv', index_col=0)
    adata.obs = adata.obs.join(meta, how='left', lsuffix='_l', rsuffix='')
    filter_num = adata.obs['x'].isna().sum()
    print(f'Filtering {filter_num} cells without spatial coordinates')
    adata = adata[~adata.obs['x'].isna()]
    
    # get spatial coordinates
    adata.obsm['spatial'] = adata.obs[['x', 'y']].values.copy()
    
    # get cluster annotation
    adata.obs = adata.obs.join(cluster_anno, on='cluster_alias')
    adata.obs = adata.obs.astype(str)
    
    # get region annotation
    ccf = pd.read_csv(data_dir / 'ccf_coordinates.csv', index_col=0)
    ccf.rename(columns={'x': 'x_ccf',
                        'y': 'y_ccf',
                        'z': 'z_ccf'},
                        inplace=True)
    adata.obs = adata.obs.join(ccf)
    filter_num = adata.obs['x_ccf'].isna().sum()
    print(f'Filtering {filter_num} cells without CCF coordinates')
    adata = adata[~adata.obs['x_ccf'].isna()]
    adata.obs = adata.obs.join(region_anno, on='parcellation_index')
    
    print(adata)
    return adata.copy()


# %% [markdown]
# ## Meta
#
# Mapping cluster to cell type, shared in all datasets
#
# https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_pivoted.csv
#

# %%
cluster_anno = pd.read_csv('./data/cluster_to_cluster_annotation_membership_pivoted.csv', index_col=0)
cluster_anno.head()

# %%
region_anno = pd.read_csv('./data/parcellation_to_parcellation_term_membership_acronym.csv')
region_anno.set_index('parcellation_index', inplace=True)
region_anno.columns = ['parcellation_%s'% x for x in  region_anno.columns]

region_anno.head()

# %% [markdown]
# ## Zeng dataset
#
# See https://github.com/AllenInstitute/abc_atlas_access/blob/main/descriptions/MERFISH-C57BL6J-638850.md

# %%
zeng = build_h5ad('./data/zeng_dataset', cluster_anno, region_anno)

# %%
zeng.write_h5ad('./data/zeng_dataset/zeng.h5ad')

# %% [markdown]
# ## Zhuang dataset
# See https://alleninstitute.github.io/abc_atlas_access/descriptions/Zhuang_dataset.html

# %% [markdown]
# ### ABCA-1

# %%
abca1 = build_h5ad('./data/zhuang_dataset/abca-1', cluster_anno, region_anno)

# %%
abca1.write_h5ad('./data/zhuang_dataset/abca1.h5ad')

# %% [markdown]
# ### ABCA-2

# %%
abca2 = build_h5ad('./data/zhuang_dataset/abca-2', cluster_anno, region_anno)

# %%
abca2.write_h5ad('./data/zhuang_dataset/abca2.h5ad')

# %% [markdown]
# ### ABCA-3

# %%
abca3 = build_h5ad('./data/zhuang_dataset/abca-3', cluster_anno, region_anno)

# %%
abca3.write_h5ad('./data/zhuang_dataset/abca3.h5ad')

# %% [markdown]
# ### ABCA-4

# %%
abca4 = build_h5ad('./data/zhuang_dataset/abca-4', cluster_anno, region_anno)

# %%
abca4.write_h5ad('./data/zhuang_dataset/abca4.h5ad')

# %% [markdown]
# ## Visualize data

# %%
abca1 = sc.read_h5ad('./data/zhuang_dataset/abca1.h5ad')
# abca2 = sc.read_h5ad('./data/zhuang_dataset/abca2.h5ad')
# abca3 = sc.read_h5ad('./data/zhuang_dataset/abca3.h5ad')
# abca4 = sc.read_h5ad('./data/zhuang_dataset/abca4.h5ad')
# zeng = sc.read_h5ad('./data/zeng_dataset/zeng.h5ad')

# %%
abca1

# %%
df_3d = abca1.obs
df_3d["size"] = 0.6

# %%
df_3d_sample = df_3d.sample(100_000)
df_3d_sample.head()
fig = px.scatter_3d(
    df_3d_sample,
    x="x_ccf",
    y="y_ccf",
    z="z_ccf",
    color="class",
    size="size",
    size_max=3,
    width=1000,
    height=700,
)
fig.show()

# %%
df_3d_sample = df_3d.sample(100_000)
df_3d_sample.head()
fig = px.scatter_3d(
    df_3d_sample,
    x="x_ccf",
    y="y_ccf",
    z="z_ccf",
    # color="parcellation_division",
    color = 'brain_section_label',
    size="size",
    size_max=3,
    width=1000,
    height=700,
)
fig.show()

# %%
df_3d_sample = df_3d.sample(100_000)
df_3d_sample['class'] = df_3d_sample['class'].astype(str)
# replace the rows' 'class' column to NA if 'class' not contain '01'
df_3d_sample.loc[~df_3d_sample['class'].str.contains('IT-ET Glut'), 'class'] = np.nan
df_3d_sample['class'].head()

fig = px.scatter_3d(
    df_3d_sample,
    x="x_ccf",
    y="y_ccf",
    z="z_ccf",
    color="class",
    size="size",
    size_max=3,
    width=1000,
    height=700,
    color_discrete_map = {'background': 'grey', 'IT-ET Glut': 'red'}
)
fig.show()

# %%
df_3d['x_ccf'].hist(bins=1000)

# %% [markdown]
# ### Build the 3D edge index across layers

# %%
df_3d.groupby('brain_section_label')['x_ccf'].mean()

# %%
ccf = df_3d[['y_ccf', 'z_ccf']].values
global_index = np.arange(df_3d.shape[0])

layer_mask_list = []
global_index_list = []

sort_layer = df_3d.groupby('brain_section_label')['x_ccf'].mean().sort_values().index
for layer in sort_layer:
    layer_mask = df_3d['brain_section_label'] == layer
    layer_mask_list.append(layer_mask)
    global_index_list.append(global_index[layer_mask])

# %%
nbr_dict = {}
# init a empty list for each cell
for i in range(df_3d.shape[0]):
    nbr_dict[i] = []

def knn_to_global_index(i, j, k):
    xi = ccf[layer_mask_list[i]]
    xj = ccf[layer_mask_list[j]]
    if i == j:
        xj = None
    nbr = knn(xi, xj, k=k, metric='euclidean')[0]
    global_index_nbr = global_index_list[j][nbr]
    global_index_self = global_index_list[i]
    for m, cell in enumerate(global_index_self):
        nbr_dict[cell] += list(global_index_nbr[m])


# %%
for i in range(len(layer_mask_list)):
    # self
    knn_to_global_index(i, i, 12)

    # last
    if i > 0:
        knn_to_global_index(i, i-1, 6)
    else:
        knn_to_global_index(i, i+2, 6)

    # next
    if i < len(layer_mask_list) - 1:
        knn_to_global_index(i, i+1, 6)
    else:
        knn_to_global_index(i, i-2, 6)

# %%
nbr = np.zeros((df_3d.shape[0], 24))
for k, v in nbr_dict.items():
    nbr[k] = np.array(v, dtype=int)

# %%
edge_index_3d = knn_to_edge_index(nbr).int()

# %%
edge_index_3d

# %%
torch.save(edge_index_3d, './results/abca-1_3d_edge_index.pt')

# %% [markdown]
# ## Single cell embedding

# %%
abca = abca1.concatenate(abca2, abca3, abca4, join='inner')

# %%
abca.X = abca.X.astype(np.int64)

# %%
abca.write_h5ad('./data/zhuang_dataset/abca.h5ad')

# %%
abca = gex_embedding(abca, rapids_after_scale=True, gpu_id=1)

# %%
abca.obs.head()

# %%
sc.pl.umap(abca, color=['feature_matrix_label', 'class',])

# %% [markdown]
# ## Prepare scaled data

# %%
abca = sc.read_h5ad('./data/zhuang_dataset/abca.h5ad')

# %%
sc.pp.normalize_total(abca, target_sum=1e3)
sc.pp.log1p(abca)
sc.pp.scale(abca, max_value=10)

# %%
# add batch
batch = LabelEncoder().fit_transform(abca.obs['feature_matrix_label'])
abca.obs['_batch'] = batch

# %%
abca.write_h5ad('./data/zhuang_dataset/abca_processed.h5ad')
