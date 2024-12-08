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
#     display_name: conda
#     language: python
#     name: python3
# ---

# %% [markdown] vscode={"languageId": "plaintext"}
# # Run LR-selection model of different niches

# %%
import scanpy as sc
import pandas as pd
import numpy as np
from addict import Dict
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

from spider import Spider, CFG
from spider.utils import scanpy_viz
from spider.utils import clip_umap

# %%
adata = sc.read_h5ad('./results/adata_analysis.h5ad')
lr = pd.read_csv('./results/spider/explain_bcell/lr.csv')


# %%
adata.obs.head(2)

# %% [markdown]
# ## Run B cell with global niches

# %%
adata.obs['cell_type_niche'] = adata.obs['cell_type'].astype(str) + '_' + adata.obs['leiden_nbr'].astype(str)
adata.obs['cell_type_niche'].value_counts()

# %%
subsets = adata.obs['cell_type_niche'].unique()
# only keep the cell types that contains 'Bcell'
subsets = [x for x in subsets if 'Bcell' in x]
subsets

# %%
model = Spider(work_dir="./results/spider", user_cfg=CFG, recover=True)

adata.X = adata.layers["counts"].copy()
adata.var_names_make_unique()
explain_cfg = Dict(epochs = 100)
model.train_gene_select(adata, cell_type='cell_type_niche', subsets=subsets, sub_dir='explain_bcell',
                        lr_mode=True, lr_data=lr, lr_radius=20, user_cfg=Dict(epochs = 100))

# %% [markdown]
# ### Analysis results

# %%
result_dir = './results/spider/explain_bcell'
top_k = 15
batch_dir = list(Path(result_dir).rglob(f'*.npy'))
batch_dir

# %%
lr.head(3)

# %%
top_k = 15
mask_index_list, gene_name_list, batch_list = [], [], []
for batch in batch_dir:
    batch_name = str(batch.parent).split('celltype_')[-1]
    batch_list.append(batch_name)
        
    mask = np.load(batch)
    mask = mask.astype(int)
    score = mask.sum(axis=0) / mask.shape[0]
    
    max_index = np.argpartition(score, -top_k)[-top_k:]
    # print(mask[max_index])
    mask_index_list.append(max_index)
    
    lr[f'score-{batch_name}'] = score
    
    gene_name = lr.iloc[max_index, 1]
    gene_name_list.append(gene_name.tolist())

    print(batch_name, mask.shape[0])


# %%
# subset lr 'interaction_name' in ['CXCL12_CXCR4', 'CXCL13_CXCR5']
region_df = lr[lr['interaction_name'].isin(['CXCL12_CXCR4', 'CXCL13_CXCR5'])]
region_df.index = region_df['interaction_name']
# only keep columns with 'score'
region_df = region_df.filter(like='score', axis=1)

# add a row ‘ratio’ which is the results of row1 / row2
region_df.loc['ratio'] = region_df.loc['CXCL12_CXCR4'] / region_df.loc['CXCL13_CXCR5']
# log2 of row 'ratio'
region_df.loc['log2_ratio'] = np.log2(region_df.loc['ratio'])
region_df = region_df.T

# sort by 'log2_ratio'
region_df = region_df.sort_values('log2_ratio', ascending=False)

region_df

# %%
# bar plot of the ratio
sns.set(style='whitegrid')
# rank the x-axis labels by value
sns.barplot(x=region_df.index, y='log2_ratio', data=region_df, palette='viridis')
# remove the grid
sns.despine()
plt.ylabel('CXCR4 / CXCR5')
plt.xlabel('')

# rotate x-axis labels
_ = plt.xticks(rotation=30, ha='right')
# rank the x-axis labels by value

# %%
region_df['cluster'] = region_df.index.str.split('_').str[1]
region_df['group'] = 'Follicular B cell'
region_df.loc[region_df['cluster'] == '2', 'group'] = 'Light Zone B cell'
region_df.loc[region_df['cluster'] == '5', 'group'] = 'Dark Zone B cell'
region_df
# group by 'group' and get the mean of 'CXCL12_CXCR4' and 'CXCL13_CXCR5' columns


# %%
group_region_df = region_df.groupby('group').mean(numeric_only=True)
# add a row ‘ratio’ which is the results of row1 / row2
group_region_df['ratio'] = group_region_df['CXCL13_CXCR5'] / group_region_df['CXCL12_CXCR4']
# log2 of row 'ratio'
group_region_df['log2_ratio'] = np.log2(group_region_df['ratio'])

# sort by 'log2_ratio'
group_region_df = group_region_df.sort_values('log2_ratio', ascending=True)
group_region_df


# %%
# bar plot of the ratio
# change figure size to 5, 5
plt.figure(figsize=(3, 4.5))
# sns.set(style='whitegrid')
# rank the x-axis labels by value
sns.barplot(x=group_region_df.index, y='log2_ratio', data=group_region_df,
            palette=['#1f77b4', '#ff7f0e', '#2ca02c'])
# remove the grid
sns.despine()
plt.ylabel('log2 CXCL13_CXCR5 score / CXCL12_CXCR4 score', fontsize=10)
plt.xlabel('')
plt.grid(False)
# rotate x-axis labels
_ = plt.xticks(rotation=30, ha='right')
# rank the x-axis labels by value

# %%
#plot the heatmap of pandas dataframe tmp

_ = group_region_df.drop(columns=['ratio', 'log2_ratio'], inplace=False)
# change the size of the heatmap
plt.figure(figsize=(3, 2.5))
sns.heatmap(_, cmap='coolwarm', annot=True)
group_region_df
plt.xlabel('')

plt.grid(False) # remove the grid
_ = plt.ylabel('')

# %%
gene_name_df = pd.DataFrame(gene_name_list, index=batch_list).T
gene_name_df

# %%
bcells = adata[adata.obs['cell_type'].str.contains('Bcell')].copy()
bcells

# %%
sc.pl.spatial(bcells, color='leiden_nbr', spot_size=5, title='Bcell')

# split by leiden_nbr
# for cluster in bcells.obs['leiden_nbr'].unique():
#     tmp = bcells[bcells.obs['leiden_nbr'] == cluster]
#     sc.pl.spatial(tmp, color='leiden_nbr', spot_size=10, title=f'Bcell_{cluster}')

# %%
# anno the spatial niches
bcells.obs['region'] = 'Follicular B cell'
bcells.obs.loc[bcells.obs['leiden_nbr'].isin(['2']), 'region'] = 'Light Zone B cell'
bcells.obs.loc[bcells.obs['leiden_nbr'].isin(['5']), 'region'] = 'Dark Zone B cell'

# %%
sc.pl.spatial(bcells, color='region', spot_size=10, title='B cell', palette=['#ff7f0e', '#1f77b4', '#2ca02c'])

# %%
sc.pl.spatial(bcells, color='region', spot_size=10, title='B cell')

# split by region
for cluster in bcells.obs['region'].unique():
    _ = bcells[bcells.obs['region'] == cluster]
    sc.pl.spatial(_, color='region', spot_size=10, title=cluster)

# %% [markdown]
# ## Re-cluster B cells

# %%
# bcells = scanpy_viz(bcells, keys=['nbr'], leiden=False)
bcells = scanpy_viz(bcells, keys=['center'], leiden=False, n_neighbors=30)

# %%
bcells.obsm['X_umap'] = bcells.obsm['X_umap_center'].copy()
bcells.obsm['X_umap'] = clip_umap(bcells.obsm['X_umap'])
sc.pl.umap(bcells, color='region', title='B cell')

# %%
bcells.obsm['X_umap'] = bcells.obsm['X_umap_nbr'].copy()
sc.pl.umap(bcells, color='region', title='UMAP of spatial context embedding')

# %%
bcells.write_h5ad('./results/bcells_analysis.h5ad')

# %%
sc.pl.spatial(bcells, color='leiden_nbr', spot_size=5, title='Bcell')

# %%
# split by leiden_nbr
for cluster in bcells.obs['leiden_nbr'].unique():
    tmp = bcells[bcells.obs['leiden_nbr'] == cluster]
    sc.pl.spatial(tmp, color='leiden_nbr', spot_size=10, title=f'Bcell_{cluster}')

# %% [markdown]
# ### Run B cell LR selection with new cluster results

# %%
# updata B cell leiden cluster info in adata.obs
adata.obs.loc[adata.obs['cell_type'].str.contains('Bcell'), 'leiden_nbr'] = bcells.obs['leiden_nbr']


# %%
adata.obs['cell_type_niche'] = adata.obs['cell_type'].astype(str) + '_' + adata.obs['leiden_nbr'].astype(str)
adata.obs['cell_type_niche'].value_counts()

# %%
subsets = adata.obs['cell_type_niche'].unique()
# only keep the cell types that contains 'Bcell'
subsets = [x for x in subsets if 'Bcell' in x]
subsets

# %%
model = Spider(work_dir="./results/spider", user_cfg=CFG, recover=True)

adata.X = adata.layers["counts"].copy()
adata.var_names_make_unique()
explain_cfg = Dict(epochs = 100)
model.train_gene_select(adata, cell_type='cell_type_niche', subsets=subsets, sub_dir='explain_bcell_recluster',
                        lr_mode=True, lr_data=lr, lr_radius=20, user_cfg=Dict(epochs = 100))

# %% [markdown]
# ### Analysis

# %%
result_dir = './results/spider/explain_bcell_recluster'
top_k = 15
batch_dir = list(Path(result_dir).rglob(f'*.npy'))
batch_dir

# %%
top_k = 15
mask_index_list, gene_name_list, batch_list = [], [], []
for batch in batch_dir:
    batch_name = str(batch.parent).split('celltype_')[-1]
    batch_list.append(batch_name)
        
    mask = np.load(batch)
    mask = mask.astype(int)
    score = mask.sum(axis=0) / mask.shape[0]
    
    max_index = np.argpartition(score, -top_k)[-top_k:]
    # print(mask[max_index])
    mask_index_list.append(max_index)
    
    lr[f'score-{batch_name}'] = score
    
    gene_name = lr.iloc[max_index, 1]
    gene_name_list.append(gene_name.tolist())

    print(batch_name, mask.shape[0])


# %%
# subset lr 'interaction_name' in ['CXCL12_CXCR4', 'CXCL13_CXCR5']
tmp = lr[lr['interaction_name'].isin(['CXCL12_CXCR4', 'CXCL13_CXCR5'])]
tmp.index = tmp['interaction_name']
# only keep columns with 'score'
tmp = tmp.filter(like='score', axis=1)

# add a row ‘ratio’ which is the results of row1 / row2
tmp.loc['ratio'] = tmp.loc['CXCL12_CXCR4'] / tmp.loc['CXCL13_CXCR5']
# log2 of row 'ratio'
tmp.loc['log2_ratio'] = np.log2(tmp.loc['ratio'])
tmp = tmp.T

# sort by 'log2_ratio'
tmp = tmp.sort_values('log2_ratio', ascending=False)

tmp

# %%
# bar plot of the ratio
sns.set(style='whitegrid')
# rank the x-axis labels by value
sns.barplot(x=tmp.index, y='log2_ratio', data=tmp, palette='viridis')
# remove the grid
sns.despine()
plt.ylabel('CXCR4 / CXCR5')
plt.xlabel('')

# rotate x-axis labels
_ = plt.xticks(rotation=30, ha='right')
# rank the x-axis labels by value

# %%
tmp

# %%
#plot the heatmap of pandas dataframe tmp
plt_df = tmp[['CXCL12_CXCR4','CXCL13_CXCR5']].T


# change the size of the heatmap
plt.figure(figsize=(9, 2.5))
sns.heatmap(plt_df, cmap='coolwarm', annot=True)

plt.grid(False) # remove the grid
_ = plt.ylabel('')

# %%
