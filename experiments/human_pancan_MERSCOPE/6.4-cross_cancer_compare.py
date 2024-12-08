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
# # Cross cancer type gene selection comparation

# %%
from pathlib import Path
from collections import Counter

import scanpy as sc
import numpy as np
import pandas as pd

import rapids_singlecell as rsc
from spider.utils import manage_gpu, scanpy_viz

# %%
adata_path = "./results/adata_analysis.h5ad"
adata = sc.read_h5ad(adata_path)
adata.obs.groupby(['Index', 'Tissue_1']).size().unstack().fillna(0)

# %%
genes = pd.read_csv('./results/genes.csv', index_col=0)
genes.head(3)

# %% [markdown]
# ## TAM

# %%
# write as a function
result_dir = './results/spider_6_10/explain/'
cell_type = 'TAM'
top_k = 15

batch_dir = list(Path(result_dir).rglob(f'*{cell_type}_batch*'))
# print(batch_dir)

counter = Counter()
mask_index_list, gene_name_list, batch_list = [], [], []
score_list = []
for batch in batch_dir:
    batch_name = batch.stem.split(':')[-1]
    batch_list.append(batch_name)
        
    mask = np.load(batch / 'gene_mask.npy')
    mask = mask.astype(int)
    n_cells = mask.shape[0]
    mask = mask.sum(axis=0)
    max_index = np.argpartition(mask, -top_k)[-top_k:]
    # print(mask[max_index])
    mask_index_list.append(max_index)
    
    score_list.append(mask / n_cells)
    
    gene_name = genes.iloc[max_index].index
    gene_name_list.append(gene_name)
    counter.update(gene_name)

    print(batch_name, n_cells)


# %%
score_df = pd.DataFrame(score_list, columns=genes.index).T
score_df.columns = batch_list
score_df

# %%
counter.most_common(10)

# %%
gene_name_df = pd.DataFrame(gene_name_list, index=batch_list).T
gene_name_df

# %%
from math import ceil, floor
import matplotlib.pyplot as plt


def plot_rss(rss, cell_type, top_n=5, max_n=None, ax=None):
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(4, 4))
    if max_n is None:
        max_n = rss.shape[0]
    data = rss[cell_type].sort_values(ascending=False)[0:max_n]
    ax.plot(np.arange(len(data)), data, ".")
    # ax.set_ylim([floor(data.min() * 100.0) / 100.0, ceil(data.max() * 100.0) / 100.0])
    ax.set_ylabel("RSS")
    ax.set_xlabel("Regulon")
    ax.set_title(cell_type)
    ax.set_xticklabels([])

    font = {
        "color": "red",
        "weight": "normal",
        "size": 4,
    }
    
    regulon_name = 'SPP1'
    # get the value of 'SPP1' row
    spp1_value = data[data.index == 'SPP1']
    # get the rank of 'SPP1' row
    spp1_rank = data.index.get_loc('SPP1')
    ax.plot([spp1_rank, spp1_rank], [spp1_value, spp1_value], "r.")
    ax.text(
        spp1_rank + (max_n / 25),
        spp1_value,
        regulon_name,
        fontdict=font,
        horizontalalignment="left",
        verticalalignment="center",
    )

    # Plot top-n
    for idx, (regulon_name, rss_val) in enumerate(
        zip(data[0:top_n].index, data[0:top_n].values)
    ):  
        if regulon_name == 'SPP1':
            continue
        # set color to green
        ax.plot([idx, idx], [rss_val, rss_val], "r.", color='green')
        ax.text(
            idx + (max_n / 25),
            rss_val,
            regulon_name,
            fontdict=font,
            horizontalalignment="left",
            verticalalignment="center",
            # set color to yellow
            color='green'
        )



# %%
from adjustText import adjust_text

cats = ['Ovarian', 'Skin', 'Melanoma', 'Uterine',
 'Lung', 'Prostate', 'Colon' ,'Breast']

fig = plt.figure(figsize=(15, 8))
for c,num in zip(cats, range(1, len(cats)+1)):
    x=score_df[c]
    ax = fig.add_subplot(2,4,num)
    plot_rss(score_df, c, top_n=5, max_n=None, ax=ax)
    # ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(14)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
# fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Gene score', ha='center', va='center', rotation='vertical', size=20)
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
# plt.savefig("PBMC10k_cellType-RSS-top5.pdf", dpi=600, bbox_inches = "tight")
plt.show()
# plot_rss(score_df, cell_type, top_n=5, max_n=10)

# %%
tam = adata[adata.obs['cell_type'] == 'TAM'].copy()
tam.X = tam.layers['counts']


# %%
# norm and log
sc.pp.normalize_total(tam, target_sum=1e4)
sc.pp.log1p(tam)

# %%
tam = scanpy_viz(tam, ['gex_0'], resolution=0.4)

# %%
tam.obsm['X_umap'] = tam.obsm['X_umap_gex_0'].copy()
sc.pl.umap(tam, color=['SPP1', 'C1QC', 'leiden_gex_0'], cmap='viridis')

# %%
tam.obs['subtype'] = 'Others'
tam.obs.loc[np.squeeze(tam.X[:, genes.index.get_loc('SPP1')].A > 4), 'subtype'] = 'SPP1+'
# tam.obs.loc[tam.obs['leiden_gex_0'].isin(['1', '2', '15']), 'subtype'] = 'SPP1+'

# %%
sc.pl.umap(tam, color=['SPP1', 'subtype'], cmap='viridis')

# %%
# change the font size
plt.rcParams.update({'font.size': 10})
sc.pl.dotplot(tam, var_names=['SPP1'], groupby='region', use_raw=False, log=True, )

# plot by subtype
for subtype in ['SPP1+', 'Others']:
    sc.pl.dotplot(tam[tam.obs['subtype'] == subtype], var_names=['SPP1'], groupby='region', use_raw=False, log=True, title=subtype)

# %%
import seaborn as sns
from scipy.stats import mannwhitneyu

# calculate the ratio of SPP1+ cells in each tissue
ratio = tam.obs.groupby('Tissue_1')['subtype'].value_counts(normalize=True).unstack().fillna(0)

ratio = pd.DataFrame(ratio.loc[:, 'SPP1+'])

# add group: high: Ovarian, Melanoma, Liver; low: Lung, Prostate, Colon, Breast, Uterine
ratio['group'] = 'Low'
ratio.loc[ratio.index.isin(['Ovarian', 'Melanoma', 'Liver']), 'group'] = 'Top'

# remove the row index 'Skin'
ratio = ratio.drop(['Skin'])
ratio.index = ratio.index.astype('str')

# plot the dot plot based on the group
plt.figure(figsize=(5, 4))
# color by the tissue 
sns.stripplot(x='group', y='SPP1+', data=ratio, jitter=True, hue='Tissue_1')
# move the legend outside
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# add statistical test
stat, p = mannwhitneyu(ratio[ratio['group'] == 'Top']['SPP1+'], ratio[ratio['group'] == 'Low']['SPP1+'])
plt.title(f'Mann-Whitney U test: p={p:.2}')
plt.ylabel('SPP1+ TAM / all TAM')
plt.xlabel('SPP1+ rank')

# %%
# 4 mins for 0.5 M cells, 8 mins for 0.8 M cells in A100-GPU
manage_gpu(2)
rsc.get.anndata_to_GPU(tam)
rsc.tl.embedding_density(tam, basis='umap', groupby='Tissue_1')
rsc.get.anndata_to_CPU(tam)

# %%
sc.pl.embedding_density(tam, basis='umap', groupby='Tissue_1', bg_dotsize=1, fg_dotsize=3)

# %% [markdown]
# ## B cell

# %%
# write as a function
result_dir = './results/spider_6_10/explain/'
cell_type = 'B_Plasma_cell'
top_k = 15

batch_dir = list(Path(result_dir).rglob(f'*{cell_type}_batch*'))
# print(batch_dir)

counter = Counter()
mask_index_list, gene_name_list, batch_list = [], [], []
score_list = []
for batch in batch_dir:
    batch_name = batch.stem.split(':')[-1]
    batch_list.append(batch_name)
        
    mask = np.load(batch / 'gene_mask.npy')
    mask = mask.astype(int)
    n_cells = mask.shape[0]
    mask = mask.sum(axis=0)
    max_index = np.argpartition(mask, -top_k)[-top_k:]
    # print(mask[max_index])
    mask_index_list.append(max_index)
    
    score_list.append(mask / n_cells)
    
    gene_name = genes.iloc[max_index].index
    gene_name_list.append(gene_name)
    counter.update(gene_name)

    print(batch_name, n_cells)


# %%
score_df = pd.DataFrame(score_list, columns=genes.index).T
score_df.columns = batch_list
score_df
