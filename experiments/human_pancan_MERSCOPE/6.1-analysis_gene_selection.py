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

# %% [markdown]
# # Pan-cancer dataset gene select analysis

# %%
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# %%
# adata = sc.read_h5ad('./data/pancancer_filter_anno.h5ad')
# adata.obs['cell_type'].value_counts()


# %%
# gene_df = pd.DataFrame(adata.var.index, columns=['gene'])
gene_df = pd.read_csv('./results/genes.csv', index_col=0)
gene_df.head()

# %%
gene_mask = np.load('./results/decipher_6_10/explain/celltype_T_cell/gene_mask.npy')
print(gene_mask.shape)

gene_df['counts'] = gene_mask.sum(0)
# sort by counts
gene_df['score'] = gene_df['counts'] / gene_mask.shape[0]
gene_df = gene_df.sort_values('counts', ascending=False)
gene_df['rank'] = range(1, gene_df.shape[0] + 1)
gene_df.head(20)

# %%
# set figure size as (4, 6)
plt.figure(figsize=(4, 2.5))
gene_df2 = gene_df.iloc[[2,3,6,9,10], :]
sns.barplot(data=gene_df2.head(5), x='score', y='gene')
plt.xlabel('Score of T cell')
plt.ylabel('Gene')
plt.show()


# %%
# find the row which gene is in ['CTAL4', 'GNLY', 'GZMB', 'IFNG', 'PDCD1', 'LAG3']
gene_df[gene_df['gene'].isin(['CTLA4', 'GNLY', 'GZMB', 'IFNG', 'PDCD1', 'LAG3'])]

# %%
gene_df.to_csv('./results/t_cell_gene_rank.csv')
