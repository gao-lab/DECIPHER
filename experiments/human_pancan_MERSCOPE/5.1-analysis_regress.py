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
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats


# %%
df = pd.read_csv('./results/spider_6_10/explain/explain_results.csv', index_col=0)
df['experiment'] = df.index
df['cell_type'] = df.index.str.split('celltype:').str[-1]
df['batch'] = df.index.str.split('_cell').str[0].str.split(':').str[1]
df

# %%
adata = sc.read_h5ad('./data/pancancer_filter_anno.h5ad')
adata

# %%
tissue_meta = adata.obs[['Tissue_1', 'Index']].value_counts().fillna(0).reset_index()
tissue_meta = pd.DataFrame(tissue_meta)
tissue_meta.head()


# %%
# merge tissue_meta with df
df = df.merge(tissue_meta, left_on="batch", right_on='Index')
df.head()

# %%
# remove rows which contrains 'reverse' in the 'experiment' column
df = df[~df['experiment'].str.contains('reverse')]
df.shape

# %%
df.to_csv('./results/explain_results_rep.csv')

# %%
rep = df[df['batch'].isin(['0', '1'])]
res_df = rep.pivot(index='cell_type', columns='batch', values='test/r2_score')
res_df.columns = ['rep1', 'rep2']
res_df


# %%
sns.lmplot(x="rep1", y="rep2", data=res_df)

# %%
stats.pearsonr(res_df['rep1'], res_df['rep2'])


# %%
stats.spearmanr(res_df['rep1'], res_df['rep2'])
