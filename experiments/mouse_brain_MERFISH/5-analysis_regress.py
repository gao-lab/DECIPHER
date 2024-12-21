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
# # Mouse braind datasert explain

# %%
import scanpy as sc
import pandas as pd
import seaborn as sns
import scipy as sp
import matplotlib.pyplot as plt
import numpy as np

# %%
adata = sc.read_h5ad('./data/zhuang_dataset/abca1.h5ad')
cell_num = adata.obs['class'].value_counts().to_frame()
cell_num['cell_type'] = cell_num.index.str.replace(' ', '_')

# %%
cell_num

# %%
df_3d = pd.read_csv('./results/decipher_abca-1_3d_0705/explain/explain_results.csv')
# sort by 'test/r2_score'
df_3d = df_3d.sort_values(by='test_r2_score', ascending=False)
colname = df_3d.columns.tolist()
colname[0] = 'index'
df_3d.columns = colname
df_3d['cell_type'] = df_3d['index'].str.split(':').str[-1]
df_3d = df_3d.merge(cell_num, left_on='cell_type', right_on='cell_type')
df_3d = df_3d[df_3d['count'] > 5000]
# df_3d.head()
df_3d

# %%
# change figure size to (3,4)
plt.figure(figsize=(5, 4))

df_3d_forward = df_3d[~df_3d['index'].str.contains('reverse')]
p = sns.barplot(df_3d_forward, x="cell_type", y="test_r2_score")
plt.xticks(rotation=90)
p.set_xticklabels(p.get_xticklabels(), size = 8)
p.set(xlabel=None, ylabel='Variance explained (R2)')
plt.tight_layout()
sns.despine()
plt.show()

# %%
plt.figure(figsize=(5, 4))
df_3d_reverse = df_3d[df_3d['index'].str.contains('reverse')]
p = sns.barplot(df_3d_reverse, x="cell_type", y="test_r2_score")
plt.xticks(rotation=90)
p.set_xticklabels(p.get_xticklabels(), size = 8)
p.set(xlabel=None, ylabel='Variance explained (R2)')
plt.tight_layout()
sns.despine()
plt.show()

# %%
plt.figure(figsize=(4, 4))

df_merge = df_3d_forward.merge(df_3d_reverse, on='cell_type', suffixes=('_forward', '_reverse'))
# sns.scatterplot(data=df_merge, x='test_r2_score_forward', y='test_r2_score_reverse')
# add diagonal line
# plt.plot([0, 1], [0, 1], transform=plt.gca().transAxes, ls="--", c="black")
p = sns.regplot(data=df_merge, x='test_r2_score_forward', y='test_r2_score_reverse')

r, p = sp.stats.pearsonr(df_merge['test_r2_score_forward'], df_merge['test_r2_score_reverse'])
ax = plt.gca()
plt.text(.05, .9, f'R={r:.2f}, p={p:.1e}', transform=ax.transAxes)
ax.set_xlabel('Variance explained (R2) in Forward', fontsize = 12)
ax.set_ylabel('Variance explained (R2) in Reverse', fontsize = 12)
sns.despine()

# add line y = x
x = np.linspace(0, 0.5, 100)
y = x
plt.plot(x, y, linestyle='--', color='black')

plt.show()

# %%
