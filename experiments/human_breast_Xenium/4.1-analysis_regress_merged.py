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
import scipy as sp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats


# %% [markdown]
# ## Different runs

# %%
df1 = pd.read_csv('./results/spider_rep1/explain/explain_results_1.csv')
df2 = pd.read_csv('./results/spider_rep1/explain/explain_results_2.csv')
df1['run'] = 'run1'
df2['run'] = 'run2'
df = pd.concat([df1, df2])
df['cell_type'] = df['Unnamed: 0'].str.split('celltype:').str[-1]


# %%
res_df = df.pivot(index='cell_type', columns='run', values='test/r2_score')
res_df.columns = ['rep1', 'rep2']

# %%
plt.figure(figsize=(4, 4))
sns.lmplot(x="rep1", y="rep2", data=res_df)
# set x lab size to 10
plt.xlim(0.2, 0.5)  
plt.ylim(0.2, 0.5)
plt.xlabel('Variance explained (R2) in Run 1', fontsize=12)
plt.ylabel('Variance explained (R2) in Run 2', fontsize=12)

# add line y = x
x = np.linspace(0, 1, 100)
y = x
plt.plot(x, y, linestyle='--', color='black')

## add r and p
r, p = stats.pearsonr(res_df['rep1'], res_df['rep2'])
ax = plt.gca()
plt.text(.05, .95, f'R={r:.2f}, p={p:.1e}', transform=ax.transAxes)

# %% [markdown]
# ## Different technology replicates

# %%
df = pd.read_csv('./results/spider_merged_6_28/explain/explain_results.csv', index_col=0)
df['experiment'] = df.index
df['cell_type'] = df.index.str.split('celltype:').str[-1]
df['batch'] = df.index.str.split('_cell').str[0].str.split(':').str[1]
df.head()

# %%
df = df[~df['experiment'].str.contains('reverse')]

res_df = df.pivot(index='cell_type', columns='batch', values='test/r2_score')
res_df.columns = ['rep1', 'rep2']
res_df


# %%
plt.figure(figsize=(4, 4))
sns.lmplot(x="rep1", y="rep2", data=res_df)
# set x lab size to 10
plt.xlim(0.0, 0.3)  
plt.ylim(0.0, 0.3)
plt.xlabel('Variance explained (R2) in Replication 1', fontsize=12)
plt.ylabel('Variance explained (R2) in Replication 2', fontsize=12)

# add line y = x
x = np.linspace(0, 1, 100)
y = x
plt.plot(x, y, linestyle='--', color='black')

## add r and p
r, p = stats.pearsonr(res_df['rep1'], res_df['rep2'])
ax = plt.gca()
plt.text(.05, .95, f'R={r:.2f}, p={p:.1e}', transform=ax.transAxes)

# %% [markdown]
# ## Forward and reverse

# %%
df = pd.read_csv('./results/spider_rep1/explain_reverse/explain_results.csv')
df['cell_type'] = df['Unnamed: 0'].str.split('celltype:').str[-1]
# create 'direct' col, if 'reverse' in Unnamed: 0, then 'direct' is 'reverse', else 'direct' is 'forward'
df['direct'] = df['Unnamed: 0'].apply(lambda x: 'reverse' if 'reverse' in x else 'forward')

# %%
res_df = df.pivot(index='cell_type', columns='direct', values='test_r2_score')
res_df.columns = ['rep1', 'rep2']
res_df = res_df[~res_df.index.str.contains('DC|Dendritic')]
res_df

# %%
plt.figure(figsize=(4, 4))
p = sns.regplot(data=res_df, x='rep1', y='rep2')
r, p = sp.stats.pearsonr(res_df['rep1'], res_df['rep2'])
ax = plt.gca()
plt.text(.05, .9, f'R={r:.2f}, p={p:.1e}', transform=ax.transAxes)
ax.set_xlabel('Variance explained (R2) in Forward', fontsize=12)
ax.set_ylabel('Variance explained (R2) in Reverse', fontsize=12)
plt.xlim(-0.03, 0.4)  
plt.ylim(-0.03, 0.4)

# add line y = x
x = np.linspace(-0.1, 1, 100)
y = x
plt.plot(x, y, linestyle='--', color='black')


sns.despine()
plt.show()

# %%
