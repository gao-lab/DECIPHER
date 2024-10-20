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
# # Ligand-receptor pairs

# %%
import pandas as pd
import scanpy as sc
import numpy as np

from spider.utils import estimate_spot_size
from spider.graphic.build import build_graph

# %%
adata = sc.read_h5ad('./data/lymph_node.h5ad')

# %%
# estimate_spot_size(adata.obsm['spatial'])
build_graph(adata.obsm['spatial'], radius = 20)

# %%
lr = pd.read_csv('../../resource/lr_data/cellchat.csv')

# %%
lr.head(3)


# %%
lr.columns

# %%
ligand_list = []
receptor_list = []
for row in lr.iterrows():
    ligand = row[1]['ligand.symbol']
    if ',' in ligand:
        ligand_list += ligand.split(',')
    else:
        ligand_list.append(ligand)
    
    receptor = row[1]['receptor.symbol']
    # print(receptor)
    if ',' in receptor:
        receptor_list += receptor.split(',')
    else:
        receptor_list.append(receptor)
ligands = list(set(ligand_list))
receptors = list(set(receptor_list))

# %%
len(ligands)

# %%
len(receptors)

# %%
valid_ligands = set(adata.var_names) & set(ligands)
valid_receptors = set(adata.var_names) & set(receptors)

# %%
valid_lr = np.zeros(len(lr), dtype=bool)
for i, row in lr.iterrows():
    ligand = row['ligand.symbol']
    if ',' in ligand:
        ligand = set(ligand.split(','))
    else:
        ligand = set([ligand])
    if len(ligand & valid_ligands) == 0:
        continue

    receptor = row['receptor.symbol']
    if ',' in receptor:
        receptor = set(receptor.split(','))
    else:
        receptor = set([receptor])
    if len(receptor & valid_receptors) == 0:
        continue
    valid_lr[i] = True

# %%
valid_lr = lr[valid_lr]

# %%
valid_lr['pathway_name'].value_counts()

# %%
