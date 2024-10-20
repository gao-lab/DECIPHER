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
import scanpy as sc

# %%
ref = sc.read_h5ad('./data/sc.h5ad')

# %%
sc.pl.dotplot(ref, ['CXCR4', 'CXCR5', 'CXCL12', 'CXCL13', 'CR2', 'FCER2'],  groupby='Subset')

# %%
sc.pl.dotplot(ref, ['STEAP4', 'CXCL9', 'TNC', 'NOTCH3', 'VEGFA'],  groupby='Subset')

# %%
