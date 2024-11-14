import scanpy as sc
import pandas as pd
from decipher import DECIPHER, CFG
from addict import Dict

CFG.omics.model.epochs = 6
CFG.omics.model.plot = True
CFG.omics.loader.batch_size = 256

# read in data
adata = sc.read_h5ad("./data/lymph_node.h5ad")
adata.X = adata.layers["counts"].copy()

# fit model
# model = DECIPHER(work_dir="./results/decipher", user_cfg=CFG)
# model.register_data(adata, cell_type="cell_type")
# model.fit_omics()

# explain model
model = DECIPHER(work_dir="./results/decipher", user_cfg=CFG, recover=True)
adata.X = adata.layers['counts'].copy()
# remove 'cell_type' == 'low_quality' cells
# adata = adata[~adata.obs['cell_type'].isin(['low_quality', 'unknown']), :]
lr = pd.read_csv('../../resource/lr_data/cellchat.csv')
explain_cfg = Dict(epochs = 100)
model.train_gene_select(adata, cell_type='cell_type', 
                        lr_mode=True, lr_data=lr, lr_radius=20, user_cfg=explain_cfg)
