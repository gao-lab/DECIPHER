import scanpy as sc

from decipher import DECIPHER, CFG


CFG.omics.model.epochs = 6
CFG.omics.model.plot = True
CFG.omics.loader.batch_size = 256

# read in data
adata = sc.read_h5ad("./data/merged_adata.h5ad")

# fit model
model = DECIPHER(work_dir="./results/decipher_merged_6_28", user_cfg=CFG)
model.register_data(adata, cell_type="cell_type", preprocess=False)
model.fit_omics()

# explain model
adata.X = adata.layers['counts'].copy()
model.train_regress_explain(adata, cell_type='cell_type', batch='batch')
