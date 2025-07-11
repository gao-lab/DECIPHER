import scanpy as sc

from decipher import DECIPHER, CFG

# read in data
adata = sc.read_h5ad("./data/merged_adata.h5ad")
adata = adata[adata.obs['batch'] == '0'].copy()

# fit model
# model = DECIPHER(work_dir="./results/decipher_rep1", user_cfg=CFG)
# model.register_data(adata, cell_type="cell_type", preprocess=False)
# model.fit_omics()

# explain model
model = DECIPHER(work_dir="./results/decipher_rep1", user_cfg=CFG, recover=True)
adata.X = adata.layers['counts'].copy()
model.train_regress_explain(adata, cell_type = 'cell_type', reverse_regress = True, explain_dir = 'explain_reverse')