import scanpy as sc

from spider import Spider, CFG


CFG.omics.model.epochs = 6
# CFG.omics.model.plot = True
CFG.omics.loader.batch_size = 256

# read in data
adata = sc.read_h5ad("./data/merged_adata.h5ad")
adata = adata[adata.obs['batch'] == '0'].copy()

# fit model
# model = Spider(work_dir="./results/spider_rep1", user_cfg=CFG)
# model.register_data(adata, cell_type="cell_type", preprocess=False)
# model.fit_omics()

# explain model
model = Spider(work_dir="./results/spider_rep1", user_cfg=CFG, recover=True)
adata.X = adata.layers['counts'].copy()
model.train_regress_explain(adata, cell_type = 'cell_type', reverse_regress = True, explain_dir = 'explain_reverse')