import scanpy as sc
import numpy as np
from spider import Spider


model = Spider(work_dir="./results/spider_6_10", recover=True)
model.nbr_emb = np.load('./results/spider_6_10/model/lightning_logs/version_3/epoch-0_nbr_emb.npy')
model.center_emb = np.load('./results/spider_6_10/model/lightning_logs/version_3/epoch-0_gex_emb.npy')

adata = sc.read_h5ad("./data/pancancer_filter_anno.h5ad")
model.train_regress_explain(adata, cell_type='cell_type', batch='Index', reverse_regress=True)
model.train_regress_explain(
    adata,
    cell_type='cell_type',
    explain_dir = 'explain_shuffle',
    batch='Index',
    user_cfg={'shuffle': True}
    )