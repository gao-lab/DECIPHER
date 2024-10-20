import scanpy as sc
import numpy as np
from decipher import DECIPHER


model = DECIPHER(work_dir="./results/decipher_6_10", recover=True)
model.nbr_emb = np.load('./results/decipher_6_10/model/lightning_logs/version_3/epoch-0_nbr_emb.npy')
model.center_emb = np.load('./results/decipher_6_10/model/lightning_logs/version_3/epoch-0_gex_emb.npy')

adata = sc.read_h5ad("./data/pancancer_filter_anno.h5ad")
adata.X = adata.layers['counts']
model.train_gene_select(
    adata,
    cell_type='cell_type',
    subsets=['T cell', 'B/Plasma cell', 'TAM', 'Mast cell', 'DC'],
    # batch='Index',
    )