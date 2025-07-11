import scanpy as sc
import torch
from decipher import CFG, DECIPHER

# set the parameters
CFG.omics.pretrain.epochs = 2
CFG.omics.model.epochs = 2
CFG.omics.model.augment.dropout_gex = 0.4
CFG.omics.loader.batch_size = 512
# disable validation plot
CFG.omics.ignore_batch = True


model = DECIPHER(work_dir="./results/decipher_abca-1_3d_0705", user_cfg=CFG)

# read in data
adata = sc.read_h5ad("./data/zhuang_dataset/abca_processed.h5ad")
adata = adata[adata.obs["feature_matrix_label"].isin(["Zhuang-ABCA-1"])].copy()
edge_index = torch.load('./results/abca-1_3d_edge_index.pt').long()
model.register_data(adata, cell_type='class', preprocess=False, edge_index=edge_index)

model.fit_omics()
