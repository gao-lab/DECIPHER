import scanpy as sc
from decipher import DECIPHER

model = DECIPHER(work_dir="./results/decipher_abca-1_3d_0705", recover=True)
adata = sc.read_h5ad("./data/zhuang_dataset/abca1.h5ad")
model.train_regress_explain(adata, cell_type='class', reverse_regress=True)