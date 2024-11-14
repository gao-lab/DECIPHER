import time

import scanpy as sc
from loguru import logger
from decipher import DECIPHER, CFG

# %%
CFG.omics.model.augment.dropout_gex = 0.4
CFG.omics.model.batch_size = 512
# disable validation plot
CFG.omics.model.max_steps = 20_000

# %%
# start_time = time.time()
# model = DECIPHER(work_dir="./results/decipher_7_15", user_cfg=CFG)

# %%
# adata_path = "./data/pancancer_filter_anno.h5ad"
# adata = sc.read_h5ad(adata_path)

# %%
# model.register_data(adata, preprocess=False, split_by='dataset')
# data_time = time.time()
# logger.info(f"Register data cost {data_time - start_time:.2f}s")

# %%
model = DECIPHER(work_dir="./results/decipher_7_15", recover=True)

# %%
model.fit_ddp(gpus=6)
# logger.info(f"Run model cost {time.time() - data_time:.2f}s")
# logger.info(f"Total cost {time.time() - start_time:.2f}s")