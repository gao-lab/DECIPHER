# Tutorials

We provide following tutorials for you to get started with `decipher`.

## Basic usage
You can check basic usages of `decipher` in following notebooks:

1. **Train model**: shows how to train `decipher` and get independent omics embedding and spatial embedding from spatial omics data.
2. **Identify localization related genes**: shows how to identity cell localization related genes via `decipher` omics embedding and spatial embedding.

```{eval-rst}
.. nbgallery::
    tutorials/train_model.ipynb
    tutorials/explain_select_genes.ipynb
```

## Advanced

### Hyperparameter setting

Default hyperparameters are robust for most cases. If you want to change the hyperparameters, you can specify them when initializing the `decipher` object. The `CFG` object is a nested dictionary that contains all the hyperparameters. You can modify the hyperparameters by changing the values in `CFG` object.

```python
from decipher import DECIPHER, CFG

# modify the hyperparameters
CFG.omics.model.batch_size = 512

# Init model with user defined hyperparameters
model = DECIPHER(work_dir='/path/to/work_dir', user_cfg=CFG)
...
```


### Multiple slices / Batch removal
You can input a list of Anndata objects, each one is a spaital slices. DECIPHER model will automatically view each object as one batch and remove the batch effects. If you do not want to remove batch effects, please set `CFG.omics.ignore_batch = True`

```python
from decipher import DECIPHER, CFG

# CFG.omics.ignore_batch = True  # uncomment this line if you do not want to remove batch effects

model = DECIPHER(work_dir='/path/to/work_dir', user_cfg=CFG)

adatas = [adata1, adata2]
model.register_data(adatas)

model.fit_omics()
```

If your input is an integrated Anndata object, you can specify `split_by` in `register_data()` function, decipher will automatically split the integrated Anndata object into multiple batches inside. If you do not want to remove batch effects, please set `CFG.omics.ignore_batch = True`

```python
from decipher import DECIPHER, CFG

# CFG.omics.ignore_batch = True  # uncomment this line if you do not want to remove batch effects

model = DECIPHER(work_dir='/path/to/work_dir', user_cfg=CFG)

adata = sc.read_h5ad('/path/to/adata.h5ad')
model.register_data(adata, split_by='batch_key')

model.fit_omics()
```


### Multi-GPU training

For spatial atlas with > 5 millons cells, use DDP (distributed data parallel) mode with multi-GPUs can save a lot of time. Here is an example showing how to run with DDP:

```python
from decipher import DECIPHER

# Init model
model = DECIPHER('/path/to/work_dir')
# Register data
adata = sc.read_h5ad('/path/to/adata.h5ad')
model.register_data(adata)
# DDP training
model.fit_ddp(gpus = 4)
```

```{note}
DDP training also consumes n times the memory
```