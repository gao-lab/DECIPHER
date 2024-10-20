# Tutorials

We provide following tutorials for you to get started with `decipher`.

## Basic usage
You can check basic usages of `decipher` in following notebooks:

1. **Train model**: shows how to train `decipher` model on spatial omics data and get high-fidelity disentangled omics and spatial embedding.
2. **Identify localization related genes**: shows how to identity cell localization related genes via `decipher`'s disentangled embeddings.

```{eval-rst}
.. nbgallery::
    tutorials/train_model.ipynb
    tutorials/explain_select_genes.ipynb
```

## Advanced topics

### Hyperparameter setting

Default hyperparameters are robust for most cases. If you want to change the hyperparameters, you can specify them when initializing the DECIPHER class. You can modify any hyperparameters by changing its values in `CFG` object, which is a nested dictionary.

```python
from decipher import DECIPHER, CFG

print(CFG)

# modify the hyperparameters
CFG.omics.model.batch_size = 512

# Init model with user defined hyperparameters
model = DECIPHER(work_dir='/path/to/work_dir', user_cfg=CFG)
...
```


### Multiple slices / Batch removal
You can input a list of Anndata objects, each one is a spaital slices. `DECIPHER` will automatically view each object as one batch and remove the batch effects. If you do not want to remove batch effects, just set `CFG.omics.ignore_batch = True`

```python
from decipher import DECIPHER, CFG

# CFG.omics.ignore_batch = True  # uncomment this line if you do not want to remove batch effects

model = DECIPHER(work_dir='/path/to/work_dir', user_cfg=CFG)

adatas = [adata1, adata2]
model.register_data(adatas)

model.fit_omics()
```

If the input is an integrated Anndata object, you can specify `split_by` in `register_data()` function, `DECIPHER` will automatically split the Anndata object into batches inside

```python
from decipher import DECIPHER, CFG

model = DECIPHER(work_dir='/path/to/work_dir', user_cfg=CFG)

adata = sc.read_h5ad('/path/to/adata.h5ad')
model.register_data(adata, split_by='batch_key')

model.fit_omics()
```


### Multi-GPU training

`decipher` support DDP (distributed data parallel on multi-GPUs) to accelerate training on big spatial atlas (especially > 1 millons cells). You just need change **1 line** of code to enable DDP:

```python
from decipher import DECIPHER

# Init model
model = DECIPHER('/path/to/work_dir')

# Register data
adata = sc.read_h5ad('/path/to/adata.h5ad')
model.register_data(adata)

# DDP training, use `fit_ddp` instead of `fit_omics`
model.fit_ddp(gpus = 4)
```

```{note}
DDP training will consumes n times the system memory (n is the number of GPUs)
```
