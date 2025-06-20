# Tutorials

We provide following tutorials for you to get started with `cell-decipher`.

## Basic usage
You can check basic usages of DECIPHER model in following notebooks:

1. **Train model**: shows how to train DECIPHER model on spatial omics data and get high-fidelity disentangled omics and spatial embedding.
2. **Remove batch effects**: shows how to remove batch effects from spatial omics data.
3. **Identify localization related LRs**: shows how to identity cell localization related ligand-receptor pairs via DECIPHER disentangled embeddings.
4. **More techniques**: shows how to use DECIPHER in other spatial techniques, such as CosMx and CODEX.

```{eval-rst}
.. nbgallery::
    tutorials/1-train_model.ipynb
    tutorials/2-remove_batch.ipynb
    tutorials/3-select_LRs.ipynb
    tutorials/4-more_techs.ipynb
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

### Multi-GPU training

```{important}
Multi-GPU training is not available in Jupyter/Colab environment. Please run it as Python script.
```

```{warning}
DDP training will also consumes n times the system memory (n is the number of GPUs)
```

`cell-decipher` uses DDP (distributed data parallel on multi-GPUs) to accelerate training on atlas-scale spatial datasets (e.g. > 1 millons cells). You just need change **1 line** of code to enable DDP:

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
