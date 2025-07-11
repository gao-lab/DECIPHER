[![stars-badge](https://img.shields.io/github/stars/gao-lab/DECIPHER?logo=GitHub&color=yellow)](https://github.com/gao-lab/DECIPHER/stargazers)
[![Downloads](https://static.pepy.tech/badge/cell-decipher)](https://pepy.tech/project/cell-decipher)
[![build-badge](https://github.com/gao-lab/DECIPHER/actions/workflows/build.yml/badge.svg)](https://github.com/gao-lab/DECIPHER/actions/workflows/build.yml)
[![docs-badge](https://readthedocs.org/projects/cell-decipher/badge/?version=latest)](https://cell-decipher.readthedocs.io/en/latest/)
[![license-badge](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![PyPI](https://img.shields.io/pypi/v/cell-decipher?label=pypi)
![Python 3.10](https://img.shields.io/badge/python->=3.10-blue.svg)

<!-- [![codecov](https://codecov.io/gh/gao-lab/DECIPHER/graph/badge.svg?token=zgwG4u9v0F)](https://codecov.io/gh/gao-lab/DECIPHER) -->

# DECIPHER
<div align="center">

[Installation](#Installation) • [Documentation](#Documentation) • [Citation](#Citation) • [FAQ](#FAQ) • [Acknowledgement](#Acknowledgement)

</div>

`DECIPHER` aims to learn cells’ disentangled intracellular molecular identity embedding and extracellular spatial context embedding from spatial omics data.

![DECIPHER](./docs/_static/Model.png)

## Installation

### PyPI

> [!IMPORTANT]
> Requires Python >= 3.10 and CUDA-enabled GPU (CPU-only device is not recommended).

We recommend to install `cell-decipher` to a new conda environment:

```sh
mamba create -n decipher python==3.11 -c conda-forge -y && conda activate decipher
pip install cell-decipher
install_pyg_dependencies
```

(Optional) You can install [RAPIDS](https://docs.rapids.ai/install) to accelerate visualization.

```sh
mamba create -n decipher -c conda-forge -c rapidsai -c nvidia python=3.11 rapids=25.06 'cuda-version>=12.0,<=12.8' -y && conda activate decipher
pip install cell-decipher
install_pyg_dependencies
```

### Docker

Build docker image from [Dockerfile](./Dockerfile) or pull image from Docker Hub directly:

```sh
docker pull huhansan666666/decipher:latest
docker run --gpus all -it --rm huhansan666666/decipher:latest
```

## Documentation

### Minimal example
Here is a minimal example for quick start:

```python
import scanpy as sc
from decipher import DECIPHER
from decipher.utils import scanpy_viz

# Init model
model = DECIPHER(work_dir='/path/to/work_dir')

# Register data (adata.X is raw counts, adata.obsm['spatial'] is spatial coordinates)
adata = sc.read_h5ad('/path/to/adata.h5ad')
model.register_data(adata)

# Fit DECIPHER model
model.fit_omics()

# Clustering disentangled embeddings
adata.obsm['X_center'] = model.center_emb  # intracellular molecular embedding
adata.obsm['X_nbr'] =  model.nbr_emb  # spatial context embedding
adata = scanpy_viz(adata, ['center', 'nbr'], rapids=False)

# Plot
adata.obsm['X_umap'] = adata.obsm['X_umap_center'].copy()
sc.pl.umap(adata, color=['cell_type'])
adata.obsm['X_umap'] = adata.obsm['X_umap_nbr'].copy()
sc.pl.umap(adata, color=['region'])
```

### Tutorials
> Please check [**documentation**](https://cell-decipher.readthedocs.io/en/latest) for all tutorials.

| Name                                    | Description                                                  | Colab                                                        |
| --------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| [Basic Model Tutorial](./docs/tutorials/1-train_model.ipynb)                | Tutorial on how to use DECIPHER                            | [![Open In Colab](https://img.shields.io/badge/Colab-PyTorch-blue?logo=googlecolab)](https://colab.research.google.com/drive/14PEtrgqlf-KbLOTfBLc9gbx0YvY6mi0S?usp=sharing) |
| [Multi-slices with Batch Effects](./docs/tutorials/2-remove_batch.ipynb)     | Tutorial on how to apply DECIPHER to multiple slices with batch effects | [![Open In Colab](https://img.shields.io/badge/Colab-PyTorch-blue?logo=googlecolab)](https://colab.research.google.com/drive/1eLJeRDZFq2tlDUWpETlSxVUdzRv9CeSD?usp=sharing) |
| [Identify Localization-related LRs](./docs/tutorials/3-select_LRs.ipynb) | Tutorial on how to identify ligand-receptors which related wtih cells’ localization based on DECIPHER embeddings | Insufficient resources |
| [Multi-GPUs Training](./docs/tutorials.md#multi-gpu-training)                        | Tutorial on how to use DECIPHER with multi-GPUs on spatial atlas | Insufficient resources |
| [More technologies](./docs/tutorials/4-more_techs.ipynb)                        | Tutorial on how to use DECIPHER with multiple spatial technologies | Insufficient resources |

## Citation

*High-fidelity disentangled cellular embeddings for large-scale heterogeneous spatial omics via DECIPHER* ([biorxiv](https://www.biorxiv.org/content/10.1101/2024.11.29.626126v1) 2024)


> If you want to repeat our benchmarks and case studies, please check the [**benchmark**](./benchmark/) and [**experiments**](./experiments/) folder.

## FAQ
> Please open a new [github issue](https://github.com/gao-lab/DECIPHER/issues/new/choose) if you meet problem.

1. Visium or ST data

DECIPHER is designed for single cell resolution data. As for Visium or ST, you can still use DECIPHER after obtaining single-cell resolution through deconvolution or spatial mapping strategies.

<!-- 2. `CUDA out of memory` error

The `model.train_gene_select()` may need a lot GPU memory. For example, it needs ~40G GPU memory in [Identify Localization-related LRs](./docs/tutorials/3-select_LRs.ipynb) tutorial (with ~700k cells and 1k LRs). If your GPU device do not have enough memory, you still can train model on GPU but set `disable_gpu=True` in `model.train_gene_select()`. -->


## Acknowledgement
We thank the following great open-source projects for their help or inspiration:

- [vit-pytorch](https://github.com/lucidrains/vit-pytorch)
- [lightly](https://github.com/lightly-ai/lightly)
- [scib](https://github.com/theislab/scib)
- [rapids_singlecell](https://github.com/scverse/rapids_singlecell/)
