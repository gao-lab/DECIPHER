# Installation guide

## PyPI

We recommend to install `cell-decipher` to a new conda environment with [RAPIDS](https://docs.rapids.ai/install) dependencies.

```sh
mamba create -n decipher -c conda-forge -c rapidsai -c nvidia python=3.11 rapids=24.04 cuda-version=11.8 cudnn cutensor cusparselt -y && conda activate decipher
pip install cell-decipher
install_pyg_dependencies
```

```{note}
RAPIDS is only avaliable at NVIDIA Volta or newer GPU with specific driver version. Please check official [website](https://docs.rapids.ai/install) for more details.
```

`RAPIDS` is not coercive, but it is highly recommanded to install it for acceleration. If you do not want to install it or meet unexpected errors, just run:

```sh
mamba create -n decipher  python==3.11 -c conda-forge -y && conda activate decipher
pip install cell-decipher
install_pyg_dependencies
```

## Docker

Build docker image from [Dockerfile](https://github.com/gao-lab/DECIPHER/blob/main/Dockerfile) or pull the latest image from Docker Hub by:
```sh
docker pull huhansan666666/decipher:latest
```
