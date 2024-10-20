# Installation guide

## PyPI

Install ``spider`` on machine with CUDA-enabled GPU is highly recommanded. We create a new conda environment with Python>=3.10 before installing ``spider``.

```sh
mamba create -p ./conda -c conda-forge -c rapidsai -c nvidia -c bioconda python=3.11 rapids=24.04 cuda-version=11.8 cudnn cutensor cusparselt -y && conda activate ./conda
pip install spider
install_pyg_dependencies
```

```{note}
RAPIDS is only avaliable at NVIDIA Volta or newer GPU with specific driver version. Please check official [website](https://docs.rapids.ai/install) for more details.
```

`RAPIDS` is not coercive, but it is highly recommanded to install it for acceleration. If you do not want to use it or meet bugs, please run:

```sh
mamba create -p ./conda python==3.11 -c conda-forge -y && conda activate ./conda
pip install spider
install_pyg_dependencies
```

## Docker

[`Dockerfile`](https://github.com/gao-lab/Spider/blob/main/Dockerfile) of `spider` package is available. You can also directly pull [image](https://hub.docker.com/repository/docker/huhansan666666/slat) by :
```
docker pull huhansan666666/spider:lastest
```
