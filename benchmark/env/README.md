# Benchmark environment setup

## Python
You can run following commands to install all dependencies for benchmarking:

```sh
mamba create -p ./conda -y -c conda-forge -c rapidsai -c nvidia -c bioconda python=3.11 \
rapids=24.04 cuda-version=11.8 cudnn cutensor cusparselt \
snakemake==7.32.4 && conda activate ./conda

# decipher & scanpy
# cd /path/to/decipher
pip install '.[dev]'
install_pyg_dependencies

# scvi-tools
pip install scvi-tools

# harmony
pip install harmonypy

# stagate
pip install git+https://github.com/QIFEIDKN/STAGATE_pyG.git

# banksy_py
# must import banksy_py from raw py file beacause it not provide a package
git clone https://github.com/prabhakarlab/Banksy_py

# spagcn
git clone https://github.com/jianhuupenn/SpaGCN
cd SpaGCN
pip install .

# scib (with anndata2ri and rpy2)
pip install scib anndata2ri rpy2

# scNiche
git clone https://github.com/ZJUFanLab/scNiche
cd scNiche
# NOTE: you should install DGL according to your torch and CUDA version
pip install  dgl -f https://data.dgl.ai/wheels/torch-2.3/cu121/repo.html
python setup.py build
python setup.py install


# export env
mamba  env export | grep -v "^prefix: " | grep -v "^name: " > conda.yaml
```

You can also check absolute version of install packages in [conda.yaml](./conda.yaml) file.

## R
`scib` use R package `kBET` to calculate kBET score. We need install it manually:

```R
install.packages("pak")
pak::pkg_install("theislab/kBET")
```
