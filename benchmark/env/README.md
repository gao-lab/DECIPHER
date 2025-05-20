# Benchmark environment setup
You can run following commands to install all dependencies for benchmarking:

## Python
We use `conda` to manage Python package dependencies. You can check absolute version of packages in [conda.yaml](./conda.yaml) file.

```sh
mamba create -p ./conda -y -c conda-forge -c rapidsai -c nvidia -c bioconda python=3.11 \
rapids=24.04 cuda-version=11.8 snakemake==7.32.4 && conda activate ./conda

# DECIPHER & scanpy
# cd /path/to/DECIPHER
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

# scniche
pip install dgl -f https://data.dgl.ai/wheels/torch-2.3/cu121/repo.html
git clone https://github.com/ZJUFanLab/scNiche
cd scNiche
python setup.py build
python setup.py install

# scib (with anndata2ri and rpy2)
pip install scib anndata2ri rpy2

## graphst
git clone https://github.com/JinmiaoChenLab/GraphST.git
cd GraphST
pip install pot
pip install .

# export env
mamba env export | grep -v "^prefix: " | grep -v "^name: " > conda.yaml
```

## R

We use `renv` to manage R package dependencies. You can also check absolute version of R packages in [renv.lock](./renv.lock) file.

```R
install.packages("renv")
renv::init(bare = TRUE)

install.packages("pak")
pak::pkg_install('IRkernel')
pak::pkg_install("theislab/kBET") # `scib` use R package `kBET` to calculate kBET score.
pak::pkg_install('zhengli09/BASS')
pak::pkg_install('yanfang-li/STADIA')

renv::snapshot(type = 'all')
```
