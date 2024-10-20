
# Developer guide
> [!WARNING]
> This section is only for developers.

## Full install
```sh
git clone git@github.com:gao-lab/Spider-dev.git && cd Spider-dev
mamba create -p ./conda -c conda-forge -c rapidsai -c nvidia -c bioconda python=3.11 rapids=24.04 cuda-version=11.8 cudnn cutensor cusparselt pandoc snakemake==7.32.4 -y && conda activate ./conda
pip install -e ".[dev, docs]"
install_pyg_dependencies
```

## Unit test
> Not finish yet
```sh
pytest --cov -n 8
```

## Build docs
```sh
cd docs
sphinx-build -b html ./ ./_build/html
```

## Enable code style checking
```sh
pre-commit install
```

## Time analysis
```sh
py-spy record -o profile.svg -- python -m <YOUR_SCRIPT.py>
```

## Build docker image
```sh
cd Spider-dev/
docker image build --no-cache . -t huhansan666666/spider:v0.1.0 --network=host
```
