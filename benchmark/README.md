# Benchmark

## Datasets

All datasets used in benchmark are public available. Here is a summary:

| **Dataset** | **Species** | **Technology** | Cells   | Genes  | **Tissue** | **Disease** | URL                                                          |
| ----------- | ----------- | -------------- | ------  | -----  | ---------- | ----------- | ------------------------------------------------------------ |
| 1           | Mouse       | MERFISH        | 378,918 | 374    | Brain      | Health      | [Link](https://www.cell.com/cell/fulltext/S0092-8674(22)01523-9) |
| 2           | Human       | Xenium         | 164,079 | 313    | Breast     | Tumor       | [Link](https://www.nature.com/articles/s41467-023-43458-x)   |
| 3           | Human       | 10X chromium   | 8,405   | 31,915 | PBMC       | Health      | [Link](https://www.10xgenomics.com/datasets)                 |


## Methods
- DECIPHER
- SLAT
- BANKSY
- STAGATE
- scVI
- Scanpy
- Hamrony


## Reproducibility

### Setup environment
Follow the [guide](./env/README.md) to set up the environment. The absolute version of packages is also provided in [conda.yaml](env/conda.yaml).

### Prepare data

Download the datasets and move them into `./input` folder.

### Run the workflow in slurm cluster

> Take about 12 hours when using 32 A100 GPUs.

We use `snakemake` to run benchmark workflows in Slurm cluster. Firstly, you need modify `partition` name in [`config`](./profiles/hpc/hpc.yaml) according to your cluster settings, then run:

```sh
snakemake -j 500 --profile profiles/hpc
```
