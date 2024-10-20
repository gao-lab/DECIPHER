# Experiments

All datasets used in experiments are publicly avaliable, please see our manuscript for more details.

Each dir contains single experiment, dir name follow: `{species}_{organ}_{dataset}` style.
Inside each dir, there are series python files, for example:

1. `1-prepare_data.py`: get the data and do necessary filter/pre-process
2. `2-run.py`: run the model
3. `3-analysis_results.py`: analysis the model output or other customized analysis

The input files are expected to be stored in `/data` subdir. And results will be stored in `/results` subdir.