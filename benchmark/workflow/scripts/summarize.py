from argparse import ArgumentParser, Namespace

import pandas as pd
import parse
import yaml
from addict import Dict
from snakemake.script import Snakemake


def parse_args() -> Namespace:
    parser = ArgumentParser()
    parser.add_argument("--input", type=str, nargs="+")
    parser.add_argument("--pattern", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()
    args.params = Namespace(pattern=args.pattern)
    args.output = [args.output]
    del args.pattern
    return args


def debug() -> Dict:
    args = Dict()
    args.input = ["results/human_pbmc_10x_mimic/seed:0/search_spider/dropout_gex:0.5-k:40-center_emb_dim:128-nbr_emb_dim:32-epochs:6-dropout:0.1-transformer_layers:3-temperature_center:0.07-temperature_nbr:0.07/metrics.yaml"]
    args.params.pattern = "results/{dataset}/seed:{seed}/search_spider/dropout_gex:{dropout_gex}-k:{k}-center_emb_dim:{center_emb_dim}-nbr_emb_dim:{nbr_emb_dim}-epochs:{epochs}-dropout:{dropout}-transformer_layers:{transformer_layers}-temperature_center:{temperature_center}-temperature_nbr:{temperature_nbr}/metrics.yaml"
    args.output = ""
    return args


def process_column(column):
    if pd.api.types.is_numeric_dtype(column):
        return column.max()
    elif pd.api.types.is_string_dtype(column):
        return column.unique().tolist()[0]
    else:
        return None


def main(snakemake: Namespace | Snakemake) -> None:
    df_list = []
    for item in set(snakemake.input):
        entry = parse.parse(snakemake.params.pattern, item)
        if entry:
            conf = entry.named
        else:
            continue
        with open(item) as f:
            result_dict = yaml.load(f, Loader=yaml.Loader)
            df = pd.concat({k: pd.DataFrame([v]) for k, v in result_dict.items()}, axis=0)
            df = df.apply(process_column)
            df = pd.DataFrame(df).T
            # df['resolution'] = df.index.get_level_values(0)
            # split resolution by : and keep -1 only
            # df['resolution'] = df['resolution'].str.split(":").str[-1]
            for k, v in conf.items():
                df[k] = v
        df_list.append(df)
    df_all = pd.concat(df_list, axis=0)
    df_all.to_csv(snakemake.output[0], index=False)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = debug()
    main(snakemake)
