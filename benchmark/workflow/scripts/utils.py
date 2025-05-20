from pathlib import Path

from snakemake.io import expand

DEBUG = False
IGNORE_FINISHED = False


def seed2range(seed: int):
    if seed > 0:
        return range(seed)
    else:
        return [0]


def filter_pool(pool) -> list[str]:
    pool_filter = []
    for file in pool:
        if 'stagate' in file and 'merfish' in file:
            continue
        
        if 'scniche_raw' in file and 'merfish' in file:
            continue

        if 'scniche_raw' in file and 'xenium' in file:
            continue

        if 'graphst' in file and 'merfish' in file:
            continue

        if 'graphst' in file and 'xenium' in file:
            continue

        if Path(file).exists():
            print(f"File {file} exists")
            if IGNORE_FINISHED:
                print("Skip!")
                continue
        pool_filter.append(file)

    if DEBUG:
        pool_filter = pool_filter[:1]
    return pool_filter


def target_benchmark_files(config: dict) -> list[str]:
    r"""
    Find benchmark results

    Parameters
    ----------
    config
        snakemake config
    """
    seeds = seed2range(config["seed"])

    pool = expand(
        "results/{dataset}/seed:{seed}/{method}/metrics.yaml",
        dataset=config["dataset"].keys(),
        seed=seeds,
        method=config["method"],
    )
    return filter_pool(pool)


# --------------------------- Hpyer-search -----------------------#
def conf_expand_pattern(conf, placeholder="null"):
    r"""
    expand the str by config, otherwise by default placeholder
    """
    expand_pattern = "-".join(f"{key}:{{{key}}}" for key in conf)
    return expand_pattern if expand_pattern else placeholder


def expand_with_default(pattern, **wildcards):
    r"""
    Extend snakemake function `expand()` to support "default" choices
    """
    has_default_choices = False
    for val in wildcards.values():  # Sanity check
        if isinstance(val, dict):
            if "default" not in val or "choices" not in val:
                print(val)
                raise ValueError("Invalid default choices!")
            has_default_choices = True

    if not has_default_choices:
        return expand(pattern, **wildcards)

    expand_set = set()
    for key, val in wildcards.items():
        if isinstance(val, dict):
            wildcards_use = {key: val["choices"]}
            for other_key, other_val in wildcards.items():
                if other_key == key:
                    continue
                if isinstance(other_val, dict):
                    wildcards_use[other_key] = other_val["default"]
                else:
                    wildcards_use[other_key] = other_val
            expand_set = expand_set.union(expand(pattern, **wildcards_use))
    return list(expand_set)


def target_hyper_search_files(config):
    r"""
    Find benchmark results

    Parameters
    ----------
    config
        snakemake config
    """
    seeds = seed2range(config["seed"])
    
    ## Search best hyper-param for all dataset
    # hyperparam_conf = config['hyperparam_conf']
    # hyperparam_conf = expand_with_default(
    #     conf_expand_pattern(hyperparam_conf, placeholder="default"), **hyperparam_conf
    # )
    # pool = expand(
    #     "results/{dataset}/seed:{seed}/search_decipher/{hyperparam_conf}/metrics.yaml",
    #     dataset=config["dataset"].keys(),
    #     seed=seeds,
    #     hyperparam_conf=hyperparam_conf,
    # )

    # Search best hyper-param for each dataset
    pool = []
    for dataset in config["dataset"].keys():
        for seed in seeds:
            config_dataset = config['hyperparam_conf'][dataset]
            print(f'config_dataset: {dataset}', config_dataset)
            hyperparam_conf_dataset = expand_with_default(
                conf_expand_pattern(config_dataset, placeholder="default"),
                **config_dataset
            )
            for conf in hyperparam_conf_dataset:
                pool.append(
                    f"results/{dataset}/seed:{seed}/search_decipher/{conf}/metrics.yaml"
                )
    return filter_pool(pool)
