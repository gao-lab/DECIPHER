
rule link_data:
    input:
        h5ad_file = lambda wildcards: config["dataset"][wildcards.dataset]["data"],
    output:
        h5ad_file = "results/{dataset}/seed:{seed}/data.h5ad",
    log:
        "results/{dataset}/seed:{seed}/link_data.log"
    threads: 1
    shell:
        """
        ln -frs {input.h5ad_file} {output.h5ad_file} > {log} 2>&1
        """
