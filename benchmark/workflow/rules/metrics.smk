
rule calc_metrics:
    input:
        center_emb = "{path}/center_emb.npy",
        nbr_emb = "{path}/nbr_emb.npy",
        time = "{path}/time.yaml",
        script = "workflow/scripts/calc_metrics.py",
    output:
        metrics = "{path}/metrics.yaml"
    params:
        notebook_in = "workflow/scripts/calc_metrics.ipynb",
        notebook_out = "{path}/calc_metrics.ipynb",
    log:
        "{path}/calc_metrics.log"
    threads:28
    shell:
        """
        jupytext --to notebook {input.script}

        timeout {config[timeout]} papermill \
            -p center_emb {input.center_emb} \
            -p nbr_emb {input.nbr_emb} \
            -p time {input.time} \
            -p metrics_file {output.metrics} \
        {params.notebook_in} {params.notebook_out} \
        > {log} 2>&1
        """
