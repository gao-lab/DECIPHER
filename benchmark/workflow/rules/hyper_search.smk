
rule search_decipher:
    input:
        data = "{path}/data.h5ad",
        script = "workflow/scripts/search_decipher.py",
    output:
        center_emb = "{path}/search_decipher/dropout_gex:{dropout_gex}-k:{k}-center_emb_dim:{center_emb_dim}-nbr_emb_dim:{nbr_emb_dim}-epochs:{epochs}-dropout:{dropout}-transformer_layers:{transformer_layers}-temperature_center:{temperature_center}-temperature_nbr:{temperature_nbr}-batch_size:{batch_size}/center_emb.npy",
        nbr_emb = "{path}/search_decipher/dropout_gex:{dropout_gex}-k:{k}-center_emb_dim:{center_emb_dim}-nbr_emb_dim:{nbr_emb_dim}-epochs:{epochs}-dropout:{dropout}-transformer_layers:{transformer_layers}-temperature_center:{temperature_center}-temperature_nbr:{temperature_nbr}-batch_size:{batch_size}/nbr_emb.npy",
        time = "{path}/search_decipher/dropout_gex:{dropout_gex}-k:{k}-center_emb_dim:{center_emb_dim}-nbr_emb_dim:{nbr_emb_dim}-epochs:{epochs}-dropout:{dropout}-transformer_layers:{transformer_layers}-temperature_center:{temperature_center}-temperature_nbr:{temperature_nbr}-batch_size:{batch_size}/time.yaml",
    params:
        notebook_in = "workflow/scripts/search_decipher.ipynb",
        notebook_out = "{path}/search_decipher/dropout_gex:{dropout_gex}-k:{k}-center_emb_dim:{center_emb_dim}-nbr_emb_dim:{nbr_emb_dim}-epochs:{epochs}-dropout:{dropout}-transformer_layers:{transformer_layers}-temperature_center:{temperature_center}-temperature_nbr:{temperature_nbr}-batch_size:{batch_size}/search_decipher.ipynb",
        work_dir = "{path}/search_decipher/dropout_gex:{dropout_gex}-k:{k}-center_emb_dim:{center_emb_dim}-nbr_emb_dim:{nbr_emb_dim}-epochs:{epochs}-dropout:{dropout}-transformer_layers:{transformer_layers}-temperature_center:{temperature_center}-temperature_nbr:{temperature_nbr}-batch_size:{batch_size}",
    log:
        "{path}/search_decipher/dropout_gex:{dropout_gex}-k:{k}-center_emb_dim:{center_emb_dim}-nbr_emb_dim:{nbr_emb_dim}-epochs:{epochs}-dropout:{dropout}-transformer_layers:{transformer_layers}-temperature_center:{temperature_center}-temperature_nbr:{temperature_nbr}-batch_size:{batch_size}/search_decipher.log"
    threads:8
    resources: gpu=1
    shell:
        """
        jupytext --to notebook {input.script}

        timeout {config[timeout]} papermill \
            -p input_file {input.data} \
            -p work_dir {params.work_dir} \
            -p time_file {output.time} \
            -p dropout_gex {wildcards.dropout_gex} \
            -p k {wildcards.k} \
            -p center_emb_dim {wildcards.center_emb_dim} \
            -p nbr_emb_dim {wildcards.nbr_emb_dim} \
            -p epochs {wildcards.epochs} \
            -p dropout {wildcards.dropout} \
            -p transformer_layers {wildcards.transformer_layers} \
            -p temperature_center {wildcards.temperature_center} \
            -p temperature_nbr {wildcards.temperature_nbr} \
            -p batch_size {wildcards.batch_size} \
        {params.notebook_in} {params.notebook_out} \
        > {log} 2>&1
        """
