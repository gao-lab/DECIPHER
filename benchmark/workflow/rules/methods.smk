
rule run_decipher:
    input:
        data = "{path}/data.h5ad",
        script = "workflow/scripts/run_decipher.py",
    output:
        center_emb = "{path}/decipher/center_emb.npy",
        nbr_emb = "{path}/decipher/nbr_emb.npy",
        time = "{path}/decipher/time.yaml",
    params:
        notebook_in = "workflow/scripts/run_decipher.ipynb",
        notebook_out = "{path}/decipher/run_decipher.ipynb",
        work_dir = "{path}/decipher",
    log:
        "{path}/decipher/run_decipher.log",
    threads:8
    resources: gpu=1
    shell:
        """
        jupytext --to notebook {input.script}

        timeout {config[timeout]} papermill \
            -p input_file {input.data} \
            -p work_dir {params.work_dir} \
            -p time_file {output.time} \
        {params.notebook_in} {params.notebook_out} \
        > {log} 2>&1
        """


rule run_banksy:
    input:
        data = "{path}/data.h5ad",
        script = "workflow/scripts/run_banksy.py",
    output:
        center_emb = "{path}/banksy/center_emb.npy",
        nbr_emb = "{path}/banksy/nbr_emb.npy",
        time = "{path}/banksy/time.yaml",
    params:
        notebook_in = "workflow/scripts/run_banksy.ipynb",
        notebook_out = "{path}/banksy/run_banksy.ipynb",
        work_dir = "{path}/banksy",
    log:
        "{path}/banksy/run_banksy.log",
    threads:8
    shell:
        """
        jupytext --to notebook {input.script}

        timeout {config[timeout]} papermill \
            -p input_file {input.data} \
            -p center_emb_file {output.center_emb} \
            -p nbr_emb_file {output.nbr_emb} \
            -p time_file {output.time} \
        {params.notebook_in} {params.notebook_out} \
        > {log} 2>&1
        """


rule run_stagate:
    input:
        data = "{path}/data.h5ad",
        script = "workflow/scripts/run_stagate.py",
    output:
        center_emb = "{path}/stagate/center_emb.npy",
        nbr_emb = "{path}/stagate/nbr_emb.npy",
        time = "{path}/stagate/time.yaml",
    params:
        notebook_in = "workflow/scripts/run_stagate.ipynb",
        notebook_out = "{path}/stagate/run_stagate.ipynb",
    log:
        "{path}/stagate/run_stagate.log"
    threads:8
    resources: gpu=1
    shell:
        """
        jupytext --to notebook {input.script}

        timeout {config[timeout]} papermill \
            -p input_file {input.data} \
            -p center_emb_file {output.center_emb} \
            -p nbr_emb_file {output.nbr_emb} \
            -p time_file {output.time} \
        {params.notebook_in} {params.notebook_out} \
        > {log} 2>&1
        """


rule run_slat:
    input:
        data = "{path}/data.h5ad",
        script = "workflow/scripts/run_slat.py",
    output:
        center_emb = "{path}/slat/center_emb.npy",
        nbr_emb = "{path}/slat/nbr_emb.npy",
        time = "{path}/slat/time.yaml",
    params:
        notebook_in = "workflow/scripts/run_slat.ipynb",
        notebook_out = "{path}/slat/run_slat.ipynb",
        work_dir = "{path}/slat",
    log:
        "{path}/slat/run_slat.log"
    threads:8
    shell:
        """
        jupytext --to notebook {input.script}

        timeout {config[timeout]} papermill \
            -p input_file {input.data} \
            -p center_emb_file {output.center_emb} \
            -p nbr_emb_file {output.nbr_emb} \
            -p time_file {output.time} \
        {params.notebook_in} {params.notebook_out} \
        > {log} 2>&1
        """


rule run_harmony:
    input:
        data = "{path}/data.h5ad",
        script = "workflow/scripts/run_harmony.py",
    output:
        center_emb = "{path}/harmony/center_emb.npy",
        nbr_emb = "{path}/harmony/nbr_emb.npy",
        time = "{path}/harmony/time.yaml",
    params:
        notebook_in = "workflow/scripts/run_harmony.ipynb",
        notebook_out = "{path}/harmony/run_harmony.ipynb",
        work_dir = "{path}/harmony",
    log:
        "{path}/harmony/run_harmony.log"
    threads:8
    shell:
        """
        jupytext --to notebook {input.script}

        timeout {config[timeout]} papermill \
            -p input_file {input.data} \
            -p center_emb_file {output.center_emb} \
            -p nbr_emb_file {output.nbr_emb} \
            -p time_file {output.time} \
        {params.notebook_in} {params.notebook_out} \
        > {log} 2>&1
        """


rule run_scvi:
    input:
        data = "{path}/data.h5ad",
        script = "workflow/scripts/run_scvi.py",
    output:
        center_emb = "{path}/scvi/center_emb.npy",
        nbr_emb = "{path}/scvi/nbr_emb.npy",
        time = "{path}/scvi/time.yaml",
    params:
        notebook_in = "workflow/scripts/run_scvi.ipynb",
        notebook_out = "{path}/scvi/run_scvi.ipynb",
        work_dir = "{path}/scvi",
    log:
        "{path}/scvi/run_scvi.log"
    threads:8
    resources: gpu=1
    shell:
        """
        jupytext --to notebook {input.script}

        timeout {config[timeout]} papermill \
            -p input_file {input.data} \
            -p center_emb_file {output.center_emb} \
            -p nbr_emb_file {output.nbr_emb} \
            -p time_file {output.time} \
        {params.notebook_in} {params.notebook_out} \
        > {log} 2>&1
        """


rule run_scanpy:
    input:
        data = "{path}/data.h5ad",
        script = "workflow/scripts/run_scanpy.py",
    output:
        center_emb = "{path}/scanpy/center_emb.npy",
        nbr_emb = "{path}/scanpy/nbr_emb.npy",
        time = "{path}/scanpy/time.yaml",
    params:
        notebook_in = "workflow/scripts/run_scanpy.ipynb",
        notebook_out = "{path}/scanpy/run_scanpy.ipynb",
        work_dir = "{path}/scanpy",
    log:
        "{path}/scanpy/run_scanpy.log"
    threads:8
    shell:
        """
        jupytext --to notebook {input.script}

        timeout {config[timeout]} papermill \
            -p input_file {input.data} \
            -p center_emb_file {output.center_emb} \
            -p nbr_emb_file {output.nbr_emb} \
            -p time_file {output.time} \
        {params.notebook_in} {params.notebook_out} \
        > {log} 2>&1
        """