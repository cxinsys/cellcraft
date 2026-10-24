import os

rule FastSCODE_Input:
    input:
        input="user/john/data/pbmc_3k.h5ad"
    output:
        expression="user/john/workflow_1/algorithm_2/results/expression.csv",
        trajectory="user/john/workflow_1/algorithm_2/results/trajectory.txt",
        cellSelect="user/john/workflow_1/algorithm_2/results/cellSelect.txt"
    params:
        geneList=lambda wildcards: "user/john/data/genes.txt" if os.path.exists("user/john/data/genes.txt") else "None",
        UMAP_lasso=lambda wildcards: "./user/john/workflow_1/ScatterPlot_4.json" if os.path.exists("./user/john/workflow_1/ScatterPlot_4.json") else "None",
        cell_group="leiden",
        pseudotime="dpt_pseudotime",
        clusters=lambda wc: ";".join(["0", "1", "B cell"])
    log:
        stdout="user/john/workflow_1/algorithm_2/logs/FastSCODE_Input.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/FastSCODE_Input.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/fastscode_input.py {input.input} {params.geneList} {params.UMAP_lasso} {output.expression} {output.trajectory} {output.cellSelect} {params.cell_group} {params.pseudotime} '{params.clusters}' > {log.stdout} 2> {log.stderr}"

rule Run:
    input:
        expression="user/john/workflow_1/algorithm_2/results/expression.csv",
        trajectory="user/john/workflow_1/algorithm_2/results/trajectory.txt",
        cellSelect="user/john/workflow_1/algorithm_2/results/cellSelect.txt"
    output:
        GRN="user/john/workflow_1/algorithm_2/results/GRN.sif",
        GRNOutdegree="user/john/workflow_1/algorithm_2/results/GRNOutdegree.txt",
        GRNTrim="user/john/workflow_1/algorithm_2/results/GRNTrim.sif",
        GRNTrimOutdegree="user/john/workflow_1/algorithm_2/results/GRNTrimOutdegree.txt"
    params:
        sampling_length="10",
        max_iteration="100",
        backend="gpu",
        number_of_devices="1",
        sampling_batch_size="1",
        gene_batch_size="0",
        num_repeat="1",
        number_of_outdegrees="1000",
        trim_threshold="0.01"
    log:
        stdout="user/john/workflow_1/algorithm_2/logs/Run.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/Run.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/scode_for_cellcraft.py {input.expression} {input.trajectory} {input.cellSelect} {params.sampling_length} {params.max_iteration} {params.backend} {params.number_of_devices} {params.sampling_batch_size} {params.gene_batch_size} {params.num_repeat} {params.number_of_outdegrees} {params.trim_threshold} {output.GRN} {output.GRNOutdegree} {output.GRNTrim} {output.GRNTrimOutdegree} > {log.stdout} 2> {log.stderr}"

