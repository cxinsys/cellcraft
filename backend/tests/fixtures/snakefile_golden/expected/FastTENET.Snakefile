import os

rule FastTENET_Input:
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
        stdout="user/john/workflow_1/algorithm_2/logs/FastTENET_Input.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/FastTENET_Input.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/fasttenet_input.py {input.input} {params.geneList} {params.UMAP_lasso} {output.expression} {output.trajectory} {output.cellSelect} {params.cell_group} {params.pseudotime} '{params.clusters}' > {log.stdout} 2> {log.stderr}"

rule Run:
    input:
        expression="user/john/workflow_1/algorithm_2/results/expression.csv",
        trajectory="user/john/workflow_1/algorithm_2/results/trajectory.txt",
        cellSelect="user/john/workflow_1/algorithm_2/results/cellSelect.txt"
    output:
        resultMatrixFDR="user/john/workflow_1/algorithm_2/results/resultMatrixFDR.sif",
        resultMatrixFdrTrim="user/john/workflow_1/algorithm_2/results/resultMatrixFdrTrim.sif",
        resultMatrixFdrOutdegree="user/john/workflow_1/algorithm_2/results/resultMatrixFdrOutdegree.txt",
        resultMatrixFdrTrimOutdegree="user/john/workflow_1/algorithm_2/results/resultMatrixFdrTrimOutdegree.txt"
    params:
        TF=lambda wildcards: "tutorials/human_TFs.txt" if os.path.exists("tutorials/human_TFs.txt") else "None",
        backend="gpu",
        number_of_devices="1",
        batch_size="32768",
        FDR="0.01",
        number_of_outdegrees="0",
        trimming_threshold="0"
    log:
        stdout="user/john/workflow_1/algorithm_2/logs/Run.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/Run.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/fasttenet_for_cellcraft.py {input.expression} {input.trajectory} {input.cellSelect} {params.TF} {output.resultMatrixFDR} {output.resultMatrixFdrTrim} {output.resultMatrixFdrOutdegree} {output.resultMatrixFdrTrimOutdegree} {params.backend} {params.number_of_devices} {params.batch_size} {params.FDR} {params.number_of_outdegrees} {params.trimming_threshold} > {log.stdout} 2> {log.stderr}"

