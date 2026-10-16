import os

rule Scribe_Input:
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
        stdout="user/john/workflow_1/algorithm_2/logs/Scribe_Input.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/Scribe_Input.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/scribe_input.py {input.input} {params.geneList} {params.UMAP_lasso} {output.expression} {output.trajectory} {output.cellSelect} {params.cell_group} {params.pseudotime} '{params.clusters}' > {log.stdout} 2> {log.stderr}"

rule Run:
    input:
        expression="user/john/workflow_1/algorithm_2/results/expression.csv",
        trajectory="user/john/workflow_1/algorithm_2/results/trajectory.txt",
        cellSelect="user/john/workflow_1/algorithm_2/results/cellSelect.txt"
    output:
        rankedEdges="user/john/workflow_1/algorithm_2/results/rankedEdges.txt",
        GRN="user/john/workflow_1/algorithm_2/results/GRN.sif",
        GRNOutdegree="user/john/workflow_1/algorithm_2/results/GRNOutdegree.txt"
    params:
        Lower_Limit="0",
        Expression_Family="uninormal",
        method="ucRDI",
        Delay="1",
        log="false",
        number_of_outdegrees="1000",
        FDR="0"
    log:
        stdout="user/john/workflow_1/algorithm_2/logs/Run.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/Run.stderr"
    shell:
        "export R_HOME=/opt/micromamba/envs/plugin_env/lib/R && export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE && /opt/micromamba/envs/plugin_env/bin/Rscript /scripts/scribe_for_cellcraft.R {input.expression} {input.trajectory} {input.cellSelect} {params.Lower_Limit} {params.Expression_Family} {params.method} {params.Delay} {params.log} {params.number_of_outdegrees} {params.FDR} {output.rankedEdges} {output.GRN} {output.GRNOutdegree} > {log.stdout} 2> {log.stderr}"

