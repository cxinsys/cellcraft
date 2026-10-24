import os

rule TENET_Input:
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
        stdout="user/john/workflow_1/algorithm_2/logs/TENET_Input.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/TENET_Input.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/tenet_input.py {input.input} {params.geneList} {params.UMAP_lasso} {output.expression} {output.trajectory} {output.cellSelect} {params.cell_group} {params.pseudotime} '{params.clusters}' > {log.stdout} 2> {log.stderr}"

rule Run:
    input:
        expression="user/john/workflow_1/algorithm_2/results/expression.csv",
        trajectory="user/john/workflow_1/algorithm_2/results/trajectory.txt",
        cellSelect="user/john/workflow_1/algorithm_2/results/cellSelect.txt"
    output:
        resultMatrix="user/john/workflow_1/algorithm_2/results/resultMatrix.txt"
    params:
        TF=lambda wildcards: "tutorials/human_TFs.txt" if os.path.exists("tutorials/human_TFs.txt") else "None",
        number_of_devices="4"
    log:
        stdout="user/john/workflow_1/algorithm_2/logs/Run.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/Run.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/tenet_for_cellcraft.py {input.expression} {input.trajectory} {input.cellSelect} {params.TF} {output.resultMatrix} {params.number_of_devices} > {log.stdout} 2> {log.stderr}"

rule GRN_Reconstruction__FDR:
    input:
        resultMatrix="user/john/workflow_1/algorithm_2/results/resultMatrix.txt"
    output:
        Fdr="user/john/workflow_1/algorithm_2/results/Fdr.sif"
    params:
        FDR="0.01"
    log:
        stdout="user/john/workflow_1/algorithm_2/logs/GRN_Reconstruction__FDR.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/GRN_Reconstruction__FDR.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/makeGRN.py {params.FDR} {input.resultMatrix} {output.Fdr} > {log.stdout} 2> {log.stderr}"

rule GRN_Reconstruction__NumLinks:
    input:
        resultMatrix="user/john/workflow_1/algorithm_2/results/resultMatrix.txt"
    output:
        Numlinks="user/john/workflow_1/algorithm_2/results/Numlinks.sif"
    params:
        number_of_outdegrees="1000"
    log:
        stdout="user/john/workflow_1/algorithm_2/logs/GRN_Reconstruction__NumLinks.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/GRN_Reconstruction__NumLinks.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/makeGRNsameNumberOfLinks.py {params.number_of_outdegrees} {input.resultMatrix} {output.Numlinks} > {log.stdout} 2> {log.stderr}"

rule Indirect_Edge_Trimming:
    input:
        Fdr="user/john/workflow_1/algorithm_2/results/Fdr.sif",
        Numlinks="user/john/workflow_1/algorithm_2/results/Numlinks.sif"
    output:
        FdrTrim="user/john/workflow_1/algorithm_2/results/FdrTrim.sif",
        NumlinksTrim="user/john/workflow_1/algorithm_2/results/NumlinksTrim.sif"
    params:
        trimming_threshold="0"
    log:
        stdout="user/john/workflow_1/algorithm_2/logs/Indirect_Edge_Trimming.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/Indirect_Edge_Trimming.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/trim_indirect.py {input.Fdr} {input.Numlinks} {output.FdrTrim} {output.NumlinksTrim} {params.trimming_threshold} > {log.stdout} 2> {log.stderr}"

rule Counting_Outdegree:
    input:
        Fdr="user/john/workflow_1/algorithm_2/results/Fdr.sif",
        Numlinks="user/john/workflow_1/algorithm_2/results/Numlinks.sif",
        FdrTrim="user/john/workflow_1/algorithm_2/results/FdrTrim.sif",
        NumlinksTrim="user/john/workflow_1/algorithm_2/results/NumlinksTrim.sif"
    output:
        FdrOutdegree="user/john/workflow_1/algorithm_2/results/FdrOutdegree.txt",
        NumlinksOutdegree="user/john/workflow_1/algorithm_2/results/NumlinksOutdegree.txt",
        FdrTrimOutdegree="user/john/workflow_1/algorithm_2/results/FdrTrimOutdegree.txt",
        NumlinksTrimOutdegree="user/john/workflow_1/algorithm_2/results/NumlinksTrimOutdegree.txt"
    log:
        stdout="user/john/workflow_1/algorithm_2/logs/Counting_Outdegree.stdout",
        stderr="user/john/workflow_1/algorithm_2/logs/Counting_Outdegree.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/countOutdegree.py {input.Fdr} {input.Numlinks} {input.FdrTrim} {input.NumlinksTrim} {output.FdrOutdegree} {output.NumlinksOutdegree} {output.FdrTrimOutdegree} {output.NumlinksTrimOutdegree} > {log.stdout} 2> {log.stderr}"

