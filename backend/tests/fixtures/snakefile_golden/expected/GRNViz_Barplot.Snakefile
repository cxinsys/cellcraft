import os
rule Barplot:
    input:
        input="user/john/workflow_1/algorithm_2/results/Numlinks.sif"
    output:
        Barplot="user/john/workflow_1/visualization_7/results/Barplot_result.json"
    params:
        Top_Genes=10
    log:
        stdout="user/john/workflow_1/visualization_7/logs/Barplot.stdout",
        stderr="user/john/workflow_1/visualization_7/logs/Barplot.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/Barplot.py {input.input} {output.Barplot} {params.Top_Genes} > {log.stdout} 2> {log.stderr}"
