import os
rule Heatmap:
    input:
        expression="user/john/workflow_1/algorithm_2/results/expression.csv",
        trajectory="user/john/workflow_1/algorithm_2/results/trajectory.txt",
        input="user/john/workflow_1/algorithm_2/results/Numlinks.sif"
    output:
        Heatmap="user/john/workflow_1/visualization_7/results/Heatmap_result.json"
    params:
        Top_Genes=10,
        Sample_Size=100
    log:
        stdout="user/john/workflow_1/visualization_7/logs/Heatmap.stdout",
        stderr="user/john/workflow_1/visualization_7/logs/Heatmap.stderr"
    shell:
        "export R_HOME=/opt/micromamba/envs/plugin_env/lib/R && export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE && /opt/micromamba/envs/plugin_env/bin/Rscript /scripts/Pseudotime_heatmap.R {input.expression} {input.trajectory} {input.input} {output.Heatmap} {params.Top_Genes} {params.Sample_Size} > {log.stdout} 2> {log.stderr}"
