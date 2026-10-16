import os
rule Network:
    input:
        input="user/john/workflow_1/algorithm_2/results/Numlinks.sif"
    output:
        Network="user/john/workflow_1/visualization_7/results/Network_result.json"
    params:
        Top_Genes=50
    log:
        stdout="user/john/workflow_1/visualization_7/logs/Network.stdout",
        stderr="user/john/workflow_1/visualization_7/logs/Network.stderr"
    shell:
        "/opt/micromamba/envs/plugin_env/bin/python /scripts/NetworkX.py {input.input} {output.Network} {params.Top_Genes} > {log.stdout} 2> {log.stderr}"