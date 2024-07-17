#!/bin/bash

#/MpirunScript.sh user/{userName}/result {input.list_jobfiles} {input.pseudotime} {input.cellSelect} {output.TE_result_all}

result_path=$1
input_list_jobfiles=$2
input_pseudotime=$3
input_cellSelect=$4
output_TE_result_all=$5

pair_jobs_path="${result_path}/pair_jobs"
outputs_path="${result_path}/outputs"

export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

if [ -d "$outputs_path" ]; then
    rm -rf "$outputs_path"
fi
mkdir "$outputs_path"

num_job=$(grep -c '^[^[:space:]]' $input_list_jobfiles)
mpirun_cmd='time mpirun'

for ((loop=1; loop<=${num_job}; loop++)); do
    input_file=$(cat $input_list_jobfiles | head -n $loop | tail -n 1)
    output_id=$(echo $input_file | awk -F'_' '{print $NF}')
    mpirun_cmd+=" -np 1 /opt/conda/envs/snakemake/bin/python /app/workflow/script/TENET/runTE_TF.py ${pair_jobs_path}/${input_file} ${outputs_path}/TE_out_${output_id}.csv $input_pseudotime $input_cellSelect {historyLength} $result_path :"
done
mpirun_cmd=$(echo "$mpirun_cmd" | sed 's/ :$//')
echo "$mpirun_cmd" > "${result_path}/mpirun_script.sh"
chmod a+x "${result_path}/mpirun_script.sh"
"${result_path}/mpirun_script.sh"
sleep 5

cat "$outputs_path/"*.csv > $output_TE_result_all
