#!/bin/bash

#/SplitPairs.sh user/{userName}/result {input.all_pairs} {output.list_jobfiles} {numOfThreads}

result_path=$1
input_all_pairs=$2
output_list_jobfiles=$3
num_jobs=$4

pair_jobs_path="${result_path}/pair_jobs"
if [ -d "$pair_jobs_path" ]; then
    rm -rf "$pair_jobs_path"
fi
mkdir "$pair_jobs_path"
mv $input_all_pairs "$pair_jobs_path/all_pairs.csv"

num_pair=$(cat "$pair_jobs_path/all_pairs.csv" | wc -l | sed -e 's/^[ \t]*//')
num_line=$(expr $(expr ${num_pair} / ${num_jobs}) + 1)
split -a 3 -l ${num_line} "$pair_jobs_path/all_pairs.csv" "$pair_jobs_path/pair_list_"
rm -f "$pair_jobs_path/all_pairs.csv"
ls -1 "$pair_jobs_path/" > $output_list_jobfiles
