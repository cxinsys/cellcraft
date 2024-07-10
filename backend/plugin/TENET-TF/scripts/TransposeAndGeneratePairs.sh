#!/bin/bash

#/TransposeAndGeneratePairs.sh user/{userName}/result {input.csv_file} {output.cell_gene_tsv} {output.gene_names}

result_path=$1
input_csv=$2
output_cell_gene_tsv=$3
output_gene_names=$4

mkdir -p $result_path
cat $input_csv | cut -d ',' -f 2- | tail -n +2 | sed 's/,/ /g' > $output_cell_gene_tsv
cat $input_csv | head -n 1 | cut -d ',' -f 2- | tr ',' '\n' > $output_gene_names
