#!/bin/bash

#/TransposeAndGeneratePairs.sh user/{userName}/result {input.expMatrix} {output.cell_gene_tsv} {output.gene_names}

input_expMatrix=$1
output_cell_gene_tsv=$2
output_gene_names=$3

# 입력 매트릭스를 cell*gene에서 gene*cell로 전치하고 모든 유전자 쌍 목록을 생성
cat $input_expMatrix | cut -d ',' -f 2- | tail -n +2 | sed 's/,/ /g' > $output_cell_gene_tsv
cat $input_expMatrix | head -n 1 | cut -d ',' -f 2- | tr ',' '\n' > $output_gene_names
