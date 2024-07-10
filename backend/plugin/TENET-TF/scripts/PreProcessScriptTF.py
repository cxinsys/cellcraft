import csv
import os
import numpy as np
import sys

# 입력으로 받은 경로를 변수에 저장
species = sys.argv[1]
cell_gene_file = sys.argv[2]
gene_names_file = sys.argv[3]
cell_gene_trsps_file = sys.argv[4]
all_pairs_file = sys.argv[5]

tf_list_file = f"/app/workflow/script/TENET/GO_symbol_{species}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"

# cell_gene.tsv 파일 읽기 및 Transpose
with open(cell_gene_file, "r") as f:
    reader = csv.reader(f, delimiter=" ")
    expression_data = np.array(list(reader)).astype("float")

expression_data = expression_data.T
with open(cell_gene_trsps_file, "w") as f:
    writer = csv.writer(f)
    writer.writerows(expression_data)

# 유전자 이름 및 전사 인자(TFs) 목록 읽기
gene_names = []
with open(os.path.join(user_result_path, "gene_names"), "r") as f:
    gene_names = [line.strip() for line in f]

# print("gene_names: ", gene_names)

TFs = []
with open(tf_list_file, "r") as f:
    TFs = [line.strip() for line in f]

# print("TFs: ", TFs)

# 전사 인자(TFs) 인덱스 생성
num_gene = len(expression_data)
TFindex = [i+1 for i in range(num_gene) if gene_names[i] in TFs]

# print("TFindex: ", TFindex)

# all_pairs.csv 파일 작성
indx = [[i, j+1] for i in TFindex for j in range(num_gene) if i != j+1]
# print("indx: ", indx)

with open(all_pairs_file, "w") as f:
    writer = csv.writer(f)
    writer.writerows(indx)
