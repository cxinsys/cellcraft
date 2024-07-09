import csv
import numpy as np
import os
import sys

# 입력으로 받은 경로를 변수에 저장
user_result_path = sys.argv[1]

cell_gene_file = sys.argv[1]
cell_gene_trsps_file = sys.argv[2]
all_pairs_file = sys.argv[3]

# 파일 경로를 지정할 때 user_result_path를 사용
cell_gene_file = os.path.join(user_result_path, "cell_gene.tsv")
cell_gene_trsps_file = os.path.join(user_result_path, "cell_gene_trsps.csv")
all_pairs_file = os.path.join(user_result_path, "all_pairs.csv")

# cell_gene.tsv 파일 읽기
with open(cell_gene_file, "r") as f:
    reader = csv.reader(f, delimiter=" ")
    expression_data = np.array(list(reader)).astype("float")

# Transpose 수행 및 cell_gene_trsps.csv 파일 작성
expression_data = expression_data.T
with open(cell_gene_trsps_file, "w") as f:
    writer = csv.writer(f)
    writer.writerows(expression_data)

# 추가 계산 및 all_pairs.csv 파일 작성
num_gene = len(expression_data)
X = np.arange(num_gene-1, 1, -1)
X_1 = np.insert(X, 0, 1)
X_2 = -np.insert(X, 0, 0) + 1
bp = X_1.cumsum() - 1
indx_1 = np.zeros(shape=(bp[-1]+1, 1))
indx_2 = np.zeros(shape=(bp[-1]+1, 1))

for i in range(len(bp)):
    indx_1[bp[i]] = 1
    indx_2[bp[i]] = X_2[i]

indx_1 = indx_1.cumsum(axis=0)
indx_2 = indx_2 + 1
indx_2 = indx_2.cumsum(axis=0)
indx = np.concatenate((indx_1, indx_2), axis=1).astype(int)

with open(all_pairs_file, "w") as f:
    writer = csv.writer(f)
    writer.writerows(indx)