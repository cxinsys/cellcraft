import os
import sys

root_path = "/app/"

# 스크립트 인자 처리
user_result_path = sys.argv[1]
output_path = root_path + sys.argv[2]

# gene_names 파일 읽기
gene_name_file = os.path.join(user_result_path, "gene_names")
with open(gene_name_file, 'r') as file:
    gene_name = [line.strip() for line in file]

# TEmatrix 초기화
TEmatrix = [[0 for _ in range(len(gene_name))] for _ in range(len(gene_name))]

# TE_result_all.csv 파일 읽기 및 TEmatrix 업데이트
te_result_file = os.path.join(user_result_path, "TE_result_all.csv")
with open(te_result_file, 'r') as file:
    for line in file:
        temp = line.strip().split(",")
        TEmatrix[int(temp[0])-1][int(temp[1])-1] = float(temp[2])
        if len(temp) > 3:
            TEmatrix[int(temp[1])-1][int(temp[0])-1] = float(temp[3])

# 결과 파일 쓰기
output_file = output_path
with open(output_file, 'w') as file:
    file.write("TE" + "\t" + "\t".join(gene_name))
    for i in range(len(gene_name)):
        file.write("\n" + gene_name[i] + "\t" + "\t".join(map(str, TEmatrix[i])))
