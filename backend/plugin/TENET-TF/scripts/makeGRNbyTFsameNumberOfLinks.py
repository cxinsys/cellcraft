import numpy as np
import sys
import os

species = sys.argv[1]
NumberOfLinks = int(sys.argv[2])
input_file_name = sys.argv[3]
output_file_name = sys.argv[4]

# 전사 인자(TFs) 목록 파일 경로 설정
tf_list_path = f"/workflow/script/TENET/GO_symbol_{species}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"
with open(tf_list_path, 'r') as file:
    TFlist = [line.strip() for line in file]

# 유전자 발현 결과 파일 경로 설정
with open(input_file_name, 'r') as file:
    gene_name = file.readline().strip().split()[1:]

    # 유전자 발현 데이터 처리
    cutOff = 0
    source, TE, target = [], [], []
    for sourceIndex, line in enumerate(file):
        if gene_name[sourceIndex] in TFlist:
            temp = line.strip().split()
            for targetIndex in range(len(temp) - 1):
                if float(temp[targetIndex + 1]) > cutOff:
                    source.append(gene_name[sourceIndex])
                    TE.append(float(temp[targetIndex + 1]))
                    target.append(gene_name[targetIndex])

# TE 값 정렬 및 상위 링크 선택
TEsortIndex = np.argsort(TE)
print(TEsortIndex)

# 결과 파일 경로 설정 및 파일 쓰기
with open(output_file_name, 'w') as file:
    if len(TEsortIndex) > 0:
        for i in range(min(NumberOfLinks, len(TEsortIndex))):
            idx = TEsortIndex[-i-1]
            file.write(f"{source[idx]}\t{TE[idx]}\t{target[idx]}\n")
    else:
        print("TEsortIndex 배열이 비어 있습니다. 빈 결과 파일이 생성됩니다.")

