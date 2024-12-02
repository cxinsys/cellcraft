import numpy as np
import statsmodels.sandbox.stats.multicomp
import scipy.stats
import sys

input_file_name = sys.argv[1]
tf_list_path = sys.argv[2]
output_file_name = sys.argv[3]
fdrCutoff = float(sys.argv[4])

try:
    with open(tf_list_path, 'r') as file:
        TFlist = [line.strip() for line in file]
except FileNotFoundError:
    print(f"경고: {tf_list_path} 파일을 찾을 수 없습니다. 빈 TF 리스트로 진행합니다.")
    TFlist = []

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

# 통계적 처리
if len(TE) > 0:
    TEzscore = (TE - np.mean(TE)) / np.std(TE)
    TEpvalue = 1 - scipy.stats.norm.cdf(TEzscore)
    # TEpvalue 배열이 비어 있지 않고, 모든 값이 같지 않은지 확인
    if len(TEpvalue) > 0 and np.std(TEpvalue) > 0:
        TEfdr = statsmodels.sandbox.stats.multicomp.multipletests(TEpvalue, alpha=0.05, method='fdr_bh')
    else:
        TEfdr = ([], [], [], [])
else:
    TEzscore = []  # 또는 적절한 기본값으로 설정
    TEfdr = ([], [], [], [])

# 결과 파일 경로 설정 및 파일 쓰기
with open(output_file_name, 'w') as file:
    # TEfdr[1]에 유효한 데이터가 있는지 확인
    if len(TEfdr[1]) > 0:
        for i in range(len(source)):
            if TEfdr[1][i] < fdrCutoff:
                file.write(f"{source[i]}\t{TE[i]}\t{target[i]}\n")
