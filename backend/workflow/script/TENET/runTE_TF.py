from jpype import *
import numpy as np
import os
import csv
import sys
import datetime

# JVM 시작 및 jar 파일 로드
jarLocation = "/app/workflow/script/TENET/infodynamics.jar"
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation, "-Xmx16G")

# 데이터 파일 로딩
trajectory_file = "/app/" + sys.argv[3]
cell_select_file = "/app/" + sys.argv[4]
trajectory1 = np.loadtxt(trajectory_file, usecols=[0])
branch = np.loadtxt(cell_select_file, usecols=[0], dtype=int)

historyLength = int(sys.argv[5])
user_result_path = sys.argv[6]

# 필터링 및 정렬
trajectory1 = trajectory1[branch == 1]
trajectory1SortIndex = np.argsort(trajectory1)

# 유전자 발현 데이터 로딩
cell_gene_all = np.genfromtxt(os.path.join(user_result_path, 'cell_gene_trsps.csv'), delimiter=',')

# 유전자 쌍 처리
input_fname = sys.argv[1]
list_pairs = np.genfromtxt(input_fname, delimiter=',', dtype=int)
TEresult = [None] * len(list_pairs)

for num_pair in range(len(list_pairs)):
    expression_data = [cell_gene_all[idx - 1][branch == 1][trajectory1SortIndex] for idx in list_pairs[num_pair]]
    
    # 전달 엔트로피 계산
    teCalcClass = JPackage("infodynamics.measures.continuous.kernel").TransferEntropyCalculatorKernel
    teCalc = teCalcClass()
    teCalc.setProperty("NORMALISE", "true")
    teCalc.initialise(historyLength, 0.5)
    teCalc.setObservations(JArray(JDouble, 1)(expression_data[0]), JArray(JDouble, 1)(expression_data[1]))
    resultTemp = [teCalc.computeAverageLocalOfObservations()]
    
    TEresult[num_pair] = np.ndarray.tolist(list_pairs[num_pair]) + resultTemp
    if (num_pair % int(len(list_pairs)/2)) == 0:
        print(datetime.datetime.now())

# 결과 출력
output_fname = sys.argv[2]
with open(output_fname, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(TEresult)
