from jpype import *
import random
import math
import numpy
import os
import csv
import sys
import datetime

# Change location of jar to match yours:
jarLocation = "/app/workflow/script/TENET/infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation,"-Xmx16G")

##ifile = open(os.environ.get('trajectory_file'))
ifile = open("/app/" + sys.argv[3])
trajectory1=[]
for line in ifile:
    trajectory1.append(float(line.split()[0]))

ifile = open("/app/" + sys.argv[4])
##ifile = open(os.environ.get('cell_select_file'))
branch=[]
for line in ifile:
    branch.append(int(line.split()[0]))
branch=numpy.array(branch)

historyLength = int(sys.argv[5])

trajectory1=numpy.array(trajectory1)
trajectory1=trajectory1[branch==1]
trajectory1SortIndex=numpy.argsort(trajectory1)

userPath = sys.argv[6]
cell_gene_all = numpy.genfromtxt(userPath + 'cell_gene_trsps.csv',delimiter=',')
#cell_gene_all = cell_gene_all.astype(int)

expression_data = []
input_fname = sys.argv[1]
#input_fname='pair_jobs/pair_list_aaa'
list_pairs = numpy.genfromtxt(input_fname,delimiter=',')
list_pairs = list_pairs.astype(int)
#cg_file = open('list_genefiles')
#lines = cg_file.readlines()
TEresult=[None] * len(list_pairs)
for num_pair in range(len(list_pairs)):
    for idx in [0, 1]:
        index = list_pairs[num_pair, idx] - 1
        if index < 0 or index >= len(cell_gene_all):
            print(f"인덱스 오류: {index}는 cell_gene_all 배열의 유효 범위를 벗어났습니다.")
            continue

        if idx == 0:
            expression_data = cell_gene_all[index][numpy.newaxis]
        else:
            expression_data = numpy.append(expression_data, cell_gene_all[index][numpy.newaxis], axis=0)

    expression_data1 = []
    for i in range(len(expression_data)):
        data_temp1 = numpy.array(expression_data[i])

        if len(data_temp1) != len(branch):
            print("불리언 인덱스 길이 오류: branch 배열의 길이가 data_temp1의 길이와 다릅니다.")
            continue

        data_temp1 = data_temp1[branch == 1]
        data_temp1 = data_temp1[trajectory1SortIndex]  # 이 부분도 동일한 인덱스 검증이 필요할 수 있음
        data_temp1 = list(data_temp1)
        expression_data1.append(data_temp1)

    expression_data = expression_data1
    del expression_data1

    if len(expression_data) < 2:
        print(f"오류: expression_data에 충분한 데이터가 없습니다. num_pair: {num_pair}")
        continue
# Create a TE calculator and run it:
    teCalcClass = JPackage("infodynamics.measures.continuous.kernel").TransferEntropyCalculatorKernel
    teCalc = teCalcClass()
    teCalc.setProperty("NORMALISE", "true") # Normalise the individual variables
    teCalc.initialise(historyLength, 0.5) # Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
    
    resultTemp=[]
    teCalc.setObservations(JArray(JDouble, 1)(expression_data[0]), JArray(JDouble, 1)(expression_data[1]))
    resultTemp.append(teCalc.computeAverageLocalOfObservations())
    
    teCalc.initialise(historyLength, 0.5) # Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
    teCalc.setObservations(JArray(JDouble, 1)(expression_data[1]), JArray(JDouble, 1)(expression_data[0]))
    resultTemp.append(teCalc.computeAverageLocalOfObservations())
    
    TEresult[num_pair] = numpy.ndarray.tolist(list_pairs[num_pair,:]) + resultTemp
    if (num_pair % int(len(list_pairs)/2)) == 0:
        print(datetime.datetime.now())

output_fname = sys.argv[2]
with open(output_fname, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(TEresult)