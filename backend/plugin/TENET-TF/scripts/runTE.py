from jpype import *
import numpy
import os
import csv
import sys
import datetime
import pandas as pd

# Change location of jar to match yours:
jarLocation = "/app/workflow/script/TENET/infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation, "-Xmx16G")

# Parse command line arguments
input_list_jobfiles = sys.argv[1]
trajectory_file = sys.argv[2]
cell_select_file = sys.argv[3]
cell_gene_trsps = sys.argv[4]
historyLength = int(sys.argv[5])
te_result_all = sys.argv[6]

# Load trajectory and cell select data
trajectory1 = []
with open(trajectory_file) as f:
    for line in f:
        trajectory1.append(float(line.split()[0]))

branch = []
with open(cell_select_file) as f:
    for line in f:
        branch.append(int(line.split()[0]))
branch = numpy.array(branch)

trajectory1 = numpy.array(trajectory1)
trajectory1 = trajectory1[branch == 1]
trajectory1SortIndex = numpy.argsort(trajectory1)

cell_gene_all = pd.read_csv(cell_gene_trsps, delimiter=',', header=None).to_numpy()

TEresult = []

# Read the list of job files and process each pair file
with open(input_list_jobfiles) as f:
    for job_file in f:
        job_file = job_file.strip()
        pair_file = os.path.join(os.path.dirname(input_list_jobfiles), "pair_jobs", job_file)
        list_pairs = pd.read_csv(pair_file, delimiter=',', header=None).to_numpy().astype(int)

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
                data_temp1 = data_temp1[trajectory1SortIndex]
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
            teCalc.setProperty("NORMALISE", "true")
            teCalc.initialise(historyLength, 0.5)

            resultTemp = []
            teCalc.setObservations(JArray(JDouble, 1)(expression_data[0]), JArray(JDouble, 1)(expression_data[1]))
            resultTemp.append(teCalc.computeAverageLocalOfObservations())

            teCalc.initialise(historyLength, 0.5)
            teCalc.setObservations(JArray(JDouble, 1)(expression_data[1]), JArray(JDouble, 1)(expression_data[0]))
            resultTemp.append(teCalc.computeAverageLocalOfObservations())

            TEresult.append(numpy.ndarray.tolist(list_pairs[num_pair, :]) + resultTemp)
            if (num_pair % int(len(list_pairs) / 2)) == 0:
                print(datetime.datetime.now())

# Write the combined results to the output file
with open(te_result_all, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(TEresult)
