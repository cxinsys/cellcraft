import numpy
import sys

# 커맨드 라인 인자 받기
input_file_1 = sys.argv[1]
input_file_2 = sys.argv[2]
input_file_3 = sys.argv[3]
input_file_4 = sys.argv[4]
output_file_1 = sys.argv[5]
output_file_2 = sys.argv[6]
output_file_3 = sys.argv[7]
output_file_4 = sys.argv[8]

def process_file(input_file, output_file):
    ifile = open(input_file)
    TFlist = []
    TFlistOutdegree = []
    
    for line in ifile:
        temp = line.split()
        if temp[0] not in TFlist:
            TFlist.append(temp[0])
            TFlistOutdegree.append(1)
        else:
            TFlistOutdegree[TFlist.index(temp[0])] = TFlistOutdegree[TFlist.index(temp[0])] + 1
    
    TFlistOutdegreeIndex = numpy.argsort(TFlistOutdegree)
    
    ofile = open(output_file, "w")
    for i in range(len(TFlist)):
        ofile.write(TFlist[TFlistOutdegreeIndex[-i-1]] + "\t" + str(TFlistOutdegree[TFlistOutdegreeIndex[-i-1]]) + "\n")
    ofile.close()
    ifile.close()

# 두 파일 각각 처리
process_file(input_file_1, output_file_1)
process_file(input_file_2, output_file_2)
process_file(input_file_3, output_file_3)
process_file(input_file_4, output_file_4)
