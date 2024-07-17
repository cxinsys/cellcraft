import numpy
import os, glob
import sys

file_name=sys.argv[1]
ifile = open(file_name)
TFlist=[];TFlistOutdegree=[]
for line in ifile:
    temp = line.split()
    if temp[0] not in TFlist:
        TFlist.append(temp[0])
        TFlistOutdegree.append(1)
    else:
        TFlistOutdegree[TFlist.index(temp[0])]=TFlistOutdegree[TFlist.index(temp[0])]+1
TFlistOutdegreeIndex=numpy.argsort(TFlistOutdegree)

output_path = sys.argv[2]
ofile = open(output_path,"w")
for i in range(len(TFlist)):
    ofile.write(TFlist[TFlistOutdegreeIndex[-i-1]]+"\t"+str(TFlistOutdegree[TFlistOutdegreeIndex[-i-1]])+"\n")
ofile.close()
