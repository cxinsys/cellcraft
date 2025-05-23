import numpy
import statsmodels.sandbox.stats.multicomp
import scipy.stats
import sys

fdrCutoff = float(sys.argv[1])
input_file_name = sys.argv[2]
output_file_name = sys.argv[3]

ifile = open(input_file_name)
line = ifile.readline()
temp = line.split()
gene_name=[]
for i in range(len(temp)-1):
    gene_name.append(temp[i+1])

cutOff=0
sourceIndex=0
TEnetwork=[]
source=[]
TE=[]
target=[]
for line in ifile:
    temp = line.split()
    for targetIndex in range(len(temp)-1):
        if float(temp[targetIndex+1])>cutOff:            
            source.append(gene_name[sourceIndex])
            TE.append(float(temp[targetIndex+1]))
            target.append(gene_name[targetIndex])
    sourceIndex=sourceIndex+1
ifile.close()

TEzscore=(TE-numpy.mean(TE))/numpy.std(TE)
TEpvalue=1-scipy.stats.norm.cdf(TEzscore)
TEfdr=statsmodels.sandbox.stats.multicomp.multipletests(TEpvalue,alpha=0.05,method='fdr_bh')

ofile = open(output_file_name,"w")
for i in range(len(source)):
    if TEfdr[1][i]<fdrCutoff:
        ofile.write(source[i]+"\t"+str(TE[i])+"\t"+target[i]+"\n")
ofile.close()