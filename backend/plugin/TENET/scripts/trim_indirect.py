import sys

def process_file(input_file, output_file, threshold):
    ifile = open(input_file)
    TF=[];TFtarget=[];TFtargetTE=[];TFtargetIndirect=[]
    
    # ... 기존 파일 처리 로직 ...
    for line in ifile:
        temp = line.split()
        if temp[0] not in TF:
            TF.append(temp[0])
            TFtarget.append([temp[2]])
            TFtargetTE.append([float(temp[1])])
            TFtargetIndirect.append([0])
        else:
            indexTemp=TF.index(temp[0])
            TFtarget[indexTemp].append(temp[2])
            TFtargetTE[indexTemp].append(float(temp[1]))
            TFtargetIndirect[indexTemp].append(0)

    # ... 간접 관계 처리 로직 ...
    for i in range(len(TF)):
        for j in range(len(TFtarget[i])):
            for k in range(len(TFtarget[i])):
                if j!=k and TFtarget[i][j] in TF:
                    indexTemp=TF.index(TFtarget[i][j])
                    if TFtarget[i][k] in TFtarget[indexTemp]:
                        if TFtargetTE[i][k]<min(TFtargetTE[i][j],TFtargetTE[indexTemp][TFtarget[indexTemp].index(TFtarget[i][k])])+threshold:
                            TFtargetIndirect[i][k]=1

    # 결과 파일 작성
    ofile = open(output_file,"w")
    for i in range(len(TF)):
        for j in range(len(TFtarget[i])):
            if TFtargetIndirect[i][j]==0:
                ofile.write(TF[i]+"\t"+str(TFtargetTE[i][j])+"\t"+TFtarget[i][j]+"\n")
    ofile.close()
    ifile.close()

# 메인 실행 부분
input_file_1 = sys.argv[1]
input_file_2 = sys.argv[2]
output_file_1 = sys.argv[3]
output_file_2 = sys.argv[4]
cutoff = float(sys.argv[5])

# 각 파일 쌍에 대해 처리 실행
process_file(input_file_1, output_file_1, cutoff)
process_file(input_file_2, output_file_2, cutoff)
