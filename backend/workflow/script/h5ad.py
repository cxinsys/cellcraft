import numpy as np
import pandas as pd
import scanpy as sc # pip install scanpy
import matplotlib.pyplot as plt

adata = sc.read_h5ad(snakemake.input[0]) # 입력받은 h5ad파일


adata.obs.to_csv(snakemake.output[0],index=False)
adata.obs.to_csv(snakemake.output[1],index=False)
adata.var.to_csv(snakemake.output[2], index=False)
np.savetxt(snakemake.output[3], adata.X, delimiter=",")
# 일단 이 세개의 csv파일만 vue에 변수로 받아오면 좋을듯 변수명은 obs, var, X 정도로 생각중