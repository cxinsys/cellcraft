import numpy as np
import pandas as pd
import scanpy as sc # pip install scanpy
import matplotlib.pyplot as plt

adata = sc.read_h5ad(snakemake.input[0]) # 입력받은 h5ad파일


adata.obs.to_csv(snakemake.output[0],index=False)
adata.obs.to_csv(snakemake.output[1],index=False)
np.savetxt(snakemake.output[2], adata.obsm['X_umap'], delimiter=",")