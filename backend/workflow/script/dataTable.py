import pandas as pd

df = pd.read_csv(snakemake.input[0])
df.to_csv(snakemake.output[0])
df.to_json(snakemake.output[1])