import pandas as pd

df = pd.read_csv(snakemake.input[1])
df.to_csv(snakemake.output[0])
df.to_html(snakemake.output[1])