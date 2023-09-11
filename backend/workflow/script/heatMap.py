import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input[0])
df.to_csv(snakemake.output[0], index=False)
df.to_csv(snakemake.output[1], index=False)
# plot = df.plot(kind='bar', rot=45)
# fig = plot.get_figure()
# fig.tight_layout()
# fig.savefig(snakemake.output[1])