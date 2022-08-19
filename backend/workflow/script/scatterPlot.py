import pandas as pd
import matplotlib.pyplot as plt

# with open(snakemake.input[0], 'rt') as f:
#     filename = f.read()
#     if filename == 'Sacramentorealestatetransactions.csv':
#         plt.scatter(df['price'], df['latitude'])
#     elif filename == 'SalesJan2009.csv':
#         plt.scatter(df['Latitude'], df['Longitude'])

#     plt.savefig(snakemake.output[1])

df = pd.read_csv(snakemake.input[1])
df.to_csv(snakemake.output[0])
plot = df.plot(kind='bar', rot=45)
fig = plot.get_figure()
fig.tight_layout()
fig.savefig(snakemake.output[1])