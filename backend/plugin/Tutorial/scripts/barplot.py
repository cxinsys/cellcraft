import sys
import pandas as pd
import plotly.express as px
import json

# Get file paths from command-line arguments
input_file_path = sys.argv[1]  # top_genes.csv
output_file_path = sys.argv[2]  # barplot.json

# Load data
df = pd.read_csv(input_file_path)

# Sort in descending order by average expression
df = df.sort_values(by='AverageExpression', ascending=False)

# Define a custom color palette
color_palette = px.colors.qualitative.Set3

# Create a horizontal bar plot using Plotly
fig = px.bar(
    df,
    x='AverageExpression',
    y='Gene',
    color='Gene',
    orientation='h',
    title="Top Genes Expression Barplot",
    color_discrete_sequence=color_palette
)

# Update layout to sort categories by total expression
fig.update_layout(yaxis={'categoryorder': 'total ascending'})

# Save the plot as JSON
fig_json = fig.to_json()
with open(output_file_path, "w") as json_file:
    json_file.write(fig_json)

print(f"Barplot JSON saved to {output_file_path}")
