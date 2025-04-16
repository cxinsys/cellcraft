import pandas as pd
import sys

exp_matrix_path = sys.argv[1]  # Path to the input expression matrix CSV file (genes as columns, cells as rows)
top_n = int(sys.argv[2])       # Number of top genes to select based on average expression
output_path = sys.argv[3]      # Output path to save the resulting top genes CSV file

# Load expression matrix data
df = pd.read_csv(exp_matrix_path, index_col=0)

# Calculate the average expression for each gene and select the top N
top_genes = df.mean(axis=0).nlargest(top_n).reset_index()
top_genes.columns = ['Gene', 'AverageExpression']

# Save the result
top_genes.to_csv(output_path, index=False)
print(f"Top {top_n} genes saved to {output_path}")
