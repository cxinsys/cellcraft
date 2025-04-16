import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np
import json
import sys

# Handle command-line arguments
h5ad = sys.argv[1]  # Path to the input .h5ad file containing single-cell data
gene_list_path = None if sys.argv[2] == "None" else sys.argv[2]  # Optional: path to a file containing a list of genes to use
selected_indices_json_path = None if sys.argv[3] == "None" else sys.argv[3]  # Optional: JSON file containing manually selected cell indices
expMatrix = sys.argv[4]  # Output path to save the expression matrix CSV
pseudotime = sys.argv[5]  # Output path to save the pseudotime values CSV
cellSelect = sys.argv[6]  # Output path to save the selected cell binary list CSV
anno_of_interest = sys.argv[7]  # Column name in .obs to use for selecting cells (e.g., cluster, cell type)
pseudo_of_interest = sys.argv[8]  # Column name in .obs containing pseudotime values
clusters_of_interest_string = sys.argv[9]  # Semicolon-separated string listing clusters (or labels) of interest (e.g., "0;1;2")

# Load h5ad file
orig_data = sc.read_h5ad(h5ad)
adata = orig_data.copy()

# Normalize data types for adata.obs columns
df = adata.obs.copy()
for column in df.columns:
    try:
        if all(df[column].astype(float) % 1 == 0):
            df[column] = df[column].astype(int)
            if column not in ['nCount_RNA', 'nFeature_RNA', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt']:
                df[column] = df[column].astype('category')
        else:
            df[column] = df[column].astype(float)
    except (TypeError, ValueError):
        df[column] = df[column].astype('category')
adata.obs = df

# Define clusters of interest
clusters_of_interest = clusters_of_interest_string.split(';')
clusters_of_interest = [int(c) if c.isdigit() else c for c in clusters_of_interest if c]

# Add 'cell_select' column based on annotation of interest
adata.obs['cell_select'] = np.where(adata.obs[anno_of_interest].isin(clusters_of_interest), '1', '0')

# If selected index JSON file is provided, overwrite cell selection
if selected_indices_json_path:
    try:
        with open(selected_indices_json_path, 'r') as f:
            selected_indices = json.load(f)
        cell_selected_bools = np.array(['0'] * adata.obs.shape[0])
        cell_selected_bools[selected_indices] = '1'
        adata.obs['cell_select'] = cell_selected_bools
    except (FileNotFoundError, Exception) as e:
        print(f"Error loading selected indices JSON file: {e}")

# Filter selected cells
adata = adata[adata.obs['cell_select'] == '1']
adata = adata[adata.obs[pseudo_of_interest] != float('inf')]
adata = adata[~adata.obs[pseudo_of_interest].isna()]

# Exit if no cells are selected
if adata.shape[0] == 0:
    print('no cell selected')
    sys.exit(1)

# Check and sync raw count layer
def check_raw_data(data):
    if isinstance(data, csr_matrix):
        return np.all(data.data % 1 == 0)
    elif isinstance(data, np.ndarray):
        return np.all(data % 1 == 0)
    return False

if hasattr(adata, 'X') and check_raw_data(adata.X):
    adata.layers['counts'] = adata.X
elif hasattr(adata, 'raw') and hasattr(adata.raw, 'X') and check_raw_data(adata.raw.X):
    adata.layers['counts'] = adata.raw.X
elif 'counts' in adata.layers and check_raw_data(adata.layers['counts']):
    pass
else:
    print("Your adata object does not contain raw counts")

# Process gene list
gene_list = []
if gene_list_path:
    try:
        with open(gene_list_path, 'r') as f:
            gene_list = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Gene list file not found: {gene_list_path}")
    except Exception as e:
        print(f"Error reading gene list file: {e}")

# If gene list is empty, use all genes
if not gene_list:
    gene_list = list(adata.var_names)
else:
    available_genes = [gene for gene in gene_list if gene in adata.var_names]
    if not available_genes:
        print("None of the genes in the input gene_list exist in the data.")
        sys.exit(1)
    missing_genes = set(gene_list) - set(available_genes)
    if missing_genes:
        print(f"The following genes are missing in the data: {', '.join(missing_genes)}")
    gene_list = available_genes

# Find indices for selected genes
gene_indices = [adata.var_names.get_loc(gene) for gene in gene_list]

# Extract expression matrix data
if isinstance(adata.layers['counts'], csr_matrix):
    matrix_data = adata.layers['counts'][:, gene_indices].todense()
else:
    matrix_data = adata.layers['counts'][:, gene_indices]

# Save expression matrix
pd.DataFrame(data=matrix_data, index=adata.obs_names, columns=gene_list).to_csv(expMatrix)

# Save pseudotime
adata.obs[pseudo_of_interest].to_csv(pseudotime, index=False, header=False)

# Save cell selection
adata.obs['cell_select'].to_csv(cellSelect, index=False, header=False)
