import os
import pandas as pd
import scanpy as sc

# def load_tab_file(file_path: str, max_rows: int = 5000):
def load_tab_file(file_path: str):
    # Check if the file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    # Extract the file extension
    file_extension = os.path.splitext(file_path)[1].lower()
    
    # Check if the file extension is either csv or tsv
    if file_extension not in ['.csv', '.tsv', '.h5ad']:
        raise ValueError(f"The file {file_path} must have a .csv or .tsv or .h5ad extension.")
    
    if file_extension == '.h5ad':
        adata = sc.read_h5ad(file_path)
        # Process and combine data to form the desired dataframe structure
        df = pd.concat([adata.obs, pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'], index=adata.obs.index)], axis=1)
        # slice the dataframe to the maximum number of rows
        # df_trimmed = df.head(max_rows)
        # return df_trimmed
        return df
    else :
        # Determine the separator based on the file extension
        sep = ',' if file_extension == '.csv' else '\t'
        
        # Load the file with pandas
        # df = pd.read_csv(file_path, sep=sep, nrows=max_rows)
        df = pd.read_csv(file_path, sep=sep)
        
        # Return the trimmed dataframe
        return df

def transform_df_to_vgt_format(df):
    columns = [
        {
            "label": col,
            "field": col,
            "type": "number" if df[col].dtype in [int, float] else "string"
        }
        for col in df.columns
    ]

    rows = df.to_dict(orient='records')
    
    return {"columns": columns, "rows": rows}