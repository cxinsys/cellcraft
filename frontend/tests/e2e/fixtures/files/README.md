# Test Fixture Files

This directory contains test files for E2E upload testing.

## Available Test Files

- **test_sample.csv**: Sample gene expression CSV file
- **test_genes.txt**: Sample gene list text file
- **test_data.h5ad**: H5AD format file (requires manual setup)

## H5AD File Setup

The `.h5ad` file format requires a valid HDF5 structure that cannot be easily created as a text fixture.

### Option 1: Download from Datasets Page (Recommended)
1. Run the application: `docker compose -f docker-compose.dev.yml up`
2. Navigate to http://localhost:8080/datasets
3. Download `PBMCLight1000.h5ad` from the tutorial datasets
4. Copy the downloaded file to this directory and rename it to `test_data.h5ad`

### Option 2: Use an Existing H5AD File
If you have access to any `.h5ad` file from your workflows:
```bash
cp /path/to/your/file.h5ad ./test_data.h5ad
```

### Option 3: Create with Python (For Advanced Users)
```python
import anndata as ad
import numpy as np
import pandas as pd

# Create minimal AnnData object
n_obs = 100  # cells
n_vars = 50  # genes

X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars))
obs = pd.DataFrame(index=[f'cell_{i}' for i in range(n_obs)])
var = pd.DataFrame(index=[f'gene_{i}' for i in range(n_vars)])

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.write('test_data.h5ad')
```

## File Naming Convention

When files are uploaded, they are automatically prefixed with the folder name (e.g., `data/` folder prefix). The test utilities handle this automatically.

## File Validation

The application validates file extensions and only accepts:
- `.h5ad` - AnnData/H5AD format for single-cell data
- `.csv` - Comma-separated values
- `.txt` - Plain text files

See `frontend/src/utils/validation.js` for validation logic.
