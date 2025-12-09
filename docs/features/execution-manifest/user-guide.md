# Execution Manifest User Guide

The execution manifest feature in CellCraft provides comprehensive tracking and documentation of your Gene Regulatory Network (GRN) analysis workflows, enabling complete reproducibility and detailed debugging capabilities.

## Table of Contents

- [Overview](#overview)
- [Getting Started](#getting-started)
- [Complete Field Reference](#complete-field-reference)
- [Understanding Execution Data](#understanding-execution-data)
- [Reproducibility Guide](#reproducibility-guide)
- [Examples and Use Cases](#examples-and-use-cases)
- [Troubleshooting](#troubleshooting)

## Overview

### What is an Execution Manifest?

An execution manifest is a comprehensive JSON file that captures every aspect of your workflow execution, including:

- **Complete workflow configuration** with all parameters and settings
- **Input/output file specifications** and data transformations
- **Step-by-step execution logs** from Snakemake workflow engine
- **Plugin metadata and dependencies** for exact environment recreation
- **Visual workflow representation** from the Drawflow interface
- **Timing and performance data** for optimization analysis

### When is it Available?

Execution manifests are **only available** for:
- ✅ Tasks with **SUCCESS** status
- ✅ **Analysis type** plugins (TENET, FastSCODE, GENIE3, etc.)
- ❌ Not available for visualization plugins or failed tasks

### File Naming Convention

Downloaded manifest files follow this pattern:
```
execution_manifest_{plugin_name}_{task_id_short}_{timestamp}.json
```

**Example**: `execution_manifest_FastSCODE_b41c3f1a_20250908_174351.json`

## Getting Started

### Step 1: Run Your Analysis

1. Create and configure your workflow in the CellCraft visual interface
2. Execute your analysis workflow normally
3. Wait for the task to complete with **SUCCESS** status

### Step 2: Download the Manifest

1. Navigate to the **Workflow** page
2. Find your completed task in the job monitoring table
3. Click the **download manifest** button (available only for successful analysis tasks)
4. The manifest JSON file will be downloaded to your computer

### Step 3: Explore the Data

Open the downloaded JSON file in any text editor or JSON viewer to explore the comprehensive execution data captured during your analysis.

## Complete Field Reference

The execution manifest contains five main sections:

### 1. Manifest Info (`manifest_info`)

Metadata about the manifest file itself:

```json
{
  "manifest_info": {
    "format_version": "1.0",           // Manifest schema version
    "generated_at": "2025-09-08T17:43:51.383934",  // ISO timestamp
    "generated_by": "username",        // User who downloaded the manifest
    "description": "CellCraft execution manifest for analysis reproducibility"
  }
}
```

### 2. Task Metadata (`task_metadata`)

Core information about the executed task:

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `task_id` | string (UUID) | Unique task identifier | `"b41c3f1a-6eb8-46aa-8b66-ced6cae0c6c4"` |
| `workflow_id` | integer | Associated workflow ID | `18` |
| `algorithm_id` | string | Algorithm execution instance ID | `"11"` |
| `plugin_name` | string | Name of the analysis plugin used | `"TENET"` |
| `task_type` | string | Type of task executed (typically "compile") | `"compile"` |
| `status` | string | Final task status (always "SUCCESS" in manifests) | `"SUCCESS"` |
| `start_time` | string (ISO-like) | Task start timestamp | `"2025-09-08 17:31:31.439645"` |
| `end_time` | string (ISO-like) | Task completion timestamp | `"2025-09-08 17:31:58.884843"` |
| `plugin_image_uri` | string (URI) | Docker image used for execution | `"ghcr.io/cxinsys/cellcraft-tenet:1.0"` |

### 3. Plugin Metadata (`plugin_metadata`)

Complete information about the analysis plugin:

| Field | Type | Description |
|-------|------|-------------|
| `id` | integer | Unique database identifier for the plugin |
| `name` | string | Plugin name |
| `description` | string | Plugin description and purpose |
| `author` | string | Plugin author |
| `version` | string | Plugin version |
| `plugin_type` | string | Always "analysis" for manifest-eligible tasks |
| `source` | string | Source type ("official" or "local") |
| `use_gpu` | boolean | Whether GPU acceleration was used |
| `dependencies` | object | Complete dependency specifications with filenames as keys and content as values |
| `drawflow` | object | Visual workflow representation data (complex nested structure) |
| `rules` | object | Snakemake rule definitions indexed by rule number |

### 4. Workflow Metadata (`workflow_metadata`)

Information about the parent workflow:

| Field | Type | Description |
|-------|------|-------------|
| `id` | integer | Workflow database ID |
| `title` | string | User-assigned workflow title |
| `workflow_info` | object | Complete workflow configuration with custom fields |
| `created_at` | string (ISO) | Workflow creation timestamp |
| `updated_at` | string (ISO) | Last workflow modification timestamp |

### 5. Execution Files (`execution_files`)

Complete capture of all execution-related files:

#### Logs (`logs`)
All log files generated during execution:

- **`run.log`**: Main Snakemake execution log with DAG building and job execution
- **`{rule_name}.stdout`**: Standard output for each workflow rule
- **`{rule_name}.stderr`**: Standard error output for each workflow rule

Each log entry contains:
```json
{
  "filename.log": {
    "content": "Complete log file content as string",
    "size": 3720,
    "modified_time": "1757320318.464426"
  }
}
```

**Field Details**:
- `content`: String containing complete log file content
- `size`: Integer file size in bytes
- `modified_time`: String containing Unix timestamp

#### Snakefile (`snakefile`)
The complete Snakefile used for execution:
```json
{
  "content": "Complete Snakefile content as string",
  "path": "./user/username/workflow_X/algorithm_Y/Snakefile"
}
```

#### Plugin Metadata (`plugin_metadata`)
The plugin's metadata.json file:
```json
{
  "content": { /* Complete parsed JSON content */ },
  "path": "/path/to/plugin/metadata.json"
}
```

#### Meta YAML (`meta_yml`)
Workflow execution metadata (if available):
```json
{
  "content": "Complete meta.yml content as string",
  "path": "./user/username/workflow_X/algorithm_Y/meta.yml"
}
```

## Understanding Execution Data

### Reading the Logs

**Main execution log (`run.log`)**:
- Shows the complete Snakemake workflow execution
- Lists all jobs and their execution order
- Includes timing information for each rule
- Shows input/output file dependencies

**Rule-specific logs**:
- `*.stdout`: Program output, results, and progress information
- `*.stderr`: Error messages, warnings, and diagnostic information

### Interpreting the Drawflow Data

The `drawflow` section contains the visual workflow representation:
- **Nodes**: Individual processing steps (rules)
- **Connections**: Data flow between steps
- **Parameters**: All input parameters and their values
- **File mappings**: Input and output file specifications

### Understanding Dependencies

The `dependencies` section shows exact package versions:
```json
{
  "requirements.txt": "fasttenet>=0.1.17\nscanpy==1.9.2\nscipy==1.10.1\n..."
}
```

This ensures you can recreate the exact computational environment.

## Reproducibility Guide

### Prerequisites

Before attempting to reproduce an analysis, ensure you have:

- **Python 3.8+** or **R 4.0+** (depending on plugin requirements)
- **Docker** (recommended) or **Snakemake 7.0+** for local execution
- **jq** command-line tool for JSON parsing
- **Git** for version control (optional but recommended)
- Access to the original input data files

### Step-by-Step Reproduction

To reproduce an analysis from an execution manifest:

#### 1. **Environment Setup**

**Option A: Using Docker (Recommended)**
```bash
# Extract Docker image from manifest
DOCKER_IMAGE=$(jq -r '.task_metadata.plugin_image_uri' manifest.json)
docker pull "$DOCKER_IMAGE"
```

**Option B: Local Python Environment**
```bash
# Extract dependencies from manifest and create requirements file
jq -r '.plugin_metadata.dependencies["requirements.txt"]' manifest.json > requirements.txt

# Create virtual environment
python -m venv cellcraft_env
source cellcraft_env/bin/activate  # On Windows: cellcraft_env\Scripts\activate

# Install exact dependencies
pip install -r requirements.txt
```

**Option C: Using Conda/Mamba**
```bash
# Create conda environment
conda create -n cellcraft_reproduction python=3.10
conda activate cellcraft_reproduction
pip install -r requirements.txt
```

#### 2. **Prepare Input Data**
- Locate the input files referenced in the manifest
- Ensure they're placed in the exact paths specified in the Snakefile
- Example: `user/testuser/data/pbmc_light_1000.h5ad`

#### 3. **Recreate Directory Structure**
```bash
# Create the workflow directory structure
mkdir -p user/{username}/workflow_{workflow_id}/algorithm_{algorithm_id}
mkdir -p user/{username}/workflow_{workflow_id}/algorithm_{algorithm_id}/results
mkdir -p user/{username}/workflow_{workflow_id}/algorithm_{algorithm_id}/logs
```

#### 4. **Extract and Run Snakefile**

**Extract Snakefile using jq:**
```bash
# Extract Snakefile content from manifest
jq -r '.execution_files.snakefile.content' manifest.json > Snakefile

# Verify extraction
head -n 10 Snakefile
```

**Execute with appropriate method:**

```bash
# Option A: Using Docker
DOCKER_IMAGE=$(jq -r '.task_metadata.plugin_image_uri' manifest.json)
docker run -v $(pwd):/workspace -w /workspace "$DOCKER_IMAGE" snakemake --cores 1

# Option B: Local Snakemake (ensure Snakemake is installed)
snakemake --cores 1 --use-conda

# Option C: Check Snakemake version compatibility
snakemake --version  # Should be 7.0+ for best compatibility
```

#### 5. **Verify Results**
- Compare output files with expected results
- Check log files match the manifest logs
- Validate execution timing and resource usage

### Parameter Extraction

**Using jq (Command Line)**
```bash
# Extract all parameters with their values
jq -r '.plugin_metadata.drawflow.drawflow.Home.data[].data.parameters[] | "\(.name): \(.defaultValue)"' manifest.json

# Extract specific parameter types
jq -r '.plugin_metadata.drawflow.drawflow.Home.data[].data.parameters[] | select(.type=="int" or .type=="float") | "\(.name): \(.defaultValue)"' manifest.json

# Get Docker image information
jq -r '.task_metadata | "Plugin: \(.plugin_name)\nDocker Image: \(.plugin_image_uri)\nExecution Time: \(.end_time as $end | .start_time as $start | $end)"' manifest.json
```

**Using Python**
```python
import json
from datetime import datetime

# Load manifest
with open('execution_manifest.json', 'r') as f:
    manifest = json.load(f)

# Extract parameters from drawflow data
def extract_parameters(manifest):
    parameters = {}
    drawflow_data = manifest['plugin_metadata']['drawflow']
    
    for node_id, node in drawflow_data['drawflow']['Home']['data'].items():
        rule_name = node['data']['title']
        rule_params = {}
        
        for param in node['data']['parameters']:
            if param['type'] not in ['inputFile', 'outputFile']:
                rule_params[param['name']] = param['defaultValue']
        
        if rule_params:  # Only add if there are non-file parameters
            parameters[rule_name] = rule_params
    
    return parameters

# Extract and display parameters
parameters = extract_parameters(manifest)
for rule_name, params in parameters.items():
    print(f"\n{rule_name}:")
    for param_name, param_value in params.items():
        print(f"  {param_name}: {param_value}")
```

### Docker Environment Recreation

Use the exact Docker image specified in the manifest:

```bash
# Extract and pull the exact image used
DOCKER_IMAGE=$(jq -r '.task_metadata.plugin_image_uri' manifest.json)
echo "Using Docker image: $DOCKER_IMAGE"
docker pull "$DOCKER_IMAGE"

# Run analysis in the same environment
docker run --rm -v "$(pwd)":/workspace -w /workspace "$DOCKER_IMAGE" snakemake --cores 1

# For GPU-enabled plugins (check manifest for use_gpu: true)
if [ $(jq -r '.plugin_metadata.use_gpu' manifest.json) = "true" ]; then
    docker run --rm --gpus all -v "$(pwd)":/workspace -w /workspace "$DOCKER_IMAGE" snakemake --cores 1
fi
```

## Examples and Use Cases

### Example 1: Debugging Failed Reproduction

**Problem**: Your reproduced analysis gives different results.

**Solution**:
1. **Check input file checksums**:
   ```bash
   # Generate checksums for your input files
   md5sum *.h5ad *.txt > input_checksums.txt
   # Compare with original checksums if available
   ```

2. **Verify parameters exactly match**:
   ```bash
   # Extract and verify specific parameters
   jq '.plugin_metadata.drawflow.drawflow.Home.data[].data.parameters[] | select(.name=="FDR") | .defaultValue' manifest.json
   # Should show: "0.01" (for TENET FDR parameter)
   
   # Check all numerical parameters
   jq -r '.plugin_metadata.drawflow.drawflow.Home.data[].data.parameters[] | select(.type=="int" or .type=="float") | "\(.name): \(.defaultValue)"' manifest.json
   ```

3. **Compare log outputs systematically**:
   ```bash
   # Extract specific log content
   jq -r '.execution_files.logs["run.log"].content' manifest.json > original_run.log
   # Compare with your reproduction
   diff original_run.log your_run.log
   ```

4. **Check Snakemake version compatibility**:
   ```bash
   snakemake --version
   # CellCraft uses Snakemake 7.14.0 - significant version differences can affect behavior
   ```

### Example 2: Performance Analysis

**Use Case**: Understanding computational requirements.

```python
import json
from datetime import datetime

with open('manifest.json', 'r') as f:
    manifest = json.load(f)

# Calculate total execution time
start = datetime.fromisoformat(manifest['task_metadata']['start_time'])
end = datetime.fromisoformat(manifest['task_metadata']['end_time'])
duration = end - start
print(f"Total execution time: {duration}")

# Analyze per-rule timing from run.log
run_log = manifest['execution_files']['logs']['run.log']['content']
# Parse timestamps to understand bottlenecks
```

### Example 3: Parameter Sensitivity Analysis

**Use Case**: Understanding which parameters affect results.

**Using jq for comparison:**
```bash
# Extract parameters from multiple manifests for comparison
jq -r '.plugin_metadata.drawflow.drawflow.Home.data[].data.parameters[] | select(.type=="int" or .type=="float") | "\(.name): \(.defaultValue)"' manifest1.json > params1.txt
jq -r '.plugin_metadata.drawflow.drawflow.Home.data[].data.parameters[] | select(.type=="int" or .type=="float") | "\(.name): \(.defaultValue)"' manifest2.json > params2.txt

# Compare parameter files
diff params1.txt params2.txt
```

**Using Python for detailed analysis:**
```python
def compare_parameters(manifest1, manifest2):
    """Compare parameters between two manifests and identify differences."""
    params1 = extract_parameters(manifest1)
    params2 = extract_parameters(manifest2)
    
    all_rules = set(params1.keys()) | set(params2.keys())
    differences = {}
    
    for rule in all_rules:
        rule_diffs = {}
        
        if rule not in params1:
            rule_diffs['missing_in_first'] = params2[rule]
        elif rule not in params2:
            rule_diffs['missing_in_second'] = params1[rule]
        else:
            # Compare parameter values
            all_params = set(params1[rule].keys()) | set(params2[rule].keys())
            for param in all_params:
                val1 = params1[rule].get(param, 'MISSING')
                val2 = params2[rule].get(param, 'MISSING')
                if val1 != val2:
                    rule_diffs[param] = {'manifest1': val1, 'manifest2': val2}
        
        if rule_diffs:
            differences[rule] = rule_diffs
    
    return differences

# Usage
differences = compare_parameters(manifest1, manifest2)
for rule, diffs in differences.items():
    print(f"\nDifferences in {rule}:")
    for param, values in diffs.items():
        print(f"  {param}: {values}")
```

### Example 4: Workflow Optimization

**Use Case**: Optimizing workflow performance based on execution data.

```python
# Analyze resource usage from stderr logs
def analyze_gpu_usage(manifest):
    for filename, log_data in manifest['execution_files']['logs'].items():
        if 'stderr' in filename:
            content = log_data['content']
            # Look for GPU utilization patterns
            if 'gpu' in content.lower():
                print(f"GPU usage in {filename}:")
                for line in content.split('\n'):
                    if 'gpu' in line.lower():
                        print(f"  {line}")
```

## Troubleshooting

### Common Issues

#### Issue 1: Manifest Download Not Available

**Symptoms**: Download button is grayed out or missing.

**Causes & Solutions**:
- ❌ **Task Status**: Only SUCCESS status tasks have manifests
  - **Solution**: Ensure your task completed successfully
- ❌ **Plugin Type**: Only Analysis plugins support manifests
  - **Solution**: Visualization plugins don't generate manifests
- ❌ **Permissions**: You can only download your own task manifests
  - **Solution**: Log in as the task owner

#### Issue 2: Incomplete Manifest Data

**Symptoms**: Missing logs or files in the manifest.

**Possible Causes**:
- File system permissions issues
- Log files were cleaned up before manifest generation
- Network interruption during download

**Solution**: Re-download the manifest or contact system administrator.

#### Issue 3: Reproduction Failures

**Symptoms**: Different results when reproducing the analysis.

**Debugging Steps**:
1. **Check Input Data**:
   ```bash
   # Verify input file checksums and sizes
   md5sum your_input.h5ad original_input.h5ad
   ls -la *.h5ad  # Check file sizes match
   
   # For H5AD files, check internal structure
   python -c "import scanpy as sc; adata=sc.read('your_input.h5ad'); print(f'Shape: {adata.shape}'); print(f'Obs columns: {list(adata.obs.columns)}')"
   ```

2. **Verify Environment Precisely**:
   ```bash
   # Extract exact versions from manifest
   jq -r '.plugin_metadata.dependencies["requirements.txt"]' manifest.json
   
   # Check your current versions match
   pip list | grep -E "(fasttenet|scanpy|pandas|numpy|scipy)"
   
   # Check Python version compatibility
   python --version  # Should be compatible with plugin requirements
   ```

3. **Compare Logs Systematically**:
   ```bash
   # Extract and compare specific log sections
   jq -r '.execution_files.logs["run.log"].content' manifest.json > original_run.log
   jq -r '.execution_files.logs["TENET_Input.stdout"].content' manifest.json > original_input.log
   
   # Compare job statistics and timing
   grep -E "Job stats:|Finished job|\[.*\]" original_run.log
   grep -E "Job stats:|Finished job|\[.*\]" your_run.log
   ```

4. **Parameter Validation with jq**:
   ```bash
   # Extract Snakefile and check parameters
   jq -r '.execution_files.snakefile.content' manifest.json > original_snakefile
   
   # Verify critical parameters
   grep -E "FDR_THRESHOLD|NUM_DEVICES|NUM_OUTDEGREES" original_snakefile
   
   # Check parameter consistency across manifest sections
   jq -r '.plugin_metadata.drawflow.drawflow.Home.data["3"].data.parameters[] | select(.name=="FDR") | .defaultValue' manifest.json
   ```

5. **Version-specific Checks**:
   ```bash
   # Check Snakemake version - critical for reproducibility
   snakemake --version
   # CellCraft uses Snakemake 7.14.0
   
   # Check if using GPU vs CPU affects results
   nvidia-smi  # Check GPU availability
   jq -r '.plugin_metadata.use_gpu' manifest.json  # Check if GPU was used originally
   ```

#### Issue 4: Large Manifest Files

**Symptoms**: Very large manifest files (>100MB).

**Causes**: Extensive log outputs from long-running analyses.

**Solutions**:
- Use JSON streaming parsers for processing
- Extract specific sections rather than loading the entire file
- Consider log rotation for future runs

**Using jq for large manifests:**
```bash
# Extract only specific sections without loading entire file
jq -r '.execution_files.logs["run.log"].content' large_manifest.json > run_log_only.txt

# Get manifest size and structure overview
ls -lh large_manifest.json
jq 'keys' large_manifest.json  # Top-level sections
jq '.execution_files | keys' large_manifest.json  # Available files

# Extract minimal information for analysis
jq '{task_id: .task_metadata.task_id, plugin: .task_metadata.plugin_name, duration: (.task_metadata.end_time | split(" ")[1] | split(":") | map(tonumber) | .[0]*3600 + .[1]*60 + .[2]) - (.task_metadata.start_time | split(" ")[1] | split(":") | map(tonumber) | .[0]*3600 + .[1]*60 + .[2])}' large_manifest.json
```

**Using Python with ijson for streaming:**
```python
import ijson
import json

def extract_logs_streaming(filename):
    """Extract logs from large manifest without loading entire file."""
    logs = {}
    
    with open(filename, 'rb') as file:
        # Stream parse only the logs section
        parser = ijson.parse(file)
        current_path = []
        current_log = None
        
        for prefix, event, value in parser:
            if prefix.startswith('execution_files.logs.'):
                if event == 'start_map' and prefix.count('.') == 2:
                    # Starting a new log file
                    current_log = prefix.split('.')[-1]
                    logs[current_log] = {}
                elif event == 'string' and current_log:
                    # Extract content and metadata
                    field = prefix.split('.')[-1]
                    logs[current_log][field] = value
    
    return logs

# Usage for large files
logs_only = extract_logs_streaming('large_manifest.json')
print(f"Extracted {len(logs_only)} log files")
print(f"Main log size: {logs_only.get('run.log', {}).get('size', 'N/A')} bytes")
```

### Getting Help

If you encounter issues not covered in this guide:

1. **Check the API logs** for error messages during manifest generation
2. **Verify file permissions** in the user data directories
3. **Contact system administrator** for infrastructure-related issues
4. **Report bugs** with specific error messages and task IDs

### Best Practices

1. **Download manifests promptly** after task completion
2. **Store manifests securely** for long-term reproducibility
3. **Document any manual modifications** made during reproduction
4. **Version control your input data** along with manifests

### Security and Data Privacy

**Important**: Execution manifests contain complete workflow information including:
- File paths and system directories
- Plugin configurations and parameters  
- Execution logs and system output

**Recommendations**:
- Review manifest contents before sharing externally
- Manifests do not contain raw data files, only metadata and logs
- Check logs for any sensitive information before publication
5. **Test reproduction workflows** in isolated environments

---

*This user guide covers the comprehensive execution manifest feature in CellCraft. For technical details about the API endpoints and schemas, see the [API Reference](api-reference.md). For a complete example of a TENET analysis manifest file, see [sample-tenet-manifest.json](examples/sample-tenet-manifest.json).*