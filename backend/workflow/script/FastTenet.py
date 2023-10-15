import subprocess
import json

import numpy as np
import fasttenet as fte
                            

def install_package(package_name, install_command):
    try:
        # Try importing the package
        __import__(package_name)
        print(f"{package_name} is already installed.")
    except ImportError:
        # Package not installed, install it
        print(f"{package_name} is not installed. Installing...")
        subprocess.check_call(install_command, shell=True)
        print(f"{package_name} has been installed successfully.")

# Check and install CuPy
# Note: Ensure that 'conda' is available in the PATH and the script has the necessary permissions.
# 12.0
# install_package("cupy", "conda install -c conda-forge cupy cuda-version=11.8")  # Replace xx.x with your CUDA version

# Check and install JAX
# install_package("jax", 'pip install --upgrade pip && pip install "jax[cuda11_cudnn86]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html')

def run_fasttenet(dpath_exp_data, dpath_trj_data, dpath_branch_data, spath_result_matrix, options):
    try:
        # Create worker
        # expression data, trajectory data, branch data path is required
        # tf data path is optional
        # save path is optional
        worker = fte.FastTENET(dpath_exp_data=dpath_exp_data, # Required
                                dpath_trj_data=dpath_trj_data, # Required
                                dpath_branch_data=dpath_branch_data, # Required
                                # dpath_tf_data=dpath_tf_data, # Optional
                                spath_result_matrix=spath_result_matrix, # Optional
                                make_binary=True) # Optional, default: False

        result_matrix = worker.run(device=options['device'], device_ids=options['device_ids'], batch_size=options['batch_size'], kp=options['kp'],
                                    percentile=options['percentile'], win_length=options['win_length'], polyorder=options['polyorder'])
        
        print(result_matrix)

        with open(snakemake.output[0], "w") as f:
            f.write("success")
    except Exception as e:
        print(e)

        with open(snakemake.output[0], "w") as f:
            f.write("fail")
        raise


options = json.load(open(snakemake.input[0]))

dpath_exp_data = snakemake.input[1]
dpath_trj_data = snakemake.input[2]
dpath_branch_data = snakemake.input[3]
spath_result_matrix = snakemake.output[1]

run_fasttenet(dpath_exp_data, dpath_trj_data, dpath_branch_data, spath_result_matrix, options)