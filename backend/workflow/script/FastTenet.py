import subprocess
import sys

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
install_package("cupy", "conda install -c conda-forge cupy cuda-version=11.8")  # Replace xx.x with your CUDA version

# Check and install JAX
install_package("jax", 'pip install --upgrade pip && pip install "jax[cuda11_cudnn86]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html')
