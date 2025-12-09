<img src="https://github.com/cxinsys/cellcraft/blob/807998fda59e15e185ea9d2835ff7b81a884460f/frontend/src/assets/cellcraft_logo_text.png"/>

[Demo website](https://cellcraft.app) • [Status page(Demo website)](https://status.cellcraft.app/) • [Docs](https://cellcraft.gitbook.io/cellcraft-docs)

### Plugin Resources

[**Official Plugins**](https://github.com/cxinsys/cellcraft-plugin) • [**Plugin Templates**](https://github.com/cxinsys/cellcraft-plugin-templates)

- **Official Plugins Repository**: Contains version information and source code for all officially supported CellCraft plugins
- **Plugin Templates Repository**: Provides development resources and templates for creating custom local plugins

### Container Images

Pre-built Docker images are available on GitHub Container Registry (GHCR) for faster deployment:

[**Frontend**](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Ffrontend) • [**Backend-CPU**](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fbackend-cpu) • [**Backend-GPU**](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fbackend-gpu) • [**Celery-CPU**](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fcelery-cpu) • [**Celery-GPU**](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fcelery-gpu)

- **Multi-platform support**: All CPU images support AMD64 and ARM64 architectures
- **GPU images**: AMD64 only, optimized for NVIDIA CUDA
- **Version**: Uses `latest` tag for most recent multi-platform builds

## Overview

**CellCraft** is a web-based application for reconstructing **gene regulatory networks (GRNs)** from single-cell RNA sequencing (scRNA-seq) data. It features an intuitive visual interface that integrates seven GRN inference tools—including **TENET** and **FastTENET**, developed by our research team—through modular plugin integration. 

Built to **lower technical barriers** in GRN analysis, CellCraft enables researchers to **configure workflows**, **run analysis**, and **explore results visually**, all without needing programming knowledge. The platform unifies access to multiple inference methods behind a consistent, guided interface and supports scalable, customizable workflows suitable for both rapid prototyping and extensive studies.


## Key Features

- **Visual Workflow Programming**: Configure and execute complex GRN analyses through an intuitive drag-and-drop interface designed for workflow-based programming.
- **Seven Integrated GRN Tools**: Built-in support for **TENET**, **FastTENET**, **FastSCODE**, **GENIE3**, **GRNBoost2**, **LEAP**, **Scribe**, and visualization tools, with continuous expansion of capabilities.
- **Interactive Visualizations**: Explore regulatory relationships and analysis results through dynamic, interactive visualizations that make complex data interpretable.
- **Reproducibility & Onboarding**: Comprehensive tutorials, example datasets, and sample input files guide users through each step from data preparation to network interpretation.
- **Scalable Workflows**: Support for both rapid proof-of-concept analyses and extensive multi-stage studies through an integrated toolchain.
- **Custom Plugin Support**: Extend functionality by creating and integrating custom analysis plugins (Windows/Linux).
- **User Data Management**: Manage datasets through user accounts, upload files, and organize projects for streamlined analysis workflows.

---

## System Requirements

Before installing CellCraft, please verify your system meets the following requirements:

| Component | Minimum | Recommended |
| :--- | :--- | :--- |
| CPU | 4 cores | 8+ cores |
| RAM | 8 GB | 16+ GB |
| Storage | 70 GB | 100+ GB |
| OS | Ubuntu 20.04 LTS, Windows 10/11, macOS 10.15 Catalina | Ubuntu 22.04 LTS |
| OS Kernel | 6.8.0 | 6.8.0+ |
| glibc | 2.39 | 2.39+ |
| Docker | 20.10.0 | 24.0.0+ |
| Docker Compose | v2.0.0 | v2.20.0+ |
| NVIDIA Driver | 535.171.04 (GPU only) | 535.171.04+ |
| CUDA | 12.1 (GPU only) | 12.2+ |

**Additional Notes**:
- For GPU-enabled installation, use `./run-gpu-mode.sh`
- For CPU-only installation, use `./run-cpu-mode.sh`
- Docker Desktop (latest version) is required
- Git with submodule support is required
- Modern web browser (Chrome, Firefox, Edge, or Safari) is required

---

## Getting Started

0. **Install Prerequisites**:

   Before starting, ensure you have the following tools installed on your system.

   **Git Installation:**
   ```bash
   # Ubuntu/Debian
   sudo apt-get update && sudo apt-get install -y git

   # macOS (using Homebrew)
   brew install git

   # Windows: Download from https://git-scm.com/download/win
   ```

   **Docker Installation:**
   ```bash
   # Ubuntu/Debian - Install Docker Engine
   sudo apt-get update
   sudo apt-get install -y ca-certificates curl gnupg
   sudo install -m 0755 -d /etc/apt/keyrings
   curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
   sudo chmod a+r /etc/apt/keyrings/docker.gpg
   echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
   sudo apt-get update
   sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

   # Add current user to docker group (logout/login required)
   sudo usermod -aG docker $USER

   # macOS/Windows: Install Docker Desktop from https://www.docker.com/products/docker-desktop
   ```

   **NVIDIA Container Toolkit (GPU users only):**
   ```bash
   # Ubuntu/Debian - Required for GPU-enabled installation
   curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
   curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
     sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
     sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
   sudo apt-get update
   sudo apt-get install -y nvidia-container-toolkit
   sudo nvidia-ctk runtime configure --runtime=docker
   sudo systemctl restart docker

   # Verify GPU access in Docker
   docker run --rm --gpus all nvidia/cuda:12.1.0-base-ubuntu22.04 nvidia-smi
   ```

   > **Note**: After installing Docker, you may need to log out and log back in for group changes to take effect. For Windows/macOS, Docker Desktop includes Docker Compose by default.

1. Clone the repository:

   ```bash
   git clone --recurse-submodules https://github.com/cxinsys/cellcraft.git
   cd cellcraft
   ```

2. Configure plugin submodule for your installation mode:

   The plugin submodule must be set to the correct branch before starting.

   **Check current submodule status:**
   ```bash
   cd backend/plugin/official
   git status
   ```

   Interpret the output:
   - `On branch release/plugins-v1.1` → Ready for GPU mode
   - `On branch release/plugins-v1.0-cpu` → Ready for CPU mode
   - `HEAD detached at ...` → Switch to your desired branch below

   Target branches:
   - `release/plugins-v1.1` for GPU-enabled installation
   - `release/plugins-v1.0-cpu` for CPU-only installation

   **For GPU-enabled installation:**
   ```bash
   git switch release/plugins-v1.1
   cd ../../..
   ```

   **For CPU-only installation:**
   ```bash
   git switch release/plugins-v1.0-cpu
   cd ../../..
   ```

3. (Optional) Manage GHCR images:

   Before installation, you can optionally **check and pre-download** required Docker images from GitHub Container Registry (GHCR). This step is recommended for faster deployment and offline environments.

   **Interactive mode (recommended for first-time users):**
   ```bash
   ./test-ghcr-check.sh
   ```

   Select from the menu:
   - Option 1 or 2: Check image availability for CPU or GPU mode
   - Option 4 or 5: Pre-download images for your preferred mode

   **Command-line mode (for automation):**
   ```bash
   ./test-ghcr-check.sh --cpu  # For CPU-only mode
   ./test-ghcr-check.sh --gpu  # For GPU-enabled mode
   ```

   **What this script does:**
   - ✅ Checks if images exist locally (instant deployment if available)
   - ⚠️ Verifies remote GHCR accessibility (download if needed)
   - 📥 Optionally pre-downloads images based on your deployment mode
   - 📊 Provides detailed status and download statistics

   **Note:** If you skip this step or if GHCR is not accessible, the installation scripts will automatically fall back to building images locally. Pre-downloading images is simply an optimization for faster deployment.

4. Start the application:

   **Option A: Using automated scripts (recommended)**

   For GPU-enabled installation:
   ```bash
   ./run-gpu-mode.sh
   ```

   For CPU-only installation:
   ```bash
   ./run-cpu-mode.sh
   ```

   **Option B: Using docker compose directly**

   For GPU-enabled setup (AMD64 with NVIDIA GPU):
   ```bash
   docker compose -f docker-compose.gpu.amd64.yml up -d --build
   ```

   For CPU-only setup (AMD64):
   ```bash
   docker compose -f docker-compose.cpu.amd64.yml up -d --build
   ```

   For CPU-only setup (ARM64):
   ```bash
   docker compose -f docker-compose.cpu.arm64.yml up -d --build
   ```

5. Access the application at [http://localhost:8080](http://localhost:8080).

6. Check the installation status:

   ```bash
   ./check-installation.sh
   ```

   This script validates the following:

   **Container Status**: Verifies all 5 required services are running (frontend, backend, db, rabbitmq, celery)
   ```bash
   # Manual check
   docker ps --format 'table {{.Names}}\t{{.Status}}' | grep cellcraft
   ```

   **HTTP Connectivity**: Tests Frontend (port 8080) and Backend API (port 8000)
   ```bash
   # Manual check
   curl -sf http://localhost:8080 && echo "Frontend OK"
   curl -sf http://localhost:8000/docs && echo "Backend OK"
   ```

   **Plugin Initialization**: On first startup, the backend downloads plugin Docker images. Monitor the progress with:
   ```bash
   # Watch plugin download progress in real-time
   docker compose -f docker-compose.gpu.amd64.yml logs -f backend   # For GPU mode
   docker compose -f docker-compose.cpu.amd64.yml logs -f backend   # For CPU mode (AMD64)
   docker compose -f docker-compose.cpu.arm64.yml logs -f backend   # For CPU mode (ARM64)
   ```

   > **Important - Plugin Download Phase:**
   > - Plugin download takes approximately **10-15 minutes** depending on your network speed
   > - During this phase, the web interface at `http://localhost:8080` **will be accessible**, but **backend API features are unavailable**
   > - Features like **sign-up, login, and workflow execution will not work** until plugin initialization completes
   > - Wait until you see `Uvicorn running on http://0.0.0.0:8000` in the logs before using the application
   > - When you see `4. Docker Images Check...` in the logs, you can monitor the real-time download progress of each plugin image (e.g., `⬇ Pulling ghcr.io/cxinsys/cellcraft-fastscode:1.0... Downloading: 58.3%`)

   For **Windows/macOS** users: You can also verify plugin images in **Docker Desktop > Images** section.

   **Using command line:**
   ```bash
   # List all plugin images
   docker images | grep -E "(tenet|fasttenet|fastscode|genie3|grnboost2|leap|scribe|grnviz)"
   ```

   Expected plugins:
   - **CPU mode (6 plugins)**: GENIE3, GRNBoost2, GRNViz, LEAP, Scribe, TENET
   - **GPU mode (8 plugins)**: Above 6 + FastSCODE, FastTENET

---

## Stopping and Cleaning Up

### Stop Services

To stop all CellCraft services while preserving your data:

```bash
# For GPU mode
docker compose -f docker-compose.gpu.amd64.yml down

# For CPU mode (AMD64)
docker compose -f docker-compose.cpu.amd64.yml down

# For CPU mode (ARM64)
docker compose -f docker-compose.cpu.arm64.yml down
```

### Stop Services and Delete All Data

To stop services and **remove all data** (database, uploaded files, etc.):

```bash
# For GPU mode (WARNING: This deletes all data!)
docker compose -f docker-compose.gpu.amd64.yml down -v

# For CPU mode (AMD64)
docker compose -f docker-compose.cpu.amd64.yml down -v

# For CPU mode (ARM64)
docker compose -f docker-compose.cpu.arm64.yml down -v
```

> **Warning**: The `-v` flag removes all Docker volumes, including your database and uploaded files. Use this only when you want a fresh start.

### Remove Docker Images

To free up disk space by removing CellCraft Docker images:

```bash
# List all CellCraft-related images
docker images | grep cellcraft

# Remove specific image (replace with actual repository:tag)
docker rmi ghcr.io/cxinsys/cellcraft/frontend:latest
docker rmi ghcr.io/cxinsys/cellcraft/backend-cpu:latest
docker rmi ghcr.io/cxinsys/cellcraft/celery-cpu:latest

# Or remove all CellCraft images at once
docker images | grep cellcraft | awk '{print $1":"$2}' | xargs docker rmi
```

To remove plugin images:
```bash
# List plugin images
docker images | grep -E "(tenet|fasttenet|fastscode|genie3|grnboost2|leap|scribe|grnviz)"

# Remove all plugin images
docker images | grep -E "(tenet|fasttenet|fastscode|genie3|grnboost2|leap|scribe|grnviz)" | awk '{print $1":"$2}' | xargs docker rmi
```

---

## Additional Resources

For more detailed information, see our documentation:

- **Installation & Deployment**
  - [GHCR Image Management](docs/installation/GHCR_IMAGE_MANAGEMENT.md) - Pre-built Docker images and management tools
  - [Troubleshooting](docs/installation/TROUBLESHOOTING.md) - Common issues and solutions

- **User Guides**
  - [Tutorials](docs/guides/TUTORIALS.md) - Step-by-step video guides
  - [Platform-Specific Guide](docs/guides/PLATFORM_SPECIFIC.md) - macOS, Windows, Linux specific notes

- **Features**
  - [Execution Manifest](docs/features/execution-manifest/README.md) - Workflow reproducibility and debugging

📺 **Video Tutorials**: [CellCraft YouTube Channel](https://www.youtube.com/@CellCraft-cislab)
