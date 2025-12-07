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
- **Seven Integrated GRN Tools**: Built-in support for **TENET**, **FastTENET**, **FastSCODE**, **GENIE3**, **GRNBOOST2**, **LEAP**, **Scribe**, and visualization tools, with continuous expansion of capabilities.
- **Interactive Visualizations**: Explore regulatory relationships and analysis results through dynamic, interactive visualizations that make complex data interpretable.
- **Reproducibility & Onboarding**: Comprehensive tutorials, example datasets, and sample input files guide users through each step from data preparation to network interpretation.
- **Scalable Workflows**: Support for both rapid proof-of-concept analyses and extensive multi-stage studies through an integrated toolchain.
- **Custom Plugin Support**: Extend functionality by creating and integrating custom analysis plugins (Windows/Linux).
- **User Data Management**: Manage datasets through user accounts, upload files, and organize projects for streamlined analysis workflows.

---

## Getting Started

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
   - `On branch release/plugins-v1.0` → Ready for GPU mode
   - `On branch release/plugins-v1.0-cpu` → Ready for CPU mode
   - `HEAD detached at ...` → Switch to your desired branch below

   Target branches:
   - `release/plugins-v1.0` for GPU-enabled installation
   - `release/plugins-v1.0-cpu` for CPU-only installation

   **For GPU-enabled installation:**
   ```bash
   git switch release/plugins-v1.0
   cd ../../..
   ```

   **For CPU-only installation:**
   ```bash
   git switch release/plugins-v1.0-cpu
   cd ../../..
   ```

   > **Note:** This step is required before first-time installation. The automated scripts (`run-gpu-mode.sh`, `run-cpu-mode.sh`) will handle this automatically, but manual configuration is needed if using docker compose directly.

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
   docker compose -f docker-compose.gpu.amd64.yml up --build
   ```

   For CPU-only setup (AMD64):
   ```bash
   docker compose -f docker-compose.cpu.amd64.yml up --build
   ```

   For CPU-only setup (ARM64):
   ```bash
   docker compose -f docker-compose.cpu.arm64.yml up --build
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

   **Plugin Registry**: Verifies all expected plugins are registered
   ```bash
   # Manual check - view backend logs
   docker compose -f [compose-file] logs backend
   ```

   If you see `Initializing plugins from CSV...` in the logs but the server startup message (e.g., `Uvicorn running on http://0.0.0.0:8000`) has not appeared yet, the backend is still loading plugins.

   For **Windows/macOS** users: You can also verify plugin images in **Docker Desktop > Images** section.

   > **Note:** On first startup, please wait **10-15 minutes** for all plugins to be fully configured. The backend service needs time to initialize and pull plugin images.

   Expected plugins:
   - **CPU mode (6 plugins)**: GENIE3, GRNBoost2, GRNViz, LEAP, Scribe, TENET
   - **GPU mode (8 plugins)**: Above 6 + FastSCODE, FastTENET

---

## GHCR Image Management

CellCraft provides **pre-built Docker images** hosted on GitHub Container Registry (GHCR) for faster deployment. The `test-ghcr-check.sh` script helps you manage these images efficiently.

### Available Container Images

| Image | Description | Platforms | GHCR Package | Latest Version |
|-------|-------------|-----------|--------------|----------------|
| **frontend** | Vue.js frontend application | AMD64, ARM64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Ffrontend) | latest |
| **backend-cpu** | FastAPI backend (CPU-only) | AMD64, ARM64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fbackend-cpu) | latest |
| **backend-gpu** | FastAPI backend (GPU-enabled) | AMD64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fbackend-gpu) | latest |
| **celery-cpu** | Celery worker (CPU-only) | AMD64, ARM64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fcelery-cpu) | latest |
| **celery-gpu** | Celery worker (GPU-enabled) | AMD64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fcelery-gpu) | latest |

**Quick Pull Commands:**

For CPU mode (all 3 images):
```bash
docker pull ghcr.io/cxinsys/cellcraft/frontend:latest
docker pull ghcr.io/cxinsys/cellcraft/backend-cpu:latest
docker pull ghcr.io/cxinsys/cellcraft/celery-cpu:latest
```

For GPU mode (all 3 images):
```bash
docker pull ghcr.io/cxinsys/cellcraft/frontend:latest
docker pull ghcr.io/cxinsys/cellcraft/backend-gpu:latest
docker pull ghcr.io/cxinsys/cellcraft/celery-gpu:latest
```

**Note:** The `test-ghcr-check.sh` script uses the `latest` tag to always pull the most recent multi-platform builds from GitHub Actions. This ensures you get the newest stable images automatically.

### Image Management Tool

The GHCR image checker provides both **interactive menu** and **command-line options** for managing CellCraft images:

**Interactive Mode:**
```bash
./test-ghcr-check.sh
```

This launches an interactive menu with the following options:
- Check CPU mode images (frontend, backend-cpu, celery-cpu)
- Check GPU mode images (frontend, backend-gpu, celery-gpu)
- Check all images (both CPU and GPU)
- Download CPU mode images
- Download GPU mode images

**Command-Line Options:**
```bash
./test-ghcr-check.sh --cpu         # Check and optionally download CPU images
./test-ghcr-check.sh --gpu         # Check and optionally download GPU images
./test-ghcr-check.sh --check-only  # Check all images without downloading
./test-ghcr-check.sh --help        # Show usage information
```

### Image Status Indicators

The script displays clear status for each image:
- ✅ **LOCAL**: Image exists locally (instant startup, no download needed)
- ⚠️ **REMOTE**: Image available remotely (will be downloaded when needed)
- ❌ **MISSING**: Image not accessible (will fall back to local build)

### Pre-downloading Images

Pre-downloading images is recommended for:
- **Faster first-time deployment** (no wait during startup)
- **Offline environments** (download once, deploy anytime)
- **Network-constrained setups** (download during off-peak hours)

The script intelligently skips images that already exist locally and provides detailed download statistics.

### Smart Image Detection

Both `run-cpu-mode.sh` and `run-gpu-mode.sh` implement **local-first checking**:
1. Check if all required images exist locally
2. If all local → **instant startup** without pulling
3. If some missing → check remote accessibility → pull only missing images
4. If remote inaccessible → automatically fall back to local build

This ensures the fastest possible deployment while maintaining reliability.

---

## Troubleshooting

### Submodule Branch Mismatch After Installation

If you have already installed CellCraft using docker compose but the plugin submodule is pointing to the wrong branch, follow these steps to fix it:

1. Switch to the correct submodule branch:

   **For GPU-enabled installation:**
   ```bash
   cd backend/plugin/official
   git switch release/plugins-v1.0
   cd ../../..
   ```

   **For CPU-only installation:**
   ```bash
   cd backend/plugin/official
   git switch release/plugins-v1.0-cpu
   cd ../../..
   ```

2. Stop the services and remove volumes to clear old plugin data:

   ```bash
   # For GPU (AMD64)
   docker compose -f docker-compose.gpu.amd64.yml down -v
   # For CPU (AMD64)
   docker compose -f docker-compose.cpu.amd64.yml down -v
   # For CPU (ARM64 - Apple Silicon/ARM Linux)
   docker compose -f docker-compose.cpu.arm64.yml down -v
   ```

3. Restart the services:

   **For GPU-enabled setup (AMD64):**
   ```bash
   docker compose -f docker-compose.gpu.amd64.yml up -d --build
   ```

   **For CPU-only setup (AMD64):**
   ```bash
   docker compose -f docker-compose.cpu.amd64.yml up -d --build
   ```

   **For CPU-only setup (ARM64):**
   ```bash
   docker compose -f docker-compose.cpu.arm64.yml up -d --build
   ```

This process ensures that:
- Plugin metadata is properly initialized in the database
- Plugin Docker images are correctly pulled from GHCR
- All plugin configurations are synchronized with the correct version

---

## Tutorials

To help you get started with CellCraft, we have prepared step-by-step tutorial videos. These tutorials cover everything from setting up your environment to performing **GRN analysis**.

📺 **Watch our tutorial series on YouTube**: [CellCraft Tutorial Playlist](https://www.youtube.com/@CellCraft-cislab)

### What You Will Learn:

1. **Exploring the Main Page** - An overview of the main page and its key contents.

   ![Exploring the Main Page](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2FHl6XamxlUoSXKnEElLbY%2Ftuto_main.gif?alt=media&token=d2fd5fb3-af62-4816-980d-57f708994087)

2. **Managing Data** - How to organize and manage data for analysis.

   ![Managing Data](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2Fe45wYiVaIBFnkeWyfSSq%2Ftuto_DataUpload.gif?alt=media&token=87adc0b1-1053-4b65-8540-a67efb5584ce)

3. **Configuring the Workflow** - Setting up the workflow before executing tasks.

   ![Configuring the Workflow](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2FKkQzRTvRyK7HkJxm2atX%2Ftuto_lasso.gif?alt=media&token=eb804547-f2fd-4e36-bdec-4ef30f3e7350)

4. **Executing the Task** - Running tasks and monitoring their progress.

   ![Executing the Task](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2Fe91usDzgphuq4hI0QQuE%2Ftuto_executeTask.gif?alt=media&token=34d65e28-8f6c-4b3d-86e4-e2f0884a2302)

5. **Analyzing the Results** - Interpreting and analyzing data based on output files.

   ![Analyzing the Results](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2FbDyVupxC3auhlGNOsWdG%2Ftuto_barplot.gif?alt=media&token=3956a66e-fb0c-418a-ab2b-91c558b4ed93)

---

## Platform-Specific Considerations

### Mac Apple Silicon Support

CellCraft provides **optimized support** for Mac Apple Silicon (M1/M2/M3/M4) with native ARM64 builds.

#### Quick Start for Mac Users

**Using the ARM64-optimized Docker Compose:**
```bash
docker compose -f docker-compose.cpu.arm64.yml up -d --build
```
---

### Plugin Compatibility

**Official Plugins**
- ✅ Fully supported on all platforms (Windows, Linux, macOS)
- Pre-configured and tested for cross-platform compatibility
- Includes: TENET, FastTENET, GENIE3, GRNBOOST2, LEAP, Scribe, and GRNViz

**Custom Local Plugins**
- ✅ Supported on Windows and Linux
- ❌ Not currently supported on macOS (planned for future updates)
- Allows users to create and integrate custom GRN analysis plugins

### System Requirements

Before installing CellCraft, please verify your system meets the following requirements:

| Component | Minimum | Recommended |
| :--- | :--- | :--- |
| CPU | 4 cores | 8+ cores |
| RAM | 8 GB | 16+ GB |
| Storage | 70 GB | 100+ GB |
| OS | Ubuntu 20.04 LTS, Window 10/11, macOS 10.15 Catalina | Ubuntu 22.04 LTS |
| OS Kernel | 6.8.0 | 6.8.0+ |
| glibc | 2.3.9 | 2.3.9+ |
| Docker | 20.10.0 | 24.0.0+ |
| Docker Compose | v2.0.0 | v2.20.0+ |
| NVIDIA Driver | 535.171.04 | 535.171.04+ |
| CUDA | 12.1 | 12.2+ |

**Additional Notes**:
- For GPU-enabled installation, use `./run-gpu-mode.sh`
- For CPU-only installation, use `./run-cpu-mode.sh`
- Docker Desktop (latest version) is required
- Git with submodule support is required
- Modern web browser (Chrome, Firefox, Edge, or Safari) is required

### Important Notes

**For macOS Users**
- Currently, only official plugins are available on macOS
- Custom local plugin development is limited to Windows and Linux environments
- We recommend using the comprehensive set of official plugins for your GRN analysis needs
- Full macOS support for custom plugins is planned for future releases

**Performance Considerations**
- GPU acceleration significantly improves performance for large-scale analyses
- CPU-only mode is suitable for smaller datasets and testing
- Always check your Docker resource allocations to ensure optimal performance

**Future Updates**
- macOS support for custom local plugins is under active development
- Additional official plugins will be added regularly
- Performance optimizations for all platforms are ongoing
