<img src="https://github.com/cxinsys/cellcraft/blob/807998fda59e15e185ea9d2835ff7b81a884460f/frontend/src/assets/cellcraft_logo_text.png"/>

[Demo Website](https://cellcraft.app) • [Docs](https://cellcraft.gitbook.io/cellcraft-docs)

### Plugin Resources

[**Official Plugins**](https://github.com/cxinsys/cellcraft-plugin) • [**Plugin Templates**](https://github.com/cxinsys/cellcraft-plugin-templates)

- **Official Plugins Repository**: Contains version information and source code for all officially supported CellCraft plugins
- **Plugin Templates Repository**: Provides development resources and templates for creating custom local plugins

## Overview

**CellCraft** is a web-based application for reconstructing **gene regulatory networks (GRNs)** from single-cell RNA sequencing (scRNA-seq) data. It features an intuitive visual interface that integrates seven GRN inference tools—including **TENET** and **FastTENET**, developed by our research team—through modular plugin integration. 

Built to **lower technical barriers** in GRN analysis, CellCraft enables researchers to **configure workflows**, **run analysis**, and **explore results visually**, all without needing programming knowledge. The platform unifies access to multiple inference methods behind a consistent, guided interface and supports scalable, customizable workflows suitable for both rapid prototyping and extensive studies.


## Key Features

- **Visual Workflow Programming**: Configure and execute complex GRN analyses through an intuitive drag-and-drop interface designed for workflow-based programming.
- **Seven Integrated GRN Tools**: Built-in support for **TENET**, **FastTENET**, **GENIE3**, **GRNBOOST2**, **LEAP**, **Scribe**, and visualization tools, with continuous expansion of capabilities.
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
   ```

2. Start the application:

   **For GPU-enabled installation:**
   ```bash
   cd cellcraft && ./run-gpu-mode.sh
   ```

   **For CPU-only installation:**
   ```bash
   cd cellcraft && ./run-cpu-mode.sh
   ```

   **If the scripts fail to execute, use these manual commands:**
   
   For GPU-enabled setup:
   ```bash
   cd backend/plugin/official && git switch release/plugins-v1.0
   cd ../../.. && docker compose -f docker-compose.gpu.yml up --build
   ```

   For CPU-only setup:
   ```bash
   cd backend/plugin/official && git switch release/plugins-v1.0-cpu
   cd ../../.. && docker compose -f docker-compose.cpu.yml up --build
   ```

3. Access the application at [http://localhost:8080](http://localhost:8080).

4. Check the installation status:

   ```bash
   ./check-installation.sh
   ```

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
