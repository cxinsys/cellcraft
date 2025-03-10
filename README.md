<img src="https://github.com/cxinsys/cellcraft/blob/807998fda59e15e185ea9d2835ff7b81a884460f/frontend/src/assets/cellcraft_logo_text.png"/>

[Demo Website](http://165.194.161.183:10001/cellcraft) • [Docs](https://cellcraft.gitbook.io/cellcraft-docs)

## Overview

CellCraft is a web-based platform designed for researchers to reconstruct Gene Regulatory Networks (GRNs) efficiently. It provides a **visual programming interface**, enabling non-programmers to conduct complex GRN analysis. Powered by **TENET** and **FastTENET**, the platform ensures high-performance GRN reconstruction. Docker-based deployment simplifies setup and facilitates seamless local server operations.

## Key Features

- **Visual Programming Interface**: Streamline data analysis with an intuitive drag-and-drop GUI.
- **Extensible GRN Tools**: Built-in support for **TENET**, **FastTENET**, and other GRN inference tools, with more plugins continuously being added to expand functionality.
- **Plugin Support**: Extend functionality by integrating custom plugins.
- **Data Management**: Users can manage their own data through their accounts, upload datasets.
- **Docker-Ready Deployment**: Deploy with ease using Docker containers for consistency and scalability.

---

## Getting Started

1. Clone the repository:

   ```bash
   git clone https://github.com/cellcraft.git
   cd cellcraft
   ```

2. Set environment variables in an `.env` file:

   ```dotenv
   POSTGRES_USER = your_user
   POSTGRES_PASSWORD = your_password
   POSTGRES_HOST = localhost
   POSTGRES_PORT = 5432
   POSTGRES_DB = cellcraft
   ```

3. Start the application:

   ```bash
   docker compose -f docker-compose.prod.yml up --build
   ```

4. Access the platform at [http://localhost:8080](http://localhost:8080).

## Tutorials

To help you get started with CellCraft, we have prepared step-by-step tutorial videos. These tutorials cover everything from setting up your environment to performing **GRN analysis**.

📺 **Watch our tutorial series on YouTube**: [CellCraft Tutorial Playlist](https://www.youtube.com/playlist?list=PLN8_i4yGKekju3EJClmRvqe4pL8xJr4hw)

### What You Will Learn:
1. **Introduction to CellCraft** - Overview of the platform and its key features.
2. **Setting Up Your Environment** - Installing dependencies and configuring your workspace.
3. **Building a GRN Workflow** - Creating and managing an analysis pipeline using the visual programming interface.
4. **Running GRN Analysis Tools** - Step-by-step guide to executing GRN inference with TENET and FastTENET.
5. **Visualizing Results** - Understanding and interpreting generated networks, heatmaps, and bar plots.
