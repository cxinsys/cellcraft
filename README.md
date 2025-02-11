<img src="https://github.com/cxinsys/cellcraft/blob/807998fda59e15e185ea9d2835ff7b81a884460f/frontend/src/assets/cellcraft_logo_text.png"/>

[Demo Website](http://165.194.161.183:10001/cellcraft) • [Docs](https://labelstud.io/guide/)

## Overview

CellCraft is a web-based platform designed for researchers to reconstruct Gene Regulatory Networks (GRNs) efficiently. It provides a **visual programming interface**, enabling non-programmers to conduct complex GRN analysis. Powered by **TENET** and **FastTENET**, the platform ensures high-performance GRN reconstruction. Docker-based deployment simplifies setup and facilitates seamless local server operations.

## Key Features
- **Visual Programming Interface**: Streamline data analysis with an intuitive drag-and-drop GUI.
- **Powerful GRN Tools**: Built-in support for TENET and FastTENET for cutting-edge GRN inference.
- **Plugin Support**: Extend functionality by integrating custom plugins.
- **Data Management**: Effortlessly manage and share data in collaborative environments.
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
