<img src="https://github.com/cxinsys/cellcraft/blob/807998fda59e15e185ea9d2835ff7b81a884460f/frontend/src/assets/cellcraft_logo_text.png"/>

[Demo Website](http://165.194.161.183:10001/cellcraft) • [Docs](https://cellcraft.gitbook.io/cellcraft-docs)

## Overview

CellCraft is a web-based application designed for researchers to efficiently reconstruct Gene Regulatory Networks (GRNs). It provides a **visual programming interface**, enabling non-programmers to conduct complex GRN analyses with ease. The platform includes built-in support for **TENET**, **FastTENET**, and various other GRN inference tools, with more plugins continuously being added to enhance its capabilities. Additionally, users can manage their own datasets through their accounts, upload files, and integrate them into projects for analysis. Docker-based deployment simplifies setup and ensures seamless local server operations.


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

1. **Exploring the Main Page** - An overview of the main page and its key contents.  
   ![Exploring the Main Page](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2F5a83O2fqqqwICe7citLT%2Ftuto_main.gif?alt=media&token=70b867bd-0674-4660-b6c1-5eba50bb75ea)

2. **Managing Data** - How to organize and manage data for analysis.  
   ![Managing Data](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2Fe45wYiVaIBFnkeWyfSSq%2Ftuto_DataUpload.gif?alt=media&token=87adc0b1-1053-4b65-8540-a67efb5584ce)

3. **Configuring the Workflow** - Setting up the workflow before executing tasks.  
   ![Configuring the Workflow](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2FKkQzRTvRyK7HkJxm2atX%2Ftuto_lasso.gif?alt=media&token=eb804547-f2fd-4e36-bdec-4ef30f3e7350)

4. **Executing the Task** - Running tasks and monitoring their progress.  
   ![Executing the Task](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2Fe91usDzgphuq4hI0QQuE%2Ftuto_executeTask.gif?alt=media&token=34d65e28-8f6c-4b3d-86e4-e2f0884a2302)

5. **Analyzing the Results** - Interpreting and analyzing data based on output files.  
   ![Analyzing the Results](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2FjRZEd1fcjAhaS66UWnMw%2Fuploads%2FbDyVupxC3auhlGNOsWdG%2Ftuto_barplot.gif?alt=media&token=3956a66e-fb0c-418a-ab2b-91c558b4ed93)
