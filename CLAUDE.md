# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Environment Setup
```bash
# Start development environment
docker compose -f docker-compose.dev.yml up --build

# Start production environment  
docker compose -f docker-compose.prod.yml up --build

# Access the application at http://localhost:8080
```

### Frontend Development
```bash
# Frontend commands (run in /frontend directory)
npm run serve    # Development server
npm run build    # Production build
npm run lint     # ESLint

# Frontend is Vue.js 2.6 with Vue Router and Vuex
```

### Backend Development
```bash
# Database migrations
alembic upgrade head              # Apply migrations
alembic revision --autogenerate  # Create new migration

# Backend runs on FastAPI with uvicorn
# Main application: backend/app/main.py
```

### Plugin Development
```bash
# Plugin structure in backend/plugin/[PLUGIN_NAME]/
# Each plugin contains:
# - metadata.json (UI flow definition)
# - Snakefile (workflow definition)
# - scripts/ (Python/R implementation)
# - dependency/ (requirements.txt or renv.lock)
```

## Architecture Overview

### Core Technologies
- **Frontend**: Vue.js 2.6, Vue Router, Vuex, Drawflow (visual workflow editor)
- **Backend**: FastAPI, SQLAlchemy, PostgreSQL, Celery (task queue)
- **Containers**: Docker + Docker Compose, GPU support via NVIDIA
- **Workflow**: Snakemake for bioinformatics pipeline execution
- **Message Queue**: RabbitMQ for async task processing

### Key Components

#### Database Models (`backend/app/database/models.py`)
- **User**: Authentication and user management
- **File**: User file storage with metadata
- **Workflow**: Visual workflow definitions (JSONB)
- **Task**: Execution tracking with status/logging
- **Plugin**: Plugin metadata and configuration

#### Plugin System
- **Visual Flow Definition**: JSON-based node configuration in metadata.json
- **Workflow Execution**: Snakemake-based pipeline with parameterized rules
- **Multi-language Support**: Python and R scripts
- **Containerized Execution**: Docker isolation for plugins
- **Visualization Components**: Plotly.js charts (heatmaps, bar plots, networks)

#### API Structure (`backend/app/routes/endpoints/`)
- **Authentication**: JWT-based user auth
- **File Management**: Upload/download with type validation
- **Workflow Management**: Visual workflow CRUD operations
- **Task Management**: Execution tracking and monitoring
- **Plugin Management**: Dynamic plugin loading and configuration

#### Frontend Structure (`frontend/src/`)
- **Visual Programming**: Drawflow-based drag-and-drop interface
- **State Management**: Vuex with persistence
- **Data Visualization**: Plotly.js integration
- **Component Architecture**: Vue.js 2.6 components with routing

### Development Workflow

1. **Frontend Changes**: Edit Vue components, run `npm run serve` for hot reload
2. **Backend Changes**: Edit Python files, Docker dev container auto-reloads
3. **Database Changes**: Create Alembic migrations, apply with `alembic upgrade head`  
4. **Plugin Development**: Create new plugin directory with metadata.json and Snakefile
5. **Testing**: No automated test framework configured - manual testing required

### Important Notes

#### Plugin System Details
- Each plugin defines a complete workflow via metadata.json
- Plugins support both computation and visualization nodes
- Visual flows are stored as JSONB in workflows table
- Snakemake handles parameter templating and rule execution
- Plugin containers inherit GPU access when configured

#### Data Management
- User files stored in backend/user/[user_id]/ directory
- Workflow definitions stored as JSONB in PostgreSQL
- Plugin execution results stored in task-specific directories
- File type validation enforced via metadata.json extensions

#### Environment Variables
Required in .env file:
- POSTGRES_USER, POSTGRES_PASSWORD, POSTGRES_HOST, POSTGRES_PORT, POSTGRES_DB
- TZ=Asia/Seoul (timezone configuration)
- GPU_COUNT for NVIDIA GPU allocation

#### Security Considerations
- JWT tokens for API authentication
- Docker container isolation for plugin execution
- Database connection pooling via SQLAlchemy
- CORS middleware configured for frontend access