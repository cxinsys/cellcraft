# CellCraft CPU-Only Mode Usage Guide

## Overview

The `run-cpu-mode.sh` script provides an automated way to launch CellCraft in CPU-only mode, automatically handling the plugin branch switching and Docker configuration.

## Quick Start

```bash
# Basic launch
./run-cpu-mode.sh

# Clean launch (recommended for first time)
./run-cpu-mode.sh --clean

# Force rebuild (when you've made changes)
./run-cpu-mode.sh --build

# Full clean and rebuild
./run-cpu-mode.sh --clean --build
```

## What the Script Does

1. **Prerequisites Check**: Verifies Docker, Git, and required files are available
2. **Git State Backup**: Saves current plugin branch state for rollback if needed
3. **Branch Switching**: Automatically switches to `release/plugins-v1.0-cpu` branch
4. **Container Management**: Optionally cleans existing containers
5. **Service Launch**: Builds and starts all required services
6. **Verification**: Runs installation checks to ensure everything works
7. **Status Report**: Provides access URLs and management commands

## Command Options

| Option | Description |
|--------|-------------|
| `--clean` | Remove existing containers before launching |
| `--build` | Force rebuild of all Docker images |
| `--skip-verify` | Skip the final installation verification |
| `--help` | Show usage information |

## CPU-Only Configuration

In CPU-only mode, the following plugins are available:
- ✅ **GENIE3** - Gene regulatory network inference
- ✅ **GRNBoost2** - Gradient boosting for GRN inference  
- ✅ **GRNViz** - Network visualization
- ✅ **LEAP** - Pseudotime trajectory analysis
- ✅ **Scribe** - RNA velocity analysis
- ✅ **TENET** - Tumor evolution network analysis

The following GPU-only plugins are excluded:
- ❌ **FastSCODE** - GPU-accelerated single-cell analysis
- ❌ **FastTENET** - GPU-accelerated tumor evolution analysis

## Access Points

After successful launch, CellCraft will be available at:

- **Frontend Application**: http://localhost:8080
- **Backend API**: http://localhost:8000
- **API Documentation**: http://localhost:8000/docs  
- **RabbitMQ Management**: http://localhost:15672 (guest/guest)

## Management Commands

```bash
# View logs
docker-compose -f docker-compose.cpu.yml logs

# View logs for specific service
docker-compose -f docker-compose.cpu.yml logs backend

# Stop services
docker-compose -f docker-compose.cpu.yml down

# Restart services  
docker-compose -f docker-compose.cpu.yml restart

# Restart specific service
docker-compose -f docker-compose.cpu.yml restart backend
```

## Troubleshooting

### Script Fails to Start
```bash
# Check Docker is running
docker info

# Check if you're in the correct directory
ls -la docker-compose.cpu.yml

# Check plugin directory exists
ls -la backend/plugin/official/
```

### Services Don't Start Properly
```bash
# Run with clean option
./run-cpu-mode.sh --clean --build

# Check detailed logs
docker-compose -f docker-compose.cpu.yml logs --tail=100
```

### Branch Switching Issues
```bash
# Check current plugin branch
cd backend/plugin/official && git status

# Manually fetch updates
cd backend/plugin/official && git fetch origin

# Check available branches
cd backend/plugin/official && git branch -a
```

### Container Health Issues
```bash
# Use the built-in check script
./check-installation.sh

# Or check individual container health
docker-compose -f docker-compose.cpu.yml ps
```

## Recovery

If something goes wrong, the script automatically restores the original git state. You can also manually restore:

```bash
cd backend/plugin/official
git status  # Check current state
git checkout main  # Return to main branch
```

## Performance Notes

**Initial Launch**: First run may take 10-15 minutes for image downloads and builds.

**Subsequent Launches**: Usually complete in 2-3 minutes using cached images.

**System Requirements**:
- Docker Desktop with at least 4GB RAM allocated
- 10GB free disk space for images and containers
- Internet connection for initial image downloads

## Support

If you encounter issues:
1. Check the troubleshooting section above
2. Run `./check-installation.sh` for detailed diagnostics
3. Review container logs for specific error messages
4. Ensure your system meets the minimum requirements

For additional help, refer to the main CellCraft documentation.