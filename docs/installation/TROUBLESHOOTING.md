# Troubleshooting

Common issues and solutions for CellCraft installation and operation.

## Submodule Branch Mismatch After Installation

If you have already installed CellCraft using docker compose but the plugin submodule is pointing to the wrong branch, follow these steps to fix it:

### 1. Switch to the correct submodule branch

**For GPU-enabled installation:**
```bash
cd backend/plugin/official
git switch release/plugins-v1.1
cd ../../..
```

**For CPU-only installation:**
```bash
cd backend/plugin/official
git switch release/plugins-v1.0-cpu
cd ../../..
```

### 2. Stop the services and remove volumes to clear old plugin data

```bash
# For GPU (AMD64)
docker compose -f docker-compose.gpu.amd64.yml down -v
# For CPU (AMD64)
docker compose -f docker-compose.cpu.amd64.yml down -v
# For CPU (ARM64 - Apple Silicon/ARM Linux)
docker compose -f docker-compose.cpu.arm64.yml down -v
```

### 3. Restart the services

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

## Plugin Download Issues

### Backend API Not Responding During Startup

**Symptom**: Web interface loads at `http://localhost:8080` but login/signup doesn't work.

**Cause**: Plugins are still being downloaded in the background.

**Solution**:
1. Monitor plugin download progress:
   ```bash
   docker compose -f docker-compose.gpu.amd64.yml logs -f backend
   ```
2. Wait until you see `Uvicorn running on http://0.0.0.0:8000`
3. Plugin download typically takes 10-15 minutes depending on network speed

### Plugin Images Not Found

**Symptom**: Error messages about missing plugin Docker images.

**Solution**:
1. Check if plugin images exist:
   ```bash
   docker images | grep -E "(tenet|fasttenet|fastscode|genie3|grnboost2|leap|scribe|grnviz)"
   ```
2. If missing, restart the backend service to trigger re-download:
   ```bash
   docker compose -f [compose-file] restart backend
   ```

## Docker Issues

### Permission Denied Errors

**Symptom**: `permission denied while trying to connect to the Docker daemon socket`

**Solution**:
```bash
# Add user to docker group
sudo usermod -aG docker $USER
# Log out and log back in, or run:
newgrp docker
```

### Port Already in Use

**Symptom**: `Bind for 0.0.0.0:8080 failed: port is already allocated`

**Solution**:
```bash
# Find process using the port
sudo lsof -i :8080
# Stop the conflicting process or use a different port
```

### Out of Disk Space

**Symptom**: Docker build fails with disk space errors.

**Solution**:
```bash
# Clean up unused Docker resources
docker system prune -a
# Remove unused volumes
docker volume prune
```

## Database Issues

### Database Connection Failed

**Symptom**: Backend fails to start with database connection errors.

**Solution**:
1. Check if the database container is running:
   ```bash
   docker ps | grep db
   ```
2. Check database logs:
   ```bash
   docker compose -f [compose-file] logs db
   ```
3. If corrupt, remove and recreate:
   ```bash
   docker compose -f [compose-file] down -v
   docker compose -f [compose-file] up -d
   ```

## Getting Help

If you encounter issues not covered here:

1. Check the [GitHub Issues](https://github.com/cxinsys/cellcraft/issues)
2. Review backend logs: `docker compose -f [compose-file] logs backend`
3. Open a new issue with:
   - Your OS and Docker version
   - Complete error messages
   - Steps to reproduce
