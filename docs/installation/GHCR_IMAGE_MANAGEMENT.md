# GHCR Image Management

CellCraft provides **pre-built Docker images** hosted on GitHub Container Registry (GHCR) for faster deployment. The `test-ghcr-check.sh` script helps you manage these images efficiently.

## Available Container Images

| Image | Description | Platforms | GHCR Package | Latest Version |
|-------|-------------|-----------|--------------|----------------|
| **frontend** | Vue.js frontend application | AMD64, ARM64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Ffrontend) | latest |
| **backend-cpu** | FastAPI backend (CPU-only) | AMD64, ARM64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fbackend-cpu) | latest |
| **backend-gpu** | FastAPI backend (GPU-enabled) | AMD64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fbackend-gpu) | latest |
| **celery-cpu** | Celery worker (CPU-only) | AMD64, ARM64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fcelery-cpu) | latest |
| **celery-gpu** | Celery worker (GPU-enabled) | AMD64 | [View Package](https://github.com/cxinsys/cellcraft/pkgs/container/cellcraft%2Fcelery-gpu) | latest |

## Quick Pull Commands

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

## Image Management Tool

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

## Image Status Indicators

The script displays clear status for each image:
- ✅ **LOCAL**: Image exists locally (instant startup, no download needed)
- ⚠️ **REMOTE**: Image available remotely (will be downloaded when needed)
- ❌ **MISSING**: Image not accessible (will fall back to local build)

## Pre-downloading Images

Pre-downloading images is recommended for:
- **Faster first-time deployment** (no wait during startup)
- **Offline environments** (download once, deploy anytime)
- **Network-constrained setups** (download during off-peak hours)

The script intelligently skips images that already exist locally and provides detailed download statistics.

## Smart Image Detection

Both `run-cpu-mode.sh` and `run-gpu-mode.sh` implement **local-first checking**:
1. Check if all required images exist locally
2. If all local → **instant startup** without pulling
3. If some missing → check remote accessibility → pull only missing images
4. If remote inaccessible → automatically fall back to local build

This ensures the fastest possible deployment while maintaining reliability.
