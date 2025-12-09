# Platform-Specific Considerations

Platform-specific guides and notes for running CellCraft on different operating systems.

## Mac Apple Silicon Support

CellCraft provides **optimized support** for Mac Apple Silicon (M1/M2/M3/M4) with native ARM64 builds.

### Quick Start for Mac Users

**Using the ARM64-optimized Docker Compose:**
```bash
docker compose -f docker-compose.cpu.arm64.yml up -d --build
```

### Notes for Mac Users

- GPU mode is not available on macOS (NVIDIA CUDA is not supported)
- Use the ARM64-optimized compose file for best performance
- Docker Desktop for Mac is required

## Windows Support

### Prerequisites

1. Install [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop)
2. Enable WSL 2 backend (recommended)
3. Allocate sufficient resources in Docker Desktop settings

### Running on Windows

```bash
# CPU mode
docker compose -f docker-compose.cpu.amd64.yml up -d --build

# GPU mode (requires NVIDIA GPU with WSL 2)
docker compose -f docker-compose.gpu.amd64.yml up -d --build
```

## Linux Support

### Ubuntu/Debian

Follow the standard installation instructions in the [main README](../../README.md).

### GPU Support on Linux

For GPU-enabled installation, ensure you have:
- NVIDIA Driver 535.171.04+
- CUDA 12.1+
- nvidia-container-toolkit installed

## Plugin Compatibility

### Official Plugins

- ✅ Fully supported on all platforms (Windows, Linux, macOS)
- Pre-configured and tested for cross-platform compatibility
- Includes: TENET, FastTENET, GENIE3, GRNBoost2, LEAP, Scribe, and GRNViz

### Custom Local Plugins

- ✅ Supported on Windows and Linux
- ❌ Not currently supported on macOS (planned for future updates)
- Allows users to create and integrate custom GRN analysis plugins

## Important Notes

### For macOS Users

- Currently, only official plugins are available on macOS
- Custom local plugin development is limited to Windows and Linux environments
- We recommend using the comprehensive set of official plugins for your GRN analysis needs
- Full macOS support for custom plugins is planned for future releases

### Performance Considerations

- GPU acceleration significantly improves performance for large-scale analyses
- CPU-only mode is suitable for smaller datasets and testing
- Always check your Docker resource allocations to ensure optimal performance

### Docker Resource Recommendations

| Platform | CPU | Memory | Disk |
|----------|-----|--------|------|
| macOS (Apple Silicon) | 4+ cores | 8+ GB | 70+ GB |
| Windows (WSL 2) | 4+ cores | 8+ GB | 70+ GB |
| Linux | 4+ cores | 8+ GB | 70+ GB |

## Future Updates

- macOS support for custom local plugins is under active development
- Additional official plugins will be added regularly
- Performance optimizations for all platforms are ongoing
