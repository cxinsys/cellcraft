# Sequential Test Execution Guide

## Overview

This guide provides step-by-step instructions for executing CellCraft deployment tests in a controlled, sequential manner. Follow this guide **before** running the automated platform-specific scripts.

**Purpose**: Validate each test category individually to identify issues early and understand system behavior at each stage.

**Estimated Total Time**: 3-4 hours (varies by platform)

---

## Pre-Execution Checklist

### 1. Environment Verification

**Required Software**:
```bash
# Check Docker
docker --version          # Should be 24.x+
docker compose version    # Should be v2+

# Check Python and pytest
python3 --version         # Should be 3.8+
pytest --version          # Should be installed

# Check available disk space
df -h .                   # Should have 50GB+ available
```

**Platform-Specific Checks**:

**Linux**:
```bash
# Check if running on native Linux
uname -s                  # Should output: Linux

# Optional: Check GPU availability
nvidia-smi                # If available, shows GPU info
```

**macOS**:
```bash
# Check architecture
uname -m                  # arm64 (Apple Silicon) or x86_64 (Intel)

# Verify Docker Desktop is running
docker info | grep "Operating System"  # Should show Docker Desktop
```

**WSL2**:
```bash
# Verify WSL2 environment
grep -i microsoft /proc/version         # Should contain 'microsoft'

# Check if project is in WSL filesystem (optimal)
pwd                                      # Should NOT start with /mnt/
```

### 2. Project Setup

```bash
# Navigate to project root
cd /home/dmshin/cellcraft

# Verify test files exist
ls tests/deployment/test_*.py
ls tests/deployment/conftest.py
ls tests/deployment/helpers.py

# Install test dependencies if needed
pip install pytest requests python-dotenv
```

### 3. Docker Environment Preparation

```bash
# Clean up any existing containers
docker compose -f docker-compose.cpu.yml down -v
docker compose -f docker-compose.gpu.yml down -v

# Verify no conflicting containers
docker ps -a | grep cellcraft

# Check available Docker resources
docker system df
```

---

## Sequential Test Execution

### Phase 1: System Installation Tests (15-20 minutes)

**Purpose**: Validate Docker Compose configuration and system installation process.

**Test Coverage**: TC-INIT-001 through TC-INIT-016
- Basic system installation
- Container orchestration
- Service initialization
- Configuration validation

**Execution**:
```bash
# Run only system installation tests
pytest tests/deployment/test_system_installation.py \
    -v \
    --tb=short \
    --log-cli-level=INFO

# Expected output:
# - 16 tests should pass
# - Containers should start and become healthy
# - Services should be accessible
```

**Validation Checkpoints**:
```bash
# Check container status
docker compose -f docker-compose.cpu.yml ps

# Verify services are accessible
curl http://localhost:8000/     # Backend API
curl http://localhost:8080/     # Frontend
curl http://localhost:15672/    # RabbitMQ UI (guest/guest)

# Check logs for errors
docker compose -f docker-compose.cpu.yml logs --tail=50
```

**Common Issues**:
- **Port conflicts**: Ensure ports 8000, 8080, 5432, 5672, 15672 are free
- **Permission errors**: Check Docker socket permissions
- **Slow startup**: Wait 30-60 seconds for services to stabilize

**Success Criteria**: ✅ All 16 tests pass, all services accessible

---

### Phase 2: CPU Mode Tests (20-25 minutes)

**Purpose**: Validate CPU-only deployment mode and workflow execution.

**Test Coverage**: TC-MODE-001 through TC-MODE-009 (CPU subset)
- CPU mode deployment
- Plugin discovery and loading
- Workflow execution without GPU
- Data processing pipelines

**Execution**:
```bash
# Run CPU mode tests only
pytest tests/deployment/test_cpu_gpu_mode.py::TestCPUMode \
    -v \
    --tb=short \
    --log-cli-level=INFO

# Expected output:
# - 9 tests should pass
# - CPU plugins should be loaded
# - Workflows should execute successfully
```

**Validation Checkpoints**:
```bash
# Check Celery worker status
docker compose -f docker-compose.cpu.yml exec celery \
    celery -A app.celery_app inspect active

# Verify plugin count
docker compose -f docker-compose.cpu.yml exec backend \
    ls -l /app/plugin/

# Check user_data directory
docker compose -f docker-compose.cpu.yml exec backend \
    ls -l /app/user_data/
```

**Platform-Specific Notes**:
- **macOS**: Expect slightly longer execution times (1.3x-1.5x multiplier)
- **WSL2**: Performance depends on filesystem location (/mnt vs native)
- **Linux**: Optimal performance, baseline timing

**Success Criteria**: ✅ All CPU mode tests pass, workflows execute correctly

---

### Phase 3: GPU Mode Tests (30-35 minutes, Optional)

**Purpose**: Validate GPU-accelerated deployment mode (if GPU available).

**Test Coverage**: TC-MODE-010 through TC-MODE-017 (GPU subset)
- GPU mode deployment
- CUDA container configuration
- GPU-accelerated workflows
- GPU plugin execution

**Prerequisites**:
```bash
# Verify GPU availability
nvidia-smi

# Check NVIDIA Container Toolkit
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi
```

**Execution**:
```bash
# Run GPU mode tests only
pytest tests/deployment/test_cpu_gpu_mode.py::TestGPUMode \
    -v \
    --tb=short \
    --log-cli-level=INFO

# Expected output:
# - 8 tests should pass
# - GPU plugins should be loaded
# - CUDA should be accessible in containers
```

**Validation Checkpoints**:
```bash
# Check GPU mode container status
docker compose -f docker-compose.gpu.yml ps

# Verify GPU access in container
docker compose -f docker-compose.gpu.yml exec backend nvidia-smi

# Check GPU plugin count
docker compose -f docker-compose.gpu.yml exec backend \
    ls -l /app/plugin/ | grep -i gpu
```

**Platform Availability**:
- **Linux**: ✅ Fully supported with NVIDIA GPU + Container Toolkit
- **WSL2**: ✅ Supported with CUDA on WSL2 setup
- **macOS**: ❌ Not supported (no NVIDIA GPU support)

**Success Criteria**: ✅ All GPU mode tests pass, GPU workflows execute

**Skip Conditions**: If no GPU available, tests will be skipped automatically

---

### Phase 4: Optimization Tests (25-30 minutes)

**Purpose**: Validate deployment optimizations and performance features.

**Test Coverage**: TC-OPT-001 through TC-OPT-006
- GHCR fallback mechanism
- Resource usage optimization
- Memory management
- Parallel container startup
- Image caching
- Startup time optimization

**Execution**:
```bash
# Run optimization tests
pytest tests/deployment/test_optimization.py \
    -v \
    --tb=short \
    --log-cli-level=INFO

# Expected output:
# - 6 tests should pass
# - GHCR fallback should work
# - Resource usage should be within limits
```

**Validation Checkpoints**:
```bash
# Monitor resource usage during tests
docker stats --no-stream

# Check Docker image layers
docker images | grep cellcraft
docker history cellcraft-backend:latest | head -10

# Verify startup timing
time ./run-cpu-mode.sh --skip-verify
```

**Performance Expectations**:
- **GHCR Fallback**: Should succeed via either GHCR or local build
- **Memory Usage**: Backend <500MB, Frontend <200MB, DB <300MB
- **Startup Time**:
  - Linux: <60s
  - macOS Intel: <78s (1.3x multiplier)
  - macOS Apple Silicon: <90s (1.5x multiplier)
  - WSL2: <72s (1.2x multiplier)

**Success Criteria**: ✅ All optimization tests pass, performance within targets

---

### Phase 5: Error Handling Tests (20-25 minutes)

**Purpose**: Validate system resilience and recovery mechanisms.

**Test Coverage**: TC-ERR-001 through TC-ERR-006
- Container failure recovery
- Database connection loss handling
- Network partition resilience
- Memory pressure handling
- Infinite loop prevention
- Disk space handling

**Execution**:
```bash
# Run error handling tests
pytest tests/deployment/test_error_handling.py \
    -v \
    --tb=short \
    --log-cli-level=INFO

# Expected output:
# - 6 tests should pass
# - System should recover from injected failures
# - No data corruption or crashes
```

**Validation Checkpoints**:
```bash
# Check container restart counts
docker compose -f docker-compose.cpu.yml ps

# Verify service recovery
curl http://localhost:8000/     # Should respond after recovery

# Check logs for error handling
docker compose -f docker-compose.cpu.yml logs --tail=100 backend
```

**Safety Notes**:
- These tests intentionally inject failures
- Containers will be stopped and restarted
- Database connections will be disrupted temporarily
- Network issues may be simulated

**Recovery Validation**:
- Backend should remain running during database disconnection
- Services should reconnect within 30 seconds after recovery
- No data loss or corruption
- Celery workers should remain responsive

**Success Criteria**: ✅ All error handling tests pass, system recovers gracefully

---

### Phase 6: Performance Tests (30-35 minutes)

**Purpose**: Validate system performance under load and stress conditions.

**Test Coverage**: TC-PERF-001 through TC-PERF-004
- Concurrent workflow execution
- Large dataset processing
- System resource limits
- Stress testing

**Execution**:
```bash
# Run performance tests
pytest tests/deployment/test_performance.py \
    -v \
    --tb=short \
    --log-cli-level=INFO

# Expected output:
# - 4 tests should pass
# - System should handle load without degradation
# - Performance metrics should be within acceptable ranges
```

**Validation Checkpoints**:
```bash
# Monitor system resources during performance tests
docker stats

# Check for memory leaks
docker compose -f docker-compose.cpu.yml exec backend \
    ps aux --sort=-%mem | head -10

# Verify database connection pool
docker compose -f docker-compose.cpu.yml exec db \
    psql -U cellcraft -c "SELECT count(*) FROM pg_stat_activity;"
```

**Performance Targets**:
- **Concurrent Workflows**: 3-5 workflows simultaneously
- **Large Datasets**: Handle 10K+ cell datasets without crashes
- **Memory Limits**: Stay within container memory limits
- **Response Time**: API responses <500ms under normal load

**Platform Adjustments**:
- Linux: Baseline performance targets
- macOS: 30-50% longer execution times expected
- WSL2: 20-30% longer if on Windows filesystem

**Success Criteria**: ✅ All performance tests pass, no degradation or crashes

---

## Post-Execution Validation

### 1. Final System Check

```bash
# Verify all services are still running
docker compose -f docker-compose.cpu.yml ps

# Check for any errors in logs
docker compose -f docker-compose.cpu.yml logs --tail=100 | grep -i error

# Verify service accessibility
curl -s http://localhost:8000/ | jq .
curl -s http://localhost:8080/ | grep -q "CellCraft"
```

### 2. Resource Cleanup

```bash
# View resource usage
docker system df

# Optional: Clean up test containers
docker compose -f docker-compose.cpu.yml down -v

# Optional: Remove test images
docker images | grep cellcraft | awk '{print $3}' | xargs docker rmi
```

### 3. Test Report Generation

```bash
# Generate HTML test report (if pytest-html installed)
pytest tests/deployment/ --html=test-report.html --self-contained-html

# Generate coverage report (if pytest-cov installed)
pytest tests/deployment/ --cov=. --cov-report=html
```

---

## Troubleshooting Guide

### Common Issues and Solutions

**Issue 1: Port Already in Use**
```bash
# Find process using port
lsof -i :8000  # or :8080, :5432, :5672, :15672

# Kill process or use different ports
docker compose -f docker-compose.cpu.yml down
```

**Issue 2: Container Fails to Start**
```bash
# Check detailed logs
docker compose -f docker-compose.cpu.yml logs backend

# Inspect container configuration
docker compose -f docker-compose.cpu.yml config

# Rebuild without cache
docker compose -f docker-compose.cpu.yml build --no-cache
```

**Issue 3: Database Connection Errors**
```bash
# Check database container status
docker compose -f docker-compose.cpu.yml exec db pg_isready

# Reset database
docker compose -f docker-compose.cpu.yml down -v
docker compose -f docker-compose.cpu.yml up -d db
```

**Issue 4: Test Timeouts**
```bash
# Increase timeout for your platform
# Edit conftest.py platform_constraints fixture
# Adjust timeout_multiplier value

# Or run tests with longer timeout
pytest tests/deployment/ --timeout=600
```

**Issue 5: WSL2 Slow Performance**
```bash
# Check if project is in Windows filesystem
pwd  # Should NOT start with /mnt/

# If in /mnt/, move to WSL filesystem
mkdir -p ~/cellcraft
cp -r /mnt/c/path/to/cellcraft ~/cellcraft/
cd ~/cellcraft
```

**Issue 6: GPU Tests Failing**
```bash
# Verify NVIDIA drivers
nvidia-smi

# Check NVIDIA Container Toolkit
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi

# Reinstall NVIDIA Container Toolkit if needed
# https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html
```

---

## Platform-Specific Execution Notes

### Linux (Native Docker Engine)

**Optimal Configuration**:
- Native Docker Engine 24.x+
- Direct hardware access
- Best performance (baseline timing)

**Expected Runtime**: 2-2.5 hours total

**Command Reference**:
```bash
# Full sequential execution
pytest tests/deployment/ -v --tb=short --log-cli-level=INFO
```

### macOS (Docker Desktop)

**Configuration Notes**:
- Docker Desktop with VM overhead
- Architecture detection (Intel vs Apple Silicon)
- No GPU support

**Expected Runtime**:
- Intel: 2.5-3 hours (1.3x multiplier)
- Apple Silicon: 3-3.5 hours (1.5x multiplier)

**Architecture Detection**:
```bash
if [ "$(uname -m)" = "arm64" ]; then
    echo "Apple Silicon detected"
else
    echo "Intel detected"
fi
```

### WSL2 (Windows Subsystem for Linux)

**Configuration Notes**:
- WSL2 with Docker Desktop or native Docker Engine
- Filesystem location critical for performance
- Optional GPU support with CUDA on WSL2

**Expected Runtime**: 2.5-3 hours (1.2x multiplier)

**Performance Optimization**:
```bash
# Verify WSL filesystem location
pwd | grep -v "^/mnt/" && echo "Optimal location" || echo "Suboptimal - move to ~/cellcraft"
```

---

## Success Criteria Summary

| Phase | Test Count | Expected Duration | Success Threshold |
|-------|-----------|-------------------|-------------------|
| System Installation | 16 | 15-20 min | 100% pass |
| CPU Mode | 9 | 20-25 min | 100% pass |
| GPU Mode (optional) | 8 | 30-35 min | 100% pass (if GPU) |
| Optimization | 6 | 25-30 min | 100% pass |
| Error Handling | 6 | 20-25 min | 100% pass |
| Performance | 4 | 30-35 min | 100% pass |
| **TOTAL** | **49** | **2-4 hours** | **100% pass** |

---

## Next Steps

After successfully completing sequential test execution:

1. **Review Test Logs**: Check for any warnings or performance issues
2. **Run Automated Scripts**: Use platform-specific scripts for full automation
   - `./tests/deployment/run-tests-linux.sh`
   - `./tests/deployment/run-tests-macos.sh`
   - `./tests/deployment/run-tests-wsl2.sh`
3. **Generate Reports**: Create comprehensive test reports for documentation
4. **CI/CD Integration**: Set up automated testing in GitHub Actions
5. **Performance Benchmarking**: Establish baseline performance metrics

---

## Additional Resources

- **Test Plan**: `docs/testing/DEPLOYMENT_TEST_PLAN_PART2.md`
- **Fixtures Reference**: `tests/deployment/conftest.py`
- **Helper Functions**: `tests/deployment/helpers.py`
- **Docker Configs**: `docker-compose.cpu.yml`, `docker-compose.gpu.yml`

---

## Appendix: Quick Reference Commands

```bash
# Single test execution
pytest tests/deployment/test_system_installation.py::TestSystemInstallation::test_basic_system_installation -v

# Run specific test class
pytest tests/deployment/test_cpu_gpu_mode.py::TestCPUMode -v

# Run with specific marker
pytest tests/deployment/ -m cpu_mode -v

# Skip slow tests
pytest tests/deployment/ -m "not slow" -v

# Run with coverage
pytest tests/deployment/ --cov=. --cov-report=term-missing

# Debug mode with print statements
pytest tests/deployment/ -v -s

# Stop on first failure
pytest tests/deployment/ -x

# Run last failed tests
pytest tests/deployment/ --lf
```
