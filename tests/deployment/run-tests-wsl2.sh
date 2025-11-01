#!/usr/bin/env bash
# WSL2 Platform Test Execution Script for CellCraft Deployment Tests
#
# Optimized for WSL2 with Docker Desktop or native Docker Engine
# Expected runtime: ~2.5-3 hours (timeout_multiplier: 1.2)
#
# Prerequisites:
#   - WSL2 with Ubuntu 20.04+ or Debian
#   - Docker Engine or Docker Desktop with WSL2 backend
#   - NVIDIA drivers + Container Toolkit (for GPU tests, if available)
#   - Python 3.8+ with pytest, requests
#   - 50GB+ disk space
#   - Execute from project root
#   - IMPORTANT: Project must be in WSL filesystem (not /mnt/c/)

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TEST_DIR="${PROJECT_ROOT}/tests/deployment"
LOG_DIR="${PROJECT_ROOT}/test-logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TIMEOUT_MULTIPLIER=1.2

# Create log directory
mkdir -p "${LOG_DIR}"

echo -e "${GREEN}=== CellCraft Deployment Tests (WSL2) ===${NC}"
echo "Timestamp: ${TIMESTAMP}"
echo "Project Root: ${PROJECT_ROOT}"
echo "Log Directory: ${LOG_DIR}"
echo "Timeout Multiplier: ${TIMEOUT_MULTIPLIER}x"
echo ""

# WSL2 Detection
if ! grep -qi microsoft /proc/version; then
    echo -e "${RED}Error: This script is for WSL2 only${NC}"
    echo "Detected system is not WSL2"
    exit 1
fi
echo "✓ WSL2 environment detected"

# Check if project is in WSL filesystem
if [[ "${PROJECT_ROOT}" == /mnt/* ]]; then
    echo -e "${YELLOW}========================================${NC}"
    echo -e "${YELLOW}WARNING: Project is in Windows filesystem (/mnt/)${NC}"
    echo -e "${YELLOW}========================================${NC}"
    echo ""
    echo "For best performance and compatibility, move project to WSL filesystem:"
    echo "  1. mkdir -p ~/cellcraft"
    echo "  2. cp -r ${PROJECT_ROOT} ~/cellcraft/"
    echo "  3. cd ~/cellcraft && ./tests/deployment/run-tests-wsl2.sh"
    echo ""
    echo "Current location may cause:"
    echo "  - 10-50x slower Docker volume I/O"
    echo "  - File permission issues"
    echo "  - Extended test execution times"
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Test execution cancelled"
        exit 0
    fi
else
    echo "✓ Project is in WSL filesystem (optimal)"
fi

# Check prerequisites
echo ""
echo -e "${YELLOW}Checking prerequisites...${NC}"

# Check Docker
if ! command -v docker &> /dev/null; then
    echo -e "${RED}Error: Docker not found${NC}"
    echo "Install Docker Engine: https://docs.docker.com/engine/install/ubuntu/"
    echo "Or use Docker Desktop with WSL2 backend"
    exit 1
fi
echo "✓ Docker $(docker --version)"

# Verify Docker is running
if ! docker info &> /dev/null; then
    echo -e "${RED}Error: Docker is not running${NC}"
    echo "If using Docker Desktop, ensure WSL2 integration is enabled"
    echo "If using Docker Engine, start it with: sudo service docker start"
    exit 1
fi
echo "✓ Docker daemon is running"

# Detect Docker engine type
DOCKER_ENGINE="Unknown"
if docker info 2>/dev/null | grep -q "Operating System.*Docker Desktop"; then
    DOCKER_ENGINE="Docker Desktop (WSL2 backend)"
elif docker info 2>/dev/null | grep -q "Server Version"; then
    DOCKER_ENGINE="Docker Engine"
fi
echo "✓ Docker engine: ${DOCKER_ENGINE}"

# Check Docker Compose
if ! docker compose version &> /dev/null; then
    echo -e "${RED}Error: Docker Compose v2 not found${NC}"
    exit 1
fi
echo "✓ Docker Compose $(docker compose version --short)"

# Check GPU availability (WSL2 GPU support)
GPU_AVAILABLE=false
if command -v nvidia-smi &> /dev/null; then
    if nvidia-smi &> /dev/null; then
        GPU_AVAILABLE=true
        echo "✓ NVIDIA GPU detected: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -n1)"
        echo "  Driver: $(nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -n1)"
    else
        echo "⚠ nvidia-smi found but failed (driver issue?)"
    fi
else
    echo "⚠ No NVIDIA GPU detected - GPU tests will be skipped"
    echo "  For WSL2 GPU support: https://docs.nvidia.com/cuda/wsl-user-guide/"
fi

# Check Python and pytest
if ! command -v pytest &> /dev/null; then
    echo -e "${RED}Error: pytest not found${NC}"
    echo "Install with: pip3 install pytest requests"
    exit 1
fi
echo "✓ pytest $(pytest --version | head -n1)"

# Check disk space
AVAILABLE_SPACE=$(df -BG "${PROJECT_ROOT}" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "${AVAILABLE_SPACE}" -lt 50 ]; then
    echo -e "${YELLOW}Warning: Only ${AVAILABLE_SPACE}GB available (recommended: 50GB+)${NC}"
fi
echo "✓ Disk space: ${AVAILABLE_SPACE}GB available"

# WSL2 specific checks
echo ""
echo -e "${YELLOW}WSL2 Platform Notes:${NC}"
echo "  - Docker Desktop backend: expect ~20% overhead vs native Linux"
echo "  - GPU tests: $(if [ "${GPU_AVAILABLE}" = true ]; then echo "Available with CUDA on WSL2"; else echo "Not available"; fi)"
echo "  - Network: May need 'host.docker.internal' for Windows host access"
echo "  - File I/O: $(if [[ "${PROJECT_ROOT}" == /mnt/* ]]; then echo "SLOW (Windows FS)"; else echo "Fast (WSL FS)"; fi)"
echo ""

echo -e "${GREEN}Starting test execution...${NC}"
echo ""

# Change to project root
cd "${PROJECT_ROOT}"

# Phase 1: System Installation Tests (~18 minutes with multiplier)
echo -e "${GREEN}[1/5] Running System Installation Tests${NC}"
pytest tests/deployment/test_system_installation.py \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    | tee "${LOG_DIR}/test_system_installation_${TIMESTAMP}.log" || {
    echo -e "${RED}System installation tests failed${NC}"
    exit 1
}

# Phase 2: CPU Mode Tests (~24 minutes with multiplier)
echo ""
echo -e "${GREEN}[2/5] Running CPU Mode Tests${NC}"
pytest tests/deployment/test_cpu_gpu_mode.py::TestCPUMode \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    | tee "${LOG_DIR}/test_cpu_mode_${TIMESTAMP}.log" || {
    echo -e "${RED}CPU mode tests failed${NC}"
    exit 1
}

# Phase 3: GPU Mode Tests (~36 minutes with multiplier, conditional)
if [ "${GPU_AVAILABLE}" = true ]; then
    echo ""
    echo -e "${GREEN}[2b/5] Running GPU Mode Tests${NC}"
    pytest tests/deployment/test_cpu_gpu_mode.py::TestGPUMode \
        -v \
        --tb=short \
        --log-cli-level=INFO \
        | tee "${LOG_DIR}/test_gpu_mode_${TIMESTAMP}.log" || {
        echo -e "${RED}GPU mode tests failed${NC}"
        exit 1
    }
else
    echo ""
    echo -e "${YELLOW}[2b/5] Skipping GPU Mode Tests (no GPU available)${NC}"
fi

# Phase 4: Optimization Tests (~30 minutes with multiplier)
echo ""
echo -e "${GREEN}[3/5] Running Optimization Tests${NC}"
pytest tests/deployment/test_optimization.py \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    | tee "${LOG_DIR}/test_optimization_${TIMESTAMP}.log" || {
    echo -e "${RED}Optimization tests failed${NC}"
    exit 1
}

# Phase 5: Error Handling Tests (~24 minutes with multiplier)
echo ""
echo -e "${GREEN}[4/5] Running Error Handling Tests${NC}"
pytest tests/deployment/test_error_handling.py \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    | tee "${LOG_DIR}/test_error_handling_${TIMESTAMP}.log" || {
    echo -e "${RED}Error handling tests failed${NC}"
    exit 1
}

# Phase 6: Performance Tests (~36 minutes with multiplier)
echo ""
echo -e "${GREEN}[5/5] Running Performance Tests${NC}"
pytest tests/deployment/test_performance.py \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    | tee "${LOG_DIR}/test_performance_${TIMESTAMP}.log" || {
    echo -e "${RED}Performance tests failed${NC}"
    exit 1
}

# Summary
echo ""
echo -e "${GREEN}=== Test Execution Complete ===${NC}"
echo "Platform: WSL2 (${DOCKER_ENGINE})"
echo "Timestamp: ${TIMESTAMP}"
echo "Logs saved to: ${LOG_DIR}"
echo ""
echo "Test Summary:"
echo "  ✓ System Installation"
echo "  ✓ CPU Mode"
if [ "${GPU_AVAILABLE}" = true ]; then
    echo "  ✓ GPU Mode"
else
    echo "  - GPU Mode (skipped)"
fi
echo "  ✓ Optimization"
echo "  ✓ Error Handling"
echo "  ✓ Performance"
echo ""

# Post-test verification
echo -e "${YELLOW}Post-test verification...${NC}"
echo ""

# Check container status
echo "Container Status:"
docker compose -f docker-compose.cpu.yml ps

# Check service accessibility
echo ""
echo "Service Accessibility:"
curl -s -o /dev/null -w "Backend API (8000): %{http_code}\n" http://localhost:8000/ || echo "Backend API: FAIL"
curl -s -o /dev/null -w "Frontend (8080): %{http_code}\n" http://localhost:8080/ || echo "Frontend: FAIL"
curl -s -o /dev/null -w "RabbitMQ UI (15672): %{http_code}\n" http://localhost:15672/ || echo "RabbitMQ UI: FAIL"

echo ""
echo -e "${GREEN}WSL2 tests completed successfully!${NC}"
echo ""

if [[ "${PROJECT_ROOT}" == /mnt/* ]]; then
    echo -e "${YELLOW}Performance Tip: ${NC}"
    echo "Consider moving project to WSL filesystem for 10-50x faster I/O:"
    echo "  mv ${PROJECT_ROOT} ~/cellcraft"
fi
