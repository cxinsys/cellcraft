#!/usr/bin/env bash
# macOS Platform Test Execution Script for CellCraft Deployment Tests
#
# Optimized for macOS with Docker Desktop
# Expected runtime: ~2.5-3 hours (timeout_multiplier: 1.3-1.5)
#
# Prerequisites:
#   - Docker Desktop for Mac 4.x+
#   - Python 3.8+ with pytest, requests
#   - 50GB+ disk space
#   - Execute from project root
#
# Note: GPU tests are not supported on macOS (no NVIDIA GPU support)
#       Custom Plugin tests may be limited due to volume mount restrictions

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

# macOS specific: Detect architecture
ARCH=$(uname -m)
if [ "${ARCH}" = "arm64" ]; then
    PLATFORM="macOS Apple Silicon"
    TIMEOUT_MULTIPLIER=1.5
else
    PLATFORM="macOS Intel"
    TIMEOUT_MULTIPLIER=1.3
fi

# Create log directory
mkdir -p "${LOG_DIR}"

echo -e "${GREEN}=== CellCraft Deployment Tests (macOS) ===${NC}"
echo "Platform: ${PLATFORM}"
echo "Timestamp: ${TIMESTAMP}"
echo "Project Root: ${PROJECT_ROOT}"
echo "Log Directory: ${LOG_DIR}"
echo "Timeout Multiplier: ${TIMEOUT_MULTIPLIER}x"
echo ""

# Check prerequisites
echo -e "${YELLOW}Checking prerequisites...${NC}"

# Check Docker Desktop
if ! command -v docker &> /dev/null; then
    echo -e "${RED}Error: Docker not found${NC}"
    echo "Please install Docker Desktop for Mac from https://www.docker.com/products/docker-desktop"
    exit 1
fi
echo "✓ Docker $(docker --version)"

# Verify Docker Desktop is running
if ! docker info &> /dev/null; then
    echo -e "${RED}Error: Docker Desktop is not running${NC}"
    echo "Please start Docker Desktop and try again"
    exit 1
fi
echo "✓ Docker Desktop is running"

# Check Docker Compose
if ! docker compose version &> /dev/null; then
    echo -e "${RED}Error: Docker Compose v2 not found${NC}"
    exit 1
fi
echo "✓ Docker Compose $(docker compose version --short)"

# Check Python and pytest
if ! command -v pytest &> /dev/null; then
    echo -e "${RED}Error: pytest not found${NC}"
    echo "Install with: pip3 install pytest requests"
    exit 1
fi
echo "✓ pytest $(pytest --version | head -n1)"

# Check disk space (macOS uses different df format)
AVAILABLE_SPACE=$(df -g "${PROJECT_ROOT}" | awk 'NR==2 {print $4}')
if [ "${AVAILABLE_SPACE}" -lt 50 ]; then
    echo -e "${YELLOW}Warning: Only ${AVAILABLE_SPACE}GB available (recommended: 50GB+)${NC}"
fi
echo "✓ Disk space: ${AVAILABLE_SPACE}GB available"

# macOS specific warnings
echo ""
echo -e "${YELLOW}macOS Platform Notes:${NC}"
echo "  - GPU tests will be skipped (no NVIDIA GPU support)"
echo "  - Custom Plugin may be excluded (expected 6 CPU plugins instead of 7)"
echo "  - Docker Desktop uses VM: expect ~30-50% longer execution times"
echo "  - Volume mounts may be slower than native Linux"
echo ""

# Confirmation prompt
read -p "Continue with macOS test execution? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Test execution cancelled"
    exit 0
fi

echo -e "${GREEN}Starting test execution...${NC}"
echo ""

# Change to project root
cd "${PROJECT_ROOT}"

# Phase 1: System Installation Tests (~20 minutes with multiplier)
echo -e "${GREEN}[1/5] Running System Installation Tests${NC}"
pytest tests/deployment/test_system_installation.py \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    -m "cpu_mode" \
    | tee "${LOG_DIR}/test_system_installation_${TIMESTAMP}.log" || {
    echo -e "${RED}System installation tests failed${NC}"
    exit 1
}

# Phase 2: CPU Mode Tests (~25 minutes with multiplier)
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

# Phase 3: GPU Mode Tests (skipped on macOS)
echo ""
echo -e "${YELLOW}[2b/5] Skipping GPU Mode Tests (not supported on macOS)${NC}"

# Phase 4: Optimization Tests (~30 minutes with multiplier)
echo ""
echo -e "${GREEN}[3/5] Running Optimization Tests${NC}"
pytest tests/deployment/test_optimization.py \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    -m "cpu_mode" \
    | tee "${LOG_DIR}/test_optimization_${TIMESTAMP}.log" || {
    echo -e "${RED}Optimization tests failed${NC}"
    exit 1
}

# Phase 5: Error Handling Tests (~25 minutes with multiplier)
echo ""
echo -e "${GREEN}[4/5] Running Error Handling Tests${NC}"
pytest tests/deployment/test_error_handling.py \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    -m "cpu_mode" \
    | tee "${LOG_DIR}/test_error_handling_${TIMESTAMP}.log" || {
    echo -e "${RED}Error handling tests failed${NC}"
    exit 1
}

# Phase 6: Performance Tests (~40 minutes with multiplier)
echo ""
echo -e "${GREEN}[5/5] Running Performance Tests${NC}"
pytest tests/deployment/test_performance.py \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    -m "cpu_mode" \
    | tee "${LOG_DIR}/test_performance_${TIMESTAMP}.log" || {
    echo -e "${RED}Performance tests failed${NC}"
    exit 1
}

# Summary
echo ""
echo -e "${GREEN}=== Test Execution Complete ===${NC}"
echo "Platform: ${PLATFORM}"
echo "Timestamp: ${TIMESTAMP}"
echo "Logs saved to: ${LOG_DIR}"
echo ""
echo "Test Summary:"
echo "  ✓ System Installation"
echo "  ✓ CPU Mode"
echo "  - GPU Mode (not supported)"
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
echo -e "${GREEN}macOS tests completed successfully!${NC}"
echo ""
echo -e "${YELLOW}Note: ${NC}If you encountered volume mount issues, ensure:"
echo "  - Docker Desktop has access to project directory"
echo "  - File Sharing is configured in Docker Desktop preferences"
