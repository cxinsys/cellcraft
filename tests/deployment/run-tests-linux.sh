#!/usr/bin/env bash
# Linux Platform Test Execution Script for CellCraft Deployment Tests
#
# Optimized for native Linux with Docker Engine
# Expected runtime: ~2-2.5 hours (with session-scoped fixtures)
#
# Prerequisites:
#   - Docker Engine 24.x+
#   - Docker Compose v2
#   - NVIDIA drivers + Container Toolkit (for GPU tests)
#   - Python 3.8+ with pytest, requests
#   - 50GB+ disk space
#   - Execute from project root

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

# Test markers
MARKERS_CPU="cpu_mode and not slow"
MARKERS_CPU_SLOW="cpu_mode and slow"
MARKERS_GPU="gpu_mode"
MARKERS_PERFORMANCE="slow"

# Create log directory
mkdir -p "${LOG_DIR}"

echo -e "${GREEN}=== CellCraft Deployment Tests (Linux) ===${NC}"
echo "Timestamp: ${TIMESTAMP}"
echo "Project Root: ${PROJECT_ROOT}"
echo "Log Directory: ${LOG_DIR}"
echo ""

# Check prerequisites
echo -e "${YELLOW}Checking prerequisites...${NC}"

# Check Docker
if ! command -v docker &> /dev/null; then
    echo -e "${RED}Error: Docker not found${NC}"
    exit 1
fi
echo "✓ Docker $(docker --version)"

# Check Docker Compose
if ! docker compose version &> /dev/null; then
    echo -e "${RED}Error: Docker Compose v2 not found${NC}"
    exit 1
fi
echo "✓ Docker Compose $(docker compose version --short)"

# Check GPU availability
GPU_AVAILABLE=false
if command -v nvidia-smi &> /dev/null && nvidia-smi &> /dev/null; then
    GPU_AVAILABLE=true
    echo "✓ NVIDIA GPU detected: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -n1)"
else
    echo "⚠ No NVIDIA GPU detected - GPU tests will be skipped"
fi

# Check Python and pytest
if ! command -v pytest &> /dev/null; then
    echo -e "${RED}Error: pytest not found${NC}"
    exit 1
fi
echo "✓ pytest $(pytest --version | head -n1)"

# Check disk space
AVAILABLE_SPACE=$(df -BG "${PROJECT_ROOT}" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "${AVAILABLE_SPACE}" -lt 50 ]; then
    echo -e "${YELLOW}Warning: Only ${AVAILABLE_SPACE}GB available (recommended: 50GB+)${NC}"
fi
echo "✓ Disk space: ${AVAILABLE_SPACE}GB available"

echo ""
echo -e "${GREEN}Starting test execution...${NC}"
echo ""

# Change to project root
cd "${PROJECT_ROOT}"

# Phase 1: System Installation Tests (~15 minutes)
echo -e "${GREEN}[1/5] Running System Installation Tests${NC}"
pytest tests/deployment/test_system_installation.py \
    -v \
    --tb=short \
    --log-cli-level=INFO \
    | tee "${LOG_DIR}/test_system_installation_${TIMESTAMP}.log" || {
    echo -e "${RED}System installation tests failed${NC}"
    exit 1
}

# Phase 2: CPU Mode Tests (~20 minutes)
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

# Phase 3: GPU Mode Tests (~30 minutes, conditional)
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

# Phase 4: Optimization Tests (~25 minutes)
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

# Phase 5: Error Handling Tests (~20 minutes)
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

# Phase 6: Performance Tests (~30 minutes)
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
echo -e "${GREEN}All tests completed successfully!${NC}"
