#!/bin/bash

# Script to verify E2E split test setup
# Checks: Jenkins jobs, Push URLs, test files, container readiness

set -e

JENKINS_URL="http://localhost:8089"
JENKINS_USER="admin"
JENKINS_TOKEN="${JENKINS_API_TOKEN:-}"

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

echo "========================================="
echo "E2E Split Test Setup Verification"
echo "========================================="

# Counters
total_checks=0
passed_checks=0
warnings=0

# Helper functions
check_pass() {
    echo -e "${GREEN}✓${NC} $1"
    ((passed_checks++))
    ((total_checks++))
}

check_fail() {
    echo -e "${RED}✗${NC} $1"
    ((total_checks++))
}

check_warn() {
    echo -e "${YELLOW}⚠${NC} $1"
    ((warnings++))
    ((total_checks++))
}

info() {
    echo -e "${BLUE}ℹ${NC} $1"
}

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo -e "\n${BLUE}=== Checking Jenkins Connection ===${NC}"

if [ -z "$JENKINS_TOKEN" ]; then
    check_warn "JENKINS_API_TOKEN not set (some checks will be skipped)"
    info "Set with: export JENKINS_API_TOKEN=your_token"
else
    if curl -s -f -u "${JENKINS_USER}:${JENKINS_TOKEN}" \
        "${JENKINS_URL}/api/json" > /dev/null 2>&1; then
        check_pass "Jenkins is accessible"
    else
        check_fail "Cannot connect to Jenkins at ${JENKINS_URL}"
    fi
fi

echo -e "\n${BLUE}=== Checking Jenkins Jobs ===${NC}"

jobs=(
    "e2e-01-file-assignment"
    "e2e-02-data-display"
    "e2e-03-scatter-plot"
    "e2e-04-algorithm-config"
    "e2e-05-workflow-execution"
    "e2e-06-result-visualization"
)

if [ -n "$JENKINS_TOKEN" ]; then
    for job in "${jobs[@]}"; do
        if curl -s -f -u "${JENKINS_USER}:${JENKINS_TOKEN}" \
            "${JENKINS_URL}/job/${job}/api/json" > /dev/null 2>&1; then
            check_pass "Job exists: ${job}"
        else
            check_fail "Job not found: ${job}"
        fi
    done
else
    info "Skipping Jenkins job checks (no API token)"
fi

echo -e "\n${BLUE}=== Checking Jenkins Job XML Files ===${NC}"

JENKINS_JOBS_DIR="${PROJECT_ROOT}/jenkins-jobs"

for job in "${jobs[@]}"; do
    xml_file="${JENKINS_JOBS_DIR}/${job}.xml"
    if [ -f "$xml_file" ]; then
        check_pass "XML file exists: ${job}.xml"
    else
        check_fail "XML file missing: ${job}.xml"
    fi
done

echo -e "\n${BLUE}=== Checking Test Files ===${NC}"

TEST_FILES_DIR="${PROJECT_ROOT}/frontend/tests/e2e/workflows"

test_files=(
    "01-file-assignment.spec.js"
    "02-data-display.spec.js"
    "03-scatter-plot.spec.js"
    "04-algorithm-config.spec.js"
    "05-workflow-execution.spec.js"
    "06-result-visualization.spec.js"
)

for test_file in "${test_files[@]}"; do
    file_path="${TEST_FILES_DIR}/${test_file}"
    if [ -f "$file_path" ]; then
        check_pass "Test file exists: ${test_file}"
    else
        check_fail "Test file missing: ${test_file}"
    fi
done

echo -e "\n${BLUE}=== Checking Push URLs Configuration ===${NC}"

PUSH_URLS_FILE="${PROJECT_ROOT}/push-urls.txt"

if [ -f "$PUSH_URLS_FILE" ]; then
    check_pass "push-urls.txt file exists"

    # Check if split E2E URLs are configured
    for i in {1..6}; do
        job_key="e2e-0${i}"
        if grep -q "^${job_key}" "$PUSH_URLS_FILE"; then
            url=$(grep "^${job_key}" "$PUSH_URLS_FILE" | cut -d'=' -f2)
            if [[ "$url" == *"YOUR_KEY"* ]]; then
                check_warn "${job_key}: URL not configured (placeholder found)"
            else
                check_pass "${job_key}: URL configured"
            fi
        else
            check_fail "${job_key}: Not found in push-urls.txt"
        fi
    done
else
    check_fail "push-urls.txt file not found"
fi

echo -e "\n${BLUE}=== Checking Docker Containers ===${NC}"

# Check Jenkins container
if docker ps --format '{{.Names}}' | grep -q "cellcraft_jenkins"; then
    check_pass "Jenkins container is running"

    # Check Node.js in Jenkins
    if docker exec cellcraft_jenkins which node > /dev/null 2>&1; then
        node_version=$(docker exec cellcraft_jenkins node --version 2>/dev/null || echo "unknown")
        check_pass "Node.js installed in Jenkins: ${node_version}"
    else
        check_warn "Node.js not found in Jenkins container"
        info "Install with: docker exec -u root cellcraft_jenkins apt-get update && apt-get install -y nodejs npm"
    fi

    # Check npm in Jenkins
    if docker exec cellcraft_jenkins which npm > /dev/null 2>&1; then
        npm_version=$(docker exec cellcraft_jenkins npm --version 2>/dev/null || echo "unknown")
        check_pass "npm installed in Jenkins: ${npm_version}"
    else
        check_warn "npm not found in Jenkins container"
    fi

    # Check Docker CLI in Jenkins
    if docker exec cellcraft_jenkins which docker > /dev/null 2>&1; then
        check_pass "Docker CLI installed in Jenkins"
    else
        check_warn "Docker CLI not found in Jenkins container"
        info "Install with: docker exec -u root cellcraft_jenkins apt-get install -y docker.io"
    fi
else
    check_fail "Jenkins container not running"
fi

# Check Frontend container
if docker ps --format '{{.Names}}' | grep -q "cellcraft-frontend-1"; then
    check_pass "Frontend container is running"
else
    check_fail "Frontend container not running"
fi

# Check Backend container
if docker ps --format '{{.Names}}' | grep -q "cellcraft-backend-1"; then
    check_pass "Backend container is running"
else
    check_fail "Backend container not running"
fi

# Check Uptime Kuma container
if docker ps --format '{{.Names}}' | grep -q "uptime-kuma"; then
    check_pass "Uptime Kuma container is running"
else
    check_fail "Uptime Kuma container not running"
fi

echo -e "\n${BLUE}=== Checking Volume Mounts ===${NC}"

# Check if Jenkins can access frontend tests
if docker exec cellcraft_jenkins test -d /workspace/frontend/tests/e2e/workflows 2>/dev/null; then
    check_pass "Jenkins can access frontend test files"
else
    check_fail "Jenkins cannot access /workspace/frontend/tests/e2e/workflows"
    info "Check docker-compose volume mount for Jenkins"
fi

# Check if Jenkins can write to /tmp
if docker exec cellcraft_jenkins test -w /tmp 2>/dev/null; then
    check_pass "Jenkins can write to /tmp"
else
    check_fail "Jenkins cannot write to /tmp"
fi

echo -e "\n${BLUE}=== Checking Network Connectivity ===${NC}"

# Check if Jenkins can reach frontend
if docker exec cellcraft_jenkins curl -s -f http://cellcraft-frontend-1:8080 > /dev/null 2>&1; then
    check_pass "Jenkins can reach frontend container"
else
    check_warn "Jenkins cannot reach frontend at http://cellcraft-frontend-1:8080"
fi

# Check if Jenkins can reach Uptime Kuma
if docker exec cellcraft_jenkins curl -s -f http://uptime-kuma:3001 > /dev/null 2>&1; then
    check_pass "Jenkins can reach Uptime Kuma"
else
    check_warn "Jenkins cannot reach Uptime Kuma at http://uptime-kuma:3001"
fi

echo -e "\n${BLUE}=== Checking Helper Scripts ===${NC}"

scripts=(
    "setup-e2e-split-jobs.sh"
    "update-e2e-push-urls.sh"
    "verify-e2e-split-setup.sh"
)

for script in "${scripts[@]}"; do
    script_path="${SCRIPT_DIR}/${script}"
    if [ -f "$script_path" ]; then
        if [ -x "$script_path" ]; then
            check_pass "Script exists and is executable: ${script}"
        else
            check_warn "Script exists but not executable: ${script}"
            info "Fix with: chmod +x ${script_path}"
        fi
    else
        check_fail "Script missing: ${script}"
    fi
done

echo -e "\n========================================="
echo -e "${GREEN}Verification Complete${NC}"
echo "========================================="
echo -e "Passed: ${GREEN}${passed_checks}${NC} / ${total_checks}"
echo -e "Warnings: ${YELLOW}${warnings}${NC}"
echo -e "Failed: ${RED}$((total_checks - passed_checks - warnings))${NC}"

if [ $passed_checks -eq $total_checks ]; then
    echo -e "\n${GREEN}✓ All checks passed! Setup is complete.${NC}"
    exit 0
elif [ $((passed_checks + warnings)) -eq $total_checks ]; then
    echo -e "\n${YELLOW}⚠ Setup complete with warnings. Review the warnings above.${NC}"
    exit 0
else
    echo -e "\n${RED}✗ Setup incomplete. Please fix the failed checks above.${NC}"
    exit 1
fi
