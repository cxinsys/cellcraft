#!/bin/bash

#######################################
# CellCraft Jenkins Setup Verification Script
#
# This script verifies that Jenkins is properly configured and
# all required jobs are created and functional.
#
# Usage:
#   ./verify-jenkins-setup.sh [jenkins-url]
#
# Example:
#   ./verify-jenkins-setup.sh http://localhost:8080
#
# Prerequisites:
#   - Jenkins must be running
#   - Docker must be accessible
#
#######################################

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

JENKINS_URL=${1:-http://localhost:8080}

echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN}CellCraft Jenkins Setup Verification${NC}"
echo -e "${GREEN}=========================================${NC}"
echo ""

# Track verification status
pass_count=0
fail_count=0
warn_count=0

# Function to print test result
print_result() {
    local status=$1
    local message=$2

    case $status in
        "PASS")
            echo -e "${GREEN}✓ PASS${NC}: $message"
            ((pass_count++))
            ;;
        "FAIL")
            echo -e "${RED}✗ FAIL${NC}: $message"
            ((fail_count++))
            ;;
        "WARN")
            echo -e "${YELLOW}⚠ WARN${NC}: $message"
            ((warn_count++))
            ;;
    esac
}

# 1. Docker Container Verification
echo -e "${BLUE}=== 1. Docker Container Status ===${NC}"
echo ""

# Check Jenkins container
if docker ps --filter "name=cellcraft_jenkins" --format "{{.Names}}" | grep -q "cellcraft_jenkins"; then
    print_result "PASS" "Jenkins container is running"
else
    print_result "FAIL" "Jenkins container is not running"
fi

# Check Uptime Kuma container
if docker ps --filter "name=cellcraft_uptime_kuma" --format "{{.Names}}" | grep -q "cellcraft_uptime_kuma"; then
    print_result "PASS" "Uptime Kuma container is running"
else
    print_result "FAIL" "Uptime Kuma container is not running"
fi

# Check other required containers
for container in cellcraft_backend cellcraft_frontend cellcraft_db cellcraft_rabbitmq; do
    if docker ps --filter "name=$container" --format "{{.Names}}" | grep -q "$container"; then
        print_result "PASS" "$container is running"
    else
        print_result "FAIL" "$container is not running"
    fi
done

echo ""

# 2. Jenkins Accessibility
echo -e "${BLUE}=== 2. Jenkins Accessibility ===${NC}"
echo ""

if curl -s -o /dev/null -w "%{http_code}" "$JENKINS_URL/login" | grep -q "200\|403"; then
    print_result "PASS" "Jenkins is accessible at $JENKINS_URL"
else
    print_result "FAIL" "Jenkins is not accessible at $JENKINS_URL"
fi

echo ""

# 3. Jenkins Job Configuration Files
echo -e "${BLUE}=== 3. Jenkins Job Configuration Files ===${NC}"
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
JOBS_DIR="$PROJECT_ROOT/jenkins-jobs"

required_jobs=(
    "backend-unit-test.xml"
    "frontend-unit-test.xml"
    "e2e-test.xml"
    "docker-health-check.xml"
)

for job in "${required_jobs[@]}"; do
    if [ -f "$JOBS_DIR/$job" ]; then
        print_result "PASS" "Job configuration exists: $job"
    else
        print_result "FAIL" "Job configuration missing: $job"
    fi
done

echo ""

# 4. Volume Verification
echo -e "${BLUE}=== 4. Docker Volume Verification ===${NC}"
echo ""

if docker volume ls | grep -q "jenkins_data"; then
    print_result "PASS" "jenkins_data volume exists"
else
    print_result "WARN" "jenkins_data volume not found (will be created on first run)"
fi

if docker volume ls | grep -q "uptime_kuma_data"; then
    print_result "PASS" "uptime_kuma_data volume exists"
else
    print_result "WARN" "uptime_kuma_data volume not found (will be created on first run)"
fi

echo ""

# 5. Docker Socket Access
echo -e "${BLUE}=== 5. Docker Socket Access ===${NC}"
echo ""

if docker exec cellcraft_jenkins test -S /var/run/docker.sock 2>/dev/null; then
    print_result "PASS" "Jenkins has access to Docker socket"
else
    print_result "WARN" "Jenkins container not running or Docker socket not mounted"
fi

echo ""

# 6. Workspace Volume Mounts
echo -e "${BLUE}=== 6. Workspace Volume Mounts ===${NC}"
echo ""

if docker exec cellcraft_jenkins test -d /workspace/backend 2>/dev/null; then
    print_result "PASS" "Backend workspace is mounted"
else
    print_result "WARN" "Backend workspace not mounted"
fi

if docker exec cellcraft_jenkins test -d /workspace/frontend 2>/dev/null; then
    print_result "PASS" "Frontend workspace is mounted"
else
    print_result "WARN" "Frontend workspace not mounted"
fi

echo ""

# 7. Network Connectivity
echo -e "${BLUE}=== 7. Network Connectivity ===${NC}"
echo ""

# Check if Jenkins can reach Uptime Kuma
if docker exec cellcraft_jenkins ping -c 1 uptime-kuma >/dev/null 2>&1; then
    print_result "PASS" "Jenkins can reach Uptime Kuma container"
else
    print_result "WARN" "Jenkins cannot ping Uptime Kuma (may not be critical)"
fi

# Check if Jenkins can reach backend
if docker exec cellcraft_jenkins ping -c 1 backend >/dev/null 2>&1; then
    print_result "PASS" "Jenkins can reach backend container"
else
    print_result "WARN" "Jenkins cannot ping backend container"
fi

echo ""

# Summary
echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN}Verification Summary${NC}"
echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN}Passed: $pass_count${NC}"
if [ $warn_count -gt 0 ]; then
    echo -e "${YELLOW}Warnings: $warn_count${NC}"
else
    echo "Warnings: 0"
fi
if [ $fail_count -gt 0 ]; then
    echo -e "${RED}Failed: $fail_count${NC}"
else
    echo "Failed: 0"
fi
echo ""

# Exit code
if [ $fail_count -gt 0 ]; then
    echo -e "${RED}✗ Verification failed with $fail_count error(s)${NC}"
    echo ""
    echo "Next steps:"
    echo "1. Start missing containers: docker-compose -f docker-compose.prod.yml up -d"
    echo "2. Check container logs: docker logs cellcraft_jenkins"
    echo "3. Verify docker-compose.prod.yml configuration"
    exit 1
else
    if [ $warn_count -gt 0 ]; then
        echo -e "${YELLOW}⚠ Verification completed with $warn_count warning(s)${NC}"
        echo ""
        echo "Review warnings above. Setup may be incomplete."
        exit 0
    else
        echo -e "${GREEN}✓ All verifications passed!${NC}"
        echo ""
        echo "Next steps:"
        echo "1. Access Jenkins at $JENKINS_URL"
        echo "2. Get initial password: docker exec cellcraft_jenkins cat /var/jenkins_home/secrets/initialAdminPassword"
        echo "3. Complete Jenkins setup wizard"
        echo "4. Run ./scripts/setup-jenkins-jobs.sh to create jobs"
        exit 0
    fi
fi
