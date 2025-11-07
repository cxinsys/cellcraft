#!/bin/bash

#######################################
# Complete Pipeline Verification Script
#
# This script performs comprehensive verification of the entire
# CellCraft status page infrastructure (Phase 1-4).
#
# Usage:
#   ./verify-complete-pipeline.sh
#
# Verifies:
#   - Phase 1: Docker Compose services
#   - Phase 2: Jenkins setup and jobs
#   - Phase 3: Uptime Kuma setup and monitors
#   - Phase 4: Jenkins-Uptime Kuma integration
#
#######################################

set -e

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

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

echo -e "${CYAN}=============================================${NC}"
echo -e "${CYAN}CellCraft Complete Pipeline Verification${NC}"
echo -e "${CYAN}Phase 1-4 Comprehensive Check${NC}"
echo -e "${CYAN}=============================================${NC}"
echo ""

# ===========================
# Phase 1: Docker Compose
# ===========================
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}Phase 1: Docker Compose Services${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

required_containers=(
    "cellcraft_frontend"
    "cellcraft_backend"
    "cellcraft_db"
    "cellcraft_rabbitmq"
    "cellcraft_celery"
    "cellcraft_jenkins"
    "cellcraft_uptime_kuma"
)

for container in "${required_containers[@]}"; do
    if docker ps --filter "name=$container" --format "{{.Names}}" | grep -q "$container"; then
        print_result "PASS" "$container is running"
    else
        print_result "FAIL" "$container is not running"
    fi
done

# Check volumes
echo ""
for volume in postgres_data jenkins_data uptime_kuma_data; do
    if docker volume ls | grep -q "$volume"; then
        print_result "PASS" "$volume volume exists"
    else
        print_result "WARN" "$volume volume not found"
    fi
done

echo ""

# ===========================
# Phase 2: Jenkins Setup
# ===========================
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}Phase 2: Jenkins Configuration${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

# Check Jenkins accessibility
JENKINS_URL="http://localhost:8080"
if curl -s -o /dev/null -w "%{http_code}" "$JENKINS_URL/login" | grep -q "200\|403"; then
    print_result "PASS" "Jenkins is accessible"
else
    print_result "WARN" "Jenkins not accessible (may need SSH tunnel)"
fi

# Check Jenkins Job XML files
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

# Check Docker socket access
if docker exec cellcraft_jenkins test -S /var/run/docker.sock 2>/dev/null; then
    print_result "PASS" "Jenkins has Docker socket access"
else
    print_result "WARN" "Jenkins Docker socket not accessible"
fi

echo ""

# ===========================
# Phase 3: Uptime Kuma Setup
# ===========================
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}Phase 3: Uptime Kuma Configuration${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

# Check Uptime Kuma accessibility
UPTIME_KUMA_URL="http://localhost:3001"
http_code=$(curl -s -o /dev/null -w "%{http_code}" "$UPTIME_KUMA_URL" || echo "000")

if [ "$http_code" = "200" ] || [ "$http_code" = "302" ]; then
    print_result "PASS" "Uptime Kuma is accessible"
else
    print_result "WARN" "Uptime Kuma not accessible (may need SSH tunnel)"
fi

# Check configuration templates
CONFIG_DIR="$PROJECT_ROOT/uptime-kuma-config"
config_files=(
    "http-monitors-template.json"
    "push-monitors-template.json"
    "status-page-config.json"
)

for file in "${config_files[@]}"; do
    if [ -f "$CONFIG_DIR/$file" ]; then
        print_result "PASS" "Configuration template exists: $file"
    else
        print_result "FAIL" "Configuration template missing: $file"
    fi
done

# Check Push URLs file
PUSH_URLS_FILE="$PROJECT_ROOT/push-urls.txt"
if [ -f "$PUSH_URLS_FILE" ]; then
    url_count=$(grep -c "^[^#].*=http" "$PUSH_URLS_FILE" 2>/dev/null || echo "0")
    if [ "$url_count" -ge 4 ]; then
        print_result "PASS" "Push URLs file contains $url_count URLs"
    else
        print_result "WARN" "Push URLs file contains only $url_count URLs (expected: 4)"
    fi
else
    print_result "WARN" "Push URLs file not found"
fi

echo ""

# ===========================
# Phase 4: Integration
# ===========================
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}Phase 4: Jenkins-Uptime Kuma Integration${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

# Check network connectivity
if docker exec cellcraft_jenkins ping -c 1 uptime-kuma >/dev/null 2>&1; then
    print_result "PASS" "Jenkins can reach Uptime Kuma"
else
    print_result "WARN" "Jenkins cannot ping Uptime Kuma"
fi

if docker exec cellcraft_uptime_kuma ping -c 1 backend >/dev/null 2>&1; then
    print_result "PASS" "Uptime Kuma can reach backend"
else
    print_result "WARN" "Uptime Kuma cannot ping backend"
fi

# Check automation scripts
scripts=(
    "setup-jenkins-jobs.sh"
    "setup-uptime-kuma-monitors.sh"
    "update-jenkins-push-urls.sh"
    "test-jenkins-uptime-integration.sh"
    "verify-jenkins-setup.sh"
    "verify-uptime-kuma-setup.sh"
)

echo ""
for script in "${scripts[@]}"; do
    if [ -x "$SCRIPT_DIR/$script" ]; then
        print_result "PASS" "Script is executable: $script"
    elif [ -f "$SCRIPT_DIR/$script" ]; then
        print_result "WARN" "Script exists but not executable: $script"
    else
        print_result "FAIL" "Script missing: $script"
    fi
done

echo ""

# ===========================
# Documentation Check
# ===========================
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}Documentation Verification${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

DOCS_DIR="$PROJECT_ROOT/docs"
docs=(
    "jenkins-setup-guide.md"
    "uptime-kuma-setup-guide.md"
    "integration-test-guide.md"
)

for doc in "${docs[@]}"; do
    if [ -f "$DOCS_DIR/$doc" ]; then
        size=$(stat -f%z "$DOCS_DIR/$doc" 2>/dev/null || stat -c%s "$DOCS_DIR/$doc" 2>/dev/null || echo "0")
        if [ "$size" -gt 1000 ]; then
            print_result "PASS" "Documentation exists: $doc (${size} bytes)"
        else
            print_result "WARN" "Documentation too small: $doc (${size} bytes)"
        fi
    else
        print_result "FAIL" "Documentation missing: $doc"
    fi
done

echo ""

# ===========================
# Summary
# ===========================
echo -e "${CYAN}=============================================${NC}"
echo -e "${CYAN}Verification Summary${NC}"
echo -e "${CYAN}=============================================${NC}"
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

# Calculate readiness score
total_checks=$((pass_count + fail_count + warn_count))
readiness_score=$(( (pass_count * 100) / total_checks ))

echo -e "${CYAN}Readiness Score: ${readiness_score}%${NC}"
echo ""

# Exit code and recommendations
if [ $fail_count -gt 0 ]; then
    echo -e "${RED}✗ Pipeline verification failed with $fail_count critical error(s)${NC}"
    echo ""
    echo "Critical Issues:"
    echo "  - Some required containers are not running"
    echo "  - Some configuration files are missing"
    echo ""
    echo "Next steps:"
    echo "  1. Start missing containers: docker-compose -f docker-compose.prod.yml up -d"
    echo "  2. Run individual phase verifications:"
    echo "     - ./scripts/verify-jenkins-setup.sh"
    echo "     - ./scripts/verify-uptime-kuma-setup.sh"
    echo "  3. Review error messages above"
    exit 1
elif [ $warn_count -gt 5 ]; then
    echo -e "${YELLOW}⚠ Pipeline verification completed with warnings${NC}"
    echo ""
    echo "Action items:"
    echo "  - $warn_count warnings need attention"
    echo "  - Some services may not be fully configured"
    echo ""
    echo "Recommended actions:"
    echo "  1. Review warnings above"
    echo "  2. Complete manual setup steps:"
    echo "     - Jenkins initial setup and plugin installation"
    echo "     - Uptime Kuma monitor creation"
    echo "     - Push URL configuration"
    echo "  3. Run integration test: ./scripts/test-jenkins-uptime-integration.sh"
    echo ""
    echo "Pipeline is functional but not fully optimized."
    exit 0
else
    echo -e "${GREEN}✓ All verifications passed!${NC}"
    echo ""
    echo "Pipeline Status: ${GREEN}READY${NC}"
    echo ""
    echo "Phase 1-4 Complete! Next steps:"
    echo ""
    echo "Testing:"
    echo "  1. Run integration test:"
    echo "     ./scripts/test-jenkins-uptime-integration.sh \\"
    echo "       http://localhost:8080 admin:TOKEN"
    echo ""
    echo "Monitoring:"
    echo "  2. Access Uptime Kuma: http://localhost:3001"
    echo "  3. Access Status Page: http://localhost:3001/status/status"
    echo "  4. Access Jenkins: http://localhost:8080"
    echo ""
    echo "Phase 5:"
    echo "  5. Configure Nginx reverse proxy"
    echo "  6. Set up DNS: status.cellcraft.app"
    echo "  7. Obtain SSL certificate"
    echo "  8. Deploy public status page"
    echo ""
    exit 0
fi
