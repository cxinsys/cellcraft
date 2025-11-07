#!/bin/bash

#######################################
# Uptime Kuma Setup Verification Script
#
# This script verifies that Uptime Kuma is properly configured
# with all required monitors and status page.
#
# Usage:
#   ./verify-uptime-kuma-setup.sh [uptime-kuma-url]
#
# Example:
#   ./verify-uptime-kuma-setup.sh http://localhost:3001
#
# Prerequisites:
#   - Uptime Kuma must be running
#   - Docker must be accessible
#
#######################################

set -e

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

UPTIME_KUMA_URL=${1:-http://localhost:3001}

echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN}Uptime Kuma Setup Verification${NC}"
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

if docker ps --filter "name=cellcraft_uptime_kuma" --format "{{.Names}}" | grep -q "cellcraft_uptime_kuma"; then
    print_result "PASS" "Uptime Kuma container is running"
else
    print_result "FAIL" "Uptime Kuma container is not running"
fi

echo ""

# 2. Uptime Kuma Accessibility
echo -e "${BLUE}=== 2. Uptime Kuma Web UI Accessibility ===${NC}"
echo ""

http_code=$(curl -s -o /dev/null -w "%{http_code}" "$UPTIME_KUMA_URL" || echo "000")

if [ "$http_code" = "200" ] || [ "$http_code" = "302" ]; then
    print_result "PASS" "Uptime Kuma web UI is accessible at $UPTIME_KUMA_URL"
else
    print_result "FAIL" "Uptime Kuma web UI not accessible (HTTP $http_code)"
fi

echo ""

# 3. Required Monitors (Manual Check)
echo -e "${BLUE}=== 3. Monitor Configuration Check ===${NC}"
echo ""

expected_monitors=(
    "CellCraft Main"
    "Backend API Health"
    "Backend Unit Tests"
    "Frontend Unit Tests"
    "E2E Integration Tests"
    "Docker Environment"
)

echo -e "${YELLOW}Expected monitors (6 total):${NC}"
for monitor in "${expected_monitors[@]}"; do
    echo "  - $monitor"
done

print_result "WARN" "Manual verification required: Check Uptime Kuma dashboard for 6 monitors"

echo ""

# 4. Status Page Check
echo -e "${BLUE}=== 4. Status Page Configuration ===${NC}"
echo ""

status_page_url="$UPTIME_KUMA_URL/status/status"
status_http_code=$(curl -s -o /dev/null -w "%{http_code}" "$status_page_url" || echo "000")

if [ "$status_http_code" = "200" ]; then
    print_result "PASS" "Status Page is accessible at $status_page_url"
elif [ "$status_http_code" = "404" ]; then
    print_result "WARN" "Status Page not found (needs manual creation)"
else
    print_result "FAIL" "Status Page check failed (HTTP $status_http_code)"
fi

echo ""

# 5. Volume Verification
echo -e "${BLUE}=== 5. Docker Volume Verification ===${NC}"
echo ""

if docker volume ls | grep -q "uptime_kuma_data"; then
    print_result "PASS" "uptime_kuma_data volume exists"

    # Check volume usage
    volume_info=$(docker volume inspect uptime_kuma_data 2>/dev/null || echo "")
    if [ -n "$volume_info" ]; then
        print_result "PASS" "Volume is properly configured"
    fi
else
    print_result "FAIL" "uptime_kuma_data volume not found"
fi

echo ""

# 6. Network Connectivity
echo -e "${BLUE}=== 6. Network Connectivity ===${NC}"
echo ""

# Check if Uptime Kuma can reach other services
if docker exec cellcraft_uptime_kuma ping -c 1 backend >/dev/null 2>&1; then
    print_result "PASS" "Uptime Kuma can reach backend container"
else
    print_result "WARN" "Cannot ping backend (may not be critical)"
fi

if docker exec cellcraft_uptime_kuma ping -c 1 frontend >/dev/null 2>&1; then
    print_result "PASS" "Uptime Kuma can reach frontend container"
else
    print_result "WARN" "Cannot ping frontend (may not be critical)"
fi

echo ""

# 7. Configuration Files Check
echo -e "${BLUE}=== 7. Configuration Files ===${NC}"
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$PROJECT_ROOT/uptime-kuma-config"

required_files=(
    "$CONFIG_DIR/http-monitors-template.json"
    "$CONFIG_DIR/push-monitors-template.json"
    "$CONFIG_DIR/status-page-config.json"
)

for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        print_result "PASS" "Configuration file exists: $(basename $file)"
    else
        print_result "FAIL" "Configuration file missing: $(basename $file)"
    fi
done

echo ""

# 8. Push URLs File Check
echo -e "${BLUE}=== 8. Push URLs Configuration ===${NC}"
echo ""

PUSH_URLS_FILE="$PROJECT_ROOT/push-urls.txt"

if [ -f "$PUSH_URLS_FILE" ]; then
    print_result "PASS" "Push URLs file exists"

    # Count URLs in file
    url_count=$(grep -c "^[^#].*=http" "$PUSH_URLS_FILE" 2>/dev/null || echo "0")
    if [ "$url_count" -ge 4 ]; then
        print_result "PASS" "Push URLs file contains $url_count URLs (expected: 4)"
    else
        print_result "WARN" "Push URLs file contains only $url_count URLs (expected: 4)"
    fi
else
    print_result "WARN" "Push URLs file not found (run setup-uptime-kuma-monitors.sh)"
fi

echo ""

# 9. Jenkins Integration Check
echo -e "${BLUE}=== 9. Jenkins Integration ===${NC}"
echo ""

if docker ps --filter "name=cellcraft_jenkins" --format "{{.Names}}" | grep -q "cellcraft_jenkins"; then
    print_result "PASS" "Jenkins container is running"

    # Check if Jenkins can reach Uptime Kuma
    if docker exec cellcraft_jenkins ping -c 1 uptime-kuma >/dev/null 2>&1; then
        print_result "PASS" "Jenkins can reach Uptime Kuma container"
    else
        print_result "WARN" "Jenkins cannot ping Uptime Kuma"
    fi
else
    print_result "WARN" "Jenkins container not running (integration not testable)"
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
    echo "1. Start Uptime Kuma container: docker-compose -f docker-compose.prod.yml up -d uptime-kuma"
    echo "2. Access Uptime Kuma UI: ssh -L 3001:localhost:3001 user@server"
    echo "3. Create monitors manually or run: ./scripts/setup-uptime-kuma-monitors.sh"
    echo "4. Create Status Page in Uptime Kuma UI"
    exit 1
else
    if [ $warn_count -gt 0 ]; then
        echo -e "${YELLOW}⚠ Verification completed with $warn_count warning(s)${NC}"
        echo ""
        echo "Manual verification needed:"
        echo "1. Login to Uptime Kuma at $UPTIME_KUMA_URL"
        echo "2. Verify all 6 monitors are created and working"
        echo "3. Check Status Page at $status_page_url"
        echo "4. Test Jenkins integration by running a job"
        echo ""
        echo "Optional automation:"
        echo "  ./scripts/setup-uptime-kuma-monitors.sh $UPTIME_KUMA_URL YOUR_API_TOKEN"
        exit 0
    else
        echo -e "${GREEN}✓ All verifications passed!${NC}"
        echo ""
        echo "Setup complete! Next steps:"
        echo "1. Access Uptime Kuma: $UPTIME_KUMA_URL"
        echo "2. Access Status Page: $status_page_url"
        echo "3. Run Jenkins jobs to test Push Monitor integration"
        echo "4. Verify Status Page shows real-time updates"
        exit 0
    fi
fi
