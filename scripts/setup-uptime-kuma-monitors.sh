#!/bin/bash

#######################################
# Uptime Kuma Monitor Setup Script
#
# This script automatically creates HTTP(s) and Push monitors
# in Uptime Kuma using the REST API.
#
# Usage:
#   ./setup-uptime-kuma-monitors.sh <uptime-kuma-url> <api-token>
#
# Example:
#   ./setup-uptime-kuma-monitors.sh http://localhost:3001 your-api-token
#
# Prerequisites:
#   - Uptime Kuma must be running
#   - API token must be generated in Settings > API Keys
#   - jq must be installed for JSON parsing
#
# Notes:
#   - This script creates monitors but NOT the Status Page
#   - Status Page must be created manually through the web UI
#   - Push Monitor URLs will be saved to push-urls.txt
#
#######################################

set -e

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check arguments
if [ $# -ne 2 ]; then
    echo -e "${RED}Error: Invalid number of arguments${NC}"
    echo "Usage: $0 <uptime-kuma-url> <api-token>"
    echo ""
    echo "Example:"
    echo "  $0 http://localhost:3001 your-api-token"
    echo ""
    echo "To generate API token:"
    echo "  1. Login to Uptime Kuma"
    echo "  2. Go to Settings > API Keys"
    echo "  3. Click 'Add API Key'"
    echo "  4. Copy the generated token"
    exit 1
fi

UPTIME_KUMA_URL=$1
API_TOKEN=$2
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
PUSH_URLS_FILE="$PROJECT_ROOT/push-urls.txt"

# Check if jq is installed
if ! command -v jq &> /dev/null; then
    echo -e "${RED}Error: jq is not installed${NC}"
    echo "Install jq: sudo apt-get install jq"
    exit 1
fi

echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN}Uptime Kuma Monitor Setup${NC}"
echo -e "${GREEN}=========================================${NC}"
echo "Uptime Kuma URL: $UPTIME_KUMA_URL"
echo ""

# Verify API connection
echo -e "${YELLOW}Verifying API connection...${NC}"
if ! curl -s -f -H "Authorization: Bearer $API_TOKEN" "$UPTIME_KUMA_URL/api/monitors" > /dev/null; then
    echo -e "${RED}✗ Failed to connect to Uptime Kuma API${NC}"
    echo ""
    echo "Troubleshooting:"
    echo "  1. Verify Uptime Kuma is running: docker ps | grep uptime"
    echo "  2. Check API token is valid"
    echo "  3. Ensure URL is correct: $UPTIME_KUMA_URL"
    exit 1
fi
echo -e "${GREEN}✓ API connection verified${NC}"
echo ""

# Function to create HTTP monitor
create_http_monitor() {
    local name=$1
    local url=$2
    local description=$3
    local tags=$4

    echo -e "${BLUE}Creating HTTP monitor: $name${NC}"

    local response=$(curl -s -X POST "$UPTIME_KUMA_URL/api/monitors" \
        -H "Authorization: Bearer $API_TOKEN" \
        -H "Content-Type: application/json" \
        -d "{
            \"type\": \"http\",
            \"name\": \"$name\",
            \"url\": \"$url\",
            \"method\": \"GET\",
            \"interval\": 60,
            \"retryInterval\": 60,
            \"maxretries\": 1,
            \"timeout\": 10,
            \"expectedStatusCodes\": [\"200-299\"],
            \"ignoreTls\": false,
            \"upsideDown\": false,
            \"description\": \"$description\",
            \"expiryNotification\": true
        }")

    if echo "$response" | jq -e '.ok' > /dev/null 2>&1; then
        echo -e "${GREEN}✓ Monitor created successfully${NC}"
        return 0
    else
        echo -e "${RED}✗ Failed to create monitor${NC}"
        echo "Response: $response"
        return 1
    fi
}

# Function to create Push monitor
create_push_monitor() {
    local name=$1
    local interval=$2
    local description=$3
    local job_key=$4

    echo -e "${BLUE}Creating Push monitor: $name${NC}"

    local response=$(curl -s -X POST "$UPTIME_KUMA_URL/api/monitors" \
        -H "Authorization: Bearer $API_TOKEN" \
        -H "Content-Type: application/json" \
        -d "{
            \"type\": \"push\",
            \"name\": \"$name\",
            \"interval\": $interval,
            \"upsideDown\": false,
            \"description\": \"$description\"
        }")

    # Extract Push URL from response
    local push_key=$(echo "$response" | jq -r '.pushToken // empty')

    if [ -n "$push_key" ]; then
        local push_url="http://uptime-kuma:3001/api/push/$push_key"
        echo -e "${GREEN}✓ Monitor created successfully${NC}"
        echo "  Push URL: $push_url"

        # Save to file
        echo "$job_key=$push_url" >> "$PUSH_URLS_FILE"
        return 0
    else
        echo -e "${RED}✗ Failed to create monitor${NC}"
        echo "Response: $response"
        return 1
    fi
}

# Create monitors
success_count=0
fail_count=0

echo -e "${YELLOW}Creating HTTP(s) monitors...${NC}"
echo ""

# HTTP Monitor 1: CellCraft Main
if create_http_monitor \
    "CellCraft Main" \
    "https://cellcraft.app" \
    "CellCraft 메인 페이지 가용성 모니터링"; then
    ((success_count++))
else
    ((fail_count++))
fi
echo ""

# HTTP Monitor 2: Backend API Health
if create_http_monitor \
    "Backend API Health" \
    "https://cellcraft.app/api/health" \
    "Backend API 헬스 체크 엔드포인트 모니터링"; then
    ((success_count++))
else
    ((fail_count++))
fi
echo ""

echo -e "${YELLOW}Creating Push monitors...${NC}"
echo ""

# Clear push URLs file
> "$PUSH_URLS_FILE"
echo "# Uptime Kuma Push URLs for Jenkins Jobs" >> "$PUSH_URLS_FILE"
echo "# Generated: $(date)" >> "$PUSH_URLS_FILE"
echo "" >> "$PUSH_URLS_FILE"

# Push Monitor 1: Backend Unit Tests
if create_push_monitor \
    "Backend Unit Tests" \
    2100 \
    "Backend pytest 유닛 테스트 결과 모니터링 (30분마다 실행)" \
    "backend-unit-test"; then
    ((success_count++))
else
    ((fail_count++))
fi
echo ""

# Push Monitor 2: Frontend Unit Tests
if create_push_monitor \
    "Frontend Unit Tests" \
    2100 \
    "Frontend Vitest 유닛 테스트 결과 모니터링 (30분마다 실행)" \
    "frontend-unit-test"; then
    ((success_count++))
else
    ((fail_count++))
fi
echo ""

# Push Monitor 3: E2E Tests
if create_push_monitor \
    "E2E Integration Tests" \
    7500 \
    "Playwright E2E 통합 테스트 결과 모니터링 (2시간마다 실행)" \
    "e2e-test"; then
    ((success_count++))
else
    ((fail_count++))
fi
echo ""

# Push Monitor 4: Docker Health
if create_push_monitor \
    "Docker Environment" \
    15000 \
    "Docker 컨테이너 헬스 체크 모니터링 (4시간마다 실행)" \
    "docker-health-check"; then
    ((success_count++))
else
    ((fail_count++))
fi
echo ""

# Summary
echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN}Setup Summary${NC}"
echo -e "${GREEN}=========================================${NC}"
echo "Monitors created: $success_count"
if [ $fail_count -gt 0 ]; then
    echo -e "${RED}Failed: $fail_count${NC}"
else
    echo "Failed: 0"
fi
echo ""

if [ $fail_count -eq 0 ]; then
    echo -e "${GREEN}✓ All monitors created successfully!${NC}"
    echo ""
    echo "Push URLs saved to: $PUSH_URLS_FILE"
    echo ""
    echo "Next steps:"
    echo "1. Review Push URLs in: $PUSH_URLS_FILE"
    echo "2. Update Jenkins jobs with Push URLs:"
    echo "   ./scripts/update-jenkins-push-urls.sh \\"
    echo "     http://localhost:8080 admin YOUR_TOKEN push-urls.txt"
    echo ""
    echo "3. Create Status Page manually in Uptime Kuma UI:"
    echo "   - Go to Status Pages > Add New Status Page"
    echo "   - Title: CellCraft System Status"
    echo "   - Slug: status"
    echo "   - Add all 6 monitors to the page"
    echo ""
    echo "4. Verify setup:"
    echo "   ./scripts/verify-uptime-kuma-setup.sh $UPTIME_KUMA_URL"
    exit 0
else
    echo -e "${RED}✗ Some monitors failed to create${NC}"
    echo ""
    echo "Troubleshooting:"
    echo "1. Check API token permissions"
    echo "2. Verify Uptime Kuma is running"
    echo "3. Check for duplicate monitor names"
    echo "4. Review error messages above"
    exit 1
fi
