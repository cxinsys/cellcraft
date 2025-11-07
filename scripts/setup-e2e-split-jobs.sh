#!/bin/bash

# Script to upload the 6 split E2E test Jenkins jobs
# Prerequisites: Jenkins must be running and accessible

set -e

JENKINS_URL="http://localhost:9090"
JENKINS_USER="admin"
JENKINS_TOKEN="${JENKINS_API_TOKEN:-}"

# Color codes for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo "========================================="
echo "E2E Test Jenkins Jobs Setup (6 jobs)"
echo "========================================="

# Check if Jenkins API token is provided
if [ -z "$JENKINS_TOKEN" ]; then
    echo -e "${YELLOW}Warning: JENKINS_API_TOKEN environment variable not set${NC}"
    echo "Please set it with: export JENKINS_API_TOKEN=your_token_here"
    echo "You can get your token from: ${JENKINS_URL}/user/admin/configure"
    exit 1
fi

# Function to create/update Jenkins job
create_job() {
    local job_name=$1
    local xml_file=$2

    echo -e "\n${YELLOW}Processing job: ${job_name}${NC}"

    # Check if job exists
    if curl -s -f -u "${JENKINS_USER}:${JENKINS_TOKEN}" \
        "${JENKINS_URL}/job/${job_name}/api/json" > /dev/null 2>&1; then
        echo "Job exists, updating..."
        curl -X POST -u "${JENKINS_USER}:${JENKINS_TOKEN}" \
            "${JENKINS_URL}/job/${job_name}/config.xml" \
            --data-binary @"${xml_file}" \
            -H "Content-Type: application/xml"
        echo -e "${GREEN}✓ Job updated: ${job_name}${NC}"
    else
        echo "Job doesn't exist, creating..."
        curl -X POST -u "${JENKINS_USER}:${JENKINS_TOKEN}" \
            "${JENKINS_URL}/createItem?name=${job_name}" \
            --data-binary @"${xml_file}" \
            -H "Content-Type: application/xml"
        echo -e "${GREEN}✓ Job created: ${job_name}${NC}"
    fi
}

# Get the script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
JENKINS_JOBS_DIR="${PROJECT_ROOT}/jenkins-jobs"

echo -e "\nProject root: ${PROJECT_ROOT}"
echo "Jenkins jobs directory: ${JENKINS_JOBS_DIR}"

# Create/update all 6 E2E test jobs
create_job "e2e-01-file-assignment" "${JENKINS_JOBS_DIR}/e2e-01-file-assignment.xml"
create_job "e2e-02-data-display" "${JENKINS_JOBS_DIR}/e2e-02-data-display.xml"
create_job "e2e-03-scatter-plot" "${JENKINS_JOBS_DIR}/e2e-03-scatter-plot.xml"
create_job "e2e-04-algorithm-config" "${JENKINS_JOBS_DIR}/e2e-04-algorithm-config.xml"
create_job "e2e-05-workflow-execution" "${JENKINS_JOBS_DIR}/e2e-05-workflow-execution.xml"
create_job "e2e-06-result-visualization" "${JENKINS_JOBS_DIR}/e2e-06-result-visualization.xml"

echo -e "\n========================================="
echo -e "${GREEN}All E2E test jobs have been set up!${NC}"
echo "========================================="

echo -e "\nExecution schedule (every hour):"
echo "  :00 - Test 01 (InputFile Node)"
echo "  :10 - Test 02 (DataTable Node)"
echo "  :20 - Test 03 (ScatterPlot Node)"
echo "  :30 - Test 04 (Algorithm Node)"
echo "  :40 - Test 05 (Workflow Execution)"
echo "  :50 - Test 06 (Result Visualization)"
echo ""
echo "Next steps:"
echo "1. Verify jobs in Jenkins UI: ${JENKINS_URL}"
echo "2. Trigger Test 05 manually first (creates SUCCESS workflow for Test 06)"
echo "3. Monitor Uptime Kuma dashboard: http://localhost:3001"
