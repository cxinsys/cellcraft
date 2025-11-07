#!/bin/bash

#######################################
# CellCraft Jenkins Job Setup Script
#
# This script automatically creates Jenkins jobs using Jenkins CLI.
#
# Usage:
#   ./setup-jenkins-jobs.sh <jenkins-url> <admin-user> <admin-token>
#
# Example:
#   ./setup-jenkins-jobs.sh http://localhost:8080 admin your-api-token
#
# Prerequisites:
#   - Jenkins must be running and accessible
#   - Jenkins CLI must be available at <jenkins-url>/jnlpJars/jenkins-cli.jar
#   - Admin user must have permission to create jobs
#   - Job XML files must exist in jenkins-jobs/ directory
#
#######################################

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check arguments
if [ $# -ne 3 ]; then
    echo -e "${RED}Error: Invalid number of arguments${NC}"
    echo "Usage: $0 <jenkins-url> <admin-user> <admin-token>"
    echo "Example: $0 http://localhost:8080 admin your-api-token"
    exit 1
fi

JENKINS_URL=$1
JENKINS_USER=$2
JENKINS_TOKEN=$3
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
JOBS_DIR="$PROJECT_ROOT/jenkins-jobs"

echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN}CellCraft Jenkins Job Setup${NC}"
echo -e "${GREEN}=========================================${NC}"
echo "Jenkins URL: $JENKINS_URL"
echo "Jenkins User: $JENKINS_USER"
echo "Jobs Directory: $JOBS_DIR"
echo ""

# Check if jenkins-jobs directory exists
if [ ! -d "$JOBS_DIR" ]; then
    echo -e "${RED}Error: jenkins-jobs directory not found at $JOBS_DIR${NC}"
    exit 1
fi

# Download Jenkins CLI
echo -e "${YELLOW}Downloading Jenkins CLI...${NC}"
CLI_JAR="$SCRIPT_DIR/jenkins-cli.jar"
wget -q -O "$CLI_JAR" "$JENKINS_URL/jnlpJars/jenkins-cli.jar" || {
    echo -e "${RED}Error: Failed to download Jenkins CLI${NC}"
    echo "Please check if Jenkins is running and accessible at $JENKINS_URL"
    exit 1
}
echo -e "${GREEN}✓ Jenkins CLI downloaded${NC}"
echo ""

# Function to create a Jenkins job
create_job() {
    local job_config=$1
    local job_name=$(basename "$job_config" .xml)

    echo -e "${YELLOW}Creating job: $job_name${NC}"

    # Check if job already exists
    if java -jar "$CLI_JAR" -s "$JENKINS_URL" -auth "$JENKINS_USER:$JENKINS_TOKEN" \
        get-job "$job_name" >/dev/null 2>&1; then
        echo -e "${YELLOW}  Job already exists. Updating...${NC}"
        java -jar "$CLI_JAR" -s "$JENKINS_URL" -auth "$JENKINS_USER:$JENKINS_TOKEN" \
            update-job "$job_name" < "$job_config" && \
            echo -e "${GREEN}  ✓ Job updated successfully${NC}" || \
            echo -e "${RED}  ✗ Failed to update job${NC}"
    else
        # Create new job
        java -jar "$CLI_JAR" -s "$JENKINS_URL" -auth "$JENKINS_USER:$JENKINS_TOKEN" \
            create-job "$job_name" < "$job_config" && \
            echo -e "${GREEN}  ✓ Job created successfully${NC}" || \
            echo -e "${RED}  ✗ Failed to create job${NC}"
    fi
    echo ""
}

# Create all jobs
echo -e "${YELLOW}Creating Jenkins jobs...${NC}"
echo ""

job_count=0
for job_config in "$JOBS_DIR"/*.xml; do
    if [ -f "$job_config" ]; then
        create_job "$job_config"
        ((job_count++))
    fi
done

# Cleanup
rm -f "$CLI_JAR"

echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN}Setup Complete!${NC}"
echo -e "${GREEN}=========================================${NC}"
echo "Total jobs created/updated: $job_count"
echo ""
echo "Next steps:"
echo "1. Access Jenkins at $JENKINS_URL"
echo "2. Verify all jobs are created"
echo "3. Update UPTIME_KUMA_PUSH_URL parameters in each job"
echo "4. Run test builds to verify functionality"
echo ""
