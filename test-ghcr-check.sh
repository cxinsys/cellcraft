#!/bin/bash

# Color codes
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly PURPLE='\033[0;35m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m'

# Logging functions
log_info() { echo -e "${BLUE}ℹ️  $1${NC}"; }
log_success() { echo -e "${GREEN}✅ $1${NC}"; }
log_step() { echo -e "${PURPLE}📋 $1${NC}"; }
log_warning() { echo -e "${YELLOW}⚠️  $1${NC}"; }
log_header() { echo -e "\n${CYAN}🚀 $1${NC}\n"; }

# The check function
check_ghcr_availability() {
    log_header "GHCR Image Availability Check"

    # Test images for CPU mode with full repository path
    local images=(
        "ghcr.io/cxinsys/cellcraft/frontend:v1.0.0"
        "ghcr.io/cxinsys/cellcraft/backend-cpu:v1.0.0"
        "ghcr.io/cxinsys/cellcraft/celery-cpu:v1.0.0"
    )

    # First, check if all images exist locally
    log_step "Checking for locally available GHCR images"
    local all_local=true
    local missing_count=0

    for image in "${images[@]}"; do
        if docker inspect "${image}" >/dev/null 2>&1; then
            log_success "Found locally: ${image}"
        else
            all_local=false
            ((missing_count++))
        fi
    done

    if [[ "$all_local" == "true" ]]; then
        log_success "All GHCR images available locally (${#images[@]}/${#images[@]})"
        return 0
    fi

    log_info "Missing ${missing_count}/${#images[@]} images locally"

    # If not all local, check remote accessibility
    local test_image="ghcr.io/cxinsys/cellcraft/frontend:v1.0.0"
    log_step "Testing GHCR remote accessibility with: $test_image"

    if timeout 30 docker manifest inspect "$test_image" >/dev/null 2>&1; then
        log_success "GHCR images are accessible (remote)"
        return 0
    else
        log_warning "GHCR access failed (network/auth issue)"
        return 1
    fi
}

# Run the test
echo "===== Testing check_ghcr_availability function ====="
if check_ghcr_availability; then
    echo "Function returned: SUCCESS (exit code 0)"
    echo "Strategy should be: GHCR images"
else
    echo "Function returned: FAILURE (exit code 1)"
    echo "Strategy should be: Local build"
fi
