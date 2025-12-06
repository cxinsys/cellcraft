#!/bin/bash

# CellCraft CPU-Only Mode Launcher
# Automatically switches to CPU-compatible plugin branch and launches CPU-only configuration
# Author: CellCraft Team
# Usage: ./run-cpu-mode.sh [--clean|--build|--help]

set -euo pipefail  # Exit on any error, undefined variable, or pipe failure

# Color codes for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m' 
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly PURPLE='\033[0;35m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m' # No Color

# Configuration
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PLUGIN_DIR="${SCRIPT_DIR}/backend/plugin/official"
readonly COMPOSE_FILE_GHCR="${SCRIPT_DIR}/docker-compose.cpu.yml"
readonly COMPOSE_FILE_LOCAL="${SCRIPT_DIR}/docker-compose.cpu.local.yml"
readonly CHECK_SCRIPT="${SCRIPT_DIR}/check-installation.sh"
readonly CPU_BRANCH="release/plugins-v1.0-cpu"
readonly VERSION_FILE="${PLUGIN_DIR}/version.json"

# Script options
CLEAN_CONTAINERS=false
FORCE_REBUILD=false
SKIP_VERIFICATION=false

# Backup variables
ORIGINAL_BRANCH=""
ORIGINAL_COMMIT=""

# Functions for colored output
log_info() { echo -e "${BLUE}ℹ️  $1${NC}" >&2; }
log_success() { echo -e "${GREEN}✅ $1${NC}" >&2; }
log_warning() { echo -e "${YELLOW}⚠️  $1${NC}" >&2; }
log_error() { echo -e "${RED}❌ $1${NC}" >&2; }
log_header() { echo -e "\n${CYAN}🚀 $1${NC}\n" >&2; }
log_step() { echo -e "${PURPLE}📋 $1${NC}" >&2; }

# Print script header
print_header() {
    echo -e "${CYAN}"
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                    CellCraft CPU-Only Mode                     ║"
    echo "║            Automated Setup and Launch Script                  ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

# Print usage information
print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Launch CellCraft in CPU-only mode with automatic plugin branch switching.

OPTIONS:
    --clean     Clean existing containers before launching
    --build     Force rebuild of all containers
    --help      Show this help message
    --skip-verify  Skip final installation verification

EXAMPLES:
    $0                    # Basic launch
    $0 --clean            # Clean containers first
    $0 --build            # Force rebuild
    $0 --clean --build    # Clean and rebuild

This script will:
1. Switch plugin submodule to CPU-compatible branch
2. Launch Docker Compose with CPU-only configuration
3. Verify installation success
4. Provide access URLs

EOF
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --clean)
                CLEAN_CONTAINERS=true
                shift
                ;;
            --build)
                FORCE_REBUILD=true
                shift
                ;;
            --skip-verify)
                SKIP_VERIFICATION=true
                shift
                ;;
            --help|-h)
                print_usage
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                print_usage
                exit 1
                ;;
        esac
    done
}

# Check GHCR image availability
# Strategy: Local-first checking for instant startup when images are cached
# 1. Check if all images exist locally → return 0 (will run without pulling)
# 2. If some missing → check remote accessibility → return 0 (will pull missing images)
# 3. If remote inaccessible → return 1 (fallback to local build)
check_ghcr_availability() {
    log_header "GHCR Image Availability Check"

    # Test images for CPU mode with full repository path
    # Use 'latest' tag to always get the most recent multi-platform builds from GitHub Actions
    local images=(
        "ghcr.io/cxinsys/cellcraft/frontend:latest"
        "ghcr.io/cxinsys/cellcraft/backend-cpu:latest"
        "ghcr.io/cxinsys/cellcraft/celery-cpu:latest"
    )

    # First, check if all images exist locally (using docker inspect)
    # This is very fast and allows instant deployment when images are cached
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

    # If all images are available locally, use them immediately without pulling
    if [[ "$all_local" == "true" ]]; then
        log_success "All GHCR images available locally (${#images[@]}/${#images[@]})"
        log_info "Deployment will use local images without pulling (instant startup)"
        return 0
    fi

    log_info "Missing ${missing_count}/${#images[@]} images locally"

    # If not all local, check remote accessibility before attempting to use GHCR strategy
    local test_image="ghcr.io/cxinsys/cellcraft/frontend:latest"
    log_step "Testing GHCR remote accessibility with: $test_image"

    if timeout 30 docker manifest inspect "$test_image" >/dev/null 2>&1; then
        log_success "GHCR images are accessible (remote)"
        log_info "Deployment will pull missing ${missing_count} images from GHCR"
        return 0
    else
        log_warning "GHCR access failed (network/auth issue)"
        return 1
    fi
}

# Determine compose strategy
determine_compose_strategy() {
    if check_ghcr_availability; then
        log_success "Strategy: GHCR images (fast deployment)"
        echo "ghcr"
    else
        log_warning "Strategy: Local build (slower but reliable)"
        echo "local"
    fi
}

# Check prerequisites
check_prerequisites() {
    log_header "Prerequisites Check"
    
    # Check if we're in the correct directory
    if [[ ! -f "${COMPOSE_FILE_GHCR}" ]]; then
        log_error "docker-compose.cpu.yml not found in current directory"
        log_error "Please run this script from the CellCraft root directory"
        exit 1
    fi
    log_success "Found docker-compose.cpu.yml"
    
    # Check local compose file
    if [[ ! -f "${COMPOSE_FILE_LOCAL}" ]]; then
        log_error "docker-compose.cpu.local.yml not found in current directory"
        log_error "Please ensure all compose files are present"
        exit 1
    fi
    log_success "Found docker-compose.cpu.local.yml"
    
    # Check if plugin directory exists
    if [[ ! -d "${PLUGIN_DIR}" ]]; then
        log_error "Plugin directory not found: ${PLUGIN_DIR}"
        log_error "Please ensure the plugin submodule is initialized"
        exit 1
    fi
    log_success "Plugin directory found"
    
    # Check Docker
    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed or not in PATH"
        exit 1
    fi
    log_success "Docker is available"
    
    # Check Docker daemon
    if ! docker info &> /dev/null; then
        log_error "Docker daemon is not running"
        log_error "Please start Docker and try again"
        exit 1
    fi
    log_success "Docker daemon is running"
    
    # Check Docker Compose
    if ! docker compose version &> /dev/null; then
        log_error "Docker Compose is not available"
        exit 1
    fi
    log_success "Docker Compose is available"
    
    # Check Git
    if ! command -v git &> /dev/null; then
        log_error "Git is not installed or not in PATH"
        exit 1
    fi
    log_success "Git is available"
}

# Backup current git state
backup_git_state() {
    log_step "Backing up current git state"
    
    cd "${PLUGIN_DIR}"
    
    # Get current branch or commit
    if ORIGINAL_BRANCH=$(git symbolic-ref --short HEAD 2>/dev/null); then
        log_info "Current branch: ${ORIGINAL_BRANCH}"
    else
        # We're in detached HEAD state
        ORIGINAL_COMMIT=$(git rev-parse HEAD)
        log_info "Currently in detached HEAD state at commit: ${ORIGINAL_COMMIT:0:8}"
    fi
    
    cd "${SCRIPT_DIR}"
}

# Switch to CPU branch
switch_to_cpu_branch() {
    log_header "Plugin Branch Configuration"
    
    cd "${PLUGIN_DIR}"
    
    # Fetch latest changes
    log_step "Fetching latest changes"
    git fetch origin || {
        log_warning "Failed to fetch from remote, continuing with local branches"
    }
    
    # Check if CPU branch exists
    if ! git show-ref --verify --quiet "refs/remotes/origin/${CPU_BRANCH}"; then
        log_error "CPU branch '${CPU_BRANCH}' not found in remote"
        log_error "Available branches:"
        git branch -r | head -10
        exit 1
    fi
    
    # Switch to CPU branch
    log_step "Switching to CPU-compatible branch: ${CPU_BRANCH}"
    if git checkout "${CPU_BRANCH}" 2>/dev/null; then
        log_success "Successfully switched to ${CPU_BRANCH}"
    elif git checkout -b "${CPU_BRANCH}" "origin/${CPU_BRANCH}" 2>/dev/null; then
        log_success "Created and switched to ${CPU_BRANCH}"
    else
        log_error "Failed to switch to ${CPU_BRANCH}"
        restore_git_state
        exit 1
    fi
    
    # Verify version.json exists
    if [[ -f "${VERSION_FILE}" ]]; then
        log_success "Version configuration file found"
        # Display version info if possible
        if command -v jq &> /dev/null && jq -e . "${VERSION_FILE}" >/dev/null 2>&1; then
            local branch_info
            branch_info=$(jq -r '.branch // "unknown"' "${VERSION_FILE}")
            local version_info
            version_info=$(jq -r '.version // "unknown"' "${VERSION_FILE}")
            log_info "Branch: ${branch_info}, Version: ${version_info}"
        fi
    else
        log_warning "Version configuration file not found, but continuing"
    fi
    
    cd "${SCRIPT_DIR}"
}

# Restore original git state
restore_git_state() {
    log_step "Restoring original git state"
    
    cd "${PLUGIN_DIR}"
    
    if [[ -n "${ORIGINAL_BRANCH}" ]]; then
        git checkout "${ORIGINAL_BRANCH}" || true
        log_info "Restored to branch: ${ORIGINAL_BRANCH}"
    elif [[ -n "${ORIGINAL_COMMIT}" ]]; then
        git checkout "${ORIGINAL_COMMIT}" || true
        log_info "Restored to commit: ${ORIGINAL_COMMIT:0:8}"
    fi
    
    cd "${SCRIPT_DIR}"
}

# Clean existing containers
clean_containers() {
    if [[ "${CLEAN_CONTAINERS}" == "true" ]]; then
        log_header "Container Cleanup"
        
        log_step "Stopping existing containers (trying both compose files)"
        docker compose -f "${COMPOSE_FILE_GHCR}" down 2>/dev/null || true
        docker compose -f "${COMPOSE_FILE_LOCAL}" down 2>/dev/null || {
            log_warning "Some containers may not have been running"
        }
        
        log_step "Removing unused images and volumes"
        docker system prune -f || {
            log_warning "System prune completed with warnings"
        }
        
        log_success "Container cleanup completed"
    fi
}

# Detect if running on Mac
detect_mac() {
    if [[ "$(uname -s)" == "Darwin" ]]; then
        # Check if it's Apple Silicon (ARM64)
        if [[ "$(uname -m)" == "arm64" ]]; then
            echo "mac-arm64"
        else
            echo "mac-intel"
        fi
    else
        echo "other"
    fi
}

# Build and launch containers with smart strategy
launch_containers() {
    log_header "Smart Container Launch"
    
    # Detect platform
    local platform_type
    platform_type=$(detect_mac)
    
    # Determine deployment strategy
    local strategy
    strategy=$(determine_compose_strategy)
    
    local compose_file
    local docker_args=""
    
    # On Mac ARM64 with GHCR strategy, use docker-compose.cpu.yml
    # Docker Compose will automatically select ARM64 variant from multi-platform manifest
    if [[ "$platform_type" == "mac-arm64" && "$strategy" == "ghcr" ]]; then
        compose_file="${COMPOSE_FILE_GHCR}"
        log_info "Using GHCR images on Mac ARM64"
        log_info "Docker will automatically pull ARM64 variant (no Rosetta needed)"
    elif [[ "$strategy" == "ghcr" ]]; then
        compose_file="${COMPOSE_FILE_GHCR}"
        log_info "Using GHCR images: docker compose up -d"
    else
        compose_file="${COMPOSE_FILE_LOCAL}"
        docker_args="--build"
        log_info "Using local build: docker compose up -d --build"
    fi
    
    log_step "Compose file: $(basename "$compose_file")"
    
    # Handle force rebuild option
    if [[ "${FORCE_REBUILD}" == "true" ]]; then
        docker_args="--build"
        log_step "Force rebuild enabled"
    fi
    
    # Launch containers
    log_step "Starting CellCraft services in CPU-only mode"
    if docker compose -f "$compose_file" up -d $docker_args; then
        log_success "All containers launched successfully"
    else
        log_error "Failed to launch containers"
        log_error "Check Docker logs for details:"
        log_error "  docker compose -f $compose_file logs"
        docker compose -f "$compose_file" logs
        restore_git_state
        exit 1
    fi
    
    # Wait for services to be ready
    log_step "Waiting for services to initialize (30 seconds)"
    sleep 30
    
    # Check container status and handle unhealthy containers
    log_step "Checking container status"

    # Get all service names
    local all_services
    all_services=$(docker compose -f "$compose_file" config --services)
    
    # Check for unhealthy containers (using dynamic compose file)
    local unhealthy_containers
    unhealthy_containers=$(docker compose -f "$compose_file" ps --format "table {{.Service}}\t{{.Status}}" | grep -E "(unhealthy|exited|restarting)" || true | awk '{print $1}' | tr '\n' ' ')
    
    if [[ -n "${unhealthy_containers// /}" ]]; then
        log_warning "Found unhealthy containers: ${unhealthy_containers}"
        log_step "Attempting to restart unhealthy containers"
        
        # Show logs for unhealthy containers
        for container in ${unhealthy_containers}; do
            if [[ -n "${container}" ]]; then
                log_info "Logs for ${container}:"
                docker compose -f "$compose_file" logs --tail=20 "${container}" || true
            fi
        done
        
        # Restart the entire stack
        log_step "Restarting all services"
        if docker compose -f "$compose_file" up -d --force-recreate; then
            log_success "Services restarted successfully"
            # Wait again for services to stabilize
            log_step "Waiting for services to stabilize (15 seconds)"
            sleep 15
        else
            log_error "Failed to restart services"
            log_error "Full service logs:"
            docker compose -f "$compose_file" logs
            restore_git_state
            exit 1
        fi
    fi
    
    # Final status check
    local running_containers
    running_containers=$(docker compose -f "$compose_file" ps --services --filter "status=running" | wc -l)
    log_info "Running containers: ${running_containers}"
    
    # Show container status
    docker compose -f "$compose_file" ps
}

# Verify installation
verify_installation() {
    if [[ "${SKIP_VERIFICATION}" == "true" ]]; then
        log_header "Skipping Installation Verification"
        return 0
    fi
    
    log_header "Installation Verification"
    
    if [[ -f "${CHECK_SCRIPT}" ]]; then
        log_step "Running installation check script"
        chmod +x "${CHECK_SCRIPT}"
        
        if "${CHECK_SCRIPT}"; then
            log_success "Installation verification passed"
        else
            log_warning "Installation verification found issues"
            log_warning "Check the output above for details"
        fi
    else
        log_warning "Installation check script not found, performing basic checks"
        
        # Basic connectivity tests
        log_step "Testing service connectivity"
        
        # Test frontend
        if curl -s --connect-timeout 10 http://localhost:8080 >/dev/null; then
            log_success "Frontend accessible at http://localhost:8080"
        else
            log_warning "Frontend not immediately accessible"
        fi
        
        # Test backend
        if curl -s --connect-timeout 10 http://localhost:8000/docs >/dev/null; then
            log_success "Backend API accessible at http://localhost:8000"
        else
            log_warning "Backend API not immediately accessible"
        fi
    fi
}

# Show final status and access information
show_final_status() {
    log_header "Launch Complete"
    
    echo -e "${GREEN}🎉 CellCraft is now running in CPU-only mode!${NC}"
    echo ""
    echo -e "${BLUE}📍 Access Points:${NC}"
    echo -e "  • ${CYAN}Frontend Application:${NC} http://localhost:8080"
    echo -e "  • ${CYAN}Backend API:${NC} http://localhost:8000"
    echo -e "  • ${CYAN}API Documentation:${NC} http://localhost:8000/docs"
    echo -e "  • ${CYAN}RabbitMQ Management:${NC} http://localhost:15672 (guest/guest)"
    echo ""
    echo -e "${BLUE}🔧 Management Commands:${NC}"
    echo -e "  • View logs: ${CYAN}docker compose -f [compose-file] logs${NC}"
    echo -e "  • Stop services: ${CYAN}docker compose -f [compose-file] down${NC}"
    echo -e "  • Restart services: ${CYAN}docker compose -f [compose-file] restart${NC}"
    echo ""
    echo -e "${BLUE}📊 Current Configuration:${NC}"
    echo -e "  • Mode: CPU-Only (GPU plugins excluded)"
    echo -e "  • Plugin Branch: ${CPU_BRANCH}"
    echo -e "  • Compatible Plugins: GENIE3, GRNBoost2, GRNViz, LEAP, Scribe, TENET"
    echo ""
    echo -e "${YELLOW}💡 Note: Keep this terminal open to see any important messages${NC}"
}

# Error handler
error_handler() {
    local exit_code=$?
    log_error "Script failed with exit code ${exit_code}"
    log_error "Attempting to restore original state"
    restore_git_state
    exit ${exit_code}
}

# Cleanup handler
cleanup_handler() {
    log_info "Script interrupted, cleaning up..."
    restore_git_state
    exit 130
}

# Set up error and interrupt handlers
trap error_handler ERR
trap cleanup_handler INT TERM

# Main execution
main() {
    print_header
    
    # Parse command line arguments
    parse_arguments "$@"
    
    # Display configuration
    echo -e "${BLUE}Configuration:${NC}"
    echo -e "  Clean containers: ${CLEAN_CONTAINERS}"
    echo -e "  Force rebuild: ${FORCE_REBUILD}"
    echo -e "  Skip verification: ${SKIP_VERIFICATION}"
    echo ""
    
    # Execute main workflow
    check_prerequisites
    backup_git_state
    switch_to_cpu_branch
    clean_containers
    launch_containers
    verify_installation
    show_final_status
    
    log_success "CPU-only mode launch completed successfully!"
}

# Execute main function with all arguments
main "$@"