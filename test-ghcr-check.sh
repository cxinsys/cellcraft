#!/usr/bin/env bash

#######################################################################
# CellCraft GHCR Image Checker and Downloader
#
# Purpose:
# - Check availability of GHCR images (local and remote)
# - Allow users to select CPU or GPU mode
# - Pre-download images for selected deployment mode
#
# Usage:
#   ./test-ghcr-check.sh              # Interactive mode
#   ./test-ghcr-check.sh --cpu        # Check/download CPU images
#   ./test-ghcr-check.sh --gpu        # Check/download GPU images
#   ./test-ghcr-check.sh --check-only # Check without download
#######################################################################

# Color codes
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly PURPLE='\033[0;35m'
readonly CYAN='\033[0;36m'
readonly RED='\033[0;31m'
readonly NC='\033[0m'

# Logging functions
log_info() { echo -e "${BLUE}ℹ️  $1${NC}"; }
log_success() { echo -e "${GREEN}✅ $1${NC}"; }
log_step() { echo -e "${PURPLE}📋 $1${NC}"; }
log_warning() { echo -e "${YELLOW}⚠️  $1${NC}"; }
log_error() { echo -e "${RED}❌ $1${NC}"; }
log_header() { echo -e "\n${CYAN}🚀 $1${NC}\n"; }

# Image definitions (using simple arrays instead of associative for compatibility)
# Use 'latest' tag to always get the most recent multi-platform builds from GitHub Actions
readonly CPU_FRONTEND="ghcr.io/cxinsys/cellcraft/frontend:latest"
readonly CPU_BACKEND="ghcr.io/cxinsys/cellcraft/backend-cpu:latest"
readonly CPU_CELERY="ghcr.io/cxinsys/cellcraft/celery-cpu:latest"
readonly GPU_FRONTEND="ghcr.io/cxinsys/cellcraft/frontend:latest"
readonly GPU_BACKEND="ghcr.io/cxinsys/cellcraft/backend-gpu:latest"
readonly GPU_CELERY="ghcr.io/cxinsys/cellcraft/celery-gpu:latest"

# Check if a single image exists locally
check_image_local() {
    local image=$1
    if docker inspect "${image}" >/dev/null 2>&1; then
        return 0
    else
        return 1
    fi
}

# Check if a single image is accessible remotely
check_image_remote() {
    local image=$1
    if timeout 30 docker manifest inspect "$image" >/dev/null 2>&1; then
        return 0
    else
        return 1
    fi
}

# Check images for a specific mode
check_mode_images() {
    local mode=$1
    local images=()

    if [[ "$mode" == "cpu" ]]; then
        images=("$CPU_FRONTEND" "$CPU_BACKEND" "$CPU_CELERY")
    elif [[ "$mode" == "gpu" ]]; then
        images=("$GPU_FRONTEND" "$GPU_BACKEND" "$GPU_CELERY")
    else
        log_error "Invalid mode: $mode"
        return 1
    fi

    log_header "Checking $(echo "$mode" | tr '[:lower:]' '[:upper:]') Mode Images"

    local local_count=0
    local remote_count=0
    local missing_count=0

    for image in "${images[@]}"; do
        local image_name=$(echo "$image" | sed 's/ghcr.io\/cxinsys\/cellcraft\///')

        if check_image_local "$image"; then
            log_success "LOCAL:  ${image_name}"
            ((local_count++))
        elif check_image_remote "$image"; then
            log_warning "REMOTE: ${image_name}"
            ((remote_count++))
        else
            log_error "MISSING: ${image_name}"
            ((missing_count++))
        fi
    done

    echo ""
    log_info "Summary for $(echo "$mode" | tr '[:lower:]' '[:upper:]') mode:"
    echo "  • Local images:  ${local_count}/${#images[@]}"
    echo "  • Remote images: ${remote_count}/${#images[@]}"
    echo "  • Missing images: ${missing_count}/${#images[@]}"

    if [[ $missing_count -gt 0 ]]; then
        return 1
    fi

    return 0
}

# Download images for a specific mode
download_mode_images() {
    local mode=$1
    local images=()

    if [[ "$mode" == "cpu" ]]; then
        images=("$CPU_FRONTEND" "$CPU_BACKEND" "$CPU_CELERY")
    elif [[ "$mode" == "gpu" ]]; then
        images=("$GPU_FRONTEND" "$GPU_BACKEND" "$GPU_CELERY")
    else
        log_error "Invalid mode: $mode"
        return 1
    fi

    log_header "Downloading $(echo "$mode" | tr '[:lower:]' '[:upper:]') Mode Images"

    local success_count=0
    local skip_count=0
    local fail_count=0

    for image in "${images[@]}"; do
        local image_name=$(echo "$image" | sed 's/ghcr.io\/cxinsys\/cellcraft\///')

        # Skip if already exists locally
        if check_image_local "$image"; then
            log_info "SKIP: ${image_name} (already exists locally)"
            ((skip_count++))
            continue
        fi

        log_step "Pulling: ${image_name}..."

        if docker pull "$image"; then
            log_success "Downloaded: ${image_name}"
            ((success_count++))
        else
            log_error "Failed: ${image_name}"
            ((fail_count++))
        fi
    done

    echo ""
    log_info "Download Summary:"
    echo "  • Successfully downloaded: ${success_count}"
    echo "  • Skipped (already local): ${skip_count}"
    echo "  • Failed: ${fail_count}"

    if [[ $fail_count -gt 0 ]]; then
        return 1
    fi

    return 0
}

# Interactive menu
show_menu() {
    echo ""
    echo "======================================"
    echo "  CellCraft GHCR Image Manager"
    echo "======================================"
    echo ""
    echo "Select an option:"
    echo ""
    echo "  1) Check CPU mode images"
    echo "  2) Check GPU mode images"
    echo "  3) Check all images (CPU + GPU)"
    echo "  4) Download CPU mode images"
    echo "  5) Download GPU mode images"
    echo "  6) Exit"
    echo ""
    read -p "Enter choice [1-6]: " choice

    case $choice in
        1)
            check_mode_images "cpu"
            ;;
        2)
            check_mode_images "gpu"
            ;;
        3)
            check_mode_images "cpu"
            check_mode_images "gpu"
            ;;
        4)
            check_mode_images "cpu"
            echo ""
            read -p "Proceed with download? [y/N]: " confirm
            if [[ "$confirm" =~ ^[Yy]$ ]]; then
                download_mode_images "cpu"
            else
                log_info "Download cancelled"
            fi
            ;;
        5)
            check_mode_images "gpu"
            echo ""
            read -p "Proceed with download? [y/N]: " confirm
            if [[ "$confirm" =~ ^[Yy]$ ]]; then
                download_mode_images "gpu"
            else
                log_info "Download cancelled"
            fi
            ;;
        6)
            log_info "Exiting..."
            exit 0
            ;;
        *)
            log_error "Invalid choice. Please enter 1-6."
            show_menu
            ;;
    esac
}

# Main execution
main() {
    # Parse command line arguments
    if [[ $# -eq 0 ]]; then
        # Interactive mode
        show_menu
    else
        case "$1" in
            --cpu)
                check_mode_images "cpu"
                echo ""
                read -p "Download missing images? [y/N]: " confirm
                if [[ "$confirm" =~ ^[Yy]$ ]]; then
                    download_mode_images "cpu"
                fi
                ;;
            --gpu)
                check_mode_images "gpu"
                echo ""
                read -p "Download missing images? [y/N]: " confirm
                if [[ "$confirm" =~ ^[Yy]$ ]]; then
                    download_mode_images "gpu"
                fi
                ;;
            --check-only)
                check_mode_images "cpu"
                check_mode_images "gpu"
                ;;
            --help|-h)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  --cpu         Check and optionally download CPU mode images"
                echo "  --gpu         Check and optionally download GPU mode images"
                echo "  --check-only  Check all images without downloading"
                echo "  --help, -h    Show this help message"
                echo ""
                echo "Without options, runs in interactive mode."
                ;;
            *)
                log_error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    fi
}

# Run main
main "$@"
