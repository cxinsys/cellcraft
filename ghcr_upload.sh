#!/bin/bash

# CellCraft GHCR Upload Script
# This script builds and uploads Docker images to GitHub Container Registry
#
# Images to be uploaded:
# 1. Frontend: frontend/Dockerfile.local -> ghcr.io/cxinsys/cellcraft/frontend
# 2. Backend CPU: backend/Dockerfile.improved -> ghcr.io/cxinsys/cellcraft/backend-cpu  
# 3. Backend GPU: backend/Dockerfile -> ghcr.io/cxinsys/cellcraft/backend-gpu
# 4. Celery CPU: backend/Dockerfile.improved.celery -> ghcr.io/cxinsys/cellcraft/celery-cpu
# 5. Celery GPU: backend/Dockerfile.celery -> ghcr.io/cxinsys/cellcraft/celery-gpu

set -e

# Configuration
GITHUB_USERNAME="${GITHUB_USERNAME:-cxinsys}"
REPO_NAME="${REPO_NAME:-cellcraft}"
REGISTRY="ghcr.io"
BASE_IMAGE_NAME="${REGISTRY}/${GITHUB_USERNAME}/${REPO_NAME}"
# Optional release tag to add alongside auto version and latest
# Defaults to initial release version; override with: RELEASE_VERSION=v1.2.3 ./ghcr_upload.sh
RELEASE_VERSION="${RELEASE_VERSION:-v1.0.0}"


# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_header() {
    echo -e "${PURPLE}[HEADER]${NC} $1"
}

log_step() {
    echo -e "${CYAN}[STEP]${NC} $1"
}

# Show image details
show_image_details() {
    echo "========================================"
    log_header "Docker Images to be built and uploaded:"
    echo "========================================"
    echo "1. Frontend (Vue.js)"
    echo "   - Dockerfile: frontend/Dockerfile.local"
    echo "   - Image: ${BASE_IMAGE_NAME}/frontend"
    echo "   - Platform: All (CPU/GPU compatible)"
    echo ""
    echo "2. Backend CPU (FastAPI)"
    echo "   - Dockerfile: backend/Dockerfile.improved"
    echo "   - Image: ${BASE_IMAGE_NAME}/backend-cpu"
    echo "   - Platform: CPU only"
    echo ""
    echo "3. Backend GPU (FastAPI)"
    echo "   - Dockerfile: backend/Dockerfile"
    echo "   - Image: ${BASE_IMAGE_NAME}/backend-gpu"
    echo "   - Platform: GPU enabled"
    echo ""
    echo "4. Celery CPU (Task Worker)"
    echo "   - Dockerfile: backend/Dockerfile.improved.celery"
    echo "   - Image: ${BASE_IMAGE_NAME}/celery-cpu"
    echo "   - Platform: CPU only"
    echo ""
    echo "5. Celery GPU (Task Worker)"
    echo "   - Dockerfile: backend/Dockerfile.celery"
    echo "   - Image: ${BASE_IMAGE_NAME}/celery-gpu"
    echo "   - Platform: GPU enabled"
    echo "========================================"
    echo ""
}

# Check prerequisites
check_prerequisites() {
    log_step "Checking prerequisites..."
    
    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed or not in PATH"
        exit 1
    fi
    
    if ! docker info &> /dev/null; then
        log_error "Docker daemon is not running or accessible"
        exit 1
    fi
    
    if [ -z "$GITHUB_TOKEN" ]; then
        log_error "GITHUB_TOKEN environment variable is not set"
        log_info "Please set your GitHub Personal Access Token:"
        log_info "export GITHUB_TOKEN=your_token_here"
        log_info ""
        log_info "Token should have 'write:packages' and 'read:packages' permissions"
        exit 1
    fi
    
    # Check if we're in the correct directory
    if [ ! -f "docker-compose.cpu.yml" ] || [ ! -f "docker-compose.gpu.yml" ]; then
        log_error "Please run this script from the CellCraft root directory"
        log_info "Expected files: docker-compose.cpu.yml, docker-compose.gpu.yml"
        exit 1
    fi
    
    log_success "Prerequisites check passed"
}

# Check buildx capabilities
check_buildx() {
    log_step "Checking Docker Buildx capabilities..."
    
    # Check if buildx is available
    if ! docker buildx version &> /dev/null; then
        log_error "Docker Buildx is not available"
        log_info "Please install Docker Buildx or use Docker Desktop"
        exit 1
    fi
    
    log_success "✅ Docker Buildx is available"
    log_info "🚀 CPU images will be built for AMD64 and ARM64"
    log_info "🎮 GPU images will be built for AMD64 only"
}

# Determine build platforms based on image type
get_build_platforms() {
    local image_name="$1"
    
    # GPU images: AMD64 only (GPU dependencies)
    if [[ "$image_name" == *"gpu"* ]]; then
        echo "linux/amd64"
    else
        # CPU images: Multi-platform
        echo "linux/amd64,linux/arm64"
    fi
}

# Login to GHCR
login_ghcr() {
    log_step "Logging in to GitHub Container Registry..."
    echo "$GITHUB_TOKEN" | docker login ghcr.io -u "$GITHUB_USERNAME" --password-stdin
    log_success "Successfully logged in to GHCR"
}

# Get current version from git or use provided version
get_version() {
    if [ -n "$1" ]; then
        echo "$1"
    elif git rev-parse --git-dir > /dev/null 2>&1; then
        # Use git commit short hash as version
        git rev-parse --short HEAD
    else
        # Fallback to timestamp
        date +%Y%m%d-%H%M%S
    fi
}

# Build and push image using buildx
build_and_push() {
    local context_dir="$1"
    local dockerfile="$2"
    local image_name="$3"
    local version="$4"
    local platform_info="$5"
    local step_num="$6"
    local total_steps="$7"
    
    local full_image_name="${BASE_IMAGE_NAME}/${image_name}"
    local versioned_tag="${full_image_name}:${version}"
    local latest_tag="${full_image_name}:latest"
    local release_tag="${full_image_name}:${RELEASE_VERSION}"
    
    # Get actual build platforms based on image type
    local build_platforms=$(get_build_platforms "$image_name")
    local platform_display=""
    
    if [[ "$build_platforms" == *","* ]]; then
        platform_display="Multi-platform ($(echo $build_platforms | tr ',' '/'))"
    else
        platform_display="Single-platform ($build_platforms)"
    fi
    
    echo ""
    echo "========================================"
    log_header "[$step_num/$total_steps] Building ${image_name}"
    echo "========================================"
    log_info "Context: ${context_dir}"
    log_info "Dockerfile: ${dockerfile}"
    log_info "Platforms: ${platform_display}"
    log_info "Tags: ${version}, ${RELEASE_VERSION}, latest"
    echo "========================================"
    
    # Check if Dockerfile exists
    if [ ! -f "${context_dir}/${dockerfile}" ]; then
        log_error "Dockerfile not found: ${context_dir}/${dockerfile}"
        return 1
    fi
    
    # Build and push using buildx
    log_step "Starting Docker buildx build and push..."
    local build_start_time=$(date +%s)
    
    # Build and push using buildx
    if docker buildx build \
        --platform "${build_platforms}" \
        --progress=plain \
        -f "${context_dir}/${dockerfile}" \
        -t "${versioned_tag}" \
        -t "${release_tag}" \
        -t "${latest_tag}" \
        --push \
        "${context_dir}"; then
        
        local build_end_time=$(date +%s)
        local build_duration=$((build_end_time - build_start_time))
        
        log_success "✅ Build and push completed in ${build_duration}s"
        log_success "Pushed: ${versioned_tag}"
        log_success "Pushed: ${release_tag}"
        log_success "Pushed: ${latest_tag}"
        log_info "Platforms: ${build_platforms}"
        
        return 0
    else
        log_error "❌ Build failed for ${image_name}"
        return 1
    fi
}

# Main execution
main() {
    echo "================================================"
    log_header "        CellCraft GHCR Upload Script"
    echo "================================================"
    
    # Parse command line arguments
    VERSION=""
    DRY_RUN=false
    IMAGES_TO_BUILD="all"
    SHOW_HELP=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            --version)
                VERSION="$2"
                shift 2
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --images)
                IMAGES_TO_BUILD="$2"
                shift 2
                ;;
            --help|-h)
                SHOW_HELP=true
                shift
                ;;
            *)
                log_error "Unknown option: $1"
                SHOW_HELP=true
                break
                ;;
        esac
    done
    
    if [ "$SHOW_HELP" = true ]; then
        echo ""
        show_image_details
        echo "Usage: $0 [options]"
        echo ""
        echo "Options:"
        echo "  --version <version>    Specify version tag (default: git short hash or timestamp)"
        echo "  --dry-run              Show what would be built without actually building"
        echo "  --images <list>        Comma-separated list of images to build:"
        echo "                         Available: frontend,backend-cpu,backend-gpu,celery-cpu,celery-gpu"
        echo "                         Use 'all' to build all images (default)"
        echo "  --help, -h             Show this help message"
        echo ""
        echo "Environment variables:"
        echo "  GITHUB_TOKEN           GitHub Personal Access Token (required)"
        echo "                         Needs 'write:packages' and 'read:packages' permissions"
        echo "  GITHUB_USERNAME        GitHub username (default: cxinsys)"
        echo "  REPO_NAME              Repository name (default: cellcraft)"
        echo "  RELEASE_VERSION        Extra semver tag to add (default: v1.0.0)"
        echo ""
        echo "Examples:"
        echo "  # Build and upload all images with auto version"
        echo "  export GITHUB_TOKEN=your_token_here"
        echo "  $0"
        echo ""
        echo "  # Build specific images with custom version"
        echo "  $0 --version v1.0.0 --images frontend,backend-cpu"
        echo ""
        echo "  # Dry run to see what would be built"
        echo "  $0 --dry-run"
        echo ""
        exit 0
    fi
    
    # Show image details
    show_image_details
    
    # Get version
    if [ -z "$VERSION" ]; then
        VERSION=$(get_version)
        log_info "Using auto-generated version: $VERSION"
    else
        log_info "Using specified version: $VERSION"
    fi
    
    # Check if this is a dry run
    if [ "$DRY_RUN" = true ]; then
        log_warning "🔍 DRY RUN MODE - No images will be built or pushed"
        echo
    fi
    
    # Define images to build (context:dockerfile:imagename:platform)
    declare -A IMAGES
    IMAGES["frontend"]="frontend:Dockerfile.local:frontend:all"
    IMAGES["backend-cpu"]="backend:Dockerfile.improved:backend-cpu:cpu"
    IMAGES["backend-gpu"]="backend:Dockerfile:backend-gpu:gpu"
    IMAGES["celery-cpu"]="backend:Dockerfile.improved.celery:celery-cpu:cpu"
    IMAGES["celery-gpu"]="backend:Dockerfile.celery:celery-gpu:gpu"
    
    # Parse images to build
    if [ "$IMAGES_TO_BUILD" = "all" ]; then
        IMAGES_LIST="frontend,backend-cpu,backend-gpu,celery-cpu,celery-gpu"
    else
        IMAGES_LIST="$IMAGES_TO_BUILD"
    fi
    
    # Convert comma-separated list to array
    IFS=',' read -ra IMAGE_ARRAY <<< "$IMAGES_LIST"
    
    # Validate image names
    for image_key in "${IMAGE_ARRAY[@]}"; do
        image_key=$(echo "$image_key" | xargs) # Trim whitespace
        if [ -z "${IMAGES[$image_key]}" ]; then
            log_error "Unknown image: ${image_key}"
            log_info "Available images: ${!IMAGES[*]}"
            exit 1
        fi
    done
    
    # Check prerequisites and login
    if [ "$DRY_RUN" = false ]; then
        check_prerequisites
        check_buildx
        login_ghcr
    else
        log_info "Dry run mode: Skipping buildx check and GHCR login"
    fi
    
    echo
    echo "================================================"
    log_info "📋 Build Configuration"
    echo "================================================"
    log_info "Images to build: ${IMAGES_LIST}"
    log_info "Version: ${VERSION}"
    log_info "Release tag: ${RELEASE_VERSION}"
    log_info "Registry: ${BASE_IMAGE_NAME}"
    log_info "Total images: ${#IMAGE_ARRAY[@]}"
    echo "================================================"
    
    # Build and push each image
    local success_count=0
    local total_count=${#IMAGE_ARRAY[@]}
    local start_time=$(date +%s)
    local step_num=1
    
    for image_key in "${IMAGE_ARRAY[@]}"; do
        image_key=$(echo "$image_key" | xargs) # Trim whitespace
        
        IFS=':' read -ra IMAGE_INFO <<< "${IMAGES[$image_key]}"
        context_dir="${IMAGE_INFO[0]}"
        dockerfile="${IMAGE_INFO[1]}"
        image_name="${IMAGE_INFO[2]}"
        platform="${IMAGE_INFO[3]}"
        
        if [ "$DRY_RUN" = true ]; then
            echo ""
            log_info "[🔍 DRY RUN $step_num/$total_count] Would build: ${BASE_IMAGE_NAME}/${image_name}:${VERSION}"
            log_info "Context: ${context_dir}, Dockerfile: ${dockerfile}, Platform: ${platform}"
            log_success "✅ Dry run completed for ${image_name}"
            success_count=$((success_count + 1))
        else
            if build_and_push "$context_dir" "$dockerfile" "$image_name" "$VERSION" "$platform" "$step_num" "$total_count"; then
                success_count=$((success_count + 1))
            else
                log_error "❌ Failed to build ${image_name}"
            fi
        fi
        
        ((step_num++))
    done
    
    local end_time=$(date +%s)
    local total_duration=$((end_time - start_time))
    
    echo ""
    echo "================================================"
    log_header "📊 BUILD SUMMARY"
    echo "================================================"
    log_info "Successfully processed: ${success_count}/${total_count} images"
    log_info "Total time: ${total_duration}s"
    
    if [ "$DRY_RUN" = false ] && [ $success_count -eq $total_count ]; then
        log_success "🎉 All images built and pushed successfully!"
        echo ""
        log_info "📦 Images available at:"
        for image_key in "${IMAGE_ARRAY[@]}"; do
            image_key=$(echo "$image_key" | xargs)
            IFS=':' read -ra IMAGE_INFO <<< "${IMAGES[$image_key]}"
            image_name="${IMAGE_INFO[2]}"
            echo "   • ${BASE_IMAGE_NAME}/${image_name}:${VERSION}"
            echo "   • ${BASE_IMAGE_NAME}/${image_name}:${RELEASE_VERSION}"
            echo "   • ${BASE_IMAGE_NAME}/${image_name}:latest"
            echo ""
        done
        
        echo "================================================"
        log_header "🚀 NEXT STEPS:"
        echo "================================================"
        log_info "1. Update docker-compose files to use GHCR images"
        log_info "2. Replace 'build' sections with 'image' references"
        log_info "3. Test image pulling:"
        echo "   docker compose -f docker-compose.cpu.yml pull"
        echo "   docker compose -f docker-compose.gpu.yml pull"
        echo ""
        log_info "4. Start services with:"
        echo "   docker compose -f docker-compose.cpu.yml up -d"
        echo "   docker compose -f docker-compose.gpu.yml up -d"
        echo "================================================"
        
    elif [ "$DRY_RUN" = false ]; then
        log_error "❌ Some images failed to build. Please check the logs above."
        exit 1
    elif [ "$DRY_RUN" = true ]; then
        log_success "🔍 Dry run completed successfully!"
        log_info "All ${total_count} images would be built and uploaded."
    fi
}

# Run main function with all arguments
main "$@"
