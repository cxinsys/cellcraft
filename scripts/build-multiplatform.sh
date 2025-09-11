#!/bin/bash
set -e

# CellCraft Multi-Platform Docker Build Script
# Builds ARM64/AMD64 compatible images for GHCR

REGISTRY="ghcr.io"
REPO_OWNER="your-username"
PROJECT_NAME="cellcraft"

# Platform support
PLATFORMS="linux/amd64,linux/arm64"

# Image tags
BACKEND_TAG="${REGISTRY}/${REPO_OWNER}/${PROJECT_NAME}/backend:latest"
CELERY_TAG="${REGISTRY}/${REPO_OWNER}/${PROJECT_NAME}/celery:latest"

echo "🚀 Building CellCraft Multi-Platform Images"
echo "Registry: ${REGISTRY}"
echo "Platforms: ${PLATFORMS}"
echo "=========================================="

# Create buildx builder if not exists
if ! docker buildx inspect multiplatform-builder &> /dev/null; then
    echo "📦 Creating buildx builder..."
    docker buildx create --name multiplatform-builder --platform ${PLATFORMS}
fi

# Use the builder
docker buildx use multiplatform-builder

# Build backend image
echo "🔨 Building Backend Image..."
docker buildx build \
    --platform ${PLATFORMS} \
    --file backend/Dockerfile \
    --tag ${BACKEND_TAG} \
    --push \
    backend/

# Build celery image  
echo "🔨 Building Celery Image..."
docker buildx build \
    --platform ${PLATFORMS} \
    --file backend/Dockerfile.celery \
    --tag ${CELERY_TAG} \
    --push \
    backend/

echo "✅ Multi-platform images built successfully!"
echo "Backend: ${BACKEND_TAG}"
echo "Celery: ${CELERY_TAG}"
echo ""
echo "To use these images, update your docker-compose.yml:"
echo "  backend:"
echo "    image: ${BACKEND_TAG}"
echo "  celery:"
echo "    image: ${CELERY_TAG}"