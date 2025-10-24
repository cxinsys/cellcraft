#!/bin/bash
# deploy-initial-version.sh

ORG_NAME=cxinsys
VERSION=1.0.0
MAJOR=1
MINOR=1.0

echo "🚀 Deploying CellCraft Plugins v${VERSION}"

# GPU 플러그인 (AMD64만 지원)
GPU_PLUGINS=("FastTENET" "FastSCODE")

# CPU 플러그인 (멀티 플랫폼 지원)
CPU_PLUGINS=("TENET" "Scribe" "LEAP" "GRNBoost2" "GENIE3" "GRNViz")

# GPU 플러그인 빌드
for plugin in "${GPU_PLUGINS[@]}"; do
  if [[ -f "${plugin}/Dockerfile" ]]; then
    echo "🎮 Building GPU plugin ${plugin} v${VERSION} (AMD64 only)..."
    plugin_lower=$(echo "$plugin" | tr '[:upper:]' '[:lower:]')
    
    docker buildx build \
      --platform linux/amd64 \
      -t ghcr.io/${ORG_NAME}/${plugin_lower}:${VERSION} \
      -t ghcr.io/${ORG_NAME}/${plugin_lower}:${MINOR} \
      -t ghcr.io/${ORG_NAME}/${plugin_lower}:${MAJOR} \
      -t ghcr.io/${ORG_NAME}/${plugin_lower}:latest \
      --push \
      "${plugin}/"
    
    echo "✅ Successfully built ${plugin}"
  else
    echo "⚠️  Skipping ${plugin}: Dockerfile not found"
  fi
done

# CPU 플러그인 빌드
for plugin in "${CPU_PLUGINS[@]}"; do
  if [[ -f "${plugin}/Dockerfile" ]]; then
    echo "📦 Building CPU plugin ${plugin} v${VERSION} (Multi-platform)..."
    plugin_lower=$(echo "$plugin" | tr '[:upper:]' '[:lower:]')
    
    docker buildx build \
      --platform linux/amd64,linux/arm64 \
      -t ghcr.io/${ORG_NAME}/${plugin_lower}:${VERSION} \
      -t ghcr.io/${ORG_NAME}/${plugin_lower}:${MINOR} \
      -t ghcr.io/${ORG_NAME}/${plugin_lower}:${MAJOR} \
      -t ghcr.io/${ORG_NAME}/${plugin_lower}:latest \
      --push \
      "${plugin}/"
    
    echo "✅ Successfully built ${plugin}"
  else
    echo "⚠️  Skipping ${plugin}: Dockerfile not found"
  fi
done

echo "🎉 All plugins deployed successfully!"