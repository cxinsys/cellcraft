#!/bin/bash
set -e

# Playwright 브라우저 자동 설치/업데이트 스크립트
# 컨테이너 시작 시 브라우저가 없거나 버전이 맞지 않으면 자동 설치

BROWSERS_PATH="${PLAYWRIGHT_BROWSERS_PATH:-/opt/playwright-browsers}"

echo "========================================"
echo "Checking Playwright browsers..."
echo "Browsers path: $BROWSERS_PATH"
echo "========================================"

# 브라우저 디렉토리가 비어있거나 chromium이 없으면 설치
if [ ! -d "$BROWSERS_PATH" ] || [ -z "$(ls -A $BROWSERS_PATH 2>/dev/null)" ] || [ -z "$(ls -d $BROWSERS_PATH/chromium-* 2>/dev/null)" ]; then
    echo "Playwright browsers not found. Installing..."
    mkdir -p "$BROWSERS_PATH"
    npx playwright install chromium
    chmod -R 755 "$BROWSERS_PATH"
    echo "Playwright browsers installed successfully."
else
    echo "Playwright browsers found:"
    ls -d "$BROWSERS_PATH"/*/ 2>/dev/null | xargs -I {} basename {}
fi

echo "========================================"
echo "Starting Jenkins..."
echo "========================================"

# Jenkins 기본 엔트리포인트 실행 (tini를 통해)
exec /usr/bin/tini -- /usr/local/bin/jenkins.sh "$@"
