#!/bin/bash

echo "========================================="
echo "Starting E2E Tests..."
echo "Base URL: $PLAYWRIGHT_BASE_URL"
echo "========================================="

# Jenkins 컨테이너 내부에서 작업 디렉토리 생성
E2E_WORK_DIR="/tmp/e2e-tests-$$"
mkdir -p "$E2E_WORK_DIR"

# 필요한 파일들 복사
cp -r /workspace/frontend/tests "$E2E_WORK_DIR/tests"
cp /workspace/frontend/playwright.config.js "$E2E_WORK_DIR/"
cp /workspace/frontend/package.json "$E2E_WORK_DIR/"

cd "$E2E_WORK_DIR"

# Playwright 설치 (첫 실행 시만 필요, 이후에는 캐시됨)
if [ ! -d "node_modules/@playwright" ]; then
  npm install @playwright/test
fi

# E2E 테스트 실행
export PLAYWRIGHT_BASE_URL="${PLAYWRIGHT_BASE_URL:-http://cellcraft-frontend-1:8080}"
export PLAYWRIGHT_SKIP_WEBSERVER='true'
export PLAYWRIGHT_HEADLESS='true'

npx playwright test \
  --reporter=junit \
  --reporter=html \
  --project=chromium

TEST_EXIT_CODE=$?

# 테스트 결과를 Jenkins workspace로 복사
mkdir -p "$WORKSPACE/test-results"
cp -r test-results/* "$WORKSPACE/test-results/" 2>/dev/null || true
cp -r playwright-report "$WORKSPACE/" 2>/dev/null || true

# 임시 디렉토리 정리
cd /
rm -rf "$E2E_WORK_DIR"

# Uptime Kuma 알림
if [ -n "$UPTIME_KUMA_PUSH_URL" ] && [ "$UPTIME_KUMA_PUSH_URL" != "http://uptime-kuma:3001/api/push/YOUR_KEY_HERE" ]; then
  if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo "========================================="
    echo "E2E Tests Completed Successfully!"
    echo "========================================="
    curl -s "${UPTIME_KUMA_PUSH_URL}?status=up&msg=E2E%20Tests%20Passed" || true
  else
    echo "========================================="
    echo "E2E Tests Failed!"
    echo "========================================="
    curl -s "${UPTIME_KUMA_PUSH_URL}?status=down&msg=E2E%20Tests%20Failed" || true
  fi
fi

exit $TEST_EXIT_CODE
