#!/bin/bash

# E2E 테스트 빠른 실행 스크립트
# 사용법: ./run-test.sh [테스트번호] [옵션]

set -e

# 색상 코드
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# 기본 설정
DEFAULT_URL="https://cellcraft.app"
DEFAULT_HEADLESS="false"

# 사용법 출력
usage() {
    echo "════════════════════════════════════════"
    echo "E2E 테스트 실행 스크립트"
    echo "════════════════════════════════════════"
    echo ""
    echo "사용법: ./run-test.sh [테스트번호] [옵션]"
    echo ""
    echo "테스트 번호:"
    echo "  01, 1    - InputFile Node Test"
    echo "  02, 2    - DataTable Node Test"
    echo "  03, 3    - ScatterPlot Node Test"
    echo "  04, 4    - Algorithm Node Test"
    echo "  05, 5    - Workflow Execution Test (긴 테스트)"
    echo "  06, 6    - Result Visualization Test"
    echo "  all      - 모든 테스트 실행"
    echo ""
    echo "옵션:"
    echo "  --headless       헤드리스 모드로 실행"
    echo "  --headed         브라우저 보면서 실행 (기본값)"
    echo "  --local          로컬 서버 사용 (http://localhost:8080)"
    echo "  --prod           프로덕션 서버 사용 (https://cellcraft.app, 기본값)"
    echo "  --debug          디버그 모드 (PWDEBUG=1)"
    echo "  --timeout=N      타임아웃 설정 (밀리초, 기본값: 자동)"
    echo ""
    echo "예시:"
    echo "  ./run-test.sh 5                  # Test 05를 브라우저 보면서 실행"
    echo "  ./run-test.sh 5 --headless       # Test 05를 헤드리스로 실행"
    echo "  ./run-test.sh 1 --local          # Test 01을 로컬 서버로 실행"
    echo "  ./run-test.sh all --headless     # 모든 테스트 헤드리스 실행"
    echo "  ./run-test.sh 5 --debug          # Test 05 디버그 모드"
    echo ""
    exit 1
}

# 인자 파싱
TEST_NUM=""
BASE_URL="$DEFAULT_URL"
HEADLESS="$DEFAULT_HEADLESS"
DEBUG_MODE=""
TIMEOUT=""

if [ $# -eq 0 ]; then
    usage
fi

TEST_NUM=$1
shift

while [ $# -gt 0 ]; do
    case "$1" in
        --headless)
            HEADLESS="true"
            ;;
        --headed)
            HEADLESS="false"
            ;;
        --local)
            BASE_URL="http://localhost:8080"
            ;;
        --prod)
            BASE_URL="https://cellcraft.app"
            ;;
        --debug)
            DEBUG_MODE="PWDEBUG=1"
            ;;
        --timeout=*)
            TIMEOUT="--timeout=${1#*=}"
            ;;
        *)
            echo -e "${RED}알 수 없는 옵션: $1${NC}"
            usage
            ;;
    esac
    shift
done

# 테스트 파일 매핑
get_test_file() {
    case "$1" in
        01|1)
            echo "tests/e2e/workflows/01-file-assignment.spec.js"
            echo "InputFile Node Test"
            ;;
        02|2)
            echo "tests/e2e/workflows/02-data-display.spec.js"
            echo "DataTable Node Test"
            ;;
        03|3)
            echo "tests/e2e/workflows/03-scatter-plot.spec.js"
            echo "ScatterPlot Node Test"
            ;;
        04|4)
            echo "tests/e2e/workflows/04-algorithm-config.spec.js"
            echo "Algorithm Node Test"
            ;;
        05|5)
            echo "tests/e2e/workflows/05-workflow-execution.spec.js"
            echo "Workflow Execution Test"
            if [ -z "$TIMEOUT" ]; then
                TIMEOUT="--timeout=600000"
            fi
            ;;
        06|6)
            echo "tests/e2e/workflows/06-result-visualization.spec.js"
            echo "Result Visualization Test"
            ;;
        all)
            echo "tests/e2e/workflows/"
            echo "전체 테스트"
            ;;
        *)
            echo ""
            echo "알 수 없는 테스트 번호: $1"
            ;;
    esac
}

# 테스트 정보 가져오기
TEST_INFO=$(get_test_file "$TEST_NUM")
TEST_FILE=$(echo "$TEST_INFO" | head -1)
TEST_NAME=$(echo "$TEST_INFO" | tail -1)

if [ -z "$TEST_FILE" ]; then
    usage
fi

# 실행 정보 출력
echo "════════════════════════════════════════"
echo -e "${BLUE}E2E 테스트 실행${NC}"
echo "════════════════════════════════════════"
echo -e "테스트: ${GREEN}${TEST_NAME}${NC}"
echo "파일: ${TEST_FILE}"
echo "URL: ${BASE_URL}"
echo "헤드리스: $([ "$HEADLESS" = "true" ] && echo "예" || echo "아니오 (브라우저 표시)")"
echo "타임아웃: ${TIMEOUT:-기본값}"
[ -n "$DEBUG_MODE" ] && echo -e "${YELLOW}디버그 모드 활성화${NC}"
echo "════════════════════════════════════════"
echo ""

# 환경 변수 설정
export PLAYWRIGHT_BASE_URL="$BASE_URL"
export PLAYWRIGHT_SKIP_WEBSERVER=true
export PLAYWRIGHT_HEADLESS="$HEADLESS"

# 디버그 모드 설정
if [ -n "$DEBUG_MODE" ]; then
    export PWDEBUG=1
fi

# 테스트 실행
echo -e "${BLUE}테스트 시작...${NC}"
echo ""

if $DEBUG_MODE npx playwright test "$TEST_FILE" --project=chromium $TIMEOUT; then
    echo ""
    echo "════════════════════════════════════════"
    echo -e "${GREEN}✅ 테스트 성공!${NC}"
    echo "════════════════════════════════════════"
    exit 0
else
    EXIT_CODE=$?
    echo ""
    echo "════════════════════════════════════════"
    echo -e "${RED}❌ 테스트 실패!${NC}"
    echo "════════════════════════════════════════"
    echo ""
    echo "결과 확인:"
    echo "  - 스크린샷: test-results/"
    echo "  - HTML 리포트: npx playwright show-report"
    echo ""
    exit $EXIT_CODE
fi
