# Priority 1 Test 실행 가이드

## ✅ 개선 사항

### 로그인 최적화 완료

**Before (비효율적):**
```javascript
// 매 테스트마다 로그인 (3개 테스트 = 3번 로그인)
test.beforeEach(async ({ page }) => {
  await loginPage.goto();
  await loginPage.login(testUser.email, testUser.password); // 느림!
  await loginPage.verifyLoginSuccess();
});
```

**After (효율적 ✅):**
```javascript
// auth fixture 사용 - Worker당 1번만 로그인, 세션 재사용
import { test, expect } from './fixtures/auth.js';

// beforeEach에서 로그인 로직 제거
test.beforeEach(async ({ page }) => {
  // 이미 인증된 상태로 시작
  await projectsPage.goto();
  await projectsPage.verifyPageLoaded();
});
```

**성능 향상:**
- 3개 Priority 1 테스트 기준
- Before: ~15-30초 로그인 시간
- After: ~5-10초 로그인 시간 (첫 실행 시 1번만)
- **약 60-70% 실행 시간 단축**

---

## 🚀 Priority 1 테스트 실행 명령어

### 1. 기본 실행 (Headless Chromium)

```bash
cd frontend

# Priority 1 테스트만 실행 (headless, chromium)
npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 1"
```

**특징:**
- ✅ Headless 모드 (백그라운드 실행)
- ✅ Chromium 브라우저만
- ✅ Priority 1 테스트만 실행 (3개 테스트)
- ✅ 가장 빠른 실행 속도

---

### 2. Headed 모드 (디버깅용)

```bash
# 브라우저를 볼 수 있는 모드
PLAYWRIGHT_HEADLESS=false npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 1"
```

**용도:**
- 테스트 실행 과정을 눈으로 확인
- 실패 원인 파악
- 셀렉터 검증

---

### 3. UI 모드 (개발/디버깅 추천)

```bash
# 대화형 UI로 테스트 실행
npm run test:e2e -- workflow-configuration.spec.js --project=chromium --ui
```

**특징:**
- ✅ 테스트 단계별 실행
- ✅ 타임 트래블 (과거 단계로 이동)
- ✅ DOM 스냅샷 확인
- ✅ 네트워크 요청 확인
- **개발 시 가장 편리한 방법**

---

### 4. 전체 Workflow Configuration 테스트

```bash
# Priority 1, 2, 3, 4, 5 모두 실행 (구현된 것만)
npm run test:e2e -- workflow-configuration.spec.js --project=chromium
```

---

### 5. Docker Compose 환경에서 실행

```bash
# 1. Docker Compose로 백엔드 실행
docker compose -f docker-compose.dev.yml up -d

# 2. 서비스가 준비될 때까지 대기
# (http://localhost:8080 접속 확인)

# 3. 테스트 실행 (webServer 자동 시작 스킵)
PLAYWRIGHT_SKIP_WEBSERVER=true npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 1"
```

---

## 🔍 인증 상태 관리

### Auth Fixture 동작 방식

1. **첫 실행 시:**
   ```
   [Worker 1] Authenticating user...
   [Worker 1] Login to http://localhost:8080/login
   [Worker 1] Authentication state saved to .auth/test-user.json
   ```

2. **이후 실행 시:**
   ```
   [Worker 1] Reusing existing authentication state
   [테스트 시작 - 이미 로그인된 상태]
   ```

3. **세션 파일 위치:**
   ```
   frontend/tests/e2e/.auth/test-user.json (gitignored)
   ```

---

## 🧪 테스트 결과 확인

### 1. 실행 중 로그

```bash
Running 3 tests using 1 worker

  ✓ [chromium] › workflow-configuration.spec.js:89:3 › Priority 1: Should assign file to InputFile node (15s)
  ✓ [chromium] › workflow-configuration.spec.js:180:3 › Priority 1: Should change assigned file to different file (12s)
  ✓ [chromium] › workflow-configuration.spec.js:220:3 › Priority 1: Should handle empty folder gracefully (8s)

  3 passed (35s)
```

### 2. HTML 리포트 보기

```bash
# 테스트 실행 후
npx playwright show-report
```

### 3. 실패 시 디버깅

```bash
# 실패한 테스트만 재실행
npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 1" --last-failed

# 트레이스 보기
npx playwright show-trace test-results/*/trace.zip
```

---

## 🎯 테스트 케이스 상세

### Priority 1-1: Should assign file to InputFile node

**검증 항목:**
- ✅ InputFile 노드가 캔버스에 존재
- ✅ 모달이 올바르게 열림
- ✅ 폴더 목록 로드
- ✅ 폴더 선택 시 하이라이트
- ✅ 파일 목록 로드
- ✅ 파일 선택 시 하이라이트
- ✅ 현재 파일 정보 표시
- ✅ Apply 버튼 상태 변경 ("Apply" → "Applied")
- ✅ 파일 할당 Vuex 전파
- ✅ 모달 재오픈 시 할당 유지 (Persistence)

### Priority 1-2: Should change assigned file to different file

**검증 항목:**
- ✅ 초기 파일 할당
- ✅ 다른 파일로 변경
- ✅ 변경된 파일 persistence

### Priority 1-3: Should handle empty folder gracefully

**검증 항목:**
- ✅ 빈 폴더 선택
- ✅ 적절한 빈 상태 처리
- ✅ "Please add data file" 메시지

---

## 🔧 환경 변수

```bash
# 기본 설정
PLAYWRIGHT_BASE_URL=http://localhost:8080
PLAYWRIGHT_HEADLESS=true
PLAYWRIGHT_USER=test1234@test.com
PLAYWRIGHT_PASS=test1234

# 커스텀 설정 예시
PLAYWRIGHT_USER=custom@test.com \
PLAYWRIGHT_PASS=custompass \
npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 1"
```

---

## 📊 성능 벤치마크

### 예상 실행 시간

| 테스트 케이스 | 예상 시간 |
|--------------|----------|
| Priority 1-1 (main) | 12-18초 |
| Priority 1-2 (change file) | 10-15초 |
| Priority 1-3 (empty folder) | 6-10초 |
| **Total (3 tests)** | **30-45초** |

**참고:**
- 첫 실행: +5-10초 (인증 세션 생성)
- 이후 실행: 세션 재사용으로 빠름

---

## 🐛 트러블슈팅

### 1. "Target closed" 에러

**원인:** 백엔드가 실행되지 않음

**해결:**
```bash
# Docker Compose 실행 확인
docker compose -f docker-compose.dev.yml ps

# 재시작
docker compose -f docker-compose.dev.yml restart
```

---

### 2. 인증 실패

**원인:** 세션 파일 손상 또는 만료

**해결:**
```bash
# 인증 상태 파일 삭제 후 재실행
rm -f frontend/tests/e2e/.auth/test-user.json

# 다시 테스트 실행 (새로운 세션 생성됨)
npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 1"
```

---

### 3. 파일/폴더를 찾을 수 없음

**원인:** 테스트 데이터 부족

**해결:**
```bash
# 1. 로그인하여 PBMCLight1000 폴더 확인
# 2. PBMCLight1000.h5ad 파일이 있는지 확인
# 3. 없다면 Datasets 페이지에서 다운로드

# 또는 테스트 코드에서 testWorkflow.folder와 expectedFile 수정
```

---

### 4. Timeout 에러

**원인:** 네트워크 느림 또는 대용량 파일

**해결:**
```bash
# playwright.config.js에서 timeout 증가
# actionTimeout: 10000 → 20000
# navigationTimeout: 30000 → 60000
```

---

## 📝 다음 단계

### Priority 2 구현 준비

```bash
# Priority 2 테스트 구현 후
npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 2"

# Priority 1 + 2 함께 실행
npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority"
```

---

## 📚 참고 문서

- **전체 테스트 플랜:** `WORKFLOW_CONFIGURATION_TEST_PLAN.md`
- **E2E 가이드:** `README.md`
- **Page Object Models:** `pages/` 디렉토리
- **Helper Utilities:** `utils/` 디렉토리
- **Playwright 공식 문서:** https://playwright.dev

---

## 🎉 요약

**빠른 실행 명령어:**
```bash
# Headless Chromium (기본, 가장 빠름)
npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 1"

# UI 모드 (개발/디버깅)
npm run test:e2e -- workflow-configuration.spec.js --project=chromium --ui
```

**핵심 개선사항:**
- ✅ Auth fixture로 로그인 최적화 (60-70% 속도 향상)
- ✅ 세션 재사용으로 반복 실행 빠름
- ✅ 명확한 실행 명령어
- ✅ 3개 Priority 1 테스트 모두 구현 완료
