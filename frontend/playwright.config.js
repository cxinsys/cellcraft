import { defineConfig, devices } from '@playwright/test';

/**
 * Playwright E2E 테스트 설정
 * @see https://playwright.dev/docs/test-configuration
 */
export default defineConfig({
  testDir: './tests/e2e',

  /* 테스트 파일 패턴 */
  testMatch: '**/*.spec.js',

  /* 병렬 실행 설정 */
  fullyParallel: true,

  /* CI 환경에서 재시도 설정 */
  retries: process.env.CI ? 2 : 0,

  /* 워커 수 설정 */
  workers: process.env.CI ? 1 : undefined,

  /* 리포터 설정 */
  reporter: [
    ['html', { outputFolder: 'playwright-report' }],
    ['json', { outputFile: 'test-results/results.json' }],
    ['list']
  ],

  /* 공통 설정 */
  use: {
    /* 베이스 URL - 환경변수로 오버라이드 가능 */
    baseURL: process.env.PLAYWRIGHT_BASE_URL || 'http://localhost:8080',

    /* 스크린샷 설정 */
    screenshot: 'only-on-failure',

    /* 비디오 설정 */
    video: 'retain-on-failure',

    /* 트레이스 설정 */
    trace: 'on-first-retry',

    /* 액션 타임아웃 */
    actionTimeout: 10000,

    /* 네비게이션 타임아웃 */
    navigationTimeout: 30000,

    /* 헤드리스 모드 설정 - 환경변수로 제어 */
    headless: process.env.PLAYWRIGHT_HEADLESS !== 'false',
  },

  /* 프로젝트별 설정 */
  projects: [
    {
      name: 'chromium',
      use: {
        ...devices['Desktop Chrome'],
        viewport: { width: 1920, height: 1080 }
      },
    },

    {
      name: 'firefox',
      use: {
        ...devices['Desktop Firefox'],
        viewport: { width: 1920, height: 1080 }
      },
    },

    /* WebKit (Safari) - Linux 환경에서 시스템 의존성 필요 */
    /* macOS 환경이나 CI/CD의 macOS runner에서만 활성화 권장 */
    // {
    //   name: 'webkit',
    //   use: {
    //     ...devices['Desktop Safari'],
    //     viewport: { width: 1920, height: 1080 }
    //   },
    // },

    /* 모바일 테스트 (선택사항) */
    // {
    //   name: 'Mobile Chrome',
    //   use: { ...devices['Pixel 5'] },
    // },
    // {
    //   name: 'Mobile Safari',
    //   use: { ...devices['iPhone 12'] },
    // },
  ],

  /* 웹서버 설정 (테스트 실행 시 자동으로 개발 서버 시작) */
  /*
   * PLAYWRIGHT_SKIP_WEBSERVER=true로 설정하면 webServer를 건너뛸 수 있습니다.
   * Docker Compose로 이미 서버가 실행 중인 경우 유용합니다.
   */
  webServer: process.env.PLAYWRIGHT_SKIP_WEBSERVER
    ? undefined
    : {
        command: 'npm run serve',
        url: process.env.PLAYWRIGHT_BASE_URL || 'http://localhost:8080',
        reuseExistingServer: !process.env.CI,
        timeout: 120000,
        stdout: 'ignore',
        stderr: 'pipe',
      },

  /*
   * 환경변수 설명:
   * - PLAYWRIGHT_BASE_URL: 테스트할 애플리케이션의 URL (기본값: http://localhost:8080)
   * - PLAYWRIGHT_HEADLESS: 헤드리스 모드 활성화/비활성화 (기본값: true)
   * - PLAYWRIGHT_USER: 테스트 사용자 이메일 (기본값: test1234@test.com)
   * - PLAYWRIGHT_PASS: 테스트 사용자 비밀번호 (기본값: test1234)
   * - PLAYWRIGHT_SKIP_WEBSERVER: webServer 자동 시작 건너뛰기 (기본값: false)
   *
   * 사용 예:
   * PLAYWRIGHT_BASE_URL=http://localhost:8081 npm run test:e2e
   * PLAYWRIGHT_HEADLESS=false npm run test:e2e -- --project=chromium
   * PLAYWRIGHT_SKIP_WEBSERVER=true npm run test:e2e
   */
});
