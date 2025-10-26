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
    /* 베이스 URL */
    baseURL: 'http://localhost:8080',

    /* 스크린샷 설정 */
    screenshot: 'only-on-failure',

    /* 비디오 설정 */
    video: 'retain-on-failure',

    /* 트레이스 설정 */
    trace: 'on-first-retry',

    /* 액션 타임아웃 */
    actionTimeout: 10000,

    /* 네비게이션 타임아웃 */
    navigationTimeout: 30000
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

    {
      name: 'webkit',
      use: {
        ...devices['Desktop Safari'],
        viewport: { width: 1920, height: 1080 }
      },
    },

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
  webServer: {
    command: 'npm run serve',
    url: 'http://localhost:8080',
    reuseExistingServer: !process.env.CI,
    timeout: 120000,
    stdout: 'ignore',
    stderr: 'pipe'
  },
});
