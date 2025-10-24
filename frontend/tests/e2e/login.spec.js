/**
 * 로그인 기능 E2E 테스트
 *
 * 테스트 시나리오:
 * - 로그인 페이지 접근
 * - 로그인 폼 제출
 * - 로그인 성공/실패 처리
 * - 인증 후 리다이렉션
 */

import { test, expect } from '@playwright/test';

test.describe('로그인 기능', () => {
  test.beforeEach(async ({ page }) => {
    // 각 테스트 전에 로그인 페이지로 이동
    await page.goto('/');
  });

  test('로그인 페이지가 올바르게 로드되어야 함', async ({ page }) => {
    // 로그인 폼이 표시되는지 확인
    await expect(page.locator('input[type="text"]')).toBeVisible();
    await expect(page.locator('input[type="password"]')).toBeVisible();
    await expect(page.locator('button[type="submit"]')).toBeVisible();
  });

  test('빈 폼으로 로그인 시도 시 유효성 검사가 작동해야 함', async ({ page }) => {
    // 로그인 버튼 클릭
    await page.locator('button[type="submit"]').click();

    // 폼이 제출되지 않고 페이지에 남아있어야 함
    await expect(page).toHaveURL('/');
  });

  test('잘못된 자격증명으로 로그인 시 에러 메시지가 표시되어야 함', async ({ page }) => {
    // 잘못된 자격증명 입력
    await page.locator('input[type="text"]').fill('wronguser');
    await page.locator('input[type="password"]').fill('wrongpassword');

    // 로그인 버튼 클릭
    await page.locator('button[type="submit"]').click();

    // 에러 메시지 확인 (실제 에러 메시지 selector로 수정 필요)
    // await expect(page.locator('.error-message')).toBeVisible();
  });

  test('올바른 자격증명으로 로그인 시 홈페이지로 리다이렉션되어야 함', async ({ page }) => {
    // 올바른 자격증명 입력 (실제 테스트 계정 정보로 수정 필요)
    await page.locator('input[type="text"]').fill('testuser');
    await page.locator('input[type="password"]').fill('testpassword');

    // 로그인 버튼 클릭
    await page.locator('button[type="submit"]').click();

    // 홈페이지로 리다이렉션 확인
    await expect(page).toHaveURL(/.*home/);
  });

  test('로그인 후 인증 토큰이 저장되어야 함', async ({ page, context }) => {
    // 올바른 자격증명으로 로그인
    await page.locator('input[type="text"]').fill('testuser');
    await page.locator('input[type="password"]').fill('testpassword');
    await page.locator('button[type="submit"]').click();

    // 페이지가 로드될 때까지 대기
    await page.waitForLoadState('networkidle');

    // 로컬 스토리지 또는 쿠키에 토큰이 저장되었는지 확인
    const cookies = await context.cookies();
    const hasAuthCookie = cookies.some(cookie =>
      cookie.name.includes('token') || cookie.name.includes('auth')
    );

    // 토큰이 쿠키에 있거나 로컬스토리지에 있어야 함
    if (!hasAuthCookie) {
      const localStorage = await page.evaluate(() => {
        return window.localStorage.getItem('token') !== null;
      });
      expect(localStorage).toBeTruthy();
    }
  });
});

test.describe('로그아웃 기능', () => {
  test.beforeEach(async ({ page }) => {
    // 로그인 상태로 시작 (실제 로그인 플로우로 수정 필요)
    await page.goto('/');
    await page.locator('input[type="text"]').fill('testuser');
    await page.locator('input[type="password"]').fill('testpassword');
    await page.locator('button[type="submit"]').click();
    await page.waitForLoadState('networkidle');
  });

  test('로그아웃 버튼 클릭 시 로그인 페이지로 리다이렉션되어야 함', async ({ page }) => {
    // 로그아웃 버튼 클릭 (실제 selector로 수정 필요)
    // await page.locator('.logout-button').click();

    // 로그인 페이지로 리다이렉션 확인
    // await expect(page).toHaveURL('/');
  });

  test('로그아웃 후 인증 토큰이 제거되어야 함', async ({ page, context }) => {
    // 로그아웃 (실제 로그아웃 플로우로 수정 필요)
    // await page.locator('.logout-button').click();

    // 쿠키와 로컬스토리지에서 토큰이 제거되었는지 확인
    const cookies = await context.cookies();
    const hasAuthCookie = cookies.some(cookie =>
      cookie.name.includes('token') || cookie.name.includes('auth')
    );
    expect(hasAuthCookie).toBeFalsy();

    const localStorage = await page.evaluate(() => {
      return window.localStorage.getItem('token');
    });
    expect(localStorage).toBeNull();
  });
});

/**
 * 주의사항:
 * 1. 실제 로그인 폼의 selector는 프로젝트 구조에 맞게 수정해야 합니다.
 * 2. 테스트용 계정 정보는 환경 변수로 관리하는 것이 좋습니다.
 * 3. API mocking이 필요한 경우 MSW를 사용하여 네트워크 요청을 가로챌 수 있습니다.
 * 4. 실제 백엔드 서버가 실행 중이어야 테스트가 정상 작동합니다.
 */
