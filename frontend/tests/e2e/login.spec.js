/**
 * 로그인 기능 E2E 테스트
 *
 * 테스트 시나리오:
 * - 로그인 페이지 로드 확인
 * - 유효성 검사
 * - 로그인 성공/실패 처리
 */

import { test, expect } from '@playwright/test';

test.describe('로그인 기능', () => {
  test.beforeEach(async ({ page }) => {
    // 로그인 페이지로 이동
    await page.goto('/login');
  });

  test('로그인 페이지가 올바르게 로드되어야 함', async ({ page }) => {
    // 로그인 폼이 표시되는지 확인
    await expect(page.getByPlaceholder('Email')).toBeVisible();
    await expect(page.getByPlaceholder('Password')).toBeVisible();
    await expect(page.getByRole('button', { name: /Sign In/i })).toBeVisible();
    await expect(page.locator('.header-text')).toHaveText('Sign in to Cellcraft');
  });

  test('빈 폼으로 로그인 시도 시 버튼이 비활성화되어야 함', async ({ page }) => {
    // 빈 상태에서 버튼이 비활성화되어 있는지 확인
    const signInButton = page.getByRole('button', { name: /Sign In/i });
    await expect(signInButton).toBeDisabled();
  });

  test('잘못된 자격증명으로 로그인 시 에러 메시지가 표시되어야 함', async ({
    page,
  }) => {
    // 잘못된 자격증명 입력
    await page.getByPlaceholder('Email').fill('wrong@test.com');
    await page.getByPlaceholder('Password').fill('wrongpassword');

    // 로그인 버튼 클릭
    await page.getByRole('button', { name: /Sign In/i }).click();

    // 에러 메시지 확인 (실제 백엔드 응답에 맞춤)
    await expect(page.locator('.login__error')).toBeVisible();
    // 에러 메시지는 백엔드 응답에 따라 다를 수 있음
    await expect(page.locator('.login__error')).toContainText(/Login failed|Invalid email or password/i);
  });

  test('올바른 자격증명으로 로그인 시 projects 페이지로 리다이렉션되어야 함', async ({
    page,
  }) => {
    const email = process.env.PLAYWRIGHT_USER ?? 'test1234@test.com';
    const password = process.env.PLAYWRIGHT_PASS ?? 'test1234';

    // 올바른 자격증명 입력
    await page.getByPlaceholder('Email').fill(email);
    await page.getByPlaceholder('Password').fill(password);

    // 로그인 버튼 클릭 및 리다이렉션 대기
    await Promise.all([
      page.waitForURL('**/projects', { timeout: 10000 }),
      page.getByRole('button', { name: /Sign In/i }).click(),
    ]);

    // Projects 페이지로 리다이렉션 확인
    await expect(page).toHaveURL(/.*\/projects/);

    // 로그인 후 헤더 메뉴 확인
    await expect(page.getByRole('link', { name: 'Workflows' })).toBeVisible();
    await expect(page.getByText('Sign Out')).toBeVisible();
  });

  test('로그인 후 인증 토큰이 저장되어야 함', async ({ page, context }) => {
    const email = process.env.PLAYWRIGHT_USER ?? 'test1234@test.com';
    const password = process.env.PLAYWRIGHT_PASS ?? 'test1234';

    // 로그인
    await page.getByPlaceholder('Email').fill(email);
    await page.getByPlaceholder('Password').fill(password);
    await page.getByRole('button', { name: /Sign In/i }).click();
    await page.waitForURL('**/projects');

    // 쿠키에 토큰이 저장되었는지 확인
    const cookies = await context.cookies();
    const hasAuthToken = cookies.some(
      (cookie) =>
        cookie.name === 'access_token' ||
        cookie.name === 'token' ||
        cookie.name.includes('auth')
    );

    expect(hasAuthToken).toBeTruthy();
  });
});

test.describe('로그아웃 기능', () => {
  test.beforeEach(async ({ page }) => {
    const email = process.env.PLAYWRIGHT_USER ?? 'test1234@test.com';
    const password = process.env.PLAYWRIGHT_PASS ?? 'test1234';

    // 로그인 상태로 시작
    await page.goto('/login');
    await page.getByPlaceholder('Email').fill(email);
    await page.getByPlaceholder('Password').fill(password);
    await page.getByRole('button', { name: /Sign In/i }).click();
    await page.waitForURL('**/projects');
  });

  test('로그아웃 버튼 클릭 시 메인 페이지로 리다이렉션되어야 함', async ({
    page,
  }) => {
    // Sign Out 버튼 클릭
    await page.getByText('Sign Out').click();

    // 메인 페이지로 리다이렉션 확인 (/ -> /main으로 리다이렉트됨)
    await page.waitForURL('**/main');
    await expect(page).toHaveURL(/.*\/main$/);

    // Sign In 링크가 다시 보이는지 확인
    await expect(page.getByRole('link', { name: 'Sign In' })).toBeVisible();
  });

  test('로그아웃 후 인증 토큰이 제거되어야 함', async ({ page, context }) => {
    // 로그아웃
    await page.getByText('Sign Out').click();
    await page.waitForURL('**/main');

    // 쿠키에서 토큰이 제거되었는지 확인
    const cookies = await context.cookies();
    const hasAuthToken = cookies.some(
      (cookie) =>
        cookie.name === 'access_token' ||
        cookie.name === 'token' ||
        cookie.name.includes('auth')
    );

    expect(hasAuthToken).toBeFalsy();
  });
});
