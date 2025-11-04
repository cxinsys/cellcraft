// frontend/tests/e2e/fixtures/auth.js
import { test as base, expect } from '@playwright/test';
import path from 'path';
import fs from 'fs/promises';

const authFile = path.resolve(__dirname, '../.auth/test-user.json');
const credentials = {
  email: process.env.PLAYWRIGHT_USER ?? 'test1234@test.com',
  password: process.env.PLAYWRIGHT_PASS ?? 'test1234',
};

/**
 * Authentication fixture that manages login session state
 * - Logs in once per worker and stores session in .auth/test-user.json
 * - Reuses stored session for subsequent tests
 * - Supports environment variable overrides for credentials
 */
export const test = base.extend({
  storageState: async ({ browser, baseURL }, use) => {
    let firstAttemptFailed = false;
    try {
      await fs.stat(authFile);
    } catch {
      firstAttemptFailed = true;
    }

    if (!firstAttemptFailed) {
      console.log('Reusing existing authentication state');
      await use(authFile);

      try {
        await fs.stat(authFile);
      } catch {
        firstAttemptFailed = true;
      }
    }

    if (firstAttemptFailed) {
      console.log('Authenticating user...');
      const context = await browser.newContext({ baseURL });
      const page = await context.newPage();

      await page.goto('/login');
      await page.getByPlaceholder('Email').fill(credentials.email);
      await page.getByPlaceholder('Password').fill(credentials.password);

      await Promise.all([
        page.waitForURL('**/projects', { timeout: 15000 }),
        page.getByRole('button', { name: /Sign In/i }).click(),
      ]);

      await expect(page).toHaveURL(/.*\/projects/);

      await context.storageState({ path: authFile });
      console.log('Authentication state saved to', authFile);

      await page.close();
      await context.close();

      await use(authFile);
    }
  },
});

export { expect };
