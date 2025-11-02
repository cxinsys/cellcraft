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
    let stateReady = true;
    try {
      await fs.stat(authFile);
    } catch {
      stateReady = false;
    }

    if (!stateReady) {
      console.log('Authenticating user...');
      // Create context with baseURL
      const context = await browser.newContext({ baseURL });
      const page = await context.newPage();

      // Navigate to login page
      await page.goto('/login');

      // Fill in credentials
      await page.getByPlaceholder('Email').fill(credentials.email);
      await page.getByPlaceholder('Password').fill(credentials.password);

      // Submit form and wait for redirect to projects page
      await Promise.all([
        page.waitForURL('**/projects', { timeout: 10000 }),
        page.getByRole('button', { name: /Sign In/i }).click(),
      ]);

      // Verify successful login by checking for projects page
      await expect(page).toHaveURL(/.*\/projects/);

      // Save authentication state
      await context.storageState({ path: authFile });
      console.log('Authentication state saved to', authFile);

      await page.close();
      await context.close();
    } else {
      console.log('Reusing existing authentication state');
    }

    await use(authFile);
  },
});

export { expect };
