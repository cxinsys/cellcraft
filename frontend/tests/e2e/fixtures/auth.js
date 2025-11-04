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
 * Check if JWT token is expired
 * @param {string} authFilePath - Path to auth file
 * @returns {Promise<boolean>} True if token is valid, false if expired or invalid
 */
async function isTokenValid(authFilePath) {
  try {
    const authData = JSON.parse(await fs.readFile(authFilePath, 'utf-8'));
    const authCookie = authData.cookies?.find(c => c.name === 'scap_auth');

    if (!authCookie) {
      console.log('No scap_auth cookie found');
      return false;
    }

    // Check if token is expired (with 10 minute buffer)
    const expiresTimestamp = authCookie.expires;
    const nowTimestamp = Date.now() / 1000;
    const bufferSeconds = 600; // 10 minutes (increased for test stability)

    if (expiresTimestamp <= nowTimestamp + bufferSeconds) {
      console.log(`Token expired or expiring soon (expires: ${new Date(expiresTimestamp * 1000).toISOString()})`);
      return false;
    }

    console.log(`Token is valid (expires: ${new Date(expiresTimestamp * 1000).toISOString()})`);
    return true;
  } catch (error) {
    console.log('Failed to validate token:', error.message);
    return false;
  }
}

/**
 * Authentication fixture that manages login session state
 * - Logs in once per worker and stores session in .auth/test-user.json
 * - Reuses stored session for subsequent tests if token is valid
 * - Automatically re-authenticates if token is expired
 * - Supports environment variable overrides for credentials
 */
export const test = base.extend({
  storageState: async ({ browser, baseURL }, use) => {
    let needsAuthentication = true;

    // Check if auth file exists and token is valid
    try {
      await fs.stat(authFile);
      if (await isTokenValid(authFile)) {
        console.log('Reusing existing authentication state');
        needsAuthentication = false;
        await use(authFile);
        return;
      } else {
        console.log('Token expired, re-authenticating...');
        await fs.unlink(authFile).catch(() => {});
      }
    } catch {
      console.log('No existing authentication state found');
    }

    if (needsAuthentication) {
      console.log('Authenticating user...');

      // Ensure .auth directory exists
      const authDir = path.dirname(authFile);
      await fs.mkdir(authDir, { recursive: true });

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
