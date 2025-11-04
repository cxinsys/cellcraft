// frontend/tests/e2e/utils/navigation.js
import { expect } from '@playwright/test';

/**
 * Navigate to a specific page using the header menu
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {string} linkText - Text of the menu link (e.g., 'Workflows', 'Files', 'Datasets', 'Plugins')
 * @param {string} expectedPath - Expected URL path after navigation (e.g., 'projects', 'files')
 */
export async function navigateTo(page, linkText, expectedPath) {
  // Click the navigation link
  await page.getByRole('link', { name: linkText, exact: true }).click();

  // Wait for URL to change
  await page.waitForURL(`**/${expectedPath}`, { timeout: 10000 });

  // Verify we're on the correct page
  await expect(page).toHaveURL(new RegExp(`/${expectedPath}$`));
}

/**
 * Navigate to Projects/Workflows page
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
export async function goToProjects(page) {
  await navigateTo(page, 'Workflows', 'projects');
}

/**
 * Navigate to Files page
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
export async function goToFiles(page) {
  await navigateTo(page, 'Files', 'files');
}

/**
 * Navigate to Datasets page
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
export async function goToDatasets(page) {
  await navigateTo(page, 'Datasets', 'datasets');
}

/**
 * Navigate to Plugins page
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
export async function goToPlugins(page) {
  await navigateTo(page, 'Plugins', 'plugins');
}

/**
 * Verify a toast/notification message appears
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {string} text - Expected toast message text
 * @param {object} options - Options for toast verification
 * @param {number} options.timeout - Maximum time to wait for toast (default: 5000ms)
 */
export async function expectToast(page, text, options = {}) {
  const { timeout = 5000 } = options;
  await expect(page.getByText(text, { exact: false })).toBeVisible({
    timeout,
  });
}

/**
 * Wait for a specific API call to complete
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {string} urlPattern - URL pattern to match (e.g., '/routes/files/upload')
 * @param {object} options - Wait options
 * @param {number} options.status - Expected HTTP status code (default: 200)
 * @param {number} options.timeout - Maximum time to wait (default: 30000ms)
 * @returns {Promise<Response>} The response object
 */
export async function waitForApiCall(page, urlPattern, options = {}) {
  const { status = 200, timeout = 30000 } = options;

  const response = await page.waitForResponse(
    (resp) => resp.url().includes(urlPattern) && resp.status() === status,
    { timeout }
  );

  return response;
}

/**
 * Verify the current page URL
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {string} expectedPath - Expected URL path
 */
export async function verifyCurrentPage(page, expectedPath) {
  await expect(page).toHaveURL(new RegExp(`/${expectedPath}$`));
}

/**
 * Wait for page to be fully loaded (networkidle + DOM content loaded)
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {number} timeout - Maximum time to wait (default: 30000ms)
 */
export async function waitForPageReady(page, timeout = 30000) {
  await Promise.all([
    page.waitForLoadState('networkidle', { timeout }),
    page.waitForLoadState('domcontentloaded', { timeout }),
  ]);
}

/**
 * Click the logo to return to main/projects page
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
export async function clickLogo(page) {
  await page.locator('.logo').click();
  await page.waitForLoadState('networkidle');
}

/**
 * Sign out the current user
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
export async function signOut(page) {
  await page.getByText('Sign Out').click();
  await page.waitForURL('**/');
}

/**
 * Verify user is authenticated by checking header menu
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
export async function verifyAuthenticated(page) {
  await expect(page.getByRole('link', { name: 'Workflows' })).toBeVisible();
  await expect(page.getByText('Sign Out')).toBeVisible();
}

/**
 * Verify user is not authenticated
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
export async function verifyNotAuthenticated(page) {
  await expect(page.getByRole('link', { name: 'Sign In' })).toBeVisible();
  await expect(page.getByText('Sign Out')).not.toBeVisible();
}
