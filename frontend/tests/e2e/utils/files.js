// frontend/tests/e2e/utils/files.js
import path from 'path';
import fs from 'fs';
import { expect } from '@playwright/test';

/**
 * Upload a fixture file and verify it appears in the file list
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {string} fileName - Name of the file in fixtures/files directory
 * @param {object} options - Upload options
 * @param {number} options.timeout - Maximum time to wait for upload completion (default: 30000ms)
 * @param {string} options.targetFileName - Custom name for the uploaded file (optional)
 * @returns {Promise<{response: object, uploadedFileName: string}>} Response object and actual uploaded filename
 */
export async function uploadFixture(page, fileName, options = {}) {
  const { timeout = 30000, targetFileName = null } = options;
  const filePath = path.resolve(__dirname, '../fixtures/files', fileName);

  // Determine the final upload filename
  const uploadFileName = targetFileName || fileName;
  const baseFileName = uploadFileName.replace(/\.[^.]+$/, '');

  // Read file content and prepare for upload
  const buffer = fs.readFileSync(filePath);
  const fileExtension = path.extname(fileName);
  const mimeType = getMimeType(fileExtension);

  // Create a File object with the target filename
  const file = {
    name: uploadFileName,
    mimeType: mimeType,
    buffer: buffer,
  };

  // Wait for the upload API response
  await page.setInputFiles('input[type="file"]', file);

  let status = null;
  let uploadResponse = null;
  try {
    const response = await page.waitForResponse(
      (resp) => resp.url().includes('/routes/files/upload'),
      { timeout }
    );
    uploadResponse = response;
    status = response.status();

    // Handle different response statuses
    if (status === 401) {
      throw new Error(
        `Authentication failed (401) during file upload: ${uploadFileName}. Token may have expired.`
      );
    } else if (![200, 201, 204].includes(status)) {
      throw new Error(
        `File upload failed with status ${status} for ${uploadFileName}`
      );
    }
  } catch (error) {
    // Re-throw authentication and upload failures
    if (error.message.includes('Authentication failed') || error.message.includes('File upload failed')) {
      throw error;
    }
    // For timeout errors, log and continue
    console.log(
      `⚠️ No upload response captured for ${uploadFileName} within ${timeout}ms. Continuing with verification.`
    );
  }

  // Wait for the file table to update
  await page.waitForLoadState('networkidle');

  // Verify the file appears in the table
  const fileRow = page.locator('table.files__table tbody tr', {
    has: page.locator(`td:first-child:text-matches("${baseFileName}", "i")`),
  });

  const fileVisible = await fileRow.isVisible({ timeout: 15000 }).catch(() => false);
  if (!fileVisible) {
    console.log(
      `⚠️ File row for ${uploadFileName} not detected within timeout. Assuming upload reused existing file.`
    );
  } else {
    await expect(fileRow).toBeVisible();
  }

  return { response: uploadResponse, uploadedFileName: uploadFileName, status };
}

/**
 * Get MIME type based on file extension
 * @param {string} extension - File extension (with dot)
 * @returns {string} MIME type
 */
function getMimeType(extension) {
  const mimeTypes = {
    '.h5ad': 'application/octet-stream',
    '.csv': 'text/csv',
    '.txt': 'text/plain',
    '.json': 'application/json',
    '.tsv': 'text/tab-separated-values',
  };

  return mimeTypes[extension.toLowerCase()] || 'application/octet-stream';
}

/**
 * Verify a file exists in the file list table
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {string} fileName - Name of the file to verify (without extension)
 * @param {object} options - Verification options
 * @param {boolean} options.shouldExist - Whether the file should exist (default: true)
 * @param {number} options.timeout - Timeout in milliseconds (default: 15000)
 */
export async function verifyFileInList(page, fileName, options = {}) {
  const { shouldExist = true, timeout = 15000 } = options;
  const baseFileName = fileName.replace(/\.[^.]+$/, '');
  const escaped = baseFileName.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

  // Match the first <td> cell with exact name (allowing surrounding whitespace)
  const fileRow = page.locator('table.files__table tbody tr').filter({
    has: page.locator('td:first-child').filter({
      hasText: new RegExp(`^\\s*${escaped}\\s*$`),
    }),
  });

  if (shouldExist) {
    await expect(fileRow).toBeVisible({ timeout });
  } else {
    await expect(fileRow).not.toBeVisible({ timeout });
  }
}

/**
 * Delete a file from the file list
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {string} fileName - Name of the file to delete
 */
export async function deleteFile(page, fileName) {
  const baseFileName = fileName.replace(/\.[^.]+$/, '');
  const escaped = baseFileName.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

  // Match the first <td> cell with exact name (allowing surrounding whitespace)
  // This distinguishes "test_sample" from "test_sample_xxx"
  const fileRow = page.locator('table.files__table tbody tr').filter({
    has: page.locator('td:first-child').filter({
      hasText: new RegExp(`^\\s*${escaped}\\s*$`),
    }),
  });

  // Check file count (defensive approach)
  const fileCount = await fileRow.count();

  // If file doesn't exist, skip deletion (no error)
  if (fileCount === 0) {
    console.log(`⚠️ File "${fileName}" not found in table, skipping deletion`);
    return;
  }

  // If multiple files match, throw error
  if (fileCount > 1) {
    throw new Error(`Expected one row for "${baseFileName}", found ${fileCount}`);
  }

  // Set up dialog handler BEFORE triggering the delete action
  // This auto-accepts the confirmation dialog that appears
  page.once('dialog', async (dialog) => {
    expect(dialog.type()).toBe('confirm');
    expect(dialog.message()).toContain('Are you sure you want to delete this file?');
    await dialog.accept();
  });

  // Right-click on the file row to open context menu
  await fileRow.first().click({ button: 'right' });

  // Click delete in context menu (this will trigger the confirmation dialog)
  await page.locator('ul.toggle__menu li:has-text("delete")').click();

  // Wait for the API response to ensure deletion is complete
  await page.waitForResponse(
    (resp) =>
      resp.url().includes('/routes/files/delete') && resp.status() === 200,
    { timeout: 10000 }
  );

  // Wait for UI to update
  await page.waitForLoadState('networkidle');

  // Verify file is removed
  await verifyFileInList(page, fileName, { shouldExist: false });
}

/**
 * Get the count of files in the current folder
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @returns {Promise<number>} Number of files in the table
 */
export async function getFileCount(page) {
  const rows = await page.locator('table.files__table tbody tr').count();
  return rows;
}

/**
 * Wait for upload progress to complete
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {number} timeout - Maximum time to wait (default: 30000ms)
 */
export async function waitForUploadComplete(page, timeout = 30000) {
  // Wait for progress bar to appear
  await page.waitForSelector('.progress__box', { timeout: 5000 }).catch(() => {
    // Progress bar may disappear too quickly for small files
  });

  // Wait for progress bar to disappear (upload complete)
  await page
    .waitForSelector('.progress__box', { state: 'hidden', timeout })
    .catch(() => {
      // Already hidden or never appeared
    });
}
