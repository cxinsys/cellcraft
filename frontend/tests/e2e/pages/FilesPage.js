// frontend/tests/e2e/pages/FilesPage.js
import { expect } from '@playwright/test';
import { uploadFixture, verifyFileInList, deleteFile } from '../utils/files.js';

/**
 * Page Object Model for the Files page
 */
export class FilesPage {
  constructor(page) {
    this.page = page;

    // Locators
    this.currentFolderText = page.locator('.files__folder');
    this.uploadButton = page.locator('label.files__button');
    this.fileInput = page.locator('input[type="file"]');
    this.filesTable = page.locator('table.files__table');
    this.progressBox = page.locator('.progress__box');
    this.contextMenu = page.locator('ul.toggle__menu');
  }

  /**
   * Navigate to the Files page
   */
  async goto() {
    await this.page.goto('/files');
    await this.page.waitForLoadState('networkidle');
  }

  /**
   * Upload a file from the fixtures directory
   * @param {string} fileName - Name of the file in fixtures/files directory
   * @param {object} options - Upload options
   * @param {string} options.targetFileName - Custom name for the uploaded file (optional)
   * @returns {Promise<{response: object, uploadedFileName: string}>} Upload response and actual filename
   */
  async uploadFile(fileName, options = {}) {
    return await uploadFixture(this.page, fileName, options);
  }

  /**
   * Verify a file exists in the file list
   * @param {string} fileName - Name of the file to verify
   * @param {number} timeout - Timeout in milliseconds (default: 15000)
   */
  async verifyFileExists(fileName, timeout = 15000) {
    await verifyFileInList(this.page, fileName, { shouldExist: true, timeout });
  }

  /**
   * Verify a file does not exist in the file list
   * @param {string} fileName - Name of the file to verify
   * @param {number} timeout - Timeout in milliseconds (default: 15000)
   */
  async verifyFileNotExists(fileName, timeout = 15000) {
    await verifyFileInList(this.page, fileName, { shouldExist: false, timeout });
  }

  /**
   * Delete a file from the file list
   * @param {string} fileName - Name of the file to delete
   */
  async deleteFile(fileName) {
    await deleteFile(this.page, fileName);
  }

  /**
   * Get the current folder name
   * @returns {Promise<string>} Current folder name
   */
  async getCurrentFolder() {
    const folderText = await this.currentFolderText.textContent();
    return folderText.trim();
  }

  /**
   * Select a folder from the navigator
   * @param {string} folderName - Name of the folder to select
   */
  async selectFolder(folderName) {
    await this.page
      .locator('.folder__item', {
        has: this.page.locator(`.folder__name:has-text("${folderName}")`),
      })
      .click();

    await this.page.waitForLoadState('networkidle');
  }

  /**
   * Get the list of files in the current folder
   * @returns {Promise<Array<{name: string, date: string, type: string, size: string}>>} Array of file objects
   */
  async getFileList() {
    const rows = await this.filesTable.locator('tbody tr').all();

    const files = [];
    for (const row of rows) {
      const cells = await row.locator('td').allTextContents();
      files.push({
        name: cells[0].trim(),
        date: cells[1].trim(),
        type: cells[2].trim(),
        size: cells[3].trim(),
      });
    }

    return files;
  }

  /**
   * Get the count of files in the current folder
   * @returns {Promise<number>} Number of files
   */
  async getFileCount() {
    return await this.filesTable.locator('tbody tr').count();
  }

  /**
   * Verify the Files page is displayed
   */
  async verifyPageLoaded() {
    await expect(this.filesTable).toBeVisible();
    await expect(this.uploadButton).toBeVisible();
  }

  /**
   * Wait for upload progress to complete
   * @param {number} timeout - Maximum time to wait (default: 30000ms)
   */
  async waitForUploadComplete(timeout = 30000) {
    // Wait for progress bar to appear (if file is large enough)
    await this.progressBox
      .waitFor({ state: 'visible', timeout: 5000 })
      .catch(() => {
        // Progress bar may not appear for small files
      });

    // Wait for progress bar to disappear (upload complete)
    await this.progressBox
      .waitFor({ state: 'hidden', timeout })
      .catch(() => {
        // Already hidden or never appeared
      });

    // Wait for table to update
    await this.page.waitForLoadState('networkidle');
  }

  /**
   * Verify file appears in the table after upload
   * @param {string} fileName - Name of the uploaded file
   */
  async verifyFileUploaded(fileName) {
    const baseFileName = fileName.replace(/\.[^.]+$/, '');
    await expect(
      this.filesTable.locator('tbody tr', {
        has: this.page.locator(
          `td:first-child:text-matches("${baseFileName}", "i")`
        ),
      })
    ).toBeVisible();
  }

  /**
   * Search for a file in the table by name
   * @param {string} fileName - Name of the file to search for
   * @returns {Promise<boolean>} True if file is found, false otherwise
   */
  async isFilePresent(fileName) {
    const baseFileName = fileName.replace(/\.[^.]+$/, '');
    const count = await this.filesTable
      .locator('tbody tr', {
        has: this.page.locator(
          `td:first-child:text-matches("${baseFileName}", "i")`
        ),
      })
      .count();

    return count > 0;
  }
}
