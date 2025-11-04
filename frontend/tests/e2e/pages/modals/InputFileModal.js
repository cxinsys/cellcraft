// frontend/tests/e2e/pages/modals/InputFileModal.js
import { expect } from '@playwright/test';

/**
 * Page Object Model for the InputFile modal
 * Handles file selection and assignment to InputFile nodes
 */
export class InputFileModal {
  constructor(page) {
    this.page = page;

    // Modal container
    this.modal = page.locator('.content-view');

    // Folder navigation
    this.folderList = page.locator('.folder__list');
    this.folderItems = page.locator('.folder__item');
    this.cloudRow = page.locator('.cloud-row');

    // File selection
    this.fileList = page.locator('.form__fileList');
    this.fileItems = page.locator('.file__item');

    // Current file display
    this.currentFileSection = page.locator('.form__info');
    this.currentFileName = page.locator('.form__info--name');
    this.currentFileSize = page.locator('.form__info--size');
    this.blankFileMessage = page.locator('.form__info--blank');

    // Action buttons
    // Use label selector to avoid interference from adjacent input element
    this.applyButton = page.locator('label.form__button--apply');
    this.fileInput = page.locator('.form__input[type="submit"]');

    // Loading indicator
    this.loadingLayout = page.locator('.loading-layout');
  }

  /**
   * Verify the InputFile modal is open and visible
   */
  async verifyModalOpen() {
    await expect(this.modal).toBeVisible();
    await expect(this.folderList).toBeVisible();
    await expect(this.applyButton).toBeVisible();
  }

  /**
   * Get the list of available folders
   * @returns {Promise<Array<string>>} Array of folder names
   */
  async getFolders() {
    const folderNamesLocator = this.folderItems.locator('.folder__name');

    await folderNamesLocator.first().waitFor({ state: 'visible', timeout: 10000 });

    const items = await folderNamesLocator.evaluateAll((elements) =>
      elements
        .map((el) => el.textContent?.trim())
        .filter((name) => name && name.length > 0)
    );

    return items;
  }

  /**
   * Select a folder to view its files
   * @param {string} folderName - Name of the folder to select
   */
  async selectFolder(folderName) {
    const folder = this.page.locator('.folder__item', {
      has: this.page.locator(`.folder__name:has-text("${folderName}")`),
    });

    await folder.click();

    // Wait for files to load
    await this.page.waitForLoadState('networkidle');

    // Wait for folder to be highlighted
    await expect(folder).toHaveClass(/toggleFolder/);
  }

  /**
   * Verify a folder is currently selected
   * @param {string} folderName - Name of the folder to verify
   */
  async verifyFolderSelected(folderName) {
    const folder = this.page.locator('.folder__item', {
      has: this.page.locator(`.folder__name:has-text("${folderName}")`),
    });

    await expect(folder).toHaveClass(/toggleFolder/);
  }

  /**
   * Get the list of files in the currently selected folder
   * @returns {Promise<Array<string>>} Array of file names
   */
  async getFiles() {
    const fileNamesLocator = this.fileItems.locator('.folder__name');

    // Wait for at least one file name to be visible before collecting all
    await fileNamesLocator.first().waitFor({ state: 'visible', timeout: 10000 });

    // Use evaluateAll for better performance and reliability
    const items = await fileNamesLocator.evaluateAll((elements) =>
      elements
        .map((el) => el.textContent?.trim())
        .filter((name) => name && name.length > 0)
    );

    return items;
  }

  /**
   * Select a file from the file list
   * @param {string} fileName - Name of the file to select
   */
  async selectFile(fileName) {
    // Select file from .form__selectFile area
    const file = this.page
      .locator('.form__selectFile .file__item')
      .filter({
        has: this.page.locator(`.folder__name:has-text("${fileName}")`),
      });

    await file.click();

    // fileClick() in InputFile.vue makes an async API call to findFile()
    // We need to wait for this network request to complete
    await this.page.waitForLoadState('networkidle', { timeout: 5000 });

    // Wait for current file info to update
    await this.page.waitForTimeout(300);

    // Verify file selection by checking Current File display
    // This ensures the API response has been received and selectFile.data is populated
    await expect(this.currentFileName).toHaveText(fileName);
  }

  /**
   * Verify a file is currently selected
   * @param {string} fileName - Name of the file to verify
   */
  async verifyFileSelected(fileName) {
    const file = this.page.locator('.file__item', {
      has: this.page.locator(`.folder__name:has-text("${fileName}")`),
    });

    await expect(file).toHaveClass(/toggleFile/);
  }

  /**
   * Get the currently selected file information
   * @returns {Promise<{name: string, size: string}|null>} File info or null if no file selected
   */
  async getCurrentFileInfo() {
    // Check if blank message is shown
    const isBlank = await this.blankFileMessage.isVisible().catch(() => false);
    if (isBlank) {
      return null;
    }

    const name = await this.currentFileName.textContent();
    const size = await this.currentFileSize.textContent();

    return {
      name: name?.trim(),
      size: size?.trim(),
    };
  }

  /**
   * Verify the current file display shows expected file
   * @param {string} fileName - Expected file name
   */
  async verifyCurrentFile(fileName) {
    await expect(this.currentFileName).toHaveText(fileName);
  }

  /**
   * Verify no file is currently selected (blank state)
   */
  async verifyNoFileSelected() {
    await expect(this.blankFileMessage).toBeVisible();
    await expect(this.blankFileMessage).toHaveText('Please add data file');
  }

  /**
   * Click the Apply button to assign the file
   */
  async clickApply() {
    // Check current button state
    const isApplied =
      (await this.applyButton.textContent())?.trim() === 'Applied';

    if (isApplied) {
      console.log('File already applied, clicking to apply again');
    }

    await this.applyButton.click();

    // applyFile() is synchronous but Vue reactivity updates DOM asynchronously
    // Wait for button text to change to "Applied" (more reliable than class check)
    await expect(this.applyButton).toContainText('Applied', { timeout: 3000 });

    // Also verify the activate class was added (green background)
    await expect(this.applyButton).toHaveClass(/activate/, { timeout: 1000 });

    // Additional wait to ensure Vuex store commit completes
    // This ensures file info is persisted before proceeding
    await this.page.waitForTimeout(500);
  }

  /**
   * Verify the Apply button state
   * @param {boolean} shouldBeApplied - Expected state (true for "Applied", false for "Apply")
   */
  async verifyApplyButtonState(shouldBeApplied) {
    const expectedText = shouldBeApplied ? 'Applied' : 'Apply';

    // Create a fresh locator to avoid stale element references after DOM re-renders
    const applyButton = this.page.locator('label.form__button--apply');

    // Check if we're still on the correct page before waiting for button
    const currentUrl = this.page.url();
    if (!currentUrl.includes('/workflow')) {
      throw new Error(
        `Unexpected navigation detected. Current URL: ${currentUrl}. ` +
          `Expected to be on /workflow page with InputFile modal open.`
      );
    }

    // Verify modal is still visible
    const modalVisible = await this.modal.isVisible().catch(() => false);
    if (!modalVisible) {
      throw new Error(
        `InputFile modal is not visible. Current URL: ${currentUrl}. ` +
          `Modal may have been closed unexpectedly.`
      );
    }

    // Wait for the element to be visible first
    await applyButton.waitFor({ state: 'visible', timeout: 10000 });

    // Use toContainText() to ignore the input value="업로드" inside the label
    await expect(applyButton).toContainText(expectedText);

    if (shouldBeApplied) {
      await expect(applyButton).toHaveClass(/activate/);
    }
  }

  /**
   * Wait for loading to complete
   * @param {number} timeout - Maximum time to wait (default: 10000ms)
   */
  async waitForLoading(timeout = 10000) {
    // Wait for loading indicator to appear
    await this.loadingLayout
      .waitFor({ state: 'visible', timeout: 2000 })
      .catch(() => {
        // Loading might be too fast to catch
      });

    // Wait for loading to disappear
    await this.loadingLayout
      .waitFor({ state: 'hidden', timeout })
      .catch(() => {
        // Already hidden
      });
  }

  /**
   * Complete workflow: select folder, select file, and apply
   * @param {string} folderName - Folder to select
   * @param {string} fileName - File to select
   */
  async assignFile(folderName, fileName) {
    await this.selectFolder(folderName);
    await this.selectFile(fileName);
    await this.verifyCurrentFile(fileName);
    await this.clickApply();
    await this.waitForLoading();
    await this.verifyApplyButtonState(true);
  }

  /**
   * Navigate to Files page (via the cloud row link)
   */
  async goToFilesPage() {
    await this.cloudRow.click();
    await this.page.waitForURL('**/files', { timeout: 5000 });
  }

  /**
   * Verify a specific file exists in the file list
   * @param {string} fileName - File name to verify
   */
  async verifyFileExists(fileName) {
    await expect(
      this.page.locator('.file__item', {
        has: this.page.locator(`.folder__name:has-text("${fileName}")`),
      })
    ).toBeVisible();
  }

  /**
   * Verify a specific folder exists in the folder list
   * @param {string} folderName - Folder name to verify
   */
  async verifyFolderExists(folderName) {
    await expect(
      this.page.locator('.folder__item', {
        has: this.page.locator(`.folder__name:has-text("${folderName}")`),
      })
    ).toBeVisible();
  }

  /**
   * Get the count of files in the current view
   * @returns {Promise<number>}
   */
  async getFileCount() {
    return await this.fileItems.count();
  }

  /**
   * Get the count of folders in the list
   * @returns {Promise<number>}
   */
  async getFolderCount() {
    return await this.folderItems.count();
  }
}
