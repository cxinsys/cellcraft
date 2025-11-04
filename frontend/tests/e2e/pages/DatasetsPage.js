// frontend/tests/e2e/pages/DatasetsPage.js
import { expect } from '@playwright/test';

/**
 * Page Object Model for the Datasets page
 */
export class DatasetsPage {
  constructor(page) {
    this.page = page;

    // Locators
    this.searchInput = page.getByPlaceholder('Search titles...');
    this.prevButton = page.getByRole('button', { name: 'Prev' });
    this.nextButton = page.getByRole('button', { name: 'Next' });
    this.pageNumber = page.locator('.pagination span');
    this.datasetsTable = page.locator('table');
  }

  /**
   * Navigate to the Datasets page
   */
  async goto() {
    await this.page.goto('/datasets');
    await this.page.waitForLoadState('networkidle');
  }

  /**
   * Search for a dataset by title
   * @param {string} searchTerm - Search term to filter datasets
   */
  async searchDataset(searchTerm) {
    await this.searchInput.fill(searchTerm);
    await this.page.waitForTimeout(300); // Wait for search filter to apply
  }

  /**
   * Download a dataset by title
   * @param {string} datasetTitle - Title of the dataset to download (e.g., 'pbmc_light_1000.h5ad')
   * @returns {Promise<Download>} Playwright Download object
   */
  async downloadDataset(datasetTitle) {
    // Wait for download event
    const downloadPromise = this.page.waitForEvent('download', {
      timeout: 60000,
    });

    // Find the row with the dataset title and click the download button
    const datasetRow = this.page.locator('tr', {
      has: this.page.locator(`.title-container:has-text("${datasetTitle}")`),
    });

    await datasetRow.locator('button.download-button').click();

    // Wait for download to start
    const download = await downloadPromise;

    return download;
  }

  /**
   * Verify a dataset download completed successfully
   * @param {Download} download - Playwright Download object
   * @param {string} expectedFilename - Expected filename of the download
   */
  async verifyDownload(download, expectedFilename) {
    const suggestedFilename = download.suggestedFilename();
    expect(suggestedFilename).toBe(expectedFilename);

    // Optionally save the file to verify it's not corrupted
    const path = await download.path();
    expect(path).toBeTruthy();
  }

  /**
   * Get the list of visible dataset titles on the current page
   * @returns {Promise<string[]>} Array of dataset titles
   */
  async getVisibleDatasets() {
    const titles = await this.page
      .locator('.title-container')
      .allTextContents();
    return titles.map((title) => title.trim());
  }

  /**
   * Navigate to the next page of datasets
   */
  async goToNextPage() {
    await this.nextButton.click();
    await this.page.waitForTimeout(300);
  }

  /**
   * Navigate to the previous page of datasets
   */
  async goToPreviousPage() {
    await this.prevButton.click();
    await this.page.waitForTimeout(300);
  }

  /**
   * Get the current page number
   * @returns {Promise<number>} Current page number
   */
  async getCurrentPage() {
    const pageText = await this.pageNumber.textContent();
    return parseInt(pageText, 10);
  }

  /**
   * Verify the Datasets page is displayed
   */
  async verifyPageLoaded() {
    await expect(this.page.locator('.header__text:has-text("Datasets")')).toBeVisible();
    await expect(this.datasetsTable).toBeVisible();
  }

  /**
   * Verify a specific dataset is visible
   * @param {string} datasetTitle - Title of the dataset
   */
  async verifyDatasetVisible(datasetTitle) {
    await expect(
      this.page.locator(`.title-container:has-text("${datasetTitle}")`)
    ).toBeVisible();
  }
}
