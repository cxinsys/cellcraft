// frontend/tests/e2e/pages/modals/DataTableModal.js
import { expect } from '@playwright/test';
import { waitForTableLoad, getTableData } from '../../utils/workflow.js';

/**
 * Page Object Model for the DataTable modal
 * Handles matrix data assertions rendered via vue-good-table
 */
export class DataTableModal {
  constructor(page) {
    this.page = page;

    // Modal root
    this.modal = page.locator('.content-view');

    // vue-good-table selectors
    this.table = page.locator('.vgt-table');
    this.headers = this.table.locator('thead th');
    this.rows = this.table.locator('tbody tr');
    this.loadingOverlay = page.locator('.vgt-loading');
    this.pagination = page.locator('.vgt-pagination');
    this.noDataMessage = page.locator('text=NO DATA FOR TABLE');
  }

  /**
   * Verify the DataTable modal is open and visible
   */
  async verifyModalOpen() {
    await expect(this.modal).toBeVisible();
    await expect(this.table).toBeVisible();
  }

  /**
   * Wait until vue-good-table has finished loading data
   * @param {number} timeout
   */
  async waitForDataLoaded(timeout = 10000) {
    await waitForTableLoad(this.page, timeout);
  }

  /**
   * Check whether the placeholder "NO DATA" message is visible
   * @returns {Promise<boolean>}
   */
  async isEmptyStateVisible() {
    return await this.noDataMessage.isVisible().catch(() => false);
  }

  /**
   * Return header labels shown in the table
   * @returns {Promise<string[]>}
   */
  async getHeaderLabels() {
    const labels = await this.headers.allTextContents();
    return labels.map((label) => label.trim()).filter(Boolean);
  }

  /**
   * Return entire matrix data as array of row objects
   * @returns {Promise<Array<object>>}
   */
  async getMatrixData() {
    return await getTableData(this.page);
  }

  /**
   * Get total number of rows currently rendered
   * @returns {Promise<number>}
   */
  async getRowCount() {
    return await this.rows.count();
  }

  /**
   * Read a specific cell value
   * @param {number} rowIndex - zero-based row index
   * @param {number} columnIndex - zero-based column index
   * @returns {Promise<string>}
   */
  async getCellValue(rowIndex, columnIndex) {
    const row = this.rows.nth(rowIndex);
    const cell = row.locator(`td:nth-child(${columnIndex + 1})`);
    const text = await cell.textContent();
    return text?.trim() ?? '';
  }

  /**
   * Determine whether pagination controls are visible
   * @returns {Promise<boolean>}
   */
  async hasPagination() {
    return await this.pagination.isVisible().catch(() => false);
  }
}


