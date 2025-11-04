// frontend/tests/e2e/pages/modals/LogsModal.js
import { expect } from '@playwright/test';

export class LogsModal {
  constructor(page) {
    this.page = page;
    this.overlay = page.locator('.logs-modal-overlay');
    this.modal = this.overlay.locator('.logs-modal');
    this.loadingIndicator = this.modal.locator('.logs-loading');
    this.content = this.modal.locator('.logs-content');
    this.emptyState = this.modal.locator('.no-logs');
    this.exportAllButton = this.modal.locator('button.export-all-btn');
    this.refreshButton = this.modal.locator('button.refresh-btn');
    this.closeButton = this.modal.locator('button.close-btn');
    this.logFiles = this.modal.locator('.log-file');
  }

  async waitForOpen() {
    await this.overlay.waitFor({ state: 'visible', timeout: 10000 });
    await expect(this.modal).toBeVisible();
  }

  async waitForLoaded() {
    await this.loadingIndicator.waitFor({ state: 'hidden', timeout: 20000 }).catch(() => null);
    if (await this.emptyState.isVisible()) {
      return;
    }
    await expect(this.content).toBeVisible();
  }

  async hasLogs() {
    return (await this.logFiles.count()) > 0;
  }

  async expectLogsAvailable() {
    await expect(this.logFiles.first()).toBeVisible();
  }

  async downloadAllLogsJson() {
    if (await this.exportAllButton.isDisabled()) {
      throw new Error('Export All (JSON) button is disabled');
    }

    const download = this.page.waitForEvent('download');
    await this.exportAllButton.click();
    return await download;
  }

  async downloadFirstLogTxt() {
    const firstLog = this.logFiles.first();
    await expect(firstLog).toBeVisible();

    const exportButton = firstLog.locator('button.export-txt-btn');
    await expect(exportButton).toBeVisible();

    const download = this.page.waitForEvent('download');
    await exportButton.click();
    return await download;
  }

  async close() {
    await this.closeButton.click();
    await this.overlay.waitFor({ state: 'hidden', timeout: 5000 });
  }
}

