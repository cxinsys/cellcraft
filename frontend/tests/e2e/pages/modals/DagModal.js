// frontend/tests/e2e/pages/modals/DagModal.js
import { expect } from '@playwright/test';

export class DagModal {
  constructor(page) {
    this.page = page;
    this.modal = page.locator('.dag-modal');
    this.content = this.modal.locator('.dag-content');
    this.loadingOverlay = this.modal.locator('.dag-plot .loading-overlay');
    this.plotContainer = this.modal.locator('.dag-plot');
    this.progressText = this.modal.locator('.dag-progress .progress-text');
    this.closeButton = this.modal.locator('.dag-controls .close-btn');
    this.errorMessage = this.modal.locator('.dag-plot .error-message');
  }

  async waitForOpen() {
    await this.modal.waitFor({ state: 'visible', timeout: 15000 });
    await expect(this.content).toBeVisible();
  }

  async waitForLoaded() {
    await this.loadingOverlay.waitFor({ state: 'hidden', timeout: 30000 }).catch(() => null);
    if (await this.errorMessage.isVisible()) {
      throw new Error(`DAG modal displayed error: ${await this.errorMessage.textContent()}`);
    }

    await expect(this.plotContainer.locator('.main-svg').first()).toBeVisible({ timeout: 20000 });
    await expect(this.progressText).toBeVisible();
  }

  async close() {
    await this.closeButton.click();
    await this.modal.waitFor({ state: 'hidden', timeout: 10000 });
  }
}

