// frontend/tests/e2e/pages/modals/CompileCheckModal.js
import { expect } from '@playwright/test';

/**
 * Page Object Model for the Compile Check modal (task execution panel)
 * Handles task summary verification, resource status checks, and execution actions
 */
export class CompileCheckModal {
  constructor(page) {
    this.page = page;

    // Core modal elements
    this.modal = page.locator('.modal-content');
    this.modalContainer = page.locator('.modal-container');
    this.loadingOverlay = this.modal.locator('.full-loading-overlay');

    // Task information
    this.taskInfoSection = this.modal.locator('.task-info');
    this.taskItems = this.modal.locator('.task-info__item');
    this.taskPlugin = this.modal.locator('.task-plugin');
    this.taskInputs = this.modal.locator('.task-input');
    this.taskOutputs = this.modal.locator('.task-output');

    // Resource information
    this.resourceSection = this.modal.locator('.resource-info');
    this.resourceBars = this.modal.locator('.resource-bar');
    this.containerInfo = this.modal.locator('.container-info');

    // Actions
    this.executeButton = this.modal.locator('.modal-actions .btn.confirm');
    this.cancelButton = this.modal.locator('.modal-actions .btn.cancel');
  }

  /**
   * Wait for the modal to be visible
   */
  async waitForOpen() {
    await this.modal.waitFor({ state: 'visible', timeout: 10000 });
    await expect(this.modalContainer).toBeVisible();
  }

  /**
   * Wait for system resource polling to finish initial load
   */
  async waitForResourcesLoaded(timeout = 20000) {
    await this.loadingOverlay.waitFor({ state: 'hidden', timeout }).catch(async () => {
      await this.loadingOverlay.waitFor({ state: 'detached', timeout });
    });

    await expect(this.modalContainer).not.toHaveClass(/loading-active/, {
      timeout,
    });
  }

  /**
   * Get all task entries including plugin, inputs, and outputs
   * @returns {Promise<Array<{plugin: string, inputs: string[], outputs: string[]}>>}
   */
  async getTaskEntries() {
    return await this.taskItems.evaluateAll((items) =>
      items.map((item) => ({
        plugin:
          item.querySelector('.task-plugin')?.textContent?.trim() ?? '',
        inputs: Array.from(item.querySelectorAll('.task-input')).map((el) =>
          el.textContent.trim()
        ),
        outputs: Array.from(item.querySelectorAll('.task-output')).map((el) =>
          el.textContent.trim()
        ),
      }))
    );
  }

  /**
   * Collect resource label texts for summary assertions
   * @returns {Promise<Array<string>>}
   */
  async getResourceLabels() {
    const labels = await this.resourceBars.locator('label').allTextContents();
    return labels.map((label) => label.trim()).filter(Boolean);
  }

  /**
   * Ensure core sections are visible to the user
   */
  async verifyCoreSectionsVisible() {
    await expect(this.taskInfoSection).toBeVisible();
    await expect(this.resourceSection).toBeVisible();
  }

  /**
   * Click the Execute button to confirm task execution
   */
  async clickExecute() {
    await this.executeButton.click();
  }

  /**
   * Cancel and close the modal
   */
  async clickCancel() {
    await this.cancelButton.click();
  }

  /**
   * Wait for the modal to fully close (detached from DOM)
   */
  async waitForClose(timeout = 15000) {
    await this.modal.waitFor({ state: 'detached', timeout });
  }
}

