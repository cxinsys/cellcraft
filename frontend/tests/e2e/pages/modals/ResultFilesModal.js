// frontend/tests/e2e/pages/modals/ResultFilesModal.js
import { expect } from '@playwright/test';

/**
 * Page Object Model for the ResultFiles modal
 * Handles section visibility, file selection, and bulk download actions
 */
export class ResultFilesModal {
  constructor(page) {
    this.page = page;

    // Modal root
    this.modal = page.locator('.result__layout');

    // Sections
    this.primarySection = page.locator('.files__section--final');
    this.intermediateSection = page.locator('.files__section--intermediate');
    this.intermediateList = this.intermediateSection.locator('.files__list--intermediate');

    // Buttons
    this.selectAllButton = page.locator('button.action__button--secondary');
    this.primarySelectAllButton = this.primarySection.locator('button.action__button--small');
    this.intermediateSelectAllButton = this.intermediateSection.locator('button.action__button--small');
    this.setFilesButton = page.locator('button.action__button--primary');
    this.downloadSelectedButton = page.locator('button.action__button--download');

    // Misc selectors
    this.selectionInfo = page.locator('.files__selection-info');
  }

  async verifyModalOpen() {
    await expect(this.modal).toBeVisible();
  }

  async waitForPrimarySection(timeout = 20000) {
    await this.primarySection.waitFor({ state: 'visible', timeout });
    await this.primarySection.locator('.file__item').first().waitFor({ state: 'visible', timeout });
  }

  async waitForIntermediateSection(timeout = 20000) {
    await this.intermediateSection.waitFor({ state: 'visible', timeout }).catch(() => null);

    const listVisible = await this.intermediateList.isVisible().catch(() => false);
    if (listVisible) {
      await this.intermediateList.locator('.file__item').first().waitFor({ state: 'visible', timeout }).catch(() => null);
    }
  }

  async isPrimarySectionVisible() {
    return await this.primarySection.isVisible().catch(() => false);
  }

  async isIntermediateSectionVisible() {
    return await this.intermediateSection.isVisible().catch(() => false);
  }

  async getPrimaryFileNames() {
    const names = await this.primarySection
      .locator('.file__name')
      .allTextContents();
    return names.map((name) => name.trim()).filter(Boolean);
  }

  async getIntermediateFileNames() {
    await this.ensureIntermediateExpanded();
    const names = await this.intermediateSection
      .locator('.file__name')
      .allTextContents();
    return names.map((name) => name.trim()).filter(Boolean);
  }

  primaryFileItem(fileName) {
    return this.primarySection.locator('.file__item').filter({
      has: this.page.locator('.file__name', { hasText: fileName }),
    });
  }

  intermediateFileItem(fileName) {
    return this.intermediateSection.locator('.file__item').filter({
      has: this.page.locator('.file__name', { hasText: fileName }),
    });
  }

  async ensureIntermediateExpanded() {
    const visible = await this.intermediateList.isVisible().catch(() => false);
    if (!visible) {
      await this.intermediateSection.locator('.section__header').click();
      await this.intermediateList.waitFor({ state: 'visible', timeout: 5000 });
    }
  }

  async selectPrimaryFile(fileName) {
    const item = this.primaryFileItem(fileName);
    await item.waitFor({ state: 'visible', timeout: 5000 });
    if (await this.isPrimaryFileSelected(fileName)) {
      return;
    }
    await item.locator('.file__checkbox-label').click();
    await expect.poll(async () => this.isPrimaryFileSelected(fileName)).toBeTruthy();
  }

  async deselectPrimaryFile(fileName) {
    const item = this.primaryFileItem(fileName);
    if (!(await this.isPrimaryFileSelected(fileName))) {
      return;
    }
    await item.locator('.file__checkbox-label').click();
    await expect.poll(async () => this.isPrimaryFileSelected(fileName)).toBeFalsy();
  }

  async selectIntermediateFile(fileName) {
    await this.ensureIntermediateExpanded();
    const item = this.intermediateFileItem(fileName);
    await item.waitFor({ state: 'visible', timeout: 5000 });
    if (await this.isIntermediateFileSelected(fileName)) {
      return;
    }
    await item.locator('.file__checkbox-label').click();
    await expect.poll(async () => this.isIntermediateFileSelected(fileName)).toBeTruthy();
  }

  async deselectIntermediateFile(fileName) {
    await this.ensureIntermediateExpanded();
    const item = this.intermediateFileItem(fileName);
    if (!(await this.isIntermediateFileSelected(fileName))) {
      return;
    }
    await item.locator('.file__checkbox-label').click();
    await expect.poll(async () => this.isIntermediateFileSelected(fileName)).toBeFalsy();
  }

  async isPrimaryFileSelected(fileName) {
    const item = this.primaryFileItem(fileName);
    const classes = await item.getAttribute('class');
    return classes?.includes('file__item--selected') ?? false;
  }

  async isIntermediateFileSelected(fileName) {
    const item = this.intermediateFileItem(fileName);
    const classes = await item.getAttribute('class');
    return classes?.includes('file__item--selected') ?? false;
  }

  async toggleSelectAll() {
    await this.selectAllButton.click();
  }

  async togglePrimarySelectAll() {
    await this.primarySelectAllButton.click();
  }

  async toggleIntermediateSelectAll() {
    await this.ensureIntermediateExpanded();
    await this.intermediateSelectAllButton.click();
  }

  async clickSetFiles() {
    await this.setFilesButton.click();
  }

  async clickDownloadSelected() {
    await this.downloadSelectedButton.click();
  }

  async isDownloadSelectedVisible() {
    return await this.downloadSelectedButton.isVisible().catch(() => false);
  }

  async isDownloadSelectedEnabled() {
    return await this.downloadSelectedButton.isEnabled().catch(() => false);
  }

  async getSelectionSummaryText() {
    const text = await this.selectionInfo.textContent();
    return text?.trim() ?? '';
  }

  async getSelectedCount() {
    const selectedItems = await this.modal.locator('.file__item--selected').count();
    return selectedItems;
  }

  async clickPrimaryFileDownloadButton(fileName) {
    const item = this.primaryFileItem(fileName);
    await item.waitFor({ state: 'visible', timeout: 5000 });
    await item.locator('button.file__download-btn').click();
  }
}


