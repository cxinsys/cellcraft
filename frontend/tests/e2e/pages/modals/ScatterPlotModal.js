// frontend/tests/e2e/pages/modals/ScatterPlotModal.js
import { expect } from '@playwright/test';
import { waitForPlotlyRender } from '../../utils/workflow.js';

/**
 * Page Object Model for the ScatterPlot modal
 */
export class ScatterPlotModal {
  constructor(page) {
    this.page = page;

    this.modal = page.locator('.content-view');
    this.plotContainer = page.locator('#plotly__scatter');
    this.blankState = page.getByText('NO DATA', { exact: false });
    this.xSelect = page.locator('.options__item').nth(0).locator('select');
    this.ySelect = page.locator('.options__item').nth(1).locator('select');
    this.groupSelect = page.locator('.options__item').nth(2).locator('select');
  }

  /**
   * Ensure the modal and Plotly container are visible
   */
  async verifyModalOpen() {
    await expect(this.modal).toBeVisible();
    await expect(this.plotContainer).toBeVisible();
  }

  /**
   * Wait until Plotly has rendered scatter traces
   * @param {number} timeout
   */
  async waitForPlotly(timeout = 15000) {
    await waitForPlotlyRender(this.page, 'plotly__scatter', timeout);
    await expect(this.plotContainer.locator('svg.main-svg').first()).toBeVisible();
  }

  /**
   * Returns whether the blank state message is visible
   * @returns {Promise<boolean>}
   */
  async isBlankStateVisible() {
    return await this.blankState.isVisible().catch(() => false);
  }

  /**
   * Get number of scatter traces currently rendered
   * @returns {Promise<number>}
   */
  async getTraceCount() {
    const svgTraces = await this.plotContainer.locator('.scatterlayer .trace').count();
    const glTraces = await this.plotContainer.locator('.gl-container canvas').count();
    return svgTraces + glTraces;
  }

  async getSelectedXAxisValue() {
    return await this.xSelect.inputValue();
  }

  async getSelectedYAxisValue() {
    return await this.ySelect.inputValue();
  }

  async getSelectedGroupValue() {
    return await this.groupSelect.inputValue();
  }

  async selectDifferentXAxis() {
    return await this.#selectDifferentOption(this.xSelect, 'X axis');
  }

  async selectDifferentYAxis() {
    return await this.#selectDifferentOption(this.ySelect, 'Y axis');
  }

  async selectDifferentGroup() {
    return await this.#selectDifferentOption(this.groupSelect, 'Group');
  }

  async #selectDifferentOption(selectLocator, label) {
    const currentValue = await selectLocator.inputValue();
    const availableValues = await selectLocator.locator('option').evaluateAll((options) =>
      options
        .filter((opt) => !opt.disabled)
        .map((opt) => opt.value)
    );

    const nextValue = availableValues.find((value) => value !== currentValue);

    if (!nextValue) {
      throw new Error(`No alternative option available for ${label} dropdown`);
    }

    await selectLocator.selectOption(nextValue);
    return { previous: currentValue, next: nextValue };
  }
}


