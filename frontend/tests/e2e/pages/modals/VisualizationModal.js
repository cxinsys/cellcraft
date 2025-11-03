// frontend/tests/e2e/pages/modals/VisualizationModal.js
import { expect } from '@playwright/test';
import { waitForPlotlyRender } from '../../utils/workflow.js';

/**
 * Page Object Model for the Visualization modal (GRNViz)
 * Handles plugin/script selection, parameter configuration, and Plotly assertions
 */
export class VisualizationModal {
  constructor(page) {
    this.page = page;

    // Layouts
    this.selectionLayout = page.locator('.plugin-selection-layout');
    this.configurationLayout = page.locator('.configuration-layout');

    // Plugin & visualization selectors
    this.pluginItems = page.locator('.plugin-item');
    this.visualizationItems = page.locator('.visualization-item');
    this.continueButton = page.locator('button.continue-button');

    // Configuration elements
    this.backButton = page.locator('.options__header .back-button');
    this.applyButton = page.locator('#apply-button');
    this.plotlyContainer = page.locator('#plotlyChart');

    // Modebar download button (available after Plotly renders)
    this.modebarDownloadButton = page.locator(
      '#plotlyChart .modebar-btn[data-title="Download plot as a png"]'
    );
  }

  async waitForPluginSelection(timeout = 10000) {
    await this.selectionLayout.waitFor({ state: 'visible', timeout });
  }

  async selectPluginByName(pluginName) {
    const plugin = this.pluginItems.filter({ hasText: pluginName });
    await plugin.first().click();
  }

  async selectVisualizationByName(scriptName) {
    const item = this.visualizationItems.filter({ hasText: scriptName });
    await item.first().click();
  }

  async proceedToConfiguration() {
    await this.continueButton.click();
    await this.configurationLayout.waitFor({ state: 'visible', timeout: 10000 });
  }

  async isConfigurationMode() {
    return await this.configurationLayout.isVisible().catch(() => false);
  }

  async returnToPluginSelection() {
    await this.backButton.click();
    await this.selectionLayout.waitFor({ state: 'visible', timeout: 10000 });
  }

  async getSelectedPluginLabel() {
    const text = await this.configurationLayout.locator('.selected-plugin').textContent();
    return text?.trim() ?? '';
  }

  async getSelectedVisualizationLabel() {
    const text = await this.configurationLayout.locator('.selected-script').textContent();
    return text?.trim() ?? '';
  }

  parameterItemByName(parameterName) {
    return this.configurationLayout.locator('.options__item').filter({
      has: this.page.locator('.options__title', { hasText: parameterName }),
    });
  }

  async getInputFileParameterNames() {
    const inputSection = this.configurationLayout.locator('.options__section').filter({
      has: this.page.locator('h4.section-subtitle', { hasText: 'Input Files' }),
    });
    const labels = await inputSection
      .locator('.options__item .options__title')
      .allTextContents()
      .catch(() => []);
    return labels.map((label) => label.trim()).filter(Boolean);
  }

  async getAvailableOptionsForParameter(parameterName) {
    const parameter = this.parameterItemByName(parameterName);
    const optionTexts = await parameter.locator('select option').allTextContents();
    return optionTexts.map((text) => text.trim());
  }

  async setInputFile(parameterName, fileName) {
    const parameter = this.parameterItemByName(parameterName);
    const select = parameter.locator('select');
    await select.selectOption(fileName);
  }

  async selectInputFileOption(parameterName, fileName) {
    const parameter = this.parameterItemByName(parameterName);
    const select = parameter.locator('select');
    await select.waitFor({ state: 'visible', timeout: 10000 });
    await select.selectOption({ label: fileName }).catch(async () => {
      await select.selectOption(fileName);
    });
  }

  async getSelectedInputFile(parameterName) {
    const parameter = this.parameterItemByName(parameterName);
    const select = parameter.locator('select');
    return await select.inputValue();
  }

  async setNumericParameter(parameterName, value) {
    const parameter = this.parameterItemByName(parameterName);
    const input = parameter.locator('input[type="number"]');
    await input.fill('');
    await input.type(String(value));
  }

  async setBooleanParameter(parameterName, checked) {
    const parameter = this.parameterItemByName(parameterName);
    const checkbox = parameter.locator('input[type="checkbox"]');
    await checkbox.setChecked(checked, { force: true });
  }

  async clickExecuteVisualization() {
    await this.applyButton.click();
  }

  async clickShowVisualization() {
    await this.applyButton.click();
  }

  async waitForPlotly(timeout = 15000) {
    await waitForPlotlyRender(this.page, 'plotlyChart', timeout);
  }

  async isApplyButtonEnabled() {
    return await this.applyButton.isEnabled().catch(() => false);
  }

  async getApplyButtonText() {
    const text = await this.applyButton.textContent();
    return text?.trim() ?? '';
  }

  async downloadPlotImage() {
    await expect(this.modebarDownloadButton).toBeVisible();
    const download = this.page.waitForEvent('download');
    await this.modebarDownloadButton.click();
    return await download;
  }

  async isPlotlyVisible() {
    return await this.plotlyContainer.locator('svg.main-svg, canvas').first().isVisible();
  }
}


