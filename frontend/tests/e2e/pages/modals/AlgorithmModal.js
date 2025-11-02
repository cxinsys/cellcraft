// frontend/tests/e2e/pages/modals/AlgorithmModal.js
import { expect } from '@playwright/test';

/**
 * Page Object Model for the Algorithm modal
 * Handles plugin selection, file assignment via dropdowns, and parameter configuration
 */
export class AlgorithmModal {
  constructor(page) {
    this.page = page;

    // Modal container
    this.modal = page.locator('.content-view');

    // Plugin selection
    this.pluginDropdown = page.locator('.plugin-dropdown');
    this.pluginOptions = page.locator('.plugin-dropdown option');

    // Input/Output file dropdowns
    this.inputDropdowns = page.locator('.plugin-input-dropdown');
    this.outputDropdowns = page.locator('.plugin-output-dropdown');

    // Parameter inputs
    this.parameterInputs = page.locator('.parameter__input');
    this.parameterLabels = page.locator('.parameter__label');

    // Checkbox parameters
    this.parameterCheckboxes = page.locator('input[type="checkbox"]');

    // Action buttons
    this.applyButton = page.locator('button:has-text("Apply")');
    this.resetButton = page.locator('button:has-text("Reset")');
  }

  /**
   * Verify the Algorithm modal is open and visible
   */
  async verifyModalOpen() {
    await expect(this.modal).toBeVisible();
  }

  /**
   * Get available plugins from the plugin dropdown
   * @returns {Promise<Array<string>>} Array of plugin names
   */
  async getAvailablePlugins() {
    const options = await this.pluginOptions.all();
    const plugins = [];

    for (const option of options) {
      const text = await option.textContent();
      if (text && text.trim() !== 'Select Plugin') {
        plugins.push(text.trim());
      }
    }

    return plugins;
  }

  /**
   * Select a plugin from the dropdown
   * @param {string} pluginName - Name of the plugin to select
   */
  async selectPlugin(pluginName) {
    await this.pluginDropdown.selectOption({ label: pluginName });
    await this.page.waitForLoadState('networkidle');
  }

  /**
   * Get the currently selected plugin
   * @returns {Promise<string>} Currently selected plugin name
   */
  async getSelectedPlugin() {
    const selectedValue = await this.pluginDropdown.inputValue();
    return selectedValue;
  }

  /**
   * Get input file dropdown by index or label
   * @param {number|string} identifier - Index (0-based) or label text
   * @returns {Promise<Locator>}
   */
  async getInputDropdown(identifier = 0) {
    if (typeof identifier === 'number') {
      const dropdowns = await this.inputDropdowns.all();
      return dropdowns[identifier];
    } else {
      // Find by label
      return this.page.locator('.plugin-input-dropdown', {
        has: this.page.locator(`label:has-text("${identifier}")`),
      });
    }
  }

  /**
   * Select a file from an input dropdown
   * @param {number|string} dropdownIdentifier - Dropdown index or label
   * @param {string} fileName - File name to select
   */
  async selectInputFile(dropdownIdentifier, fileName) {
    const dropdown = await this.getInputDropdown(dropdownIdentifier);
    await dropdown.selectOption({ label: fileName });
    await this.page.waitForTimeout(300);
  }

  /**
   * Get the selected file from an input dropdown
   * @param {number|string} dropdownIdentifier - Dropdown index or label
   * @returns {Promise<string>} Selected file name
   */
  async getSelectedInputFile(dropdownIdentifier) {
    const dropdown = await this.getInputDropdown(dropdownIdentifier);
    const selectedValue = await dropdown.inputValue();
    return selectedValue;
  }

  /**
   * Verify a file is available in the input dropdown
   * @param {number|string} dropdownIdentifier - Dropdown index or label
   * @param {string} fileName - File name to verify
   */
  async verifyFileInDropdown(dropdownIdentifier, fileName) {
    const dropdown = await this.getInputDropdown(dropdownIdentifier);
    const options = await dropdown.locator('option').allTextContents();
    const fileExists = options.some((opt) =>
      opt.trim().includes(fileName.trim())
    );

    expect(fileExists).toBeTruthy();
  }

  /**
   * Get all available files in an input dropdown
   * @param {number|string} dropdownIdentifier - Dropdown index or label
   * @returns {Promise<Array<string>>} Array of available file names
   */
  async getAvailableFiles(dropdownIdentifier) {
    const dropdown = await this.getInputDropdown(dropdownIdentifier);
    const options = await dropdown.locator('option').all();
    const files = [];

    for (const option of options) {
      const text = await option.textContent();
      if (text && text.trim() !== 'Select File' && text.trim() !== '') {
        files.push(text.trim());
      }
    }

    return files;
  }

  /**
   * Get parameter input by name or index
   * @param {string|number} identifier - Parameter name or index
   * @returns {Promise<Locator>}
   */
  async getParameterInput(identifier) {
    if (typeof identifier === 'number') {
      const inputs = await this.parameterInputs.all();
      return inputs[identifier];
    } else {
      // Find by label
      return this.page.locator('.parameter__input', {
        has: this.page.locator(`label:has-text("${identifier}")`),
      });
    }
  }

  /**
   * Set a parameter value
   * @param {string|number} parameterIdentifier - Parameter name or index
   * @param {string|number|boolean} value - Value to set
   */
  async setParameter(parameterIdentifier, value) {
    const input = await this.getParameterInput(parameterIdentifier);
    const inputType = await input.getAttribute('type');

    if (inputType === 'checkbox') {
      const isChecked = await input.isChecked();
      if ((value === true && !isChecked) || (value === false && isChecked)) {
        await input.click();
      }
    } else {
      await input.fill(String(value));
      await input.blur(); // Trigger change event
    }

    await this.page.waitForTimeout(300);
  }

  /**
   * Get a parameter value
   * @param {string|number} parameterIdentifier - Parameter name or index
   * @returns {Promise<string|boolean>} Parameter value
   */
  async getParameterValue(parameterIdentifier) {
    const input = await this.getParameterInput(parameterIdentifier);
    const inputType = await input.getAttribute('type');

    if (inputType === 'checkbox') {
      return await input.isChecked();
    } else {
      return await input.inputValue();
    }
  }

  /**
   * Verify a parameter has a specific default value
   * @param {string|number} parameterIdentifier - Parameter name or index
   * @param {string|number|boolean} expectedValue - Expected default value
   */
  async verifyDefaultValue(parameterIdentifier, expectedValue) {
    const actualValue = await this.getParameterValue(parameterIdentifier);

    if (typeof expectedValue === 'boolean') {
      expect(actualValue).toBe(expectedValue);
    } else {
      expect(actualValue).toBe(String(expectedValue));
    }
  }

  /**
   * Get all parameter names
   * @returns {Promise<Array<string>>} Array of parameter names
   */
  async getParameterNames() {
    const labels = await this.parameterLabels.all();
    const names = [];

    for (const label of labels) {
      const text = await label.textContent();
      names.push(text?.trim());
    }

    return names;
  }

  /**
   * Click the Apply button
   */
  async clickApply() {
    await this.applyButton.click();
    await this.page.waitForLoadState('networkidle');
  }

  /**
   * Click the Reset button
   */
  async clickReset() {
    await this.resetButton.click();
    await this.page.waitForTimeout(300);
  }

  /**
   * Verify the dropdown is enabled (not disabled)
   * @param {number|string} dropdownIdentifier - Dropdown identifier
   */
  async verifyDropdownEnabled(dropdownIdentifier) {
    const dropdown = await this.getInputDropdown(dropdownIdentifier);
    await expect(dropdown).toBeEnabled();
  }

  /**
   * Verify the dropdown is disabled
   * @param {number|string} dropdownIdentifier - Dropdown identifier
   */
  async verifyDropdownDisabled(dropdownIdentifier) {
    const dropdown = await this.getInputDropdown(dropdownIdentifier);
    await expect(dropdown).toBeDisabled();
  }

  /**
   * Get the count of input file dropdowns
   * @returns {Promise<number>}
   */
  async getInputDropdownCount() {
    return await this.inputDropdowns.count();
  }

  /**
   * Get the count of parameters
   * @returns {Promise<number>}
   */
  async getParameterCount() {
    return await this.parameterInputs.count();
  }
}
