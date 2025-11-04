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
    this.pluginLogo = page.locator('.algorithm-logo');

    // Input/Output file dropdowns
    this.inputDropdowns = page.locator('.plugin-input-dropdown');
    this.outputDropdowns = page.locator('.plugin-output-dropdown');

    // Parameter inputs
    this.parameterInputs = page.locator('.parameter__input');
    this.parameterLabels = page.locator('.parameter__label');
    this.parameterContainers = page.locator('.parameters');

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
   * Get the plugin logo text rendered in the modal header
   * @returns {Promise<string>}
   */
  async getPluginLogoText() {
    const text = await this.pluginLogo.textContent();
    return text?.trim() ?? '';
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
    const dropdownCount = await this.pluginDropdown.count();
    if (dropdownCount === 0) {
      return;
    }

    const dropdown = this.pluginDropdown.first();
    await dropdown.waitFor({ state: 'visible', timeout: 10000 });
    await dropdown.selectOption({ label: pluginName }).catch(async () => {
      // fallback: try value-based selection if label is missing
      await dropdown.selectOption(pluginName).catch(() => null);
    });
    await this.page.waitForTimeout(300);
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
    const selectedOption = dropdown.locator('option:checked');
    const text = await selectedOption.textContent();
    return text?.trim() ?? '';
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
   * Get a parameter container by name
   * @param {string} parameterName
   * @returns {Locator}
   */
  _getParameterContainer(parameterName) {
    return this.parameterContainers.filter({
      has: this.page.locator('.parameter-id', { hasText: parameterName }),
    }).first();
  }

  /**
   * Internal helper to set cluster selections using the checkbox dropdown
   * @param {Locator} container - Parameter container locator
   * @param {Array<string>} clusterNames - Cluster names to select
   */
  async _setClusterSelection(container, clusterNames) {
    const dropdown = container.locator('.parameter__dropdown--checkbox');

    if (!(await dropdown.count())) {
      throw new Error('Cluster dropdown not found for parameter');
    }

    // Ensure dropdown is active so that checkboxes are interactable
    await dropdown.click();

    const menu = dropdown.locator('.parameter__dropdown--menu');
    const checkboxes = menu.locator('input[type="checkbox"]');

    await checkboxes.first().waitFor({ state: 'visible', timeout: 10000 });

    const desiredSet = new Set((clusterNames || []).map((name) => name.trim()));
    const checkboxCount = await checkboxes.count();

    for (let index = 0; index < checkboxCount; index += 1) {
      const isActive = await dropdown.evaluate((node) => node.classList.contains('isactive'));
      if (!isActive) {
        await dropdown.click();
        await this.page.waitForTimeout(100);
      }

      const checkbox = checkboxes.nth(index);
      const optionName = await checkbox.evaluate(
        (node) => node.name || node.getAttribute('name') || node.value || ''
      );
      const shouldSelect = desiredSet.has(optionName.trim());
      await checkbox.setChecked(shouldSelect);
    }

    // Collapse dropdown to avoid covering other controls
    await dropdown.click();
    await this.page.waitForTimeout(300);
  }

  /**
   * Retrieve the current parameter value from the UI
   * @param {string} parameterName
   * @returns {Promise<string|boolean|null>}
   */
  async getParameterValueByName(parameterName) {
    const container = this._getParameterContainer(parameterName);
    await expect(container, `Parameter container for "${parameterName}" not found`).toBeVisible();

    const numberInput = container.locator('input[type="number"]');
    if (await numberInput.count()) {
      return await numberInput.first().inputValue();
    }

    const textInput = container.locator('input[type="text"]:not(.umap-status-input)');
    if (await textInput.count()) {
      return await textInput.first().inputValue();
    }

    const clusterDropdown = container.locator('.parameter__dropdown--checkbox');
    if (await clusterDropdown.count()) {
      const selectedClusters = await clusterDropdown
        .locator('input[type="checkbox"]')
        .evaluateAll((inputs) =>
          inputs
            .filter((input) => input.checked)
            .map((input) => input.name || input.getAttribute('name') || input.value)
        );
      return selectedClusters;
    }

    const checkbox = container.locator('input[type="checkbox"]');
    if (await checkbox.count()) {
      return await checkbox.first().isChecked();
    }

    const select = container.locator('select');
    if (await select.count()) {
      return await select.first().inputValue();
    }

    return null;
  }

  /**
   * Set a parameter value by name, handling different input types
   * @param {string} parameterName
   * @param {string|number|boolean} value
   */
  async setParameterValueByName(parameterName, value) {
    const container = this._getParameterContainer(parameterName);
    await expect(container, `Parameter container for "${parameterName}" not found`).toBeVisible();

    const numberInput = container.locator('input[type="number"]');
    if (await numberInput.count()) {
      await numberInput.first().fill(String(value));
      await numberInput.first().blur();
      await this.page.waitForTimeout(300);
      return;
    }

    const textInput = container.locator('input[type="text"]:not(.umap-status-input)');
    if (await textInput.count()) {
      await textInput.first().fill(String(value));
      await textInput.first().blur();
      await this.page.waitForTimeout(300);
      return;
    }

    const select = container.locator('select');
    if (await select.count()) {
      await select.first().selectOption(String(value));
      await this.page.waitForTimeout(300);
      return;
    }

    const clusterDropdown = container.locator('.parameter__dropdown--checkbox');
    if (await clusterDropdown.count()) {
      const values = Array.isArray(value) ? value : [value];
      await this._setClusterSelection(container, values.filter(Boolean));
      return;
    }

    const checkbox = container.locator('input[type="checkbox"]');
    if (await checkbox.count()) {
      const shouldCheck = Boolean(value);
      const target = checkbox.first();
      const isChecked = await target.isChecked();
      if (shouldCheck !== isChecked) {
        await target.click();
        await this.page.waitForTimeout(200);
      }
      return;
    }

    throw new Error(`Unsupported parameter input type for "${parameterName}"`);
  }

  /**
   * Retrieve available dropdown options for a parameter
   * @param {string} parameterName
   * @returns {Promise<Array<string>>}
   */
  async getParameterDropdownOptions(parameterName) {
    const container = this._getParameterContainer(parameterName);
    const select = container.locator('select');
    if (!(await select.count())) {
      return [];
    }

    return await select
      .first()
      .locator('option')
      .evaluateAll((options) =>
        options
          .filter((opt) => !opt.disabled)
          .map((opt) => opt.value)
      );
  }

  /**
   * Select specific option from parameter dropdown, or the first different option if none provided
   * @param {string} parameterName
   * @param {string|null} optionValue
   * @returns {Promise<{ previous: string|null, next: string }>} previous and selected values
   */
  async selectParameterDropdownOption(parameterName, optionValue = null) {
    const container = this._getParameterContainer(parameterName);
    const select = container.locator('select');
    await expect(select, `Dropdown for parameter "${parameterName}" not found`).toBeVisible();

    const currentValue = await select.first().inputValue();
    const options = await this.getParameterDropdownOptions(parameterName);

    let nextValue = optionValue;
    if (!nextValue) {
      nextValue = options.find((opt) => opt !== currentValue && opt !== '') ?? currentValue;
    }

    if (!nextValue) {
      throw new Error(`No selectable option available for "${parameterName}" dropdown`);
    }

    if (nextValue === currentValue) {
      // Nothing to change
      return { previous: currentValue, next: currentValue };
    }

    await select.first().selectOption(nextValue);
    await this.page.waitForTimeout(300);
    return { previous: currentValue, next: nextValue };
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
