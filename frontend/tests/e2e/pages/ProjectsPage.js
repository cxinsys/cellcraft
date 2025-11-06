// frontend/tests/e2e/pages/ProjectsPage.js
import { expect } from '@playwright/test';

/**
 * Page Object Model for the Projects page
 */
export class ProjectsPage {
  constructor(page) {
    this.page = page;

    // Locators
    this.newWorkflowCard = page.locator('.project-component', {
      has: page.locator('.project__title:has-text("New Workflow")'),
    });
    this.projectsList = page.locator('.project-view__list');
    this.pluginSelectModal = page.locator('.plugin-select-modal');
    this.deleteModal = page.locator('.delete-modal');
    this.contextMenu = page.locator('ul.toggle__menu');
    this.messageBox = page.locator('.message');
  }

  /**
   * Navigate to the Projects page
   */
  async goto() {
    await this.page.goto('/projects');
    await this.page.waitForLoadState('networkidle');
  }

  /**
   * Click the "New Workflow" card to start creating a workflow
   */
  async clickNewWorkflow() {
    await this.newWorkflowCard.click();
    await this.pluginSelectModal.waitFor({ state: 'visible', timeout: 5000 });
  }

  /**
   * Select a plugin template from the modal
   * @param {string} pluginName - Name of the plugin (e.g., 'TENET', 'GENIE3')
   */
  async selectPluginTemplate(pluginName) {
    // Wait for plugin modal to be visible
    await expect(this.pluginSelectModal).toBeVisible();

    // Scope the plugin item search to within the modal only
    // This prevents accidentally clicking elements behind the modal (e.g., header links)
    // Use exact text matching to avoid partial matches (e.g., "TENET" vs "FastTENET")
    const pluginItem = this.pluginSelectModal
      .locator('.plugin-item')
      .filter({
        has: this.page.locator('.plugin-name').getByText(pluginName, { exact: true })
      });

    await pluginItem.click();

    // Wait for redirect to workflow page
    await this.page.waitForURL('**/workflow**', { timeout: 10000 });
  }

  /**
   * Create a workflow with default (blank) template
   */
  async selectDefaultTemplate() {
    await expect(this.pluginSelectModal).toBeVisible();
    await this.page.locator('.plugin-item.default-option').click();
    await this.page.waitForURL('**/workflow**', { timeout: 10000 });
  }

  /**
   * Close the plugin selection modal
   */
  async closePluginModal() {
    await this.page.locator('.plugin-select-modal .close-button').click();
    await this.pluginSelectModal.waitFor({ state: 'hidden', timeout: 3000 });
  }

  /**
   * Get the list of available plugin templates
   * @returns {Promise<Array<{name: string, source: string, description: string}>>} Array of plugin objects
   */
  async getAvailablePlugins() {
    await expect(this.pluginSelectModal).toBeVisible();

    const pluginItems = await this.page.locator('.plugin-item:not(.default-option)').all();

    const plugins = [];
    for (const item of pluginItems) {
      const name = await item.locator('.plugin-name').textContent();
      const badge = await item.locator('.plugin-badge').textContent();
      const description = await item.locator('.plugin-desc').textContent();

      plugins.push({
        name: name.trim(),
        source: badge.trim(),
        description: description.trim(),
      });
    }

    return plugins;
  }

  /**
   * Verify a specific plugin template is available
   * @param {string} pluginName - Name of the plugin
   */
  async verifyPluginAvailable(pluginName) {
    await expect(this.pluginSelectModal).toBeVisible();
    await expect(
      this.page.locator('.plugin-name', { hasText: pluginName })
    ).toBeVisible();
  }

  /**
   * Get the list of existing workflows
   * @returns {Promise<Array<{title: string, date: string}>>} Array of workflow objects
   */
  async getWorkflows() {
    const workflowCards = await this.page
      .locator('.project-component:not(:has(.project__title:has-text("New Workflow")))')
      .all();

    const workflows = [];
    for (const card of workflowCards) {
      const title = await card.locator('.project__title').textContent();
      const date = await card.locator('.project__date').textContent();

      workflows.push({
        title: title.trim(),
        date: date.trim(),
      });
    }

    return workflows;
  }

  /**
   * Open an existing workflow by title
   * @param {string} workflowTitle - Title of the workflow to open
   */
  async openWorkflow(workflowTitle) {
    const workflowCard = this.page.locator('.project-component', {
      has: this.page.locator(`.project__title:has-text("${workflowTitle}")`),
    });

    await workflowCard.click();
    await this.page.waitForURL('**/workflow**', { timeout: 10000 });
  }

  /**
   * Open the most recently updated workflow (first workflow card after "New Workflow")
   * @returns {Promise<{ title: string }>} Recently opened workflow info
   */
  async openMostRecentWorkflow() {
    const workflowCards = this.page.locator('.project-component');
    const count = await workflowCards.count();

    if (count <= 1) {
      throw new Error('No existing workflows available to open.');
    }

    const mostRecentCard = workflowCards.nth(1);
    const titleText = await mostRecentCard.locator('.project__title').textContent();

    await mostRecentCard.click();
    await this.page.waitForURL('**/workflow**', { timeout: 10000 });

    return {
      title: titleText?.trim() ?? '',
    };
  }

  /**
   * Delete a workflow by title
   * @param {string} workflowTitle - Title of the workflow to delete
   */
  async deleteWorkflow(workflowTitle) {
    const workflowCard = this.page.locator('.project-component', {
      has: this.page.locator(`.project__title:has-text("${workflowTitle}")`),
    });

    // Right-click to open context menu
    await workflowCard.click({ button: 'right' });

    // Wait for context menu to appear
    await expect(this.contextMenu).toBeVisible();

    // Click delete option
    await this.contextMenu.locator('li:has-text("Delete")').click();

    // Confirm deletion in modal
    await expect(this.deleteModal).toBeVisible();
    await this.page
      .locator('.delete-modal__content button:has-text("Yes")')
      .click();

    // Wait for deletion to complete
    await this.page.waitForLoadState('networkidle');
  }

  /**
   * Verify the Projects page is displayed
   */
  async verifyPageLoaded() {
    await expect(
      this.page.locator('.profile__title:has-text("Projects")')
    ).toBeVisible();
    await expect(this.projectsList).toBeVisible();
    await expect(this.newWorkflowCard).toBeVisible();
  }

  /**
   * Verify a workflow exists in the list
   * @param {string} workflowTitle - Title of the workflow
   */
  async verifyWorkflowExists(workflowTitle) {
    await expect(
      this.page.locator('.project__title', { hasText: workflowTitle })
    ).toBeVisible();
  }

  /**
   * Verify a workflow does not exist in the list
   * @param {string} workflowTitle - Title of the workflow
   */
  async verifyWorkflowNotExists(workflowTitle) {
    await expect(
      this.page.locator('.project__title', { hasText: workflowTitle })
    ).not.toBeVisible();
  }

  /**
   * Get the user profile information displayed on the page
   * @returns {Promise<{username: string, email: string}>} User profile object
   */
  async getUserProfile() {
    const username = await this.page.locator('.profile__username').textContent();
    const email = await this.page.locator('.profile__email').textContent();

    return {
      username: username.trim(),
      email: email.trim(),
    };
  }
}
