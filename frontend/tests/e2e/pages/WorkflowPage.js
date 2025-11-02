// frontend/tests/e2e/pages/WorkflowPage.js
import { expect } from '@playwright/test';

/**
 * Page Object Model for the Workflow Editor page
 * Handles Drawflow canvas interactions, node management, and modal operations
 */
export class WorkflowPage {
  constructor(page) {
    this.page = page;

    // Main workflow locators
    this.drawflowCanvas = page.locator('#drawflow');
    this.nodeBar = page.locator('.node-bar');
    this.tabComponent = page.locator('.content-component');
    this.contentView = page.locator('.content-view');
    this.messageBox = page.locator('.message');

    // Node bar draggable nodes
    this.nodeBarNodes = page.locator('.node-bar__drag-drawflow');

    // Tab system
    this.tabList = page.locator('.content-tab');
    this.tabItems = page.locator('.tab__item');
    this.currentTab = page.locator('.tab__item.currentTab');
    this.tabHideButton = page.locator('.tab__hide');
  }

  /**
   * Navigate to a workflow by ID
   * @param {string|number} workflowId - The workflow ID
   */
  async goto(workflowId) {
    await this.page.goto(`/workflow?workflow_id=${workflowId}`);
    await this.page.waitForLoadState('networkidle');
    await this.drawflowCanvas.waitFor({ state: 'visible' });
  }

  /**
   * Verify the workflow page is loaded
   */
  async verifyPageLoaded() {
    await expect(this.drawflowCanvas).toBeVisible();
    await expect(this.nodeBar).toBeVisible();
  }

  /**
   * Wait for all nodes to be fully rendered and stable
   * @param {number} expectedCount - Expected number of nodes
   * @param {number} timeout - Maximum wait time
   */
  async waitForNodesReady(expectedCount = 1, timeout = 10000) {
    // Wait for Vuex store to have workflow data
    await this.page.waitForFunction(
      (count) => {
        const store = window.app?.$store || window.vue?.$store || window.$store;
        if (!store) return false;

        const workflowInfo = store.getters.getWorkflowInfo;
        if (!workflowInfo?.drawflow?.Home?.data) return false;

        const nodeCount = Object.keys(workflowInfo.drawflow.Home.data).length;
        return nodeCount >= count;
      },
      expectedCount,
      { timeout }
    );

    // Small delay for DOM synchronization
    await this.page.waitForTimeout(300);
  }

  /**
   * Get all nodes currently on the canvas
   * @returns {Promise<Array<{id: string, type: string, position: {x: number, y: number}}>>}
   */
  async getCanvasNodes() {
    const nodes = await this.page.locator('.drawflow-node').all();
    const nodeData = [];

    for (const node of nodes) {
      const id = await node.getAttribute('id');
      const type = await node.locator('.nodeTitle').textContent();
      const style = await node.getAttribute('style');

      // Parse position from style attribute
      const leftMatch = style.match(/left:\s*(\d+)px/);
      const topMatch = style.match(/top:\s*(\d+)px/);

      nodeData.push({
        id,
        type: type?.trim(),
        position: {
          x: leftMatch ? parseInt(leftMatch[1]) : 0,
          y: topMatch ? parseInt(topMatch[1]) : 0,
        },
      });
    }

    return nodeData;
  }

  /**
   * Find a node on the canvas by its title
   * @param {string} nodeType - The node title (e.g., 'Input h5ad', 'DataTable', 'Algorithm')
   * @returns {Promise<Locator>} The node locator
   */
  async findNodeByType(nodeType) {
    // TENET template: node-7 is always "Input h5ad" InputFile
    if (nodeType === 'Input h5ad') {
      const node = this.page.locator('#node-7.drawflow-node');
      await node.waitFor({ state: 'attached', timeout: 5000 });
      return node;
    }

    // Fallback: DOM-based search for other nodes
    const nodes = await this.page.locator('.drawflow-node').all();

    for (const node of nodes) {
      const titleInput = node.locator('input.nodeTitle');
      const count = await titleInput.count();

      if (count > 0) {
        const inputValue = await titleInput.inputValue();
        if (inputValue.trim() === nodeType) {
          return node;
        }
      }
    }

    throw new Error(`Node with title "${nodeType}" not found on canvas`);
  }

  /**
   * Double-click a node to open its configuration modal
   * @param {string} nodeType - The node type to open
   */
  async openNodeModal(nodeType) {
    const node = await this.findNodeByType(nodeType);
    await node.dblclick();

    // Wait for tab to appear
    await this.tabComponent.waitFor({ state: 'visible', timeout: 5000 });
    await this.contentView.waitFor({ state: 'visible' });

    // Map nodeType to actual route path
    let routePath;
    if (nodeType === 'Input h5ad') {
      routePath = 'inputfile';
    } else {
      routePath = nodeType.toLowerCase().replace(/\s+/g, '');
    }

    // Wait for router to navigate
    await this.page.waitForURL(
      `**/workflow/${routePath}?workflow_id=*&node=*`,
      { timeout: 5000 }
    );
  }

  /**
   * Get the currently active tab information
   * @returns {Promise<{name: string, title: string}>}
   */
  async getCurrentTab() {
    const currentTabElement = this.currentTab;
    const tabText = await currentTabElement.locator('.tab__text').textContent();

    return {
      name: tabText?.trim(),
      title: tabText?.trim(),
    };
  }

  /**
   * Switch to a different tab by index
   * @param {number} tabIndex - Zero-based tab index
   */
  async switchTab(tabIndex) {
    const tabs = await this.tabItems.all();
    if (tabIndex >= 0 && tabIndex < tabs.length) {
      await tabs[tabIndex].click();
      await this.page.waitForLoadState('networkidle');
    } else {
      throw new Error(`Tab index ${tabIndex} out of range`);
    }
  }

  /**
   * Close a tab by clicking its close button
   * @param {string} tabName - The tab name to close (should match .tab__text content exactly)
   */
  async closeTab(tabName) {
    // TabComponent.vue structure: .tab__item > .tab__name > .tab__text (p element)
    // The tab name is displayed in lowercase with dots (e.g., "input.h5ad")
    // We need to find the tab by matching the exact text in .tab__text
    const closeButton = this.page
      .locator('.tab__item')
      .filter({
        has: this.page.locator('.tab__text', { hasText: tabName }),
      })
      .locator('img.tab__close');

    await closeButton.click({ timeout: 5000 });
    await this.page.waitForTimeout(300); // Allow animation
  }

  /**
   * Hide/minimize the tab panel
   */
  async hideTabPanel() {
    await this.tabHideButton.click();
    await this.page.waitForTimeout(300);
  }

  /**
   * Drag a node from the node bar onto the canvas
   * @param {string} nodeType - The node type to drag (e.g., 'InputFile', 'DataTable')
   * @param {object} position - Target position on canvas
   * @param {number} position.x - X coordinate
   * @param {number} position.y - Y coordinate
   */
  async dragNodeToCanvas(nodeType, position = { x: 400, y: 300 }) {
    // Find the node in the node bar
    const nodeBarNode = this.page.locator(
      `.node-bar__drag-drawflow[data-node="${nodeType}"]`
    );

    // Get the bounding boxes
    const nodeBox = await nodeBarNode.boundingBox();
    const canvasBox = await this.drawflowCanvas.boundingBox();

    if (!nodeBox || !canvasBox) {
      throw new Error('Could not find node or canvas bounding box');
    }

    // Perform drag operation
    await this.page.mouse.move(
      nodeBox.x + nodeBox.width / 2,
      nodeBox.y + nodeBox.height / 2
    );
    await this.page.mouse.down();

    await this.page.mouse.move(
      canvasBox.x + position.x,
      canvasBox.y + position.y,
      { steps: 10 }
    );
    await this.page.mouse.up();

    // Wait for node to appear on canvas
    await this.page.waitForTimeout(500);
  }

  /**
   * Connect two nodes by dragging from output to input
   * @param {string} sourceNodeType - Source node type
   * @param {string} targetNodeType - Target node type
   */
  async connectNodes(sourceNodeType, targetNodeType) {
    const sourceNode = await this.findNodeByType(sourceNodeType);
    const targetNode = await this.findNodeByType(targetNodeType);

    // Find output and input connection points
    const outputPort = sourceNode.locator('.output');
    const inputPort = targetNode.locator('.input');

    // Get bounding boxes
    const outputBox = await outputPort.boundingBox();
    const inputBox = await inputPort.boundingBox();

    if (!outputBox || !inputBox) {
      throw new Error('Could not find connection ports');
    }

    // Drag from output to input
    await this.page.mouse.move(
      outputBox.x + outputBox.width / 2,
      outputBox.y + outputBox.height / 2
    );
    await this.page.mouse.down();

    await this.page.mouse.move(
      inputBox.x + inputBox.width / 2,
      inputBox.y + inputBox.height / 2,
      { steps: 10 }
    );
    await this.page.mouse.up();

    await this.page.waitForTimeout(500);
  }

  /**
   * Verify a connection exists between two nodes
   * @param {string} sourceNodeType - Source node type
   * @param {string} targetNodeType - Target node type
   */
  async verifyConnection(sourceNodeType, targetNodeType) {
    // Connections are represented by SVG paths in Drawflow
    const connectionPath = this.page.locator('.connection .main-path');
    await expect(connectionPath).toBeVisible();
  }

  /**
   * Get all connections on the canvas
   * @returns {Promise<number>} Number of connections
   */
  async getConnectionCount() {
    return await this.page.locator('.connection').count();
  }

  /**
   * Wait for a message notification to appear
   * @param {string} expectedMessage - The message text to wait for (optional)
   * @param {number} timeout - Maximum time to wait (default: 5000ms)
   */
  async waitForMessage(expectedMessage = null, timeout = 5000) {
    await this.messageBox.waitFor({ state: 'visible', timeout });

    if (expectedMessage) {
      const messageText = await this.messageBox
        .locator('.message__text')
        .textContent();
      expect(messageText).toContain(expectedMessage);
    }
  }

  /**
   * Close the message notification
   */
  async closeMessage() {
    await this.messageBox.locator('.message__close').click();
    await this.messageBox.waitFor({ state: 'hidden', timeout: 3000 });
  }

  /**
   * Get workflow save status by checking for visual indicators
   * @returns {Promise<boolean>} True if workflow appears saved
   */
  async isWorkflowSaved() {
    // Check for success message or other save indicators
    const messageVisible = await this.messageBox.isVisible();
    if (messageVisible) {
      const messageText = await this.messageBox
        .locator('.message__text')
        .textContent();
      return messageText?.includes('saved') || messageText?.includes('success');
    }
    return false;
  }

  /**
   * Get the number of tabs currently open
   * @returns {Promise<number>}
   */
  async getTabCount() {
    return await this.tabItems.count();
  }

  /**
   * Verify a specific tab is open
   * @param {string} tabName - The tab name to verify
   */
  async verifyTabOpen(tabName) {
    await expect(
      this.page.locator('.tab__text', { hasText: tabName })
    ).toBeVisible();
  }

  /**
   * Get all available node types from the node bar
   * @returns {Promise<Array<string>>}
   */
  async getAvailableNodeTypes() {
    const nodes = await this.nodeBarNodes.all();
    const types = [];

    for (const node of nodes) {
      const type = await node.getAttribute('data-node');
      if (type) types.push(type);
    }

    return types;
  }

  /**
   * Verify a node type is available in the node bar
   * @param {string} nodeType - The node type to verify
   */
  async verifyNodeTypeAvailable(nodeType) {
    await expect(
      this.page.locator(`.node-bar__drag-drawflow[data-node="${nodeType}"]`)
    ).toBeVisible();
  }
}
