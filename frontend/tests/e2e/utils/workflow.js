// frontend/tests/e2e/utils/workflow.js

/**
 * Utility functions for workflow and Drawflow interactions
 */

/**
 * Wait for Vuex store to update with specific condition
 * @param {Page} page - Playwright page object
 * @param {Function} condition - Condition function to check store state
 * @param {number} timeout - Maximum time to wait (default: 5000ms)
 * @returns {Promise<any>} Store state when condition is met
 */
export async function waitForVuexState(page, condition, timeout = 5000) {
  return await page.waitForFunction(
    (conditionString) => {
      const conditionFn = new Function('store', `return ${conditionString}`);
      const store =
        window.app?.$store || window.vue?.$store || window.$store;
      return store && conditionFn(store.state);
    },
    condition.toString(),
    { timeout }
  );
}

/**
 * Get Drawflow node data from the page
 * @param {Page} page - Playwright page object
 * @param {string} nodeId - Node ID to retrieve
 * @returns {Promise<object>} Node data object
 */
export async function getDrawflowNodeData(page, nodeId) {
  return await page.evaluate((id) => {
    const editor = window.editor || window.app?.editor;
    if (!editor) return null;

    const drawflow = editor.drawflow.drawflow;
    const homeData = drawflow?.Home?.data || drawflow?.home?.data;

    return homeData ? homeData[id] : null;
  }, nodeId);
}

/**
 * Get all Drawflow connections
 * @param {Page} page - Playwright page object
 * @returns {Promise<Array>} Array of connection objects
 */
export async function getDrawflowConnections(page) {
  return await page.evaluate(() => {
    const editor = window.editor || window.app?.editor;
    if (!editor) return [];

    const drawflow = editor.drawflow.drawflow;
    const homeData = drawflow?.Home?.data || drawflow?.home?.data;

    const connections = [];
    if (homeData) {
      Object.values(homeData).forEach((node) => {
        if (node.outputs) {
          Object.entries(node.outputs).forEach(([outputKey, output]) => {
            output.connections.forEach((conn) => {
              connections.push({
                from: node.id,
                to: conn.node,
                outputKey,
                inputKey: conn.output,
              });
            });
          });
        }
      });
    }

    return connections;
  });
}

/**
 * Programmatically create a node on the Drawflow canvas
 * @param {Page} page - Playwright page object
 * @param {string} nodeType - Type of node to create
 * @param {object} position - Position on canvas
 * @param {number} position.x - X coordinate
 * @param {number} position.y - Y coordinate
 * @returns {Promise<string>} Created node ID
 */
export async function createDrawflowNode(page, nodeType, position = {}) {
  return await page.evaluate(
    ({ type, pos }) => {
      const editor = window.editor || window.app?.editor;
      if (!editor) throw new Error('Drawflow editor not found');

      const x = pos.x || 200;
      const y = pos.y || 200;

      // This would need to match your actual node creation logic
      // Adjust based on your WorkFlowPage.vue implementation
      const nodeId = editor.addNode(
        type,
        1,
        1,
        x,
        y,
        type.toLowerCase(),
        {},
        type
      );

      return nodeId;
    },
    { type: nodeType, pos: position }
  );
}

/**
 * Programmatically create a connection between two nodes
 * @param {Page} page - Playwright page object
 * @param {string} fromNodeId - Source node ID
 * @param {string} toNodeId - Target node ID
 * @returns {Promise<boolean>} Success status
 */
export async function connectDrawflowNodes(page, fromNodeId, toNodeId) {
  return await page.evaluate(
    ({ from, to }) => {
      const editor = window.editor || window.app?.editor;
      if (!editor) throw new Error('Drawflow editor not found');

      // Assuming single output/input for simplicity
      // Adjust based on your actual connection logic
      try {
        editor.addConnection(from, to, 'output_1', 'input_1');
        return true;
      } catch (error) {
        console.error('Connection failed:', error);
        return false;
      }
    },
    { from: fromNodeId, to: toNodeId }
  );
}

/**
 * Wait for Plotly chart to render
 * @param {Page} page - Playwright page object
 * @param {string} containerId - Plotly container ID (default: 'plotly__scatter')
 * @param {number} timeout - Maximum time to wait (default: 10000ms)
 */
export async function waitForPlotlyRender(
  page,
  containerId = 'plotly__scatter',
  timeout = 10000
) {
  await page.waitForFunction(
    (id) => {
      const container = document.getElementById(id);
      if (!container) return false;

      const svg = container.querySelector('svg.main-svg');
      if (!svg) return false;

      const scatterLayer = container.querySelector('.scatterlayer');
      if (scatterLayer && scatterLayer.children.length > 0) {
        return true;
      }

      const heatmapLayer = container.querySelector('.heatmaplayer, .hm');
      if (heatmapLayer) {
        const hasImage = heatmapLayer.querySelector('image');
        const hasRects = heatmapLayer.querySelector('rect');
        if (hasImage || hasRects || heatmapLayer.children.length > 0) {
          return true;
        }
      }

      const imageLayer = container.querySelector('.layer-above image, .layer-above .imagelayer image');
      if (imageLayer) {
        return true;
      }

      const webglCanvas = container.querySelector('.gl-container canvas');
      if (webglCanvas) {
        return true;
      }

      return false;
    },
    containerId,
    { timeout }
  );
}

/**
 * Wait for vue-good-table to load data
 * @param {Page} page - Playwright page object
 * @param {number} timeout - Maximum time to wait (default: 10000ms)
 */
export async function waitForTableLoad(page, timeout = 10000) {
  // Wait for loading overlay to disappear
  await page
    .locator('.vgt-loading')
    .waitFor({ state: 'hidden', timeout })
    .catch(() => {
      // Loading might not be visible if data loads quickly
    });

  // Wait for table rows to appear
  await page.waitForSelector('.vgt-table tbody tr', { timeout });
}

/**
 * Get vue-good-table data
 * @param {Page} page - Playwright page object
 * @returns {Promise<Array<object>>} Table data rows
 */
export async function getTableData(page) {
  return await page.evaluate(() => {
    const rows = document.querySelectorAll('.vgt-table tbody tr');
    const data = [];

    rows.forEach((row) => {
      const cells = row.querySelectorAll('td');
      const rowData = {};

      cells.forEach((cell, index) => {
        const header = document.querySelector(
          `.vgt-table thead th:nth-child(${index + 1})`
        );
        const headerText = header?.textContent?.trim() || `column_${index}`;
        rowData[headerText] = cell.textContent?.trim();
      });

      data.push(rowData);
    });

    return data;
  });
}

/**
 * Clear localStorage to reset Vuex persistence
 * @param {Page} page - Playwright page object
 */
export async function clearWorkflowStorage(page) {
  await page.evaluate(() => {
    // Clear workflow-related localStorage items
    const keysToRemove = [];
    for (let i = 0; i < localStorage.length; i++) {
      const key = localStorage.key(i);
      if (key && (key.includes('workflow') || key.includes('vuex'))) {
        keysToRemove.push(key);
      }
    }

    keysToRemove.forEach((key) => localStorage.removeItem(key));
  });
}

/**
 * Get workflow metadata from Vuex store
 * @param {Page} page - Playwright page object
 * @returns {Promise<object>} Workflow metadata
 */
export async function getWorkflowMetadata(page) {
  return await page.evaluate(() => {
    const store = window.app?.$store || window.vue?.$store || window.$store;
    if (!store) return null;

    return {
      workflowId: store.state.workflow?.currentWorkflowId,
      workflowInfo: store.state.workflow?.workflow_info,
      nodes: store.state.workflow?.nodes,
    };
  });
}

/**
 * Wait for file propagation to complete (file shared to connected nodes)
 * @param {Page} page - Playwright page object
 * @param {string} fileName - File name to wait for
 * @param {number} timeout - Maximum time to wait (default: 5000ms)
 */
/**
 * Get node file assignment from Vuex store
 * @param {Page} page - Playwright page object
 * @param {string} nodeId - Node ID to check
 * @returns {Promise<string|object|null>} Assigned file(s)
 */
export async function getNodeFileAssignment(page, nodeId) {
  return await page.evaluate((id) => {
    const store = window.app?.$store || window.vue?.$store || window.$store;
    if (!store) return null;

    const workflowInfo = store.state.workflow?.workflow_info;
    if (!workflowInfo) return null;

    const drawflowData =
      workflowInfo.drawflow?.Home?.data ||
      workflowInfo.drawflow?.home?.data ||
      workflowInfo.data ||
      null;

    if (!drawflowData) return null;

    const node = drawflowData[id];
    if (!node) return null;

    // Return file or files object (single string, array, or mapping)
    if (node.data?.file) {
      return node.data.file;
    }

    if (node.data?.files) {
      return node.data.files;
    }

    return null;
  }, nodeId);
}

/**
 * Wait for modal animation to complete
 * @param {Page} page - Playwright page object
 * @param {number} duration - Animation duration (default: 300ms)
 */
export async function waitForModalAnimation(page, duration = 300) {
  await page.waitForTimeout(duration);
}

/**
 * Check if a node exists on the canvas
 * @param {Page} page - Playwright page object
 * @param {string} nodeType - Node title to check (e.g., 'Input h5ad', 'DataTable')
 * @returns {Promise<boolean>} True if node exists
 */
export async function nodeExistsOnCanvas(page, nodeType) {
  // Wait for nodes to be in DOM
  await page.waitForSelector('.drawflow-node', { timeout: 5000 });

  // Additional wait for Vue reactivity
  await page.waitForTimeout(500);

  const nodes = await page.locator('.drawflow-node').all();

  for (const node of nodes) {
    const titleInput = node.locator('input.nodeTitle');
    const count = await titleInput.count();

    if (count > 0) {
      // Wait for input to stabilize
      await titleInput.waitFor({ state: 'attached', timeout: 3000 });

      const inputValue = await titleInput.inputValue();

      if (inputValue.trim() === nodeType) {
        return true;
      }
    }
  }

  return false;
}

/**
 * Get the number of nodes of a specific type on canvas
 * @param {Page} page - Playwright page object
 * @param {string} nodeType - Node title to count (e.g., 'Input h5ad', 'DataTable')
 * @returns {Promise<number>} Count of nodes
 */
export async function countNodesOnCanvas(page, nodeType = null) {
  if (nodeType) {
    // Wait for nodes to be in DOM
    await page.waitForSelector('.drawflow-node', { timeout: 5000 });

    // Additional wait for Vue reactivity
    await page.waitForTimeout(500);

    const nodes = await page.locator('.drawflow-node').all();
    let count = 0;

    for (const node of nodes) {
      const titleInput = node.locator('input.nodeTitle');
      const inputCount = await titleInput.count();

      if (inputCount > 0) {
        // Wait for input to stabilize
        await titleInput.waitFor({ state: 'attached', timeout: 3000 });

        const inputValue = await titleInput.inputValue();

        if (inputValue.trim() === nodeType) {
          count++;
        }
      }
    }

    return count;
  } else {
    return await page.locator('.drawflow-node').count();
  }
}

/**
 * Check if InputFile node exists on canvas (hardcoded for TENET template)
 * @param {Page} page - Playwright page object
 * @returns {Promise<boolean>} True if node exists
 */
export async function inputFileNodeExists(page) {
  try {
    const node = page.locator('#node-7.drawflow-node');
    await node.waitFor({ state: 'attached', timeout: 3000 });
    return true;
  } catch (error) {
    return false;
  }
}
