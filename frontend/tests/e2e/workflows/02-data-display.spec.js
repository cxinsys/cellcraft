// frontend/tests/e2e/workflows/02-data-display.spec.js
import { test, expect } from '../fixtures/auth.js';
import { testWorkflow } from './support/workflow-constants.js';
import { setupPageObjects, setupTestFiles, cleanupTestFiles, createWorkflowWithUniqueTitle, cleanupWorkflows } from './support/workflow-setup.js';
import { getNodeFileAssignment, getWorkflowMetadata } from '../utils/workflow.js';

/**
 * Test Suite: DataTable Node - Matrix Display
 *
 * This test suite verifies the DataTable node functionality:
 * - File propagation from InputFile to DataTable via Drawflow connections
 * - Backend API response for matrix data loading
 * - Table rendering with proper columns and rows
 * - Empty state handling
 *
 * Success criteria:
 * - File assignment propagates through connections to DataTable node
 * - API call succeeds with valid payload (columns and rows)
 * - Table renders without empty state
 * - Row count is greater than 0
 */
test.describe('DataTable Node - Matrix Display', () => {
  test.describe.configure({ mode: 'serial' });

  let pageObjects;
  const uploadedFiles = [];
  const createdWorkflows = [];
  let currentTestFileName = null;
  let currentWorkflowTitle = null;

  test.beforeEach(async ({ page }) => {
    pageObjects = setupPageObjects(page);
    const { uploadedFileName } = await setupTestFiles(pageObjects.filesPage, testWorkflow, uploadedFiles);
    currentTestFileName = uploadedFileName;
    await pageObjects.projectsPage.goto();
    await pageObjects.projectsPage.verifyPageLoaded();
  });

  test.afterEach(async ({ page }) => {
    await cleanupTestFiles(pageObjects.filesPage, uploadedFiles);
    await cleanupWorkflows(pageObjects.projectsPage, createdWorkflows);
  });

  /**
   * Test: Should display matrix data in DataTable node
   *
   * Verifies DataTable functionality:
   * - File assignment to InputFile node
   * - File propagation to DataTable via connections
   * - API call for data loading
   * - Table rendering with data
   */
  test('Should display matrix data in DataTable node', async ({ page }) => {
    test.setTimeout(60000);

    // Create workflow with unique title
    currentWorkflowTitle = await createWorkflowWithUniqueTitle(
      pageObjects.projectsPage,
      pageObjects.workflowPage,
      testWorkflow,
      createdWorkflows
    );

    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    // Assign file to InputFile node
    await pageObjects.workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await pageObjects.inputFileModal.assignFile(testWorkflow.folder, currentTestFileName);

    // Close InputFile tab to focus on DataTable later
    await pageObjects.workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(300);

    // Determine DataTable node id from DOM (fallback to known default if missing)
    const dataTableLocator = await pageObjects.workflowPage.findNodeByType('DataTable');
    const dataTableNodeIdAttr = await dataTableLocator.getAttribute('id');
    const dataTableNodeId = dataTableNodeIdAttr?.replace('node-', '') ?? '8';

    // Debug: log workflow drawflow data before polling
    const metadataBefore = await getWorkflowMetadata(page);
    console.log(
      '📦 Vuex drawflow data (before polling):',
      JSON.stringify(metadataBefore?.workflowInfo?.drawflow?.Home?.data ?? null, null, 2)
    );

    // Verify file propagation via Vuex store (ensures connection + assignment)
    await expect
      .poll(async () => await getNodeFileAssignment(page, dataTableNodeId), {
        message: `Expected DataTable node (${dataTableNodeId}) to receive file ${currentTestFileName}`,
        timeout: 10000,
      })
      .toBe(currentTestFileName);

    const metadataAfter = await getWorkflowMetadata(page);
    console.log(
      '📦 Vuex drawflow data (after polling):',
      JSON.stringify(metadataAfter?.workflowInfo?.drawflow?.Home?.data ?? null, null, 2)
    );

    // Observe DataTable API call to ensure backend responds successfully
    const dataRequestPromise = page.waitForResponse(
      (resp) =>
        resp.url().includes('/routes/datatable/load_data') &&
        resp.request().method() === 'POST',
      { timeout: 15000 }
    );

    await pageObjects.workflowPage.openNodeModal('DataTable');
    await pageObjects.dataTableModal.verifyModalOpen();

    const dataResponse = await dataRequestPromise;
    const payload = await dataResponse.json();
    console.log('📡 DataTable API payload:', payload);

    if (Object.prototype.hasOwnProperty.call(payload, 'success')) {
      expect(payload.success).toBeTruthy();
    }

    const payloadColumns = Array.isArray(payload.columns) ? payload.columns.length : 0;
    const payloadRows = Array.isArray(payload.rows) ? payload.rows.length : 0;

    expect(payloadColumns).toBeGreaterThan(0);
    expect(payloadRows).toBeGreaterThan(0);

    // Wait for table to render
    await pageObjects.dataTableModal.waitForDataLoaded();

    const emptyStateVisible = await pageObjects.dataTableModal.isEmptyStateVisible();
    expect(emptyStateVisible).toBeFalsy();

    await expect(page.locator('.table-layout')).toBeVisible();

    const rowCount = await pageObjects.dataTableModal.getRowCount();
    expect(rowCount).toBeGreaterThan(0);
  });
});
