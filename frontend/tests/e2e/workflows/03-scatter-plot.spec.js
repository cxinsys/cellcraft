// frontend/tests/e2e/workflows/03-scatter-plot.spec.js
import { test, expect } from '../fixtures/auth.js';
import { testWorkflow } from './support/workflow-constants.js';
import { setupPageObjects, setupTestFiles, cleanupTestFiles, createWorkflowWithUniqueTitle, cleanupWorkflows } from './support/workflow-setup.js';
import { getNodeFileAssignment, getWorkflowMetadata } from '../utils/workflow.js';

/**
 * Test Suite: ScatterPlot Node - UMAP Visualization
 *
 * This test suite verifies the ScatterPlot node functionality:
 * - File propagation from InputFile to ScatterPlot
 * - Plotly rendering with UMAP data
 * - Interactive dropdown controls (X axis, Y axis, Group)
 * - Plot re-rendering on parameter changes
 *
 * Success criteria:
 * - File assignment propagates to ScatterPlot node
 * - API call succeeds and returns data
 * - Plotly chart renders without blank state
 * - Trace count is greater than 0
 * - Dropdown changes trigger re-rendering
 */
test.describe('ScatterPlot Node - UMAP Visualization', () => {
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
   * Test: Should render UMAP scatter plot in ScatterPlot node
   *
   * Verifies ScatterPlot functionality:
   * - File assignment to InputFile node
   * - File propagation to ScatterPlot via connections
   * - API call for UMAP data loading
   * - Plotly chart rendering
   * - Interactive dropdown controls
   */
  test('Should render UMAP scatter plot in ScatterPlot node', async ({ page }) => {
    test.setTimeout(60000);

    // Create workflow with unique title
    currentWorkflowTitle = await createWorkflowWithUniqueTitle(
      pageObjects.projectsPage,
      pageObjects.workflowPage,
      testWorkflow,
      createdWorkflows
    );

    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    await pageObjects.workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await pageObjects.inputFileModal.assignFile(testWorkflow.folder, currentTestFileName);

    await pageObjects.workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(300);

    const scatterLocator = await pageObjects.workflowPage.findNodeByType('ScatterPlot');
    const scatterNodeIdAttr = await scatterLocator.getAttribute('id');
    const scatterNodeId = scatterNodeIdAttr?.replace('node-', '') ?? '9';

    const scatterMetadataBefore = await getWorkflowMetadata(page);
    console.log(
      '📦 Vuex drawflow data before ScatterPlot polling:',
      JSON.stringify(scatterMetadataBefore?.workflowInfo?.drawflow?.Home?.data ?? null, null, 2)
    );

    await expect
      .poll(async () => {
        const assignment = await getNodeFileAssignment(page, scatterNodeId);
        if (!assignment) return null;
        if (typeof assignment === 'string') return assignment;
        if (Array.isArray(assignment)) {
          return assignment.find((item) => item === currentTestFileName) ?? null;
        }
        if (typeof assignment === 'object') {
          const values = Object.values(assignment);
          return values.find((item) => item === currentTestFileName) ?? null;
        }
        return null;
      }, {
        message: `Expected ScatterPlot node (${scatterNodeId}) to receive file ${currentTestFileName}`,
        timeout: 10000,
      })
      .toBe(currentTestFileName);

    const scatterMetadataAfter = await getWorkflowMetadata(page);
    console.log(
      '📦 Vuex drawflow data after ScatterPlot polling:',
      JSON.stringify(scatterMetadataAfter?.workflowInfo?.drawflow?.Home?.data ?? null, null, 2)
    );

    const scatterDataResponsePromise = page.waitForResponse(
      (resp) =>
        resp.url().includes('/routes/files/data/') &&
        resp.request().method() === 'GET',
      { timeout: 15000 }
    );

    await pageObjects.workflowPage.openNodeModal('ScatterPlot');
    await pageObjects.scatterPlotModal.verifyModalOpen();

    const scatterDataResponse = await scatterDataResponsePromise;
    const scatterPayload = await scatterDataResponse.json();
    console.log('📡 ScatterPlot API payload:', scatterPayload);

    if (Object.prototype.hasOwnProperty.call(scatterPayload, 'success')) {
      expect(scatterPayload.success).toBeTruthy();
    }

    await pageObjects.scatterPlotModal.waitForPlotly();

    const blankStateVisible = await pageObjects.scatterPlotModal.isBlankStateVisible();
    expect(blankStateVisible).toBeFalsy();

    await expect(page.locator('#plotly__scatter')).toBeVisible();

    const traceCount = await pageObjects.scatterPlotModal.getTraceCount();
    expect(traceCount).toBeGreaterThan(0);

    await test.step('Update ScatterPlot X axis dropdown', async () => {
      const { previous, next } = await pageObjects.scatterPlotModal.selectDifferentXAxis();
      console.log(`🔁 ScatterPlot X axis changed: ${previous} → ${next}`);
      expect(next).not.toBe(previous);
      await pageObjects.scatterPlotModal.waitForPlotly();
      const current = await pageObjects.scatterPlotModal.getSelectedXAxisValue();
      expect(current).toBe(next);
    });

    await test.step('Update ScatterPlot Y axis dropdown', async () => {
      const { previous, next } = await pageObjects.scatterPlotModal.selectDifferentYAxis();
      console.log(`🔁 ScatterPlot Y axis changed: ${previous} → ${next}`);
      expect(next).not.toBe(previous);
      await pageObjects.scatterPlotModal.waitForPlotly();
      const current = await pageObjects.scatterPlotModal.getSelectedYAxisValue();
      expect(current).toBe(next);
    });

    await test.step('Update ScatterPlot group dropdown', async () => {
      const { previous, next } = await pageObjects.scatterPlotModal.selectDifferentGroup();
      console.log(`🔁 ScatterPlot group changed: ${previous} → ${next}`);
      expect(next).not.toBe(previous);
      await pageObjects.scatterPlotModal.waitForPlotly();
      const current = await pageObjects.scatterPlotModal.getSelectedGroupValue();
      expect(current).toBe(next);
      const tracesAfterGroup = await pageObjects.scatterPlotModal.getTraceCount();
      expect(tracesAfterGroup).toBeGreaterThan(0);
    });
  });
});
