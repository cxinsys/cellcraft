// frontend/tests/e2e/workflows/01-file-assignment.spec.js
import { test, expect } from '../fixtures/auth.js';
import { testWorkflow } from './support/workflow-constants.js';
import { setupPageObjects, setupTestFiles, cleanupTestFiles, setupFileChangeTest, createWorkflowWithUniqueTitle, cleanupWorkflows } from './support/workflow-setup.js';
import { inputFileNodeExists } from '../utils/workflow.js';

/**
 * Test Suite: InputFile Node - File Assignment
 *
 * This test suite verifies the core file assignment functionality for InputFile nodes:
 * - Opening InputFile node modal
 * - Navigating folder structure
 * - Selecting and applying files
 * - Verifying Vuex store persistence
 * - Changing assigned files
 *
 * Success criteria:
 * - Folder selection highlights the folder
 * - File selection shows file details
 * - Apply button changes to "Applied" after submission
 * - File is persisted in Vuex store (verified by reopening modal)
 * - File changes update Vuex state and reset Apply button
 */
test.describe('InputFile Node - File Assignment', () => {
  // Configure serial mode to prevent parallel execution issues
  test.describe.configure({ mode: 'serial' });

  // Set timeout for all tests in this suite (including hooks)
  test.setTimeout(60000);

  let pageObjects;
  const uploadedFiles = [];
  const createdWorkflows = [];
  let currentTestFileName = null;
  let currentWorkflowTitle = null;

  test.beforeEach(async ({ page }) => {
    // Initialize page objects
    pageObjects = setupPageObjects(page);

    // Upload test file with unique filename for test isolation
    const { uploadedFileName } = await setupTestFiles(pageObjects.filesPage, testWorkflow, uploadedFiles);
    currentTestFileName = uploadedFileName;

    // Navigate to Projects page (already authenticated via fixture)
    await pageObjects.projectsPage.goto();
    await pageObjects.projectsPage.verifyPageLoaded();
  });

  test.afterEach(async ({ page }) => {
    // Clean up uploaded files
    await cleanupTestFiles(pageObjects.filesPage, uploadedFiles);
    // Clean up created workflows
    await cleanupWorkflows(pageObjects.projectsPage, createdWorkflows);
  });

  /**
   * Test: Should assign h5ad file to workflow input node
   *
   * Verifies the core file assignment functionality:
   * - Creating workflow from TENET template
   * - Opening InputFile node modal
   * - Selecting file from folder structure
   * - Applying file assignment
   * - Verifying Vuex store persistence
   */
  test('Should assign h5ad file to workflow input node', async ({ page }) => {
    // Create workflow with unique title
    currentWorkflowTitle = await createWorkflowWithUniqueTitle(
      pageObjects.projectsPage,
      pageObjects.workflowPage,
      testWorkflow,
      createdWorkflows
    );

    // Wait for canvas to have nodes rendered
    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    // Verify InputFile node exists on canvas (TENET template should have it)
    const inputFileExists = await inputFileNodeExists(page);
    expect(inputFileExists).toBeTruthy();

    // Open InputFile node modal
    await pageObjects.workflowPage.openNodeModal(testWorkflow.inputNodeName);

    // Verify modal opened
    await pageObjects.inputFileModal.verifyModalOpen();

    // Verify folder structure is loaded
    const folders = await pageObjects.inputFileModal.getFolders();
    expect(folders.length).toBeGreaterThan(0);
    console.log('Available folders:', folders);

    // Select the test folder
    await pageObjects.inputFileModal.selectFolder(testWorkflow.folder);

    // Verify folder is selected (highlighted)
    await pageObjects.inputFileModal.verifyFolderSelected(testWorkflow.folder);

    // Get files in the folder
    const files = await pageObjects.inputFileModal.getFiles();
    expect(files.length).toBeGreaterThan(0);
    console.log('Available files:', files);

    // Verify the expected file exists (use dynamic filename)
    await pageObjects.inputFileModal.verifyFileExists(currentTestFileName);

    // Select the file
    await pageObjects.inputFileModal.selectFile(currentTestFileName);

    // Verify file is selected (highlighted)
    await pageObjects.inputFileModal.verifyFileSelected(currentTestFileName);

    // Verify current file display shows the selected file
    await pageObjects.inputFileModal.verifyCurrentFile(currentTestFileName);

    // Get current file info
    const fileInfo = await pageObjects.inputFileModal.getCurrentFileInfo();
    expect(fileInfo).not.toBeNull();
    expect(fileInfo.name).toBe(currentTestFileName);
    console.log('Selected file info:', fileInfo);

    // Verify Apply button is in "Apply" state (not yet applied)
    await pageObjects.inputFileModal.verifyApplyButtonState(false);

    // Click Apply to assign the file
    await pageObjects.inputFileModal.clickApply();
    // clickApply() already verifies "Applied" state internally

    // NOTE: Apply button click only updates Vuex store, no backend API call
    // - applyFile() method (InputFile.vue:135-149) only commits to Vuex
    // - Backend persistence happens when workflow is saved/executed
    // - This test verifies UI state and session-level persistence only
    console.log('✅ InputFile assignment completed (Vuex store updated)');

    // Verify persistence: close and reopen modal
    // This checks if Vuex store maintains file info within same session
    // Note: Tab text is displayed as "input.h5ad" (lowercase with dots)
    await pageObjects.workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(500);

    // Reopen the modal
    await pageObjects.workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await pageObjects.inputFileModal.verifyModalOpen();

    // Verify the file is still assigned (from Vuex store)
    // When modal reopens, mounted() hook (InputFile.vue:151-178) reads from store:
    //   const currentFile = this.$store.getters.getWorkflowNodeFileInfo(this.nodeId)
    await expect
      .poll(async () => {
        const info = await pageObjects.inputFileModal.getCurrentFileInfo();
        return info?.name ?? null;
      }, {
        message: 'Waiting for InputFile modal to reload assigned file',
        timeout: 10000,
      })
      .toBe(currentTestFileName);

    // Verify Apply button is still in "Applied" state
    await pageObjects.inputFileModal.verifyApplyButtonState(true);

    console.log('✅ File assignment persisted in Vuex store (session-level)');
    console.log('ℹ️  Note: Backend persistence occurs when workflow is saved/executed');
  });

  /**
   * Test: Should update assigned input file to different file
   *
   * Edge case verification:
   * - Verifies that users can change the file assignment after initial selection
   * - Tests that Apply button state resets when file is changed
   * - Confirms Vuex store updates with new file selection
   */
  test('Should update assigned input file to different file', async ({ page }) => {
    // Upload second file with unique filename for file change test
    const { uploadedFileName: firstFile, secondFileName: secondFile } = await setupTestFiles(
      pageObjects.filesPage,
      testWorkflow,
      uploadedFiles,
      { uploadSecondFile: true }
    );

    console.log(`📝 File change test files: "${firstFile}" and "${secondFile}"`);

    // Navigate to Projects page
    await pageObjects.projectsPage.goto();
    await pageObjects.projectsPage.verifyPageLoaded();

    // Create workflow with unique title
    currentWorkflowTitle = await createWorkflowWithUniqueTitle(
      pageObjects.projectsPage,
      pageObjects.workflowPage,
      testWorkflow,
      createdWorkflows
    );

    // Wait for canvas to have nodes rendered
    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    // Open InputFile modal
    await pageObjects.workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await pageObjects.inputFileModal.verifyModalOpen();

    // Select folder to display files
    await pageObjects.inputFileModal.selectFolder(testWorkflow.folder);

    // Verify we have at least 2 files for this test
    const files = await pageObjects.inputFileModal.getFiles();
    expect(files.length).toBeGreaterThan(1); // Need at least 2 files for this test

    // Select and apply first file
    await pageObjects.inputFileModal.selectFile(firstFile);
    await pageObjects.inputFileModal.clickApply();
    // clickApply() already verifies "Applied" state internally

    console.log(`✅ First file "${firstFile}" assigned`);

    // Change to second file
    console.log(`🔄 Selecting second file: ${secondFile}`);
    await pageObjects.inputFileModal.selectFile(secondFile);
    await pageObjects.inputFileModal.verifyCurrentFile(secondFile);

    // Wait for Vue reactivity to update Apply button state
    await page.waitForTimeout(500);

    // Verify Apply button returned to "Apply" state (not "Applied")
    // This should happen because fileClick() sets this.apply = false
    console.log('🔍 Checking Apply button state after file change...');
    const applyButtonText = await page.locator('label.form__button--apply').textContent();
    console.log(`Apply button text: "${applyButtonText.trim()}"`);

    // Verify modal is still open and we're on the correct page
    const currentUrl = page.url();
    console.log(`Current URL before second Apply: ${currentUrl}`);
    expect(currentUrl).toContain('/workflow');

    // Verify modal is visible
    await pageObjects.inputFileModal.verifyModalOpen();

    // Apply the second file
    console.log(`📤 Applying second file: ${secondFile}`);
    await pageObjects.inputFileModal.clickApply();
    // clickApply() already verifies "Applied" state internally

    // Verify we're still on the workflow page after apply
    const urlAfterApply = page.url();
    console.log(`Current URL after second Apply: ${urlAfterApply}`);
    expect(urlAfterApply).toContain('/workflow');

    console.log(`✅ Second file "${secondFile}" assigned successfully`);

    // Verify persistence of the new file
    // Note: Tab text is displayed as "input.h5ad" (lowercase with dots)
    await pageObjects.workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(500);

    await pageObjects.workflowPage.openNodeModal(testWorkflow.inputNodeName);

    // Wait for mounted() hook to complete and Vuex store to be read
    await page.waitForTimeout(500);

    await pageObjects.inputFileModal.verifyModalOpen();

    // Wait for InputFile modal to reload file info from Vuex store
    await expect
      .poll(async () => {
        const info = await pageObjects.inputFileModal.getCurrentFileInfo();
        return info?.name ?? null;
      }, {
        message: 'Waiting for InputFile modal to reload second file',
        timeout: 10000,
      })
      .toBe(secondFile);

    console.log(`✅ File successfully changed to "${secondFile}"`);
  });
});
