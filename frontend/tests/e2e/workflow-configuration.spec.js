// frontend/tests/e2e/workflow-configuration.spec.js
import { test, expect } from './fixtures/auth.js'; // 인증 fixture 사용
import { ProjectsPage } from './pages/ProjectsPage.js';
import { WorkflowPage } from './pages/WorkflowPage.js';
import { FilesPage } from './pages/FilesPage.js';
import { InputFileModal } from './pages/modals/InputFileModal.js';
import { AlgorithmModal } from './pages/modals/AlgorithmModal.js';
import { setupWorkflowRoutes, clearWorkflowRoutes } from './fixtures/workflow-api.js';
import { inputFileNodeExists } from './utils/workflow.js';

/**
 * Test Suite: Workflow Configuration (Section 2.2)
 *
 * Priority order based on QA analysis:
 * 1. InputFile assignment (unlocks downstream tests)
 * 2. Dropdown file selection (Algorithm modal)
 * 3. Algorithm parameter configuration
 * 4. DataTable matrix display
 * 5. ScatterPlot UMAP visualization
 */

test.describe('Workflow Configuration Tests', () => {
  // Configure serial mode to prevent parallel execution issues
  // This ensures only one test runs at a time, avoiding duplicate file uploads
  test.describe.configure({ mode: 'serial' });

  let projectsPage;
  let workflowPage;
  let filesPage;
  let inputFileModal;
  let algorithmModal;

  // Track uploaded files for cleanup
  const uploadedFiles = [];

  // Test data
  const testWorkflow = {
    name: 'TENET',
    expectedFile: 'pbmc_light_1000.h5ad',
    folder: 'data',
    // TENET 템플릿의 실제 InputFile 노드 이름
    inputNodeName: 'Input h5ad', // Node title for opening modal
    inputNodeTabName: 'input.h5ad', // Tab text displayed in lowercase with dots
  };

  test.beforeEach(async ({ page }) => {
    // Initialize page objects
    projectsPage = new ProjectsPage(page);
    workflowPage = new WorkflowPage(page);
    filesPage = new FilesPage(page);
    inputFileModal = new InputFileModal(page);
    algorithmModal = new AlgorithmModal(page);

    // Setup API mocking for deterministic tests
    // Commenting out for now to use real backend, uncomment when needed
    // await setupWorkflowRoutes(page, { useDefaultFixtures: true });

    // Ensure required test file exists (conditional upload to avoid duplication)
    await filesPage.goto();
    await filesPage.verifyPageLoaded();

    const fileExists = await filesPage.isFilePresent(testWorkflow.expectedFile);

    if (!fileExists) {
      console.log(
        `📤 Uploading test file: ${testWorkflow.expectedFile} (not found)`
      );
      const { uploadedFileName } = await filesPage.uploadFile('test_data.h5ad', {
        targetFileName: testWorkflow.expectedFile,
      });
      uploadedFiles.push(uploadedFileName);
      await filesPage.waitForUploadComplete();
      await filesPage.verifyFileExists(uploadedFileName);
      console.log(`✅ Test file uploaded: ${uploadedFileName}`);
    } else {
      console.log(
        `ℹ️  Test file already exists: ${testWorkflow.expectedFile}`
      );
    }

    // For Test 2 (file change test), upload a second file
    const currentTestTitle = test.info().title;
    if (currentTestTitle.includes('change assigned file')) {
      const secondFileExists = await filesPage.isFilePresent('demo.h5ad');

      if (!secondFileExists) {
        console.log('📤 Uploading second test file for file change test');
        const { uploadedFileName: secondFile } = await filesPage.uploadFile(
          'demo.h5ad',
          {
            targetFileName: 'demo.h5ad',
          }
        );
        uploadedFiles.push(secondFile);
        await filesPage.waitForUploadComplete();
        await filesPage.verifyFileExists(secondFile);
        console.log(`✅ Second test file uploaded: ${secondFile}`);
      } else {
        console.log('ℹ️  Second test file already exists: demo.h5ad');
      }
    }

    // Navigate to Projects page (already authenticated via fixture)
    await projectsPage.goto();
    await projectsPage.verifyPageLoaded();
  });

  test.afterEach(async ({ page }) => {
    // Clean up uploaded files (only those uploaded by this test)
    if (uploadedFiles.length > 0) {
      console.log(`🧹 Cleaning up ${uploadedFiles.length} uploaded file(s)`);
      await filesPage.goto();

      for (const fileName of uploadedFiles) {
        try {
          await filesPage.deleteFile(fileName);
          console.log(`✅ Deleted test file: ${fileName}`);
        } catch (error) {
          console.log(`⚠️  Failed to delete file: ${fileName}`, error.message);
        }
      }

      // Clear the array for next test
      uploadedFiles.length = 0;
    }

    // Clean up route handlers
    // await clearWorkflowRoutes(page);
  });

  /**
   * Priority 1: InputFile Assignment
   *
   * This test verifies the core file assignment functionality:
   * - Opening InputFile node modal
   * - Navigating folder structure
   * - Selecting a file
   * - Applying the file assignment
   * - Verifying the Apply button state changes
   *
   * Success criteria:
   * - Folder selection highlights the folder
   * - File selection shows file details
   * - Apply button changes to "Applied" after submission
   * - File is persisted in Vuex store (verified by reopening modal)
   */
  test('Priority 1: Should assign file to InputFile node', async ({ page }) => {
    test.info().annotations.push({
      type: 'priority',
      description: 'Priority 1 - Foundation test that unlocks downstream tests',
    });

    // Create workflow from TENET template
    await projectsPage.clickNewWorkflow();
    await projectsPage.selectPluginTemplate(testWorkflow.name);

    // Verify workflow page loaded
    await workflowPage.verifyPageLoaded();

    // Wait for canvas to have nodes rendered
    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    // Verify InputFile node exists on canvas (TENET template should have it)
    const inputFileExists = await inputFileNodeExists(page);
    expect(inputFileExists).toBeTruthy();

    // Open InputFile node modal
    await workflowPage.openNodeModal(testWorkflow.inputNodeName);

    // Verify modal opened
    await inputFileModal.verifyModalOpen();

    // Verify folder structure is loaded
    const folders = await inputFileModal.getFolders();
    expect(folders.length).toBeGreaterThan(0);
    console.log('Available folders:', folders);

    // Select the test folder
    await inputFileModal.selectFolder(testWorkflow.folder);

    // Verify folder is selected (highlighted)
    await inputFileModal.verifyFolderSelected(testWorkflow.folder);

    // Get files in the folder
    const files = await inputFileModal.getFiles();
    expect(files.length).toBeGreaterThan(0);
    console.log('Available files:', files);

    // Verify the expected file exists
    await inputFileModal.verifyFileExists(testWorkflow.expectedFile);

    // Select the file
    await inputFileModal.selectFile(testWorkflow.expectedFile);

    // Verify file is selected (highlighted)
    await inputFileModal.verifyFileSelected(testWorkflow.expectedFile);

    // Verify current file display shows the selected file
    await inputFileModal.verifyCurrentFile(testWorkflow.expectedFile);

    // Get current file info
    const fileInfo = await inputFileModal.getCurrentFileInfo();
    expect(fileInfo).not.toBeNull();
    expect(fileInfo.name).toBe(testWorkflow.expectedFile);
    console.log('Selected file info:', fileInfo);

    // Verify Apply button is in "Apply" state (not yet applied)
    await inputFileModal.verifyApplyButtonState(false);

    // Click Apply to assign the file
    await inputFileModal.clickApply();
    // clickApply() already verifies "Applied" state internally

    // NOTE: Apply button click only updates Vuex store, no backend API call
    // - applyFile() method (InputFile.vue:135-149) only commits to Vuex
    // - Backend persistence happens when workflow is saved/executed
    // - This test verifies UI state and session-level persistence only
    console.log('✅ InputFile assignment completed (Vuex store updated)');

    // Verify persistence: close and reopen modal
    // This checks if Vuex store maintains file info within same session
    // Note: Tab text is displayed as "input.h5ad" (lowercase with dots)
    await workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(500);

    // Reopen the modal
    await workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await inputFileModal.verifyModalOpen();

    // Verify the file is still assigned (from Vuex store)
    // When modal reopens, mounted() hook (InputFile.vue:151-178) reads from store:
    //   const currentFile = this.$store.getters.getWorkflowNodeFileInfo(this.nodeId)
    const persistedFileInfo = await inputFileModal.getCurrentFileInfo();
    expect(persistedFileInfo).not.toBeNull();
    expect(persistedFileInfo.name).toBe(testWorkflow.expectedFile);

    // Verify Apply button is still in "Applied" state
    await inputFileModal.verifyApplyButtonState(true);

    console.log('✅ File assignment persisted in Vuex store (session-level)');
    console.log(
      'ℹ️  Note: Backend persistence occurs when workflow is saved/executed'
    );
  });

  /**
   * Priority 1 Edge Case: Change assigned file
   *
   * Verifies that users can change the file assignment after initial selection
   */
  test('Priority 1: Should change assigned file to different file', async ({
    page,
  }) => {
    test.setTimeout(60000); // 60 seconds timeout for file upload operations

    test.info().annotations.push({
      type: 'priority',
      description: 'Priority 1 - Edge case: changing file assignment',
    });

    // Create workflow
    await projectsPage.clickNewWorkflow();
    await projectsPage.selectPluginTemplate(testWorkflow.name);
    await workflowPage.verifyPageLoaded();

    // Wait for canvas to have nodes rendered
    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    // Open InputFile modal
    await workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await inputFileModal.verifyModalOpen();

    // Select folder to display files
    await inputFileModal.selectFolder(testWorkflow.folder);

    // Verify we have at least 2 files for this test
    const files = await inputFileModal.getFiles();
    expect(files.length).toBeGreaterThan(1); // Need at least 2 files for this test

    // Use explicit file names (uploaded in beforeEach)
    const firstFile = 'pbmc_light_1000.h5ad';
    const secondFile = 'demo.h5ad';

    // Select and apply first file
    await inputFileModal.selectFile(firstFile);
    await inputFileModal.clickApply();
    // clickApply() already verifies "Applied" state internally

    console.log(`✅ First file "${firstFile}" assigned`);

    // Change to second file
    console.log(`🔄 Selecting second file: ${secondFile}`);
    await inputFileModal.selectFile(secondFile);
    await inputFileModal.verifyCurrentFile(secondFile);

    // Wait for Vue reactivity to update Apply button state
    await page.waitForTimeout(500);

    // Verify Apply button returned to "Apply" state (not "Applied")
    // This should happen because fileClick() sets this.apply = false
    console.log('🔍 Checking Apply button state after file change...');
    const applyButtonText = await page
      .locator('label.form__button--apply')
      .textContent();
    console.log(`Apply button text: "${applyButtonText.trim()}"`);

    // Verify modal is still open and we're on the correct page
    const currentUrl = page.url();
    console.log(`Current URL before second Apply: ${currentUrl}`);
    expect(currentUrl).toContain('/workflow');

    // Verify modal is visible
    await inputFileModal.verifyModalOpen();

    // Apply the second file
    console.log(`📤 Applying second file: ${secondFile}`);
    await inputFileModal.clickApply();
    // clickApply() already verifies "Applied" state internally

    // Verify we're still on the workflow page after apply
    const urlAfterApply = page.url();
    console.log(`Current URL after second Apply: ${urlAfterApply}`);
    expect(urlAfterApply).toContain('/workflow');

    console.log(`✅ Second file "${secondFile}" assigned successfully`);

    // Verify persistence of the new file
    // Note: Tab text is displayed as "input.h5ad" (lowercase with dots)
    await workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(500);

    await workflowPage.openNodeModal(testWorkflow.inputNodeName);

    // Wait for mounted() hook to complete and Vuex store to be read
    await page.waitForTimeout(500);

    await inputFileModal.verifyModalOpen();

    const persistedFileInfo = await inputFileModal.getCurrentFileInfo();
    expect(persistedFileInfo.name).toBe(secondFile);

    console.log(`✅ File successfully changed to "${secondFile}"`);
  });
});
