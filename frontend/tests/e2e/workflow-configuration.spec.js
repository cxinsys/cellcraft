// frontend/tests/e2e/workflow-configuration.spec.js
import { test, expect } from './fixtures/auth.js'; // 인증 fixture 사용
import { readFile } from 'fs/promises';
import { ProjectsPage } from './pages/ProjectsPage.js';
import { WorkflowPage } from './pages/WorkflowPage.js';
import { FilesPage } from './pages/FilesPage.js';
import { InputFileModal } from './pages/modals/InputFileModal.js';
import { AlgorithmModal } from './pages/modals/AlgorithmModal.js';
import { DataTableModal } from './pages/modals/DataTableModal.js';
import { ScatterPlotModal } from './pages/modals/ScatterPlotModal.js';
import { LogsModal } from './pages/modals/LogsModal.js';
import { DagModal } from './pages/modals/DagModal.js';
import { CompileCheckModal } from './pages/modals/CompileCheckModal.js';
import { ResultFilesModal } from './pages/modals/ResultFilesModal.js';
import { VisualizationModal } from './pages/modals/VisualizationModal.js';
import {
  inputFileNodeExists,
  getNodeFileAssignment,
  getWorkflowMetadata,
  getDrawflowConnections,
} from './utils/workflow.js';

function generateUniqueFileName(baseFileName) {
  const timestamp = Date.now();
  const randomId = Math.random().toString(36).substring(2, 8);
  const lastDot = baseFileName.lastIndexOf('.');
  if (lastDot === -1) {
    return `${baseFileName}_${timestamp}_${randomId}`;
  }
  const name = baseFileName.substring(0, lastDot);
  const extension = baseFileName.substring(lastDot);
  return `${name}_${timestamp}_${randomId}${extension}`;
}

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
  let dataTableModal;
  let scatterPlotModal;
  let logsModal;
  let dagModal;
  let compileCheckModal;
  let resultFilesModal;
  let visualizationModal;

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
    dataTableModal = new DataTableModal(page);
    scatterPlotModal = new ScatterPlotModal(page);
    logsModal = new LogsModal(page);
    dagModal = new DagModal(page);
    compileCheckModal = new CompileCheckModal(page);
    resultFilesModal = new ResultFilesModal(page);
    visualizationModal = new VisualizationModal(page);

    // Setup API mocking for deterministic tests
    // Commenting out for now to use real backend, uncomment when needed
    // await setupWorkflowRoutes(page, { useDefaultFixtures: true });

    // Ensure required test file exists (conditional upload to avoid duplication)
    await filesPage.goto();
    await filesPage.verifyPageLoaded();

    try {
      await filesPage.selectFolder(testWorkflow.folder);
    } catch (error) {
      console.log(`⚠️ Unable to select folder "${testWorkflow.folder}":`, error.message);
    }

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
    const ensureDemoUpload = async () => {
      const secondFileExists = await filesPage.isFilePresent('demo.h5ad');

      if (!secondFileExists) {
        console.log('📤 Uploading second test file for file change test');
        const { uploadedFileName: secondFile } = await filesPage.uploadFile(
          'demo.h5ad',
          {
            targetFileName: 'demo.h5ad',
            timeout: 60000,
          }
        );
        uploadedFiles.push(secondFile);
        await filesPage.waitForUploadComplete();
        await filesPage.verifyFileExists(secondFile);
        console.log(`✅ Second test file uploaded: ${secondFile}`);
      } else {
        console.log('ℹ️  Second test file already exists: demo.h5ad');
      }
    };

    const currentTestTitle = test.info().title;
    if (currentTestTitle.includes('change assigned file')) {
      await ensureDemoUpload();

    try {
      await filesPage.deleteFile(testWorkflow.expectedFile);
      console.log(`🧹 Removed existing ${testWorkflow.expectedFile} before test`);
    } catch (error) {
      console.log(`ℹ️ No existing ${testWorkflow.expectedFile} to delete (expected in some runs)`);
    }

    await filesPage.goto();
    await filesPage.selectFolder(testWorkflow.folder);
    await filesPage.uploadFile('test_data.h5ad', {
      targetFileName: testWorkflow.expectedFile,
      timeout: 60000,
    });
    await filesPage.waitForUploadComplete();
    await filesPage.verifyFileExists(testWorkflow.expectedFile);
    }

    // Navigate to Projects page (already authenticated via fixture)
    await projectsPage.goto();
    await projectsPage.verifyPageLoaded();
  });

  test('Priority 9: Should visualize TENET results with GRNViz', async ({ page }) => {
    test.setTimeout(300000);
    test.info().annotations.push({
      type: 'priority',
      description: 'Priority 9 - Result visualization validation with GRNViz',
    });

    await test.step('Open prepared SUCCESS workflow', async () => {
      await projectsPage.verifyWorkflowExists('SUCCESS');
      await projectsPage.openWorkflow('SUCCESS');
      await workflowPage.verifyPageLoaded();
      await workflowPage.waitForNodesReady(5, 15000);
    });

    const workflowTitle = await workflowPage.getWorkflowTitle();
    expect(workflowTitle).toContain('SUCCESS');

    await test.step('Validate existing SUCCESS job in monitoring panel', async () => {
      await workflowPage.openJobTable();
      await workflowPage.waitForJobRows(1, 60000);
      const entries = await workflowPage.getJobTableEntries();
      const pluginMatcher = 'TENET';
      const typeMatcher = 'Analysis';
      const matchingJobs = entries
        .filter((entry) => entry.name === workflowTitle)
        .filter((entry) => (entry.plugin ?? '').includes(pluginMatcher))
        .filter((entry) => (entry.type ?? '').includes(typeMatcher));

      expect(matchingJobs.length).toBeGreaterThan(0);

      matchingJobs.sort((a, b) => b.startTimestamp - a.startTimestamp);
      const latestMatchingJob = matchingJobs[0];
      expect(latestMatchingJob.status?.toUpperCase?.()).toBe('SUCCESS');

      await page.evaluate(({ job, pluginText, typeText }) => {
        window.__latestSuccessJob = {
          name: job.name,
          plugin: job.plugin,
          type: job.type,
          pluginMatcher: pluginText,
          typeMatcher: typeText,
        };
      }, {
        job: latestMatchingJob,
        pluginText: pluginMatcher,
        typeText: typeMatcher,
      });

      await workflowPage.closeJobTable();
    });

    await test.step('Inspect ResultFiles modal, select downloads, and verify sections', async () => {
      await workflowPage.openNodeModal('ResultFiles');
      await resultFilesModal.verifyModalOpen();
      await resultFilesModal.waitForPrimarySection(20000);
      await resultFilesModal.waitForIntermediateSection(20000);
      expect(await resultFilesModal.isPrimarySectionVisible()).toBeTruthy();
      expect(await resultFilesModal.isIntermediateSectionVisible()).toBeTruthy();

      const primaryFiles = await resultFilesModal.getPrimaryFileNames();
      const intermediateFiles = await resultFilesModal.getIntermediateFileNames();
      expect(primaryFiles.length + intermediateFiles.length).toBeGreaterThan(1);

      const selectionTargets = primaryFiles.slice(0, 2).map((name) => ({ name, type: 'primary' }));
      if (selectionTargets.length < 2) {
        throw new Error('Expected at least two primary files to download');
      }

      console.log('Primary files found:', primaryFiles);
      console.log('Intermediate files found:', intermediateFiles);
      console.log('Selection targets:', selectionTargets);

      const downloadRequests = [];
      const requestListener = async (request) => {
        if (
          request.url().includes('/routes/workflow/result') &&
          request.method() === 'POST'
        ) {
          try {
            const payload = await request.postDataJSON();
            downloadRequests.push({ url: request.url(), payload });
          } catch (error) {
            console.log('Failed to read download request payload:', error);
          }
        }
      };
      page.on('request', requestListener);

      const downloads = [];
      for (const target of selectionTargets) {
        const downloadPromise = page.waitForEvent('download', { timeout: 10000 });
        await resultFilesModal.clickPrimaryFileDownloadButton(target.name);
        downloads.push(await downloadPromise);
      }

      const downloadedNames = await Promise.all(
        downloads.map(async (dl) => dl.suggestedFilename())
      );
      const normalizeFilename = (name) => {
        if (!name) return '';
        return name.replace(/ \(\d+\)\./, '.').trim();
      };
      const normalizedDownloadedNames = downloadedNames.map((name) => normalizeFilename(name));
      console.log('Downloaded filenames:', downloadedNames);
      console.log('Normalized filenames:', normalizedDownloadedNames);
      console.log('Download request payloads:', downloadRequests);

      for (const target of selectionTargets) {
        expect(
          normalizedDownloadedNames.some((filename) => filename.endsWith(target.name))
        ).toBeTruthy();
      }

      for (const dl of downloads) {
        try {
          const filePath = await dl.path();
          let preview = '<binary>'; // default placeholder
          if (filePath) {
            const content = await readFile(filePath, 'utf-8').catch(() => null);
            if (content !== null) {
              preview = content.substring(0, 200);
            }
          }
          console.log('Downloaded file detail:', {
            suggestedFilename: dl.suggestedFilename(),
            path: await dl.path(),
            preview,
          });
        } catch (error) {
          console.log('Failed to inspect downloaded file:', error);
        }
      }

      page.off('request', requestListener);

      await Promise.all(downloads.map((dl) => dl.delete().catch(() => {})));
    });

    await test.step('Download execution manifest and validate sections', async () => {
      await workflowPage.openJobTable();
      await workflowPage.waitForJobRows(1, 60000);

      const latestJobContext = await page.evaluate(() => window.__latestSuccessJob || null);
      if (!latestJobContext) {
        throw new Error('Latest success job context not found for manifest download');
      }

      const manifestDownloadPromise = page.waitForEvent('download');
      await workflowPage.openJobContextMenuForTitle(latestJobContext.name, {
        pluginSubstring: latestJobContext.pluginMatcher,
        typeSubstring: latestJobContext.typeMatcher,
      });
      await workflowPage.selectJobContextOption('Download manifest');

      const manifestDownload = await manifestDownloadPromise;
      const manifestPath = await manifestDownload.path();
      const manifestContent = await readFile(manifestPath, 'utf-8');
      const manifestJson = JSON.parse(manifestContent);

      expect(manifestJson).toHaveProperty('manifest_info');
      expect(manifestJson).toHaveProperty('task_metadata');
      expect(manifestJson).toHaveProperty('plugin_metadata');
      expect(manifestJson).toHaveProperty('workflow_metadata');
      expect(manifestJson).toHaveProperty('execution_files');

      await manifestDownload.delete().catch(() => {});
      await workflowPage.closeJobTable();
    });

    await test.step('Execute GRNViz visualization and validate Plotly rendering', async () => {
      await workflowPage.openNodeModal('Visualization');

      if (!(await visualizationModal.isConfigurationMode())) {
        await visualizationModal.waitForPluginSelection();
        await visualizationModal.selectPluginByName('GRNViz');
        await expect(visualizationModal.visualizationItems.first()).toBeVisible();

        let visualizationName = 'Bar plot';
        const barPlotCount = await visualizationModal.visualizationItems
          .filter({ hasText: 'Bar plot' })
          .count();

        if (barPlotCount === 0) {
          const firstVisualizationText = (await visualizationModal.visualizationItems
            .first()
            .textContent())?.trim();

          if (!firstVisualizationText) {
            throw new Error('No visualization scripts available for GRNViz plugin');
          }

          visualizationName = firstVisualizationText;
        }

        await visualizationModal.selectVisualizationByName(visualizationName);
        await visualizationModal.proceedToConfiguration();
      }

      expect(await visualizationModal.getSelectedPluginLabel()).toContain('GRNViz');
      const selectedVisualizationLabel = await visualizationModal.getSelectedVisualizationLabel();
      expect(selectedVisualizationLabel.length).toBeGreaterThan(0);

      const inputParameters = await visualizationModal.getInputFileParameterNames();
      for (const parameterName of inputParameters) {
        const normalizedName = parameterName.trim().toLowerCase();

        if (normalizedName === 'input') {
          await visualizationModal.selectInputFileOption(parameterName, 'FdrOutdegree.txt');
          await expect
            .poll(async () => await visualizationModal.getSelectedInputFile(parameterName))
            .toBe('FdrOutdegree.txt');
        } else if (normalizedName.includes('expression')) {
          await visualizationModal.selectInputFileOption(parameterName, 'expression.csv');
          await expect
            .poll(async () => await visualizationModal.getSelectedInputFile(parameterName))
            .toBe('expression.csv');
        } else if (normalizedName.includes('trajectory')) {
          await visualizationModal.selectInputFileOption(parameterName, 'trajectory.txt');
          await expect
            .poll(async () => await visualizationModal.getSelectedInputFile(parameterName))
            .toBe('trajectory.txt');
        } else {
          await expect
            .poll(async () => (await visualizationModal.getAvailableOptionsForParameter(parameterName)).length > 0, {
              message: `Waiting for options to load for parameter ${parameterName}`,
            })
            .toBeTruthy();

          const options = await visualizationModal.getAvailableOptionsForParameter(parameterName);
          const firstNonEmpty = options.find((opt) => opt.trim() !== '' && opt.trim() !== 'Select File');
          if (!firstNonEmpty) {
            throw new Error(`No selectable option available for parameter ${parameterName}`);
          }
          await visualizationModal.selectInputFileOption(parameterName, firstNonEmpty.trim());
          await expect
            .poll(async () => await visualizationModal.getSelectedInputFile(parameterName))
            .toBe(firstNonEmpty.trim());
        }
      }

      expect(await visualizationModal.isApplyButtonEnabled()).toBeTruthy();

      const runResponsePromise = page.waitForResponse(
        (resp) =>
          resp.url().includes('/routes/workflow/visualization') &&
          resp.request().method() === 'POST',
        { timeout: 60000 }
      );

      const resultResponsePromise = page.waitForResponse(
        (resp) =>
          resp.url().includes('/routes/workflow/visualization/result') &&
          resp.request().method() === 'POST',
        { timeout: 60000 }
      );

      await visualizationModal.clickExecuteVisualization();

      const runResponse = await runResponsePromise;
      const runPayload = await runResponse.json().catch(() => ({}));
      if (Object.prototype.hasOwnProperty.call(runPayload, 'success')) {
        expect(runPayload.success).toBeTruthy();
      }

      const resultResponse = await resultResponsePromise;
      const resultPayload = await resultResponse.json().catch(() => ({}));
      if (Object.prototype.hasOwnProperty.call(resultPayload, 'success')) {
        expect(resultPayload.success).toBeTruthy();
      }
      if (Array.isArray(resultPayload.data)) {
        expect(resultPayload.data.length).toBeGreaterThan(0);
      }

      await visualizationModal.waitForPlotly();
      if ((await visualizationModal.getApplyButtonText()).includes('Show Visualization')) {
        await visualizationModal.clickShowVisualization();
        await visualizationModal.waitForPlotly();
      }
      expect(await visualizationModal.isPlotlyVisible()).toBeTruthy();

      const plotDownload = await visualizationModal.downloadPlotImage();
      const plotFilename = plotDownload.suggestedFilename();
      expect(plotFilename).toMatch(/\.png$/i);
      await plotDownload.delete().catch(() => {});
    });
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

    // Clean up route handlers (enable if API mocking is used in future tests)
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
    await expect
      .poll(async () => {
        const info = await inputFileModal.getCurrentFileInfo();
        return info?.name ?? null;
      }, {
        message: 'Waiting for InputFile modal to reload assigned file',
        timeout: 10000,
      })
      .toBe(testWorkflow.expectedFile);

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

  test('Priority 4: Should display matrix data in DataTable node', async ({ page }) => {
    test.setTimeout(60000);
    test.info().annotations.push({
      type: 'priority',
      description: 'Priority 4 - DataTable matrix display verification',
    });

    // Create workflow from TENET template
    await projectsPage.clickNewWorkflow();
    await projectsPage.selectPluginTemplate(testWorkflow.name);
    await workflowPage.verifyPageLoaded();
    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    // Assign file to InputFile node
    await workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await inputFileModal.assignFile(testWorkflow.folder, testWorkflow.expectedFile);

    // Close InputFile tab to focus on DataTable later
    await workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(300);

    // Determine DataTable node id from DOM (fallback to known default if missing)
    const dataTableLocator = await workflowPage.findNodeByType('DataTable');
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
        message: `Expected DataTable node (${dataTableNodeId}) to receive file ${testWorkflow.expectedFile}`,
        timeout: 10000,
      })
      .toBe(testWorkflow.expectedFile);

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

    await workflowPage.openNodeModal('DataTable');
    await dataTableModal.verifyModalOpen();

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
    await dataTableModal.waitForDataLoaded();

    const emptyStateVisible = await dataTableModal.isEmptyStateVisible();
    expect(emptyStateVisible).toBeFalsy();

    await expect(page.locator('.table-layout')).toBeVisible();

    const rowCount = await dataTableModal.getRowCount();
    expect(rowCount).toBeGreaterThan(0);
  });

  test('Priority 5: Should render UMAP scatter plot in ScatterPlot node', async ({ page }) => {
    test.setTimeout(60000);
    test.info().annotations.push({
      type: 'priority',
      description: 'Priority 5 - ScatterPlot UMAP visualization verification',
    });

    await projectsPage.clickNewWorkflow();
    await projectsPage.selectPluginTemplate(testWorkflow.name);
    await workflowPage.verifyPageLoaded();
    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    await workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await inputFileModal.assignFile(testWorkflow.folder, testWorkflow.expectedFile);

    await workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(300);

    const scatterLocator = await workflowPage.findNodeByType('ScatterPlot');
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
          return assignment.find((item) => item === testWorkflow.expectedFile) ?? null;
        }
        if (typeof assignment === 'object') {
          const values = Object.values(assignment);
          return values.find((item) => item === testWorkflow.expectedFile) ?? null;
        }
        return null;
      }, {
        message: `Expected ScatterPlot node (${scatterNodeId}) to receive file ${testWorkflow.expectedFile}`,
        timeout: 10000,
      })
      .toBe(testWorkflow.expectedFile);

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

    await workflowPage.openNodeModal('ScatterPlot');
    await scatterPlotModal.verifyModalOpen();

    const scatterDataResponse = await scatterDataResponsePromise;
    const scatterPayload = await scatterDataResponse.json();
    console.log('📡 ScatterPlot API payload:', scatterPayload);

    if (Object.prototype.hasOwnProperty.call(scatterPayload, 'success')) {
      expect(scatterPayload.success).toBeTruthy();
    }

    await scatterPlotModal.waitForPlotly();

    const blankStateVisible = await scatterPlotModal.isBlankStateVisible();
    expect(blankStateVisible).toBeFalsy();

    await expect(page.locator('#plotly__scatter')).toBeVisible();

    const traceCount = await scatterPlotModal.getTraceCount();
    expect(traceCount).toBeGreaterThan(0);

    await test.step('Update ScatterPlot X axis dropdown', async () => {
      const { previous, next } = await scatterPlotModal.selectDifferentXAxis();
      console.log(`🔁 ScatterPlot X axis changed: ${previous} → ${next}`);
      expect(next).not.toBe(previous);
      await scatterPlotModal.waitForPlotly();
      const current = await scatterPlotModal.getSelectedXAxisValue();
      expect(current).toBe(next);
    });

    await test.step('Update ScatterPlot Y axis dropdown', async () => {
      const { previous, next } = await scatterPlotModal.selectDifferentYAxis();
      console.log(`🔁 ScatterPlot Y axis changed: ${previous} → ${next}`);
      expect(next).not.toBe(previous);
      await scatterPlotModal.waitForPlotly();
      const current = await scatterPlotModal.getSelectedYAxisValue();
      expect(current).toBe(next);
    });

    await test.step('Update ScatterPlot group dropdown', async () => {
      const { previous, next } = await scatterPlotModal.selectDifferentGroup();
      console.log(`🔁 ScatterPlot group changed: ${previous} → ${next}`);
      expect(next).not.toBe(previous);
      await scatterPlotModal.waitForPlotly();
      const current = await scatterPlotModal.getSelectedGroupValue();
      expect(current).toBe(next);
      const tracesAfterGroup = await scatterPlotModal.getTraceCount();
      expect(tracesAfterGroup).toBeGreaterThan(0);
    });
  });

  test('Priority 6: Should configure Algorithm node parameters', async ({ page }) => {
    test.setTimeout(60000);
    test.info().annotations.push({
      type: 'priority',
      description: 'Priority 6 - Algorithm parameter configuration verification',
    });

    await projectsPage.clickNewWorkflow();
    await projectsPage.selectPluginTemplate(testWorkflow.name);
    await workflowPage.verifyPageLoaded();
    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    await workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await inputFileModal.assignFile(testWorkflow.folder, testWorkflow.expectedFile);

    await workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(300);

    const algorithmLocator = await workflowPage.findNodeByType('Algorithm');
    const algorithmNodeIdAttr = await algorithmLocator.getAttribute('id');
    const algorithmNodeId = algorithmNodeIdAttr?.replace('node-', '') ?? '12';

    await expect
      .poll(async () => {
        const assignment = await getNodeFileAssignment(page, algorithmNodeId);
        if (!assignment) return null;

        if (typeof assignment === 'string') {
          return assignment.includes(testWorkflow.expectedFile) ? assignment : null;
        }

        if (Array.isArray(assignment)) {
          return assignment.find((value) =>
            typeof value === 'string' && value.includes(testWorkflow.expectedFile)
          ) ?? null;
        }

        if (typeof assignment === 'object') {
          return (
            Object.values(assignment).find(
              (value) => typeof value === 'string' && value.includes(testWorkflow.expectedFile)
            ) ?? null
          );
        }

        return null;
      }, {
        message: `Expected Algorithm node (${algorithmNodeId}) to receive file ${testWorkflow.expectedFile}`,
        timeout: 10000,
      })
      .not.toBeNull();

    await workflowPage.openNodeModal('Algorithm');
    await algorithmModal.verifyModalOpen();

    await expect
      .poll(async () => {
        const logoText = await algorithmModal.getPluginLogoText();
        return logoText;
      }, {
        message: `Waiting for algorithm logo to display plugin name ${testWorkflow.name}`,
        timeout: 15000,
      })
      .toContain(testWorkflow.name);

    const pluginLogoText = await algorithmModal.getPluginLogoText();
    console.log('🔖 Algorithm logo text:', pluginLogoText);

    await algorithmModal.verifyFileInDropdown(0, testWorkflow.expectedFile);
    const selectedInputFile = await algorithmModal.getSelectedInputFile(0);
    expect(selectedInputFile).toContain(testWorkflow.expectedFile);

    const metadataBefore = await getWorkflowMetadata(page);
    console.log(
      '📦 Vuex drawflow data before Algorithm parameter changes:',
      JSON.stringify(metadataBefore?.workflowInfo?.drawflow?.Home?.data ?? null, null, 2)
    );

    const algorithmNodeDataBefore = metadataBefore?.workflowInfo?.drawflow?.Home?.data?.[algorithmNodeId];
    const pluginRulesBefore = algorithmNodeDataBefore?.data?.selectedPluginRules ?? [];
    console.log('🧮 Algorithm selectedPluginRules (before):', JSON.stringify(pluginRulesBefore, null, 2));

    const flattenedParamsBefore = pluginRulesBefore.flatMap((rule) => rule.parameters ?? []);

    const numericParam = flattenedParamsBefore.find((param) =>
      ['int', 'float', 'number'].includes(param?.type)
    );
    let numericParamName = null;
    let numericNewValue = null;
    if (numericParam) {
      numericParamName = numericParam.name;
      const numericInitialValue = Number(numericParam.defaultValue ?? 0);
      numericNewValue = String(numericInitialValue + 1);
      await algorithmModal.setParameterValueByName(numericParamName, numericNewValue);
      const numericUiValue = await algorithmModal.getParameterValueByName(numericParamName);
      expect(numericUiValue).toBe(numericNewValue);
    }

    const stringParam = flattenedParamsBefore.find((param) => param?.type === 'string');
    let stringParamName = null;
    let stringNewValue = null;
    if (!numericParam && stringParam) {
      stringParamName = stringParam.name;
      const stringInitialValue = stringParam.defaultValue ?? '';
      stringNewValue = stringInitialValue === '' ? 'test-value' : `${stringInitialValue}-updated`;
      await algorithmModal.setParameterValueByName(stringParamName, stringNewValue);
      const stringUiValue = await algorithmModal.getParameterValueByName(stringParamName);
      expect(stringUiValue).toBe(stringNewValue);
    }

    const booleanParam = flattenedParamsBefore.find((param) => param?.type === 'boolean');
    let booleanParamName = null;
    let booleanNewValue = null;
    if (booleanParam) {
      booleanParamName = booleanParam.name;
      const booleanInitialValue = booleanParam.defaultValue === true || booleanParam.defaultValue === 'true';
      booleanNewValue = !booleanInitialValue;
      await algorithmModal.setParameterValueByName(booleanParamName, booleanNewValue);
      const booleanUiValue = await algorithmModal.getParameterValueByName(booleanParamName);
      expect(booleanUiValue).toBe(booleanNewValue);
    }

    const dropdownParam = flattenedParamsBefore.find(
      (param) => param?.type === 'h5adParameter' && param?.name !== 'clusters'
    );
    let dropdownParamName = null;
    let dropdownNewValue = null;
    if (dropdownParam) {
      dropdownParamName = dropdownParam.name;
      const dropdownOptions = await algorithmModal.getParameterDropdownOptions(dropdownParamName);
      console.log(`📝 Dropdown options for ${dropdownParamName}:`, dropdownOptions);
      const preferredOption = dropdownOptions.find((opt) => opt && opt !== dropdownParam.defaultValue) ?? dropdownOptions[0];

      if (preferredOption && preferredOption !== dropdownParam.defaultValue) {
        const { next } = await algorithmModal.selectParameterDropdownOption(
          dropdownParamName,
          preferredOption
        );
        dropdownNewValue = next;
        const dropdownUiValue = await algorithmModal.getParameterValueByName(dropdownParamName);
        expect(dropdownUiValue).toBe(dropdownNewValue);
      }
    }

    await page.waitForTimeout(500);

    const metadataAfter = await getWorkflowMetadata(page);
    console.log(
      '📦 Vuex drawflow data after Algorithm polling:',
      JSON.stringify(metadataAfter?.workflowInfo?.drawflow?.Home?.data ?? null, null, 2)
    );

    const algorithmNodeDataAfter = metadataAfter?.workflowInfo?.drawflow?.Home?.data?.[algorithmNodeId];
    const pluginRulesAfter = algorithmNodeDataAfter?.data?.selectedPluginRules ?? [];
    const flattenedParamsAfter = pluginRulesAfter.flatMap((rule) => rule.parameters ?? []);

    if (numericParamName) {
      const updatedNumericParam = flattenedParamsAfter.find((param) => param?.name === numericParamName);
      expect(updatedNumericParam?.defaultValue).toBe(numericNewValue);
    }

    if (stringParamName) {
      const updatedStringParam = flattenedParamsAfter.find((param) => param?.name === stringParamName);
      expect(updatedStringParam?.defaultValue).toBe(stringNewValue);
    }

    if (booleanParamName) {
      const updatedBooleanParam = flattenedParamsAfter.find((param) => param?.name === booleanParamName);
      const booleanAfterValue = updatedBooleanParam?.defaultValue;
      const booleanAfterNormalized = booleanAfterValue === true || booleanAfterValue === 'true';
      expect(booleanAfterNormalized).toBe(booleanNewValue);
    }

    if (dropdownParamName && dropdownNewValue) {
      const updatedDropdownParam = flattenedParamsAfter.find((param) => param?.name === dropdownParamName);
      expect(updatedDropdownParam?.defaultValue).toBe(dropdownNewValue);
    }
  });

  test('Priority 7-8: Should execute TENET workflow, monitor metadata, and clean up resources', async ({ page }) => {
    test.setTimeout(600000);
    test.info().annotations.push({
      type: 'priority',
      description: 'Priority 7 & 8 - Combined execution, monitoring, and cleanup flow',
    });

    const desiredClusters = ['CD4+ T', 'CD14+ Mono', 'NK'];
    const cleanupUploads = [];
    const newTitle = `Playwright Workflow ${Date.now()}`;
    const fixturesToUpload = ['test_data.h5ad'];
    let workflowInputFileName = testWorkflow.expectedFile;

    await test.step('Upload additional files for workflow run', async () => {
      await filesPage.goto();
      await filesPage.verifyPageLoaded();

      try {
        await filesPage.selectFolder(testWorkflow.folder);
      } catch (error) {
        console.warn(`⚠️ Unable to select folder "${testWorkflow.folder}":`, error.message);
      }

      for (const fixtureName of fixturesToUpload) {
        const uniqueFileName = generateUniqueFileName(fixtureName);
        const { uploadedFileName } = await filesPage.uploadFile(fixtureName, {
          targetFileName: uniqueFileName,
        });
        cleanupUploads.push(uploadedFileName);
        uploadedFiles.push(uploadedFileName);
        await filesPage.waitForUploadComplete();
        await filesPage.verifyFileExists(uploadedFileName);
        console.log(`✅ Uploaded fixture ${fixtureName} as ${uploadedFileName}`);

        // Use the newly uploaded PBMC file as workflow input
        workflowInputFileName = uploadedFileName;
      }

      await projectsPage.goto();
      await projectsPage.verifyPageLoaded();
    });

    let originalTitle;
    let currentJobTitle;

    await test.step('Create workflow from TENET template', async () => {
      await projectsPage.clickNewWorkflow();
      await projectsPage.selectPluginTemplate(testWorkflow.name);
      await workflowPage.verifyPageLoaded();
      await page.waitForSelector('.drawflow-node', { timeout: 10000 });
      originalTitle = await workflowPage.getWorkflowTitle();
      currentJobTitle = originalTitle;
    });

    await test.step('Assign InputFile node to pbmc dataset', async () => {
      await workflowPage.openNodeModal(testWorkflow.inputNodeName);
      await inputFileModal.assignFile(testWorkflow.folder, workflowInputFileName);
      await workflowPage.closeTab(testWorkflow.inputNodeTabName);
      await page.waitForTimeout(300);
    });

    await test.step('Configure Algorithm node parameters for TENET', async () => {
      await workflowPage.openNodeModal('Algorithm');
      await algorithmModal.verifyModalOpen();

      await expect
        .poll(async () => {
          const options = await algorithmModal.getParameterDropdownOptions('cell group');
          return options.includes('seurat_annotation') ? 'ready' : null;
        }, {
          message: 'Waiting for Cell group options to load',
          timeout: 20000,
        })
        .toBe('ready');

      await algorithmModal.setParameterValueByName('cell group', 'seurat_annotation');
      expect(await algorithmModal.getParameterValueByName('cell group')).toBe('seurat_annotation');

      await expect
        .poll(async () => {
          const options = await algorithmModal.getParameterDropdownOptions('pseudotime');
          return options.includes('Pseudotime') ? 'ready' : null;
        }, {
          message: 'Waiting for pseudotime options to load',
          timeout: 20000,
        })
        .toBe('ready');

      await algorithmModal.setParameterValueByName('pseudotime', 'Pseudotime');
      expect(await algorithmModal.getParameterValueByName('pseudotime')).toBe('Pseudotime');

      await algorithmModal.setParameterValueByName('clusters', desiredClusters);
      const selectedClusters = await algorithmModal.getParameterValueByName('clusters');
      expect(selectedClusters).toEqual(expect.arrayContaining(desiredClusters));

      const pluginLogo = await algorithmModal.getPluginLogoText();
      console.log('Algorithm modal plugin:', pluginLogo);
      expect(pluginLogo).toContain(testWorkflow.name);
    });

    await test.step('Open compile check panel and verify task summary', async () => {
      await workflowPage.openCompileCheck();
      await compileCheckModal.waitForOpen();
      await compileCheckModal.waitForResourcesLoaded();
      await compileCheckModal.verifyCoreSectionsVisible();

      const taskEntries = await compileCheckModal.getTaskEntries();
      console.log('Compile check task entries:', taskEntries);
      expect(taskEntries.length).toBeGreaterThan(0);
      expect(
        taskEntries.some((entry) => entry.plugin && entry.plugin.includes(testWorkflow.name))
      ).toBeTruthy();

      const resourceLabels = await compileCheckModal.getResourceLabels();
      console.log('Resource summary labels:', resourceLabels);
      expect(resourceLabels.length).toBeGreaterThan(0);
    });

    await test.step('Execute workflow and wait for RUNNING status', async () => {
      page.once('dialog', async (dialog) => {
        console.log('Execution confirmation dialog:', dialog.message());
        await dialog.accept();
      });

      await compileCheckModal.clickExecute();
      await compileCheckModal.waitForClose();

      await page.waitForTimeout(1000);
      await workflowPage.openJobTable();
      await workflowPage.waitForJobRows(1, 180000);

      await workflowPage.waitForLatestJobStatus(currentJobTitle, 'RUNNING', 240000);
      const latestJob = await workflowPage.getLatestJobEntryByTitle(currentJobTitle);

      console.log('Latest job entry:', latestJob);
      expect(latestJob).not.toBeNull();
      expect(latestJob?.status?.toUpperCase()).toContain('RUNNING');
      expect(latestJob?.plugin ?? '').toContain(testWorkflow.name);

      const pluginFormatRegex = /^[^/]+\/[^ :]+ : v\d+(?:\.\d+)*$/;
      expect(latestJob?.plugin ?? '').toMatch(pluginFormatRegex);
    });

    await test.step('Wait and rename workflow title', async () => {
      await page.waitForTimeout(30000);
      await workflowPage.closeJobTable();
      await workflowPage.updateWorkflowTitle(newTitle);
      await workflowPage.saveWorkflow();
      await workflowPage.closeMessage().catch(() => {});
      currentJobTitle = newTitle;
    });

    await test.step('Validate RUNNING job and inspect logs', async () => {
      await workflowPage.openJobTable();
      await workflowPage.waitForJobRows(1, 180000);

      await expect
        .poll(async () => {
          const entry = await workflowPage.getLatestJobEntryByTitle(currentJobTitle);
          return entry?.status?.toUpperCase() ?? null;
        }, {
          timeout: 240000,
          message: `Waiting for job "${currentJobTitle}" to remain RUNNING`,
        })
        .toBe('RUNNING');

      const latestJob = await workflowPage.getLatestJobEntryByTitle(currentJobTitle);
      console.log('Latest job entry:', latestJob);
      expect(latestJob).not.toBeNull();
      expect(latestJob?.name).toBe(currentJobTitle);
      expect(latestJob?.plugin ?? '').toContain(testWorkflow.name);
      const pluginFormatRegex = /^[^/]+\/[^ :]+ : v\d+(?:\.\d+)*$/;
      expect(latestJob?.plugin ?? '').toMatch(pluginFormatRegex);

      await workflowPage.openJobContextMenuForTitle(currentJobTitle);
      await workflowPage.selectJobContextOption('View logs');
      await logsModal.waitForOpen();
      await logsModal.waitForLoaded();
      await logsModal.expectLogsAvailable();

      try {
        const jsonDownload = await logsModal.downloadAllLogsJson();
        await jsonDownload.delete().catch(() => {});
      } catch (error) {
        console.warn('⚠️ Failed to download JSON logs:', error.message);
      }

      try {
        const txtDownload = await logsModal.downloadFirstLogTxt();
        await txtDownload.delete().catch(() => {});
      } catch (error) {
        console.warn('⚠️ Failed to download TXT log:', error.message);
      }

      await logsModal.close();
      await workflowPage.closeMessage().catch(() => {});
    });

    await test.step('Inspect workflow progress visualization', async () => {
      const dagStructurePromise = page.waitForResponse(
        (resp) => resp.url().includes('/dag-structure') && resp.request().method() === 'GET',
        { timeout: 20000 }
      );

      const ruleStatusPromise = page.waitForResponse(
        (resp) => resp.url().includes('/rule-status') && resp.request().method() === 'GET',
        { timeout: 20000 }
      );

      await workflowPage.openJobContextMenuForTitle(currentJobTitle);
      await workflowPage.selectJobContextOption('View progress');

      const dagStructureResponse = await dagStructurePromise;
      expect(dagStructureResponse.ok()).toBeTruthy();

      const ruleStatusResponse = await ruleStatusPromise;
      expect(ruleStatusResponse.ok()).toBeTruthy();

      await dagModal.waitForOpen();
      await dagModal.waitForLoaded();
      await dagModal.close();
    });

    await test.step('Cancel running job and verify message', async () => {
      await workflowPage.cancelJobByTitle(currentJobTitle);
      await workflowPage.waitForMessage('Cancel task successfully!', 15000);
      await workflowPage.closeMessage().catch(() => {});

      await workflowPage.openJobTable();
      await workflowPage.waitForJobRows(1, 60000);
      await workflowPage.waitForLatestJobStatus(currentJobTitle, 'REVOKED', 240000);

      const cancelledJob = await workflowPage.getLatestJobEntryByTitle(currentJobTitle);
      console.log('Cancelled job entry:', cancelledJob);
      expect(cancelledJob?.status?.toUpperCase()).toBe('REVOKED');

      await workflowPage.closeJobTable();
    });

    await test.step('Clean up uploaded test files', async () => {
      if (cleanupUploads.length === 0) {
        return;
      }

      await filesPage.goto();
      await filesPage.verifyPageLoaded();

      try {
        await filesPage.selectFolder(testWorkflow.folder);
      } catch (error) {
        console.warn(`⚠️ Unable to re-open folder "${testWorkflow.folder}":`, error.message);
      }

      for (const fileName of cleanupUploads) {
        try {
          await filesPage.deleteFile(fileName);
          await filesPage.verifyFileNotExists(fileName);
          console.log(`🧹 Deleted uploaded file: ${fileName}`);
        } catch (error) {
          console.warn(`⚠️ Failed to delete uploaded file ${fileName}:`, error.message);
        } finally {
          const idx = uploadedFiles.indexOf(fileName);
          if (idx !== -1) {
            uploadedFiles.splice(idx, 1);
          }
        }
      }

      cleanupUploads.length = 0;
    });
  });
});
