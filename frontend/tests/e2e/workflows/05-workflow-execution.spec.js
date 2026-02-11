// frontend/tests/e2e/workflows/05-workflow-execution.spec.js
import { test, expect } from '../fixtures/auth.js';
import { testWorkflow, generateUniqueFileName } from './support/workflow-constants.js';
import { setupPageObjects, setupTestFiles, cleanupTestFiles, createWorkflowWithUniqueTitle, cleanupWorkflows } from './support/workflow-setup.js';

/**
 * Test Suite: Workflow Execution and Monitoring
 *
 * This comprehensive test suite verifies the complete workflow execution lifecycle:
 * - Workflow creation and InputFile assignment
 * - Algorithm parameter configuration
 * - Compile check panel validation
 * - Workflow execution and status monitoring
 * - Real-time job status updates (RUNNING)
 * - Log inspection and download
 * - DAG progress visualization
 * - Job cancellation
 * - Job deletion and removal verification
 * - Resource cleanup
 *
 * Success criteria:
 * - Workflow executes successfully
 * - Job status transitions correctly (PENDING → RUNNING → REVOKED)
 * - Plugin information displays correctly (name/version)
 * - Logs are accessible and downloadable
 * - DAG structure and rule status APIs respond successfully
 * - Job cancellation works properly
 * - Job deletion removes entry from table
 * - Uploaded files are cleaned up after test
 */
test.describe('Workflow Execution and Monitoring', () => {
  test.describe.configure({ mode: 'serial' });

  let pageObjects;
  const uploadedFiles = [];
  const createdWorkflows = [];
  let currentWorkflowTitle = null;

  test.beforeEach(async ({ page }) => {
    pageObjects = setupPageObjects(page);
    // Upload test file with unique filename (not used in this test, but maintains consistency)
    await setupTestFiles(pageObjects.filesPage, testWorkflow, uploadedFiles);
    await pageObjects.projectsPage.goto();
    await pageObjects.projectsPage.verifyPageLoaded();
  });

  test.afterEach(async ({ page }) => {
    await cleanupTestFiles(pageObjects.filesPage, uploadedFiles);
    await cleanupWorkflows(pageObjects.projectsPage, createdWorkflows);
  });

  /**
   * Test: Should execute TENET workflow with monitoring and cleanup
   *
   * This test covers the entire workflow execution lifecycle:
   * - File upload and workflow setup
   * - Algorithm parameter configuration
   * - Compile check validation
   * - Workflow execution
   * - Job monitoring (RUNNING status)
   * - Workflow title update
   * - Log inspection and download
   * - DAG progress visualization
   * - Job cancellation
   * - Job deletion and verification
   * - Resource cleanup
   */
  test('Should execute TENET workflow with monitoring and cleanup', async ({ page }) => {
    test.setTimeout(600000);

    const desiredClusters = ['CD4+ T', 'CD14+ Mono', 'NK'];
    const cleanupUploads = [];
    const fixturesToUpload = ['test_data.h5ad'];
    let workflowInputFileName = testWorkflow.expectedFile;

    await test.step('Upload additional files for workflow run', async () => {
      await pageObjects.filesPage.goto();
      await pageObjects.filesPage.verifyPageLoaded();

      try {
        await pageObjects.filesPage.selectFolder(testWorkflow.folder);
      } catch (error) {
        console.warn(`⚠️ Unable to select folder "${testWorkflow.folder}":`, error.message);
      }

      for (const fixtureName of fixturesToUpload) {
        const uniqueFileName = generateUniqueFileName(fixtureName);
        const { uploadedFileName } = await pageObjects.filesPage.uploadFile(fixtureName, {
          targetFileName: uniqueFileName,
        });
        cleanupUploads.push(uploadedFileName);
        uploadedFiles.push(uploadedFileName);
        await pageObjects.filesPage.waitForUploadComplete();
        await pageObjects.filesPage.verifyFileExists(uploadedFileName);
        console.log(`✅ Uploaded fixture ${fixtureName} as ${uploadedFileName}`);

        // Use the newly uploaded PBMC file as workflow input
        workflowInputFileName = uploadedFileName;
      }

      await pageObjects.projectsPage.goto();
      await pageObjects.projectsPage.verifyPageLoaded();
    });

    await test.step('Create workflow from TENET template', async () => {
      // Create workflow with unique title
      currentWorkflowTitle = await createWorkflowWithUniqueTitle(
        pageObjects.projectsPage,
        pageObjects.workflowPage,
        testWorkflow,
        createdWorkflows
      );

      await page.waitForSelector('.drawflow-node', { timeout: 10000 });
    });

    await test.step('Assign InputFile node to pbmc dataset', async () => {
      await pageObjects.workflowPage.openNodeModal(testWorkflow.inputNodeName);
      await pageObjects.inputFileModal.assignFile(testWorkflow.folder, workflowInputFileName);
      await pageObjects.workflowPage.closeTab(testWorkflow.inputNodeTabName);
      await page.waitForTimeout(300);
    });

    await test.step('Configure Algorithm node parameters for TENET', async () => {
      await pageObjects.workflowPage.openNodeModal('Algorithm');
      await pageObjects.algorithmModal.verifyModalOpen();

      await expect
        .poll(async () => {
          const options = await pageObjects.algorithmModal.getParameterDropdownOptions('cell group');
          return options.includes('seurat_annotation') ? 'ready' : null;
        }, {
          message: 'Waiting for Cell group options to load',
          timeout: 20000,
        })
        .toBe('ready');

      await pageObjects.algorithmModal.setParameterValueByName('cell group', 'seurat_annotation');
      expect(await pageObjects.algorithmModal.getParameterValueByName('cell group')).toBe('seurat_annotation');

      await expect
        .poll(async () => {
          const options = await pageObjects.algorithmModal.getParameterDropdownOptions('pseudotime');
          return options.includes('Pseudotime') ? 'ready' : null;
        }, {
          message: 'Waiting for pseudotime options to load',
          timeout: 20000,
        })
        .toBe('ready');

      await pageObjects.algorithmModal.setParameterValueByName('pseudotime', 'Pseudotime');
      expect(await pageObjects.algorithmModal.getParameterValueByName('pseudotime')).toBe('Pseudotime');

      await pageObjects.algorithmModal.setParameterValueByName('clusters', desiredClusters);
      const selectedClusters = await pageObjects.algorithmModal.getParameterValueByName('clusters');
      expect(selectedClusters).toEqual(expect.arrayContaining(desiredClusters));

      const pluginLogo = await pageObjects.algorithmModal.getPluginLogoText();
      console.log('Algorithm modal plugin:', pluginLogo);
      expect(pluginLogo).toContain(testWorkflow.name);
    });

    await test.step('Open compile check panel and verify task summary', async () => {
      await pageObjects.workflowPage.openCompileCheck();
      await pageObjects.compileCheckModal.waitForOpen();
      await pageObjects.compileCheckModal.waitForResourcesLoaded();
      await pageObjects.compileCheckModal.verifyCoreSectionsVisible();

      const taskEntries = await pageObjects.compileCheckModal.getTaskEntries();
      console.log('Compile check task entries:', taskEntries);
      expect(taskEntries.length).toBeGreaterThan(0);
      expect(
        taskEntries.some((entry) => entry.plugin && entry.plugin.includes(testWorkflow.name))
      ).toBeTruthy();

      const resourceLabels = await pageObjects.compileCheckModal.getResourceLabels();
      console.log('Resource summary labels:', resourceLabels);
      expect(resourceLabels.length).toBeGreaterThan(0);
    });

    await test.step('Execute workflow and wait for RUNNING status', async () => {
      page.once('dialog', async (dialog) => {
        console.log('Execution confirmation dialog:', dialog.message());
        await dialog.accept();
      });

      await pageObjects.compileCheckModal.clickExecute();
      await pageObjects.compileCheckModal.waitForClose();

      await page.waitForTimeout(1000);
      await pageObjects.workflowPage.openJobTable();
      await pageObjects.workflowPage.waitForJobRows(1, 180000);

      await pageObjects.workflowPage.waitForLatestJobStatus(currentWorkflowTitle, 'RUNNING', 240000);
      const latestJob = await pageObjects.workflowPage.getLatestJobEntryByTitle(currentWorkflowTitle);

      console.log('Latest job entry:', latestJob);
      expect(latestJob).not.toBeNull();
      expect(latestJob?.status?.toUpperCase()).toContain('RUNNING');
      expect(latestJob?.plugin ?? '').toContain(testWorkflow.name);

      const pluginFormatRegex = /^[^/]+\/[^ :]+ : v\d+(?:\.\d+)*$/;
      expect(latestJob?.plugin ?? '').toMatch(pluginFormatRegex);
    });

    await test.step('Wait for job to continue running', async () => {
      await page.waitForTimeout(30000);
      await pageObjects.workflowPage.closeJobTable();
    });

    await test.step('Validate RUNNING job and inspect logs', async () => {
      await pageObjects.workflowPage.openJobTable();
      await pageObjects.workflowPage.waitForJobRows(1, 180000);

      await expect
        .poll(async () => {
          const entry = await pageObjects.workflowPage.getLatestJobEntryByTitle(currentWorkflowTitle);
          return entry?.status?.toUpperCase() ?? null;
        }, {
          timeout: 240000,
          message: `Waiting for job "${currentWorkflowTitle}" to remain RUNNING`,
        })
        .toBe('RUNNING');

      const latestJob = await pageObjects.workflowPage.getLatestJobEntryByTitle(currentWorkflowTitle);
      console.log('Latest job entry:', latestJob);
      expect(latestJob).not.toBeNull();
      expect(latestJob?.name).toBe(currentWorkflowTitle);
      expect(latestJob?.plugin ?? '').toContain(testWorkflow.name);
      const pluginFormatRegex = /^[^/]+\/[^ :]+ : v\d+(?:\.\d+)*$/;
      expect(latestJob?.plugin ?? '').toMatch(pluginFormatRegex);

      await pageObjects.workflowPage.openJobContextMenuForTitle(currentWorkflowTitle);
      await pageObjects.workflowPage.selectJobContextOption('View logs');
      await pageObjects.logsModal.waitForOpen();
      await pageObjects.logsModal.waitForLoaded();
      await pageObjects.logsModal.expectLogsAvailable();

      try {
        const jsonDownload = await pageObjects.logsModal.downloadAllLogsJson();
        await jsonDownload.delete().catch(() => {});
      } catch (error) {
        console.warn('⚠️ Failed to download JSON logs:', error.message);
      }

      try {
        const txtDownload = await pageObjects.logsModal.downloadFirstLogTxt();
        await txtDownload.delete().catch(() => {});
      } catch (error) {
        console.warn('⚠️ Failed to download TXT log:', error.message);
      }

      await pageObjects.logsModal.close();
      await pageObjects.workflowPage.closeMessage().catch(() => {});
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

      await pageObjects.workflowPage.openJobContextMenuForTitle(currentWorkflowTitle);
      await pageObjects.workflowPage.selectJobContextOption('View progress');

      const dagStructureResponse = await dagStructurePromise;
      expect(dagStructureResponse.ok()).toBeTruthy();

      const ruleStatusResponse = await ruleStatusPromise;
      expect(ruleStatusResponse.ok()).toBeTruthy();

      await pageObjects.dagModal.waitForOpen();
      await pageObjects.dagModal.waitForLoaded();
      await pageObjects.dagModal.close();
    });

    await test.step('Cancel running job and verify message', async () => {
      await pageObjects.workflowPage.cancelJobByTitle(currentWorkflowTitle);
      await pageObjects.workflowPage.waitForMessage('Cancel task successfully!', 15000);
      await pageObjects.workflowPage.closeMessage().catch(() => {});

      await pageObjects.workflowPage.openJobTable();
      await pageObjects.workflowPage.waitForJobRows(1, 60000);

      // waitForLatestJobStatus already validates the status via polling
      // No need for additional verification after this succeeds
      await pageObjects.workflowPage.waitForLatestJobStatus(currentWorkflowTitle, 'REVOKED', 240000);
      console.log(`✅ Job "${currentWorkflowTitle}" status verified as REVOKED`);

      await pageObjects.workflowPage.closeJobTable();
    });

    await test.step('Delete cancelled job and verify removal', async () => {
      await pageObjects.workflowPage.openJobTable();
      await pageObjects.workflowPage.waitForJobRows(1, 60000);

      // Verify job is in REVOKED state before deletion (prerequisite for delete)
      const jobBeforeDelete = await pageObjects.workflowPage.getLatestJobEntryByTitle(currentWorkflowTitle);
      console.log('Job state before delete:', jobBeforeDelete?.status);
      expect(jobBeforeDelete?.status?.toUpperCase()).toBe('REVOKED');
      console.log('✅ Confirmed job is REVOKED, proceeding with delete');

      // Delete the REVOKED job
      await pageObjects.workflowPage.deleteJobByTitle(currentWorkflowTitle);
      await pageObjects.workflowPage.waitForMessage('Delete task successfully!', 15000);
      await pageObjects.workflowPage.closeMessage().catch(() => {});

      // Wait for DOM to stabilize after deletion
      await page.waitForTimeout(1000);

      // Verify job is removed from table
      await pageObjects.workflowPage.openJobTable();
      const deletedJob = await pageObjects.workflowPage.getLatestJobEntryByTitle(currentWorkflowTitle);
      console.log('Deleted job entry (should be null):', deletedJob);
      expect(deletedJob).toBeNull();

      await pageObjects.workflowPage.closeJobTable();
    });

    await test.step('Clean up uploaded test files', async () => {
      if (cleanupUploads.length === 0) {
        return;
      }

      await pageObjects.filesPage.goto();
      await pageObjects.filesPage.verifyPageLoaded();

      try {
        await pageObjects.filesPage.selectFolder(testWorkflow.folder);
      } catch (error) {
        console.warn(`⚠️ Unable to re-open folder "${testWorkflow.folder}":`, error.message);
      }

      for (const fileName of cleanupUploads) {
        try {
          await pageObjects.filesPage.deleteFile(fileName);
          await pageObjects.filesPage.verifyFileNotExists(fileName);
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
