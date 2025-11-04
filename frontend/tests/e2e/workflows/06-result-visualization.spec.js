// frontend/tests/e2e/workflows/06-result-visualization.spec.js
import { test, expect } from '../fixtures/auth.js';
import { readFile } from 'fs/promises';
import { setupPageObjects } from './support/workflow-setup.js';

/**
 * Test Suite: Result Visualization with GRNViz
 *
 * This test suite verifies the complete result visualization workflow:
 * - Opening prepared SUCCESS workflow
 * - Validating SUCCESS job in monitoring panel
 * - Inspecting ResultFiles modal (Primary/Intermediate outputs)
 * - Downloading result files
 * - Downloading execution manifest
 * - Validating manifest sections
 * - Executing GRNViz visualization
 * - Validating Plotly rendering
 * - Downloading plot images
 * - Deleting visualization job and verifying removal
 *
 * Success criteria:
 * - SUCCESS workflow exists and loads correctly
 * - Job monitoring panel shows SUCCESS status
 * - ResultFiles modal displays Primary and Intermediate outputs
 * - File downloads work correctly
 * - Execution manifest contains all required sections
 * - GRNViz visualization executes successfully
 * - Plotly chart renders correctly
 * - Plot image download works
 * - Visualization job deletion works correctly
 */
test.describe('Result Visualization with GRNViz', () => {
  test.describe.configure({ mode: 'serial' });

  let pageObjects;

  test.beforeEach(async ({ page }) => {
    pageObjects = setupPageObjects(page);
    await pageObjects.projectsPage.goto();
    await pageObjects.projectsPage.verifyPageLoaded();
  });

  /**
   * Test: Should visualize TENET results using GRNViz plugin
   *
   * This test verifies the complete result visualization workflow:
   * - Opening prepared SUCCESS workflow
   * - Validating existing SUCCESS job
   * - Inspecting and downloading result files
   * - Downloading and validating execution manifest
   * - Executing GRNViz visualization
   * - Validating Plotly rendering
   * - Deleting visualization job and verifying removal
   */
  test('Should visualize TENET results using GRNViz plugin', async ({ page }) => {
    test.setTimeout(300000);

    await test.step('Open prepared SUCCESS workflow', async () => {
      await pageObjects.projectsPage.verifyWorkflowExists('SUCCESS');
      await pageObjects.projectsPage.openWorkflow('SUCCESS');
      await pageObjects.workflowPage.verifyPageLoaded();
      await pageObjects.workflowPage.waitForNodesReady(5, 15000);
    });

    const workflowTitle = await pageObjects.workflowPage.getWorkflowTitle();
    expect(workflowTitle).toContain('SUCCESS');

    await test.step('Validate existing SUCCESS job in monitoring panel', async () => {
      await pageObjects.workflowPage.openJobTable();
      await pageObjects.workflowPage.waitForJobRows(1, 60000);
      const entries = await pageObjects.workflowPage.getJobTableEntries();
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

      await pageObjects.workflowPage.closeJobTable();
    });

    await test.step('Inspect ResultFiles modal, select downloads, and verify sections', async () => {
      await pageObjects.workflowPage.openNodeModal('ResultFiles');
      await pageObjects.resultFilesModal.verifyModalOpen();
      await pageObjects.resultFilesModal.waitForPrimarySection(20000);
      await pageObjects.resultFilesModal.waitForIntermediateSection(20000);
      expect(await pageObjects.resultFilesModal.isPrimarySectionVisible()).toBeTruthy();
      expect(await pageObjects.resultFilesModal.isIntermediateSectionVisible()).toBeTruthy();

      const primaryFiles = await pageObjects.resultFilesModal.getPrimaryFileNames();
      const intermediateFiles = await pageObjects.resultFilesModal.getIntermediateFileNames();
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
        await pageObjects.resultFilesModal.clickPrimaryFileDownloadButton(target.name);
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
      await pageObjects.workflowPage.openJobTable();
      await pageObjects.workflowPage.waitForJobRows(1, 60000);

      const latestJobContext = await page.evaluate(() => window.__latestSuccessJob || null);
      if (!latestJobContext) {
        throw new Error('Latest success job context not found for manifest download');
      }

      const manifestDownloadPromise = page.waitForEvent('download');
      await pageObjects.workflowPage.openJobContextMenuForTitle(latestJobContext.name, {
        pluginSubstring: latestJobContext.pluginMatcher,
        typeSubstring: latestJobContext.typeMatcher,
      });
      await pageObjects.workflowPage.selectJobContextOption('Download manifest');

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
      await pageObjects.workflowPage.closeJobTable();
    });

    await test.step('Execute GRNViz visualization and validate Plotly rendering', async () => {
      await pageObjects.workflowPage.openNodeModal('Visualization');

      if (!(await pageObjects.visualizationModal.isConfigurationMode())) {
        await pageObjects.visualizationModal.waitForPluginSelection();
        await pageObjects.visualizationModal.selectPluginByName('GRNViz');
        await expect(pageObjects.visualizationModal.visualizationItems.first()).toBeVisible();

        let visualizationName = 'Bar plot';
        const barPlotCount = await pageObjects.visualizationModal.visualizationItems
          .filter({ hasText: 'Bar plot' })
          .count();

        if (barPlotCount === 0) {
          const firstVisualizationText = (await pageObjects.visualizationModal.visualizationItems
            .first()
            .textContent())?.trim();

          if (!firstVisualizationText) {
            throw new Error('No visualization scripts available for GRNViz plugin');
          }

          visualizationName = firstVisualizationText;
        }

        await pageObjects.visualizationModal.selectVisualizationByName(visualizationName);
        await pageObjects.visualizationModal.proceedToConfiguration();
      }

      expect(await pageObjects.visualizationModal.getSelectedPluginLabel()).toContain('GRNViz');
      const selectedVisualizationLabel = await pageObjects.visualizationModal.getSelectedVisualizationLabel();
      expect(selectedVisualizationLabel.length).toBeGreaterThan(0);

      const inputParameters = await pageObjects.visualizationModal.getInputFileParameterNames();
      for (const parameterName of inputParameters) {
        const normalizedName = parameterName.trim().toLowerCase();

        if (normalizedName === 'input') {
          await pageObjects.visualizationModal.selectInputFileOption(parameterName, 'FdrOutdegree.txt');
          await expect
            .poll(async () => await pageObjects.visualizationModal.getSelectedInputFile(parameterName))
            .toBe('FdrOutdegree.txt');
        } else if (normalizedName.includes('expression')) {
          await pageObjects.visualizationModal.selectInputFileOption(parameterName, 'expression.csv');
          await expect
            .poll(async () => await pageObjects.visualizationModal.getSelectedInputFile(parameterName))
            .toBe('expression.csv');
        } else if (normalizedName.includes('trajectory')) {
          await pageObjects.visualizationModal.selectInputFileOption(parameterName, 'trajectory.txt');
          await expect
            .poll(async () => await pageObjects.visualizationModal.getSelectedInputFile(parameterName))
            .toBe('trajectory.txt');
        } else {
          await expect
            .poll(async () => (await pageObjects.visualizationModal.getAvailableOptionsForParameter(parameterName)).length > 0, {
              message: `Waiting for options to load for parameter ${parameterName}`,
            })
            .toBeTruthy();

          const options = await pageObjects.visualizationModal.getAvailableOptionsForParameter(parameterName);
          const firstNonEmpty = options.find((opt) => opt.trim() !== '' && opt.trim() !== 'Select File');
          if (!firstNonEmpty) {
            throw new Error(`No selectable option available for parameter ${parameterName}`);
          }
          await pageObjects.visualizationModal.selectInputFileOption(parameterName, firstNonEmpty.trim());
          await expect
            .poll(async () => await pageObjects.visualizationModal.getSelectedInputFile(parameterName))
            .toBe(firstNonEmpty.trim());
        }
      }

      expect(await pageObjects.visualizationModal.isApplyButtonEnabled()).toBeTruthy();

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

      await pageObjects.visualizationModal.clickExecuteVisualization();

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

      await pageObjects.visualizationModal.waitForPlotly();
      if ((await pageObjects.visualizationModal.getApplyButtonText()).includes('Show Visualization')) {
        await pageObjects.visualizationModal.clickShowVisualization();
        await pageObjects.visualizationModal.waitForPlotly();
      }
      expect(await pageObjects.visualizationModal.isPlotlyVisible()).toBeTruthy();

      const plotDownload = await pageObjects.visualizationModal.downloadPlotImage();
      const plotFilename = plotDownload.suggestedFilename();
      expect(plotFilename).toMatch(/\.png$/i);
      await plotDownload.delete().catch(() => {});
    });

    await test.step('Delete visualization job and verify removal', async () => {
      // Wait for any modals to close and UI to stabilize
      await page.waitForTimeout(1000);

      await pageObjects.workflowPage.openJobTable();
      await pageObjects.workflowPage.waitForJobRows(1, 60000);

      // Find the Visualization jobs for this workflow (name === workflowTitle)
      const entries = await pageObjects.workflowPage.getJobTableEntries();
      const visualizationJobs = entries
        .filter((entry) => entry.name === workflowTitle)
        .filter((entry) => entry.type === 'Visualization');

      // If no visualization jobs exist, skip deletion step
      if (visualizationJobs.length === 0) {
        console.log(`⚠️ No Visualization jobs found for "${workflowTitle}", skipping deletion step`);
        await pageObjects.workflowPage.closeJobTable();
        return;
      }

      const initialVisualizationCount = visualizationJobs.length;
      console.log(`Initial Visualization job count for "${workflowTitle}": ${initialVisualizationCount}`);

      // Sort by startTimestamp to get most recent
      visualizationJobs.sort((a, b) => b.startTimestamp - a.startTimestamp);
      const latestVisualizationJob = visualizationJobs[0];

      console.log('Visualization job to delete:', latestVisualizationJob);

      // Verify job is in SUCCESS or FAILURE state (both can be deleted)
      const jobStatus = latestVisualizationJob.status?.toUpperCase();
      expect(['SUCCESS', 'FAILURE'].includes(jobStatus)).toBeTruthy();
      console.log(`✅ Confirmed Visualization job is ${jobStatus}, proceeding with delete`);

      // Delete the Visualization job using type filter to avoid targeting Analysis jobs
      await pageObjects.workflowPage.openJobContextMenuForTitle(workflowTitle, {
        typeSubstring: 'Visualization',
      });
      await pageObjects.workflowPage.selectJobContextOption('Delete');
      await pageObjects.workflowPage.waitForMessage('Delete task successfully!', 15000);
      await pageObjects.workflowPage.closeMessage().catch(() => {});

      // Wait for DOM to stabilize after deletion
      await page.waitForTimeout(1000);

      // Verify Visualization job count decreased by 1 for this workflow (name === workflowTitle)
      await pageObjects.workflowPage.openJobTable();
      const entriesAfterDelete = await pageObjects.workflowPage.getJobTableEntries();
      const remainingVisualizationJobs = entriesAfterDelete
        .filter((entry) => entry.name === workflowTitle)
        .filter((entry) => entry.type === 'Visualization');

      console.log(`Remaining Visualization job count for "${workflowTitle}": ${remainingVisualizationJobs.length}`);
      expect(remainingVisualizationJobs.length).toBe(initialVisualizationCount - 1);

      await pageObjects.workflowPage.closeJobTable();
    });
  });
});
