// frontend/tests/e2e/workflows/support/workflow-setup.js

import { ProjectsPage } from '../../pages/ProjectsPage.js';
import { WorkflowPage } from '../../pages/WorkflowPage.js';
import { FilesPage } from '../../pages/FilesPage.js';
import { InputFileModal } from '../../pages/modals/InputFileModal.js';
import { AlgorithmModal } from '../../pages/modals/AlgorithmModal.js';
import { DataTableModal } from '../../pages/modals/DataTableModal.js';
import { ScatterPlotModal } from '../../pages/modals/ScatterPlotModal.js';
import { LogsModal } from '../../pages/modals/LogsModal.js';
import { DagModal } from '../../pages/modals/DagModal.js';
import { CompileCheckModal } from '../../pages/modals/CompileCheckModal.js';
import { ResultFilesModal } from '../../pages/modals/ResultFilesModal.js';
import { VisualizationModal } from '../../pages/modals/VisualizationModal.js';
import { generateUniqueFileName, generateUniqueWorkflowTitle } from './workflow-constants.js';

/**
 * Initialize all page objects for workflow tests
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @returns {Object} Object containing all initialized page objects
 */
export function setupPageObjects(page) {
  return {
    projectsPage: new ProjectsPage(page),
    workflowPage: new WorkflowPage(page),
    filesPage: new FilesPage(page),
    inputFileModal: new InputFileModal(page),
    algorithmModal: new AlgorithmModal(page),
    dataTableModal: new DataTableModal(page),
    scatterPlotModal: new ScatterPlotModal(page),
    logsModal: new LogsModal(page),
    dagModal: new DagModal(page),
    compileCheckModal: new CompileCheckModal(page),
    resultFilesModal: new ResultFilesModal(page),
    visualizationModal: new VisualizationModal(page),
  };
}

/**
 * Setup test files: always upload files with unique names for test isolation
 * @param {FilesPage} filesPage - FilesPage instance
 * @param {Object} testWorkflow - Test workflow configuration
 * @param {Array<string>} uploadedFiles - Array to track uploaded files for cleanup
 * @param {Object} options - Additional options
 * @param {boolean} options.uploadSecondFile - Whether to upload second file (for file change test)
 * @returns {Promise<Object>} Returns { uploadedFileName, secondFileName } with uploaded filenames
 */
export async function setupTestFiles(filesPage, testWorkflow, uploadedFiles, options = {}) {
  await filesPage.goto();
  await filesPage.verifyPageLoaded();

  // Wait for page to stabilize
  await filesPage.page.waitForTimeout(1000);

  // Note: Folder selection removed as it's not necessary for file upload
  // Files can be uploaded directly without explicitly selecting folder first

  // Always upload with unique filename to ensure test isolation
  const uniqueFileName = generateUniqueFileName(testWorkflow.expectedFile);
  console.log(`📤 Uploading test file with unique name: ${uniqueFileName}`);

  const { uploadedFileName } = await filesPage.uploadFile('test_data.h5ad', {
    targetFileName: uniqueFileName,
  });
  uploadedFiles.push(uploadedFileName);
  await filesPage.waitForUploadComplete();

  // Verify file appears in table with extended timeout (15s)
  await filesPage.verifyFileExists(uploadedFileName, 15000);
  console.log(`✅ Test file uploaded: ${uploadedFileName}`);

  let secondFileName = null;

  // Upload second file if requested (for file change tests)
  if (options.uploadSecondFile) {
    const uniqueSecondFileName = generateUniqueFileName('demo.h5ad');
    console.log(`📤 Uploading second test file with unique name: ${uniqueSecondFileName}`);

    const { uploadedFileName: secondFile } = await filesPage.uploadFile('demo.h5ad', {
      targetFileName: uniqueSecondFileName,
      timeout: 60000,
    });
    uploadedFiles.push(secondFile);
    await filesPage.waitForUploadComplete();

    // Verify second file appears in table with extended timeout (15s)
    await filesPage.verifyFileExists(secondFile, 15000);
    console.log(`✅ Second test file uploaded: ${secondFile}`);
    secondFileName = secondFile;
  }

  return { uploadedFileName, secondFileName };
}

/**
 * Cleanup uploaded test files - always attempt cleanup even if array is empty
 * @param {FilesPage} filesPage - FilesPage instance
 * @param {Array<string>} uploadedFiles - Array of uploaded files to delete
 * @param {string} folder - Folder name to select before cleanup (default: 'data')
 * @returns {Promise<void>}
 */
export async function cleanupTestFiles(filesPage, uploadedFiles, folder = 'data') {
  if (uploadedFiles.length === 0) {
    console.log('🧹 No files to cleanup (uploadedFiles array is empty)');
    return;
  }

  console.log(`🧹 [Cleanup] Starting file cleanup at ${new Date().toISOString()}`);
  console.log(`📋 Files to cleanup: ${uploadedFiles.length}`);
  uploadedFiles.forEach((file, index) => {
    console.log(`   ${index + 1}. ${file}`);
  });

  await filesPage.goto();
  await filesPage.verifyPageLoaded();

  // Wait for page to stabilize after navigation
  await filesPage.page.waitForTimeout(1000);

  // Check if already in correct folder
  let currentFolder = '';
  try {
    currentFolder = await filesPage.getCurrentFolder();
    console.log(`Current folder: ${currentFolder}`);
  } catch (error) {
    console.log(`⚠️ Unable to get current folder:`, error.message);
  }

  // Select folder before cleanup with retry (skip if already in correct folder)
  if (currentFolder !== folder) {
    let folderSelected = false;
    for (let attempt = 1; attempt <= 3; attempt++) {
      try {
        await filesPage.selectFolder(folder);
        console.log(`✓ Selected folder: ${folder} (attempt ${attempt})`);
        folderSelected = true;
        break;
      } catch (error) {
        if (attempt === 3) {
          console.warn(`⚠️ Unable to select folder "${folder}" after ${attempt} attempts:`, error.message);
          console.warn(`⚠️ Skipping file cleanup as folder cannot be selected`);
          uploadedFiles.length = 0;
          return; // Early return if folder selection fails
        }
        console.log(`⚠️ Folder selection attempt ${attempt} failed, retrying...`);
        await filesPage.page.waitForTimeout(1000);
      }
    }
  } else {
    console.log(`✓ Already in folder: ${folder}`);
  }

  let successCount = 0;
  let failCount = 0;

  for (const fileName of uploadedFiles) {
    const startTime = Date.now();
    try {
      // Delete the file
      await filesPage.deleteFile(fileName);

      // Verify deletion
      await filesPage.verifyFileNotExists(fileName);

      const elapsed = Date.now() - startTime;
      console.log(`✅ [${elapsed}ms] Deleted and verified: ${fileName}`);
      successCount++;
    } catch (error) {
      const elapsed = Date.now() - startTime;
      console.log(`❌ Failed to delete: ${fileName}`);
      console.log(`   Error: ${error.message}`);
      console.log(`   Type: ${error.constructor.name}`);
      console.log(`   Time elapsed: ${elapsed}ms`);
      failCount++;
    }
  }

  console.log(`🧹 [Cleanup] File cleanup complete: ${successCount} succeeded, ${failCount} failed (${uploadedFiles.length} total)`);

  // Warn if any deletions failed
  if (failCount > 0) {
    console.warn(`⚠️ WARNING: ${failCount} file(s) failed to cleanup. May cause test interference.`);
  }

  // Clear the array for next test
  uploadedFiles.length = 0;
}

/**
 * Setup test file for file change test - uploads with unique filename
 * @param {FilesPage} filesPage - FilesPage instance
 * @param {Object} testWorkflow - Test workflow configuration
 * @param {Array<string>} uploadedFiles - Array to track uploaded files for cleanup
 * @returns {Promise<string>} Returns uploaded filename
 */
export async function setupFileChangeTest(filesPage, testWorkflow, uploadedFiles) {
  await filesPage.goto();

  // Note: Folder selection removed as it's not necessary for file upload
  // Files can be uploaded directly without explicitly selecting folder first

  // Always use unique filename for test isolation
  const uniqueFileName = generateUniqueFileName(testWorkflow.expectedFile);
  console.log(`📤 Uploading file for file change test with unique name: ${uniqueFileName}`);

  const { uploadedFileName } = await filesPage.uploadFile('test_data.h5ad', {
    targetFileName: uniqueFileName,
    timeout: 60000,
  });

  // Track uploaded file for cleanup
  uploadedFiles.push(uploadedFileName);

  await filesPage.waitForUploadComplete();

  // Verify file appears in table with extended timeout (15s)
  await filesPage.verifyFileExists(uploadedFileName, 15000);
  console.log(`✅ File change test file uploaded: ${uploadedFileName}`);

  return uploadedFileName;
}

/**
 * Create workflow with unique title to ensure test isolation
 * @param {ProjectsPage} projectsPage - ProjectsPage instance
 * @param {WorkflowPage} workflowPage - WorkflowPage instance
 * @param {Object} testWorkflow - Test workflow configuration
 * @param {Array<string>} createdWorkflows - Array to track created workflows for cleanup
 * @returns {Promise<string>} Returns unique workflow title
 */
export async function createWorkflowWithUniqueTitle(projectsPage, workflowPage, testWorkflow, createdWorkflows) {
  // 1. Create workflow from template (starts as "Untitled" and auto-saved)
  await projectsPage.clickNewWorkflow();
  await projectsPage.selectPluginTemplate(testWorkflow.name);
  await workflowPage.verifyPageLoaded();

  // 2. Generate unique title
  const uniqueTitle = generateUniqueWorkflowTitle(testWorkflow.name);

  // 3. Update workflow title
  await workflowPage.updateWorkflowTitle(uniqueTitle);

  // 4. Save workflow (ensures title is persisted to database)
  await workflowPage.saveWorkflow();
  await workflowPage.closeMessage().catch(() => {}); // Close success message

  // 5. Track for cleanup
  createdWorkflows.push(uniqueTitle);

  console.log(`✅ Workflow created with unique title: ${uniqueTitle}`);

  return uniqueTitle;
}

/**
 * Cleanup created workflows - delete from database via Projects page
 * @param {ProjectsPage} projectsPage - ProjectsPage instance
 * @param {Array<string>} createdWorkflows - Array of created workflow titles to delete
 * @returns {Promise<void>}
 */
export async function cleanupWorkflows(projectsPage, createdWorkflows) {
  if (createdWorkflows.length === 0) {
    console.log('🧹 No workflows to cleanup (createdWorkflows array is empty)');
    return;
  }

  console.log(`🧹 [Cleanup] Starting workflow cleanup at ${new Date().toISOString()}`);
  console.log(`📋 Workflows to cleanup: ${createdWorkflows.length}`);
  createdWorkflows.forEach((title, index) => {
    console.log(`   ${index + 1}. ${title}`);
  });

  // Navigate to Projects page for deletion
  await projectsPage.goto();
  await projectsPage.verifyPageLoaded();

  let successCount = 0;
  let failCount = 0;

  for (const workflowTitle of createdWorkflows) {
    const startTime = Date.now();
    try {
      // Delete the workflow
      await projectsPage.deleteWorkflow(workflowTitle);

      // Verify deletion
      await projectsPage.verifyWorkflowNotExists(workflowTitle);

      const elapsed = Date.now() - startTime;
      console.log(`✅ [${elapsed}ms] Deleted and verified: ${workflowTitle}`);
      successCount++;
    } catch (error) {
      const elapsed = Date.now() - startTime;
      console.log(`❌ Failed to delete: ${workflowTitle}`);
      console.log(`   Error: ${error.message}`);
      console.log(`   Type: ${error.constructor.name}`);
      console.log(`   Time elapsed: ${elapsed}ms`);
      failCount++;
    }
  }

  console.log(`🧹 [Cleanup] Workflow cleanup complete: ${successCount} succeeded, ${failCount} failed (${createdWorkflows.length} total)`);

  // Warn if any deletions failed
  if (failCount > 0) {
    console.warn(`⚠️ WARNING: ${failCount} workflow(s) failed to cleanup. May cause test interference.`);
  }

  // Clear the array for next test
  createdWorkflows.length = 0;
}
