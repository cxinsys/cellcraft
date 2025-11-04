// frontend/tests/e2e/authenticated-workflow.spec.js
import { test, expect } from './fixtures/auth.js';
import { DatasetsPage } from './pages/DatasetsPage.js';
import { FilesPage } from './pages/FilesPage.js';
import { ProjectsPage } from './pages/ProjectsPage.js';
import {
  goToDatasets,
  goToFiles,
  goToProjects,
  verifyAuthenticated,
} from './utils/navigation.js';

/**
 * Generate a unique filename with timestamp and random suffix
 * @param {string} baseFileName - Base filename (e.g., "test_data.h5ad")
 * @returns {string} Unique filename (e.g., "test_data_1730419234567_abc123.h5ad")
 */
function generateUniqueFileName(baseFileName) {
  const timestamp = Date.now();
  const randomId = Math.random().toString(36).substring(2, 8);
  const extension = baseFileName.substring(baseFileName.lastIndexOf('.'));
  const baseName = baseFileName.substring(0, baseFileName.lastIndexOf('.'));

  return `${baseName}_${timestamp}_${randomId}${extension}`;
}

/**
 * Clean up all test files from the files page
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
async function cleanupTestFiles(page) {
  const filesPage = new FilesPage(page);

  try {
    // Navigate directly to Files page (works from any page)
    await page.goto('/files');
    await page.waitForLoadState('networkidle');

    // Get all files in the current folder
    const allFiles = await filesPage.getFileList();
    const testFiles = allFiles.filter((file) => file.name.startsWith('test_'));

    if (testFiles.length > 0) {
      console.log(`🧹 Cleanup: Found ${testFiles.length} test files to remove`);

      // Delete all files that start with 'test_'
      for (const file of testFiles) {
        try {
          await filesPage.deleteFile(file.name);
          console.log(`  ✓ Removed ${file.name}`);
        } catch (error) {
          console.warn(`  ⚠ Could not delete ${file.name}:`, error.message);
        }
      }

      // Wait for deletions to complete
      await page.waitForLoadState('networkidle');
      console.log(`✅ Cleanup complete`);
    }
  } catch (error) {
    console.warn('Cleanup failed:', error.message);
  }
}

/**
 * Test Suite: Project Initialization Workflow
 *
 * This test suite validates the complete workflow for initializing a new project
 * using tutorial datasets and creating a TENET-based workflow.
 *
 * Test Scenario:
 * 1. Download tutorial dataset (PBMCLight1000.h5ad) from Datasets page
 * 2. Upload various file types (H5AD, CSV, TXT) via Files page
 * 3. Verify all uploaded files appear in the file list
 * 4. Create a new workflow project using TENET template
 *
 * NOTE: This test suite runs serially to avoid file conflicts in shared folders
 */
test.describe.serial('Project Initialization Workflow', () => {
  // Track uploaded files for cleanup
  let uploadedFiles = [];

  test.beforeEach(async ({ page }) => {
    // Verify user is authenticated before each test
    await page.goto('/projects');
    await verifyAuthenticated(page);

    // Reset uploaded files tracking
    uploadedFiles = [];
  });

  test.afterEach(async ({ page }) => {
    // Cleanup: Delete all uploaded test files
    if (uploadedFiles.length > 0) {
      const filesPage = new FilesPage(page);

      try {
        // Navigate directly to Files page (works from any page including workflow)
        await page.goto('/files');
        await page.waitForLoadState('networkidle');

        for (const fileName of uploadedFiles) {
          try {
            await filesPage.deleteFile(fileName);
            console.log(`✓ Cleaned up test file: ${fileName}`);
          } catch (error) {
            console.warn(`Could not delete ${fileName}:`, error.message);
          }
        }
      } catch (error) {
        console.warn('Cleanup failed:', error.message);
      }
    }
  });

  test('should complete full project initialization flow', async ({ page }) => {
    // ============================================
    // Pre-cleanup: Remove old test files
    // ============================================
    await cleanupTestFiles(page);

    // ============================================
    // Step 1: Download Tutorial Dataset
    // ============================================
    await test.step('Download PBMCLight1000 dataset', async () => {
      const datasetsPage = new DatasetsPage(page);

      // Navigate to Datasets page
      await goToDatasets(page);
      await datasetsPage.verifyPageLoaded();

      // Search for PBMC dataset
      await datasetsPage.searchDataset('pbmc_light_1000');
      await datasetsPage.verifyDatasetVisible('pbmc_light_1000.h5ad');

      // Download the dataset
      const download = await datasetsPage.downloadDataset(
        'pbmc_light_1000.h5ad'
      );

      // Verify download completed successfully
      await datasetsPage.verifyDownload(download, 'pbmc_light_1000.h5ad');

      console.log('✓ Dataset download completed successfully');
    });

    // ============================================
    // Step 2: Upload Test Files with Unique Names
    // ============================================
    await test.step('Upload H5AD file', async () => {
      const filesPage = new FilesPage(page);

      // Navigate to Files page
      await goToFiles(page);
      await filesPage.verifyPageLoaded();

      // Verify we're in the 'data' folder
      const currentFolder = await filesPage.getCurrentFolder();
      expect(currentFolder).toBe('data');

      // Generate unique filename and upload
      const uniqueFileName = generateUniqueFileName('test_data.h5ad');
      const { uploadedFileName } = await filesPage.uploadFile('test_data.h5ad', {
        targetFileName: uniqueFileName,
      });
      await filesPage.waitForUploadComplete();

      // Verify file was uploaded
      await filesPage.verifyFileExists(uploadedFileName);
      uploadedFiles.push(uploadedFileName); // Track for cleanup

      console.log(`✓ H5AD file uploaded successfully as ${uploadedFileName}`);
    });

    await test.step('Upload CSV file', async () => {
      const filesPage = new FilesPage(page);

      // Generate unique filename and upload
      const uniqueFileName = generateUniqueFileName('test_sample.csv');
      const { uploadedFileName } = await filesPage.uploadFile('test_sample.csv', {
        targetFileName: uniqueFileName,
      });
      await filesPage.waitForUploadComplete();

      // Verify file was uploaded
      await filesPage.verifyFileExists(uploadedFileName);
      uploadedFiles.push(uploadedFileName); // Track for cleanup

      console.log(`✓ CSV file uploaded successfully as ${uploadedFileName}`);
    });

    await test.step('Upload TXT file', async () => {
      const filesPage = new FilesPage(page);

      // Generate unique filename and upload
      const uniqueFileName = generateUniqueFileName('test_genes.txt');
      const { uploadedFileName } = await filesPage.uploadFile('test_genes.txt', {
        targetFileName: uniqueFileName,
      });
      await filesPage.waitForUploadComplete();

      // Verify file was uploaded
      await filesPage.verifyFileExists(uploadedFileName);
      uploadedFiles.push(uploadedFileName); // Track for cleanup

      console.log(`✓ TXT file uploaded successfully as ${uploadedFileName}`);
    });

    // ============================================
    // Step 3: Verify File List
    // ============================================
    await test.step('Verify all uploaded files in file list', async () => {
      const filesPage = new FilesPage(page);

      // Get the complete file list
      const fileList = await filesPage.getFileList();
      console.log('Current file list:', fileList);

      // Verify minimum number of files
      const fileCount = await filesPage.getFileCount();
      expect(fileCount).toBeGreaterThanOrEqual(3);

      // Verify each uploaded file is present
      for (const fileName of uploadedFiles) {
        const isPresent = await filesPage.isFilePresent(fileName);
        expect(isPresent).toBe(true);
      }

      console.log('✓ All files verified in file list');
    });

    // ============================================
    // Step 4: Create Workflow with TENET Template
    // ============================================
    await test.step('Create workflow project with TENET template', async () => {
      const projectsPage = new ProjectsPage(page);

      // Navigate to Projects page
      await goToProjects(page);
      await projectsPage.verifyPageLoaded();

      // Click New Workflow to open plugin selection modal
      await projectsPage.clickNewWorkflow();

      // Verify TENET plugin is available
      await projectsPage.verifyPluginAvailable('TENET');

      // Get list of all available plugins for debugging
      const availablePlugins = await projectsPage.getAvailablePlugins();
      console.log('Available plugin templates:', availablePlugins);

      // Select TENET template
      await projectsPage.selectPluginTemplate('TENET');

      // Verify redirect to workflow page
      await expect(page).toHaveURL(/.*\/workflow.*/);

      console.log('✓ TENET workflow created successfully');
    });
  });

  /**
   * Additional test: Verify file operations
   */
  test('should upload and delete files successfully', async ({ page }) => {
    // Pre-cleanup: Remove old test files
    await cleanupTestFiles(page);

    const filesPage = new FilesPage(page);
    let uploadedFileName;

    await test.step('Navigate to Files page', async () => {
      await goToFiles(page);
      await filesPage.verifyPageLoaded();
    });

    await test.step('Upload a test file', async () => {
      const initialCount = await filesPage.getFileCount();

      // Generate unique filename and upload
      const uniqueFileName = generateUniqueFileName('test_sample.csv');
      const result = await filesPage.uploadFile('test_sample.csv', {
        targetFileName: uniqueFileName,
      });
      uploadedFileName = result.uploadedFileName;
      await filesPage.waitForUploadComplete();

      const newCount = await filesPage.getFileCount();
      expect(newCount).toBe(initialCount + 1);

      console.log(`✓ File uploaded successfully as ${uploadedFileName}`);
    });

    await test.step('Delete the uploaded file', async () => {
      await filesPage.verifyFileExists(uploadedFileName);
      await filesPage.deleteFile(uploadedFileName);
      await filesPage.verifyFileNotExists(uploadedFileName);

      console.log(`✓ File ${uploadedFileName} deleted successfully`);
    });
  });

  /**
   * Additional test: Verify dataset search functionality
   */
  test('should search and filter datasets correctly', async ({ page }) => {
    const datasetsPage = new DatasetsPage(page);

    await test.step('Navigate to Datasets page', async () => {
      await goToDatasets(page);
      await datasetsPage.verifyPageLoaded();
    });

    await test.step('Search for PBMC dataset', async () => {
      await datasetsPage.searchDataset('pbmc');

      const visibleDatasets = await datasetsPage.getVisibleDatasets();
      console.log('Filtered datasets:', visibleDatasets);

      // Verify PBMC dataset is in results
      const hasPBMC = visibleDatasets.some((title) =>
        title.toLowerCase().includes('pbmc')
      );
      expect(hasPBMC).toBe(true);
    });

    await test.step('Clear search and verify all datasets shown', async () => {
      await datasetsPage.searchDataset('');

      const allDatasets = await datasetsPage.getVisibleDatasets();
      expect(allDatasets.length).toBeGreaterThan(0);
    });
  });

  /**
   * Additional test: Verify plugin template selection
   */
  test('should display and select plugin templates', async ({ page }) => {
    const projectsPage = new ProjectsPage(page);

    await test.step('Navigate to Projects page', async () => {
      await goToProjects(page);
      await projectsPage.verifyPageLoaded();
    });

    await test.step('Open plugin selection modal', async () => {
      await projectsPage.clickNewWorkflow();

      await expect
        .poll(async () => {
          const plugins = await projectsPage.getAvailablePlugins();
          console.log('Available plugins:', plugins);
          return plugins.some((p) => p.name && p.name.includes('TENET'));
        }, {
          message: 'Waiting for TENET plugin to appear in plugin list',
          timeout: 20000,
        })
        .toBe(true);
    });

    await test.step('Select default template', async () => {
      await projectsPage.selectDefaultTemplate();
      await expect(page).toHaveURL(/.*\/workflow.*/);
    });
  });
});
