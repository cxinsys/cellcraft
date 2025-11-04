// frontend/tests/e2e/workflows/support/workflow-constants.js

/**
 * Test data configuration for TENET workflow tests
 */
export const testWorkflow = {
  name: 'TENET',
  expectedFile: 'pbmc_light_1000.h5ad',
  folder: 'data',
  // TENET 템플릿의 실제 InputFile 노드 이름
  inputNodeName: 'Input h5ad', // Node title for opening modal
  inputNodeTabName: 'input.h5ad', // Tab text displayed in lowercase with dots
};

/**
 * Generate unique filename with timestamp and random ID
 * @param {string} baseFileName - Base filename (e.g., 'test_data.h5ad')
 * @returns {string} Unique filename (e.g., 'test_data_1234567890_abc123.h5ad')
 */
export function generateUniqueFileName(baseFileName) {
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
 * Generate unique workflow title with timestamp and random ID
 * @param {string} baseTitle - Base title (default: 'Test Workflow')
 * @returns {string} Unique title (e.g., 'Test Workflow_1234567890_abc123')
 */
export function generateUniqueWorkflowTitle(baseTitle = 'Test Workflow') {
  const timestamp = Date.now();
  const randomId = Math.random().toString(36).substring(2, 8);
  return `${baseTitle}_${timestamp}_${randomId}`;
}
