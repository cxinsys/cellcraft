// API functions for Plugin Architecture Separation features
import { instance } from './index';

/**
 * Get plugins filtered by category (algorithm or visualization)
 * @param {string} category - Plugin category to filter by
 * @returns {Promise} API response with filtered plugins
 */
function getPluginsByCategory(category) {
  return instance.get('/routes/plugin/list', {
    params: { category }
  });
}

/**
 * Get plugin recommendations based on workflow context
 * @param {Object} context - Workflow context for recommendations
 * @returns {Promise} API response with recommendations
 */
function getVisualizationRecommendations(context) {
  return instance.post('/routes/plugin/recommendations', context);
}

/**
 * Validate workflow with separated plugin architecture
 * @param {Object} workflowData - Workflow data to validate
 * @returns {Promise} API response with validation results
 */
function validateSeparatedWorkflow(workflowData) {
  return instance.post('/routes/workflow/validate/separated', workflowData);
}

/**
 * Save file mappings for visualization nodes
 * @param {Object} mappingData - File mapping configuration
 * @returns {Promise} API response
 */
function saveFileMappings(mappingData) {
  return instance.post('/routes/workflow/file-mappings', mappingData);
}

/**
 * Get file mappings for a workflow
 * @param {string} workflowId - Workflow ID
 * @returns {Promise} API response with file mappings
 */
function getFileMappings(workflowId) {
  return instance.get(`/routes/workflow/file-mappings/${workflowId}`);
}

/**
 * Get result files from algorithm execution
 * @param {string} taskId - Algorithm task ID
 * @returns {Promise} API response with result files
 */
function getAlgorithmResultFiles(taskId) {
  return instance.get(`/routes/workflow/result-files/${taskId}`);
}

/**
 * Migration-related API functions
 */

/**
 * Analyze workflow for migration compatibility
 * @param {Object} workflowData - Workflow to analyze
 * @returns {Promise} API response with migration analysis
 */
function analyzeWorkflowForMigration(workflowData) {
  return instance.post('/routes/workflow/migration/analyze', workflowData);
}

/**
 * Apply migration to convert legacy workflow to separated architecture
 * @param {Object} migrationData - Migration configuration and workflow data
 * @returns {Promise} API response with migrated workflow
 */
function applyWorkflowMigration(migrationData) {
  return instance.post('/routes/workflow/migration/apply', migrationData);
}

/**
 * Create backup of workflow before migration
 * @param {Object} workflowData - Workflow to backup
 * @returns {Promise} API response with backup ID
 */
function createWorkflowBackup(workflowData) {
  return instance.post('/routes/workflow/backup', workflowData);
}

/**
 * Get migration status for a workflow
 * @param {string} workflowId - Workflow ID
 * @returns {Promise} API response with migration status
 */
function getMigrationStatus(workflowId) {
  return instance.get(`/routes/workflow/migration/status/${workflowId}`);
}

/**
 * Plugin template and configuration functions
 */

/**
 * Get plugin template with input/output specifications
 * @param {string} pluginId - Plugin ID
 * @param {string} pluginType - Plugin type (algorithm or visualization)
 * @returns {Promise} API response with plugin template
 */
function getPluginTemplate(pluginId, pluginType) {
  return instance.get(`/routes/plugin/template/${pluginId}`, {
    params: { type: pluginType }
  });
}

/**
 * Get plugin configuration options
 * @param {string} pluginId - Plugin ID
 * @returns {Promise} API response with configuration options
 */
function getPluginConfiguration(pluginId) {
  return instance.get(`/routes/plugin/configuration/${pluginId}`);
}

/**
 * Save plugin configuration for a workflow node
 * @param {Object} configData - Plugin configuration data
 * @returns {Promise} API response
 */
function savePluginConfiguration(configData) {
  return instance.post('/routes/plugin/configuration/save', configData);
}

/**
 * Workflow execution functions for separated architecture
 */

/**
 * Execute algorithm plugin separately
 * @param {Object} executionData - Algorithm execution configuration
 * @returns {Promise} API response with task ID
 */
function executeAlgorithmPlugin(executionData) {
  return instance.post('/routes/workflow/execute/algorithm', executionData);
}

/**
 * Execute visualization plugin with mapped files
 * @param {Object} visualizationData - Visualization execution configuration
 * @returns {Promise} API response with visualization results
 */
function executeVisualizationPlugin(visualizationData) {
  return instance.post('/routes/workflow/execute/visualization', visualizationData);
}

/**
 * Get execution status for separated plugin execution
 * @param {string} taskId - Task ID
 * @param {string} pluginType - Plugin type (algorithm or visualization)
 * @returns {Promise} API response with execution status
 */
function getPluginExecutionStatus(taskId, pluginType) {
  return instance.get(`/routes/workflow/execution/status/${taskId}`, {
    params: { type: pluginType }
  });
}

/**
 * File management functions
 */

/**
 * Get available files for file mapping
 * @param {Object} filters - File filters (type, source, etc.)
 * @returns {Promise} API response with available files
 */
function getAvailableFilesForMapping(filters = {}) {
  return instance.get('/routes/files/available-for-mapping', {
    params: filters
  });
}

/**
 * Validate file compatibility with plugin input requirements
 * @param {Object} validationData - File and plugin compatibility data
 * @returns {Promise} API response with compatibility results
 */
function validateFileCompatibility(validationData) {
  return instance.post('/routes/files/validate-compatibility', validationData);
}

/**
 * Preview file content for mapping decisions
 * @param {string} fileId - File ID to preview
 * @param {number} lines - Number of lines to preview (default: 10)
 * @returns {Promise} API response with file preview
 */
function previewFileForMapping(fileId, lines = 10) {
  return instance.get(`/routes/files/preview/${fileId}`, {
    params: { lines }
  });
}

/**
 * Statistics and analytics functions
 */

/**
 * Get plugin usage statistics
 * @param {Object} filters - Usage statistics filters
 * @returns {Promise} API response with usage statistics
 */
function getPluginUsageStatistics(filters = {}) {
  return instance.get('/routes/plugin/statistics/usage', {
    params: filters
  });
}

/**
 * Get recommendation effectiveness metrics
 * @param {string} userId - User ID (optional)
 * @returns {Promise} API response with recommendation metrics
 */
function getRecommendationMetrics(userId = null) {
  return instance.get('/routes/plugin/recommendations/metrics', {
    params: userId ? { user_id: userId } : {}
  });
}

/**
 * Export functions for use in components
 */
export {
  // Plugin category functions
  getPluginsByCategory,
  getVisualizationRecommendations,
  
  // Workflow validation and file mapping
  validateSeparatedWorkflow,
  saveFileMappings,
  getFileMappings,
  getAlgorithmResultFiles,
  
  // Migration functions
  analyzeWorkflowForMigration,
  applyWorkflowMigration,
  createWorkflowBackup,
  getMigrationStatus,
  
  // Plugin configuration
  getPluginTemplate,
  getPluginConfiguration,
  savePluginConfiguration,
  
  // Execution functions
  executeAlgorithmPlugin,
  executeVisualizationPlugin,
  getPluginExecutionStatus,
  
  // File management
  getAvailableFilesForMapping,
  validateFileCompatibility,
  previewFileForMapping,
  
  // Analytics
  getPluginUsageStatistics,
  getRecommendationMetrics
};