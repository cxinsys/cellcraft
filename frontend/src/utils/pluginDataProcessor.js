/**
 * Plugin Data Processing Utilities
 *
 * Pure functions for enriching and transforming plugin data
 * with user associations and build status information.
 */

/**
 * Enrich a single plugin with user association and build status
 *
 * @param {Object} plugin - Raw plugin data from API
 * @param {number} plugin.id - Plugin ID
 * @param {string} plugin.name - Plugin name
 * @param {Array<Object>} plugin.users - Array of associated users
 * @param {Object} plugin.latest_build - Latest build information
 * @param {string} plugin.version - Plugin version
 * @param {string} currentUsername - Current logged-in username
 * @returns {Object} Enriched plugin object with computed properties
 *
 * @example
 * const plugin = {
 *   id: 1,
 *   name: 'TENET',
 *   users: [{ username: 'user1' }],
 *   latest_build: { task_id: '123', status: 'RUNNING' },
 *   version: '1.0.0'
 * };
 * const enriched = enrichPluginData(plugin, 'user1');
 * // enriched.checked === true
 * // enriched.building === true
 */
export function enrichPluginData(plugin, currentUsername) {
  // Check if current user is associated with this plugin
  const userIncluded = plugin.users && Array.isArray(plugin.users)
    ? plugin.users.some(user => user.username === currentUsername)
    : false;

  // Extract build information
  const buildInfo = plugin.latest_build || {};
  const isBuilding = buildInfo.status === 'RUNNING' || buildInfo.status === 'PENDING';

  return {
    ...plugin,
    checked: userIncluded,
    building: isBuilding,
    imageExists: false, // Will be updated by separate image check
    buildTaskId: buildInfo.task_id || null,
    buildStatus: buildInfo.status || null,
    current_version: plugin.version || plugin.current_version || 'latest',
  };
}

/**
 * Enrich an array of plugins with user association and build status
 *
 * @param {Array<Object>} plugins - Array of raw plugin data from API
 * @param {string} currentUsername - Current logged-in username
 * @returns {Array<Object>} Array of enriched plugin objects
 *
 * @example
 * const plugins = [
 *   { id: 1, name: 'TENET', users: [{ username: 'user1' }] },
 *   { id: 2, name: 'GENIE3', users: [] }
 * ];
 * const enriched = enrichPlugins(plugins, 'user1');
 * // enriched[0].checked === true
 * // enriched[1].checked === false
 */
export function enrichPlugins(plugins, currentUsername) {
  if (!Array.isArray(plugins)) {
    return [];
  }

  return plugins.map(plugin => enrichPluginData(plugin, currentUsername));
}

/**
 * Check if a plugin's build status indicates it is currently building
 *
 * @param {string} status - Build status string
 * @returns {boolean} True if plugin is currently building
 */
export function isBuildingStatus(status) {
  return status === 'RUNNING' || status === 'PENDING';
}

/**
 * Check if a plugin's build status indicates completion (success or failure)
 *
 * @param {string} status - Build status string
 * @returns {boolean} True if build is complete
 */
export function isTerminalBuildStatus(status) {
  return ['SUCCESS', 'FAILURE', 'REVOKED', 'ERROR'].includes(status);
}
