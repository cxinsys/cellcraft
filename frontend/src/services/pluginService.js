/**
 * Plugin Service
 *
 * Service layer for plugin-related API operations.
 * Encapsulates all plugin API calls and provides a clean interface
 * for components to interact with plugin data.
 */

import {
  getUser,
  getPlugins,
  checkPluginImage,
  associatePlugin as associatePluginAPI,
  dissociatePlugin as dissociatePluginAPI,
  getBuildTasks as getBuildTasksAPI
} from '@/api/index';
import { enrichPlugins } from '@/utils/pluginDataProcessor';

export class PluginService {
  /**
   * Get the current user's profile information
   *
   * @returns {Promise<Object>} User profile data
   * @throws {Error} If API call fails
   *
   * @example
   * const service = new PluginService();
   * const profile = await service.getUserProfile();
   * console.log(profile.username);
   */
  async getUserProfile() {
    const response = await getUser();
    return response.data;
  }

  /**
   * Get the list of all available plugins
   *
   * @returns {Promise<Array<Object>>} Array of plugin objects
   * @throws {Error} If API call fails
   *
   * @example
   * const service = new PluginService();
   * const plugins = await service.getPluginsList();
   */
  async getPluginsList() {
    const response = await getPlugins();
    return response.data.plugins;
  }

  /**
   * Check if a Docker image exists for a specific plugin
   *
   * @param {string} pluginName - Name of the plugin to check
   * @returns {Promise<boolean>} True if image exists, false otherwise
   * @throws {Error} If API call fails
   *
   * @example
   * const service = new PluginService();
   * const exists = await service.checkPluginImage('TENET');
   */
  async checkPluginImage(pluginName) {
    const response = await checkPluginImage(pluginName);
    return response.data.image_exists;
  }

  /**
   * Check Docker image existence for multiple plugins
   * Updates each plugin object with imageExists property
   *
   * @param {Array<Object>} plugins - Array of plugin objects
   * @returns {Promise<Array<Object>>} Array of plugins with imageExists property
   *
   * @example
   * const service = new PluginService();
   * const plugins = [{ name: 'TENET' }, { name: 'GENIE3' }];
   * const withImages = await service.checkMultiplePluginImages(plugins);
   * // withImages[0].imageExists === true/false
   */
  async checkMultiplePluginImages(plugins) {
    const promises = plugins.map(async (plugin) => {
      try {
        const imageExists = await this.checkPluginImage(plugin.name);
        return { ...plugin, imageExists };
      } catch (error) {
        console.error(`Error checking image for plugin ${plugin.name}:`, error);
        return { ...plugin, imageExists: false };
      }
    });

    return Promise.all(promises);
  }

  /**
   * Get complete plugin data with user associations and image status
   * This is a convenience method that combines multiple API calls
   *
   * @param {string} currentUsername - Current logged-in username
   * @returns {Promise<Object>} Object containing profile and enriched plugins
   * @returns {Object} result.profile - User profile data
   * @returns {Array<Object>} result.plugins - Array of enriched plugin objects
   *
   * @example
   * const service = new PluginService();
   * const { profile, plugins } = await service.getCompletePluginData('user1');
   * // plugins[0].checked === true/false (user association)
   * // plugins[0].building === true/false (build status)
   */
  async getCompletePluginData(currentUsername) {
    // Get user profile and plugins list in parallel
    const [profile, plugins] = await Promise.all([
      this.getUserProfile(),
      this.getPluginsList()
    ]);

    // Enrich plugins with user associations and build status
    const enrichedPlugins = enrichPlugins(plugins, currentUsername);

    return {
      profile,
      plugins: enrichedPlugins
    };
  }

  /**
   * Associate a plugin with the current user
   *
   * @param {number} pluginId - ID of the plugin to associate
   * @returns {Promise<Object>} API response data
   * @throws {Error} If API call fails
   *
   * @example
   * const service = new PluginService();
   * await service.associatePlugin(123);
   */
  async associatePlugin(pluginId) {
    const response = await associatePluginAPI(pluginId);
    return response.data;
  }

  /**
   * Dissociate a plugin from the current user
   *
   * @param {number} pluginId - ID of the plugin to dissociate
   * @returns {Promise<Object>} API response data
   * @throws {Error} If API call fails
   *
   * @example
   * const service = new PluginService();
   * await service.dissociatePlugin(123);
   */
  async dissociatePlugin(pluginId) {
    const response = await dissociatePluginAPI(pluginId);
    return response.data;
  }

  /**
   * Get build tasks list, optionally filtered by plugin name
   *
   * @param {string|null} pluginName - Optional plugin name to filter tasks
   * @returns {Promise<Array<Object>>} Array of build task objects
   * @throws {Error} If API call fails
   *
   * @example
   * const service = new PluginService();
   * const allTasks = await service.getBuildTasks();
   * const tenetTasks = await service.getBuildTasks('TENET');
   */
  async getBuildTasks(pluginName = null) {
    const response = await getBuildTasksAPI(pluginName);
    return response.data.tasks || [];
  }
}
