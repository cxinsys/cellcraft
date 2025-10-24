/**
 * Build Service
 *
 * Service for managing plugin build operations including
 * single and batch builds, log retrieval, and result processing.
 */

import {
  buildPluginDocker,
  getBuildLogs,
  checkPluginImage,
  cancelBuildTask as cancelBuildTaskAPI
} from '@/api/index';
import { analyzeError } from '@/utils/errorAnalyzer';

export class BuildService {
  /**
   * Create a new BuildService instance
   *
   * @param {Object} [apiClient] - API client for testing (optional)
   */
  constructor(apiClient = null) {
    this.apiClient = apiClient || {
      buildPluginDocker,
      getBuildLogs,
      checkPluginImage,
      cancelBuildTask: cancelBuildTaskAPI
    };
  }

  /**
   * Build a single plugin
   *
   * @param {string} pluginName - Name of the plugin to build
   * @param {boolean} [useGpu=false] - Whether to use GPU for build
   * @returns {Promise<Object>} Build result
   * @returns {boolean} result.success - Whether build was successful
   * @returns {string} result.plugin - Plugin name
   * @returns {string} [result.taskId] - Build task ID if successful
   * @returns {Error} [result.error] - Error object if failed
   *
   * @example
   * const service = new BuildService();
   * const result = await service.buildPlugin('TENET', false);
   * if (result.success) {
   *   console.log('Build started:', result.taskId);
   * }
   */
  async buildPlugin(pluginName, useGpu = false) {
    try {
      const response = await this.apiClient.buildPluginDocker(pluginName, useGpu);
      return {
        success: true,
        plugin: pluginName,
        taskId: response.data.task_id
      };
    } catch (error) {
      console.error(`Error building plugin ${pluginName}:`, error);
      return {
        success: false,
        plugin: pluginName,
        error
      };
    }
  }

  /**
   * Build multiple plugins in parallel
   *
   * @param {Array<Object>} plugins - Array of plugin objects
   * @param {string} plugins[].name - Plugin name
   * @param {boolean} [useGpu=false] - Whether to use GPU for builds
   * @param {Function} [onProgress] - Progress callback (currentIndex, total) => void
   * @returns {Promise<Array<Object>>} Array of build results
   *
   * @example
   * const service = new BuildService();
   * const plugins = [{ name: 'TENET' }, { name: 'GENIE3' }];
   * const results = await service.buildMultiplePlugins(plugins, false, (current, total) => {
   *   console.log(`Building ${current}/${total}`);
   * });
   */
  async buildMultiplePlugins(plugins, useGpu = false, onProgress = null) {
    if (!Array.isArray(plugins) || plugins.length === 0) {
      return [];
    }

    const buildPromises = plugins.map(async (plugin, index) => {
      try {
        const response = await this.apiClient.buildPluginDocker(plugin.name, useGpu);

        if (onProgress) {
          onProgress(index + 1, plugins.length);
        }

        // Wait a bit and check image existence
        await this.delay(1000);

        let imageExists = false;
        try {
          const checkResult = await this.apiClient.checkPluginImage(plugin.name);
          imageExists = checkResult.data.image_exists;
        } catch (error) {
          console.error(`Error checking image for plugin ${plugin.name}:`, error);
        }

        return {
          success: true,
          plugin: plugin.name,
          taskId: response.data.task_id,
          imageExists
        };
      } catch (error) {
        console.error(`Error building plugin ${plugin.name}:`, error);

        if (onProgress) {
          onProgress(index + 1, plugins.length);
        }

        return {
          success: false,
          plugin: plugin.name,
          error
        };
      }
    });

    return Promise.all(buildPromises);
  }

  /**
   * Get build logs for a plugin
   *
   * @param {string} pluginName - Name of the plugin
   * @returns {Promise<Object>} Build logs data
   * @returns {string} result.logs - Log content
   * @returns {string} [result.error] - Error message if retrieval failed
   *
   * @example
   * const service = new BuildService();
   * const logs = await service.getBuildLogs('TENET');
   * console.log(logs.logs);
   */
  async getBuildLogs(pluginName) {
    try {
      const response = await this.apiClient.getBuildLogs(pluginName);
      return {
        logs: response.data.log_content || 'No logs available',
        success: true
      };
    } catch (error) {
      console.error(`Error fetching build logs for plugin ${pluginName}:`, error);
      return {
        logs: null,
        error: error.message || 'Failed to load build logs',
        success: false
      };
    }
  }

  /**
   * Process build results and generate summary statistics
   *
   * @param {Array<Object>} results - Array of build results
   * @returns {Object} Build summary
   * @returns {number} summary.total - Total number of builds
   * @returns {number} summary.successful - Number of successful builds
   * @returns {number} summary.failed - Number of failed builds
   * @returns {Array<string>} summary.successfulPlugins - Names of successfully built plugins
   * @returns {Array<string>} summary.failedPlugins - Names of failed plugins
   *
   * @example
   * const summary = service.processBuildResults(results);
   * console.log(`${summary.successful}/${summary.total} builds succeeded`);
   */
  processBuildResults(results) {
    if (!Array.isArray(results)) {
      return {
        total: 0,
        successful: 0,
        failed: 0,
        successfulPlugins: [],
        failedPlugins: []
      };
    }

    const successful = results.filter(r => r.success);
    const failed = results.filter(r => !r.success);

    return {
      total: results.length,
      successful: successful.length,
      failed: failed.length,
      successfulPlugins: successful.map(r => r.plugin),
      failedPlugins: failed.map(r => r.plugin)
    };
  }

  /**
   * Get a summary message for build results
   *
   * @param {Object} summary - Build summary from processBuildResults
   * @returns {Object} Message object
   * @returns {string} message.text - Summary message text
   * @returns {string} message.type - Message type ('success', 'warning', 'error')
   *
   * @example
   * const summary = service.processBuildResults(results);
   * const message = service.getBuildSummaryMessage(summary);
   * // message.text: "All 5 plugin(s) were built successfully!"
   * // message.type: "success"
   */
  getBuildSummaryMessage(summary) {
    const { total, successful, failed } = summary;

    if (total === 0) {
      return {
        text: 'No plugins were built.',
        type: 'info'
      };
    }

    if (failed === 0) {
      return {
        text: `All ${successful} plugin(s) were built successfully!`,
        type: 'success'
      };
    }

    if (successful === 0) {
      return {
        text: `All ${failed} plugin(s) failed to build.`,
        type: 'error'
      };
    }

    return {
      text: `${successful} plugin(s) built successfully, ${failed} plugin(s) failed.`,
      type: 'warning'
    };
  }

  /**
   * Helper function to delay execution
   *
   * @param {number} ms - Milliseconds to delay
   * @returns {Promise<void>}
   * @private
   */
  delay(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
  }

  /**
   * Filter plugins that need to be built
   *
   * @param {Array<Object>} plugins - Array of plugin objects
   * @param {boolean} plugins[].building - Whether plugin is currently building
   * @param {boolean} plugins[].imageExists - Whether Docker image exists
   * @returns {Array<Object>} Filtered array of plugins that need building
   *
   * @example
   * const service = new BuildService();
   * const toBuild = service.filterPluginsToBuild(allPlugins);
   */
  filterPluginsToBuild(plugins) {
    if (!Array.isArray(plugins)) {
      return [];
    }

    return plugins.filter(plugin => !plugin.building && !plugin.imageExists);
  }

  /**
   * Cancel a build task
   *
   * @param {string} taskId - ID of the build task to cancel
   * @returns {Promise<Object>} Result object
   * @returns {boolean} result.success - Whether cancellation was successful
   * @returns {string} [result.message] - Success message if successful
   * @returns {string} [result.error] - Error message if failed
   *
   * @example
   * const service = new BuildService();
   * const result = await service.cancelBuildTask('task-123');
   * if (result.success) {
   *   console.log('Task cancelled');
   * }
   */
  async cancelBuildTask(taskId) {
    try {
      await this.apiClient.cancelBuildTask(taskId);
      return {
        success: true,
        message: 'Build task has been cancelled.'
      };
    } catch (error) {
      console.error('Failed to cancel build task:', error);
      const errorAnalysis = analyzeError(error);
      return {
        success: false,
        error: `Failed to cancel build task: ${errorAnalysis.message}`
      };
    }
  }
}
