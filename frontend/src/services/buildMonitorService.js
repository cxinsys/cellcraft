/**
 * Build Monitor Service
 *
 * Service for monitoring plugin build status with automatic polling.
 * Handles build status checks and notifies when builds complete or fail.
 */

import { getBuildStatus } from '@/api/index';

/**
 * Build status constants
 */
export const BuildStatus = {
  PENDING: 'PENDING',
  RUNNING: 'RUNNING',
  SUCCESS: 'SUCCESS',
  FAILURE: 'FAILURE',
  REVOKED: 'REVOKED',
  ERROR: 'ERROR'
};

export class BuildMonitor {
  /**
   * Create a new BuildMonitor instance
   *
   * @param {IntervalManager} intervalManager - Interval manager for scheduling checks
   * @param {Object} apiClient - API client object (for testing)
   * @param {Function} apiClient.getBuildStatus - Function to get build status
   */
  constructor(intervalManager, apiClient = { getBuildStatus }) {
    this.intervalManager = intervalManager;
    this.apiClient = apiClient;

    /**
     * Map of task IDs to monitoring metadata
     * @type {Map<string, Object>}
     */
    this.monitors = new Map();
  }

  /**
   * Check the current build status for a task
   *
   * @param {string} taskId - The task ID to check
   * @returns {Promise<Object>} Status check result
   * @returns {string} result.status - Current build status
   * @returns {boolean} result.success - Whether the API call succeeded
   * @returns {Error} [result.error] - Error object if call failed
   *
   * @example
   * const monitor = new BuildMonitor(intervalManager);
   * const result = await monitor.checkBuildStatus('task-123');
   * console.log(result.status); // 'RUNNING', 'SUCCESS', etc.
   */
  async checkBuildStatus(taskId) {
    try {
      const result = await this.apiClient.getBuildStatus(taskId);
      return {
        status: result.data.state,
        success: true
      };
    } catch (error) {
      console.error(`Error checking build status for task ${taskId}:`, error);
      return {
        status: BuildStatus.ERROR,
        success: false,
        error
      };
    }
  }

  /**
   * Check if a status indicates build completion (terminal state)
   *
   * @param {string} status - Build status to check
   * @returns {boolean} True if status is terminal (SUCCESS, FAILURE, REVOKED, ERROR)
   *
   * @example
   * monitor.isTerminalStatus('SUCCESS') // true
   * monitor.isTerminalStatus('RUNNING') // false
   */
  isTerminalStatus(status) {
    return [
      BuildStatus.SUCCESS,
      BuildStatus.FAILURE,
      BuildStatus.REVOKED,
      BuildStatus.ERROR
    ].includes(status);
  }

  /**
   * Start monitoring a build task
   * Polls the build status at regular intervals and calls the callback with updates
   * Automatically stops when a terminal status is reached
   *
   * @param {string} pluginName - Name of the plugin being built
   * @param {string} taskId - Task ID to monitor
   * @param {Function} onStatusUpdate - Callback for status updates
   * @param {Object} onStatusUpdate.update - Status update object
   * @param {string} onStatusUpdate.update.pluginName - Plugin name
   * @param {string} onStatusUpdate.update.taskId - Task ID
   * @param {string} onStatusUpdate.update.status - Current status
   * @param {Error} [onStatusUpdate.update.error] - Error if status check failed
   * @param {Object} [options] - Monitoring options
   * @param {number} [options.interval=2000] - Polling interval in milliseconds
   * @returns {string} Monitor key for this task
   *
   * @example
   * const monitor = new BuildMonitor(intervalManager);
   * monitor.startMonitoring('TENET', 'task-123', (update) => {
   *   console.log(`${update.pluginName}: ${update.status}`);
   *   if (update.status === BuildStatus.SUCCESS) {
   *     console.log('Build completed!');
   *   }
   * });
   */
  startMonitoring(pluginName, taskId, onStatusUpdate, options = {}) {
    const { interval = 2000 } = options;
    const monitorKey = `build-${taskId}`;

    // Define the status check function
    const checkStatus = async () => {
      const result = await this.checkBuildStatus(taskId);

      // Call the update callback
      onStatusUpdate({
        pluginName,
        taskId,
        status: result.status,
        error: result.error
      });

      // Stop monitoring if terminal status reached
      if (this.isTerminalStatus(result.status)) {
        this.stopMonitoring(taskId);
      }
    };

    // Start the interval
    this.intervalManager.start(monitorKey, checkStatus, interval);

    // Store monitoring metadata
    this.monitors.set(taskId, {
      pluginName,
      startTime: Date.now(),
      monitorKey
    });

    // Perform initial check immediately
    checkStatus();

    return monitorKey;
  }

  /**
   * Stop monitoring a specific build task
   *
   * @param {string} taskId - Task ID to stop monitoring
   * @returns {boolean} True if monitoring was stopped, false if not found
   *
   * @example
   * monitor.stopMonitoring('task-123');
   */
  stopMonitoring(taskId) {
    const metadata = this.monitors.get(taskId);

    if (metadata) {
      this.intervalManager.stop(metadata.monitorKey);
      this.monitors.delete(taskId);
      return true;
    }

    return false;
  }

  /**
   * Stop all active build monitors
   *
   * @example
   * // Clean up on component destroy
   * beforeDestroy() {
   *   this.buildMonitor.stopAll();
   * }
   */
  stopAll() {
    this.monitors.forEach((metadata, taskId) => {
      this.stopMonitoring(taskId);
    });
  }

  /**
   * Get information about all active monitors
   *
   * @returns {Array<Object>} Array of active monitor info
   */
  getActiveMonitors() {
    return Array.from(this.monitors.entries()).map(([taskId, metadata]) => ({
      taskId,
      ...metadata
    }));
  }

  /**
   * Check if a task is currently being monitored
   *
   * @param {string} taskId - Task ID to check
   * @returns {boolean} True if task is being monitored
   */
  isMonitoring(taskId) {
    return this.monitors.has(taskId);
  }
}
