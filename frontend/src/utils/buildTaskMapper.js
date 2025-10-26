/**
 * Build Task Data Transformation Utilities
 *
 * Pure functions for transforming build task data from API responses
 * into the format needed by the UI components.
 */

/**
 * Transform a single build task from API format to UI format
 *
 * @param {Object} task - Raw task data from API
 * @param {string} task.task_id - Unique task identifier
 * @param {string} task.plugin_name - Name of the plugin being built
 * @param {string} task.start_time - ISO timestamp when task started
 * @param {string} task.end_time - ISO timestamp when task ended
 * @param {string} task.state - Current state of the task
 * @param {string} task.error - Error message if task failed
 * @param {Object} task.info - Additional task information
 * @param {Function} durationCalculator - Function to calculate duration string
 * @returns {Object} Transformed task object for UI
 */
export function transformBuildTask(task, durationCalculator) {
  return {
    task_id: task.task_id,
    plugin_name: task.plugin_name,
    start_time: task.start_time,
    end_time: task.end_time,
    running_time: durationCalculator(task),
    status: task.state,
    error: task.error,
    info: task.info
  };
}

/**
 * Transform an array of build tasks from API format to UI format
 *
 * @param {Array<Object>} tasks - Array of raw task data from API
 * @param {Function} durationCalculator - Function to calculate duration string
 * @returns {Array<Object>} Array of transformed task objects for UI
 *
 * @example
 * const tasks = [
 *   { task_id: '123', plugin_name: 'TENET', state: 'SUCCESS', ... }
 * ];
 * const transformedTasks = transformBuildTasks(tasks, calculateDuration);
 */
export function transformBuildTasks(tasks, durationCalculator) {
  if (!Array.isArray(tasks)) {
    return [];
  }

  return tasks.map(task => transformBuildTask(task, durationCalculator));
}
