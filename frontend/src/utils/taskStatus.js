/**
 * Task status constants
 */
export const TASK_STATUS = {
  SUCCESS: 'SUCCESS',
  FAILURE: 'FAILURE',
  REVOKED: 'REVOKED',
  RETRY: 'RETRY',
  RUNNING: 'RUNNING',
  PENDING: 'PENDING',
  INSTALLING: 'INSTALLING'
};

/**
 * Get CSS class name for task status indicator
 * @param {string} status - Task status
 * @returns {string} CSS class name
 */
export function getStatusClass(status) {
  if (status === TASK_STATUS.SUCCESS) return "status-success";
  if (status === TASK_STATUS.FAILURE ||
      status === TASK_STATUS.REVOKED) return "status-failure";
  if (status === TASK_STATUS.RETRY) return "status-retry";
  if (status === TASK_STATUS.RUNNING ||
      status === TASK_STATUS.PENDING ||
      status === TASK_STATUS.INSTALLING) return "status-running";
  return "";
}

/**
 * Check if task is completed (finished execution)
 * @param {string} status - Task status
 * @returns {boolean} True if task is completed
 */
export function isTaskCompleted(status) {
  return [
    TASK_STATUS.SUCCESS,
    TASK_STATUS.FAILURE,
    TASK_STATUS.REVOKED
  ].includes(status);
}

/**
 * Check if task is currently running or pending
 * @param {string} status - Task status
 * @returns {boolean} True if task is running
 */
export function isTaskRunning(status) {
  return [
    TASK_STATUS.RUNNING,
    TASK_STATUS.PENDING,
    TASK_STATUS.INSTALLING,
    TASK_STATUS.RETRY
  ].includes(status);
}
