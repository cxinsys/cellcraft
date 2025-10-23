/**
 * Filename generation utilities for task-related file exports
 * Centralizes file naming logic for consistency across the application
 */

/**
 * Format timestamp for filenames
 * Converts Date object to ISO string format without colons and dots
 * @param {Date} date - Date object to format (defaults to current date)
 * @returns {string} Formatted timestamp (e.g., "2025-10-23T12-34-56")
 */
export function formatTimestamp(date = new Date()) {
  return date.toISOString().replace(/[:.]/g, '-').slice(0, -5);
}

/**
 * Extract short identifier from task ID
 * @param {string} taskId - Full task ID
 * @param {number} length - Number of characters to extract (default: 8)
 * @returns {string} Short task ID (e.g., "abc12345")
 */
export function getTaskShortId(taskId, length = 8) {
  if (!taskId || typeof taskId !== 'string') {
    return '';
  }
  return taskId.substring(0, length);
}

/**
 * Extract base name from filename by removing extension
 * @param {string} filename - Filename with extension
 * @returns {string} Filename without extension
 */
export function getBaseName(filename) {
  if (!filename || typeof filename !== 'string') {
    return '';
  }
  const baseName = filename.split('.').slice(0, -1).join('.') || filename;
  return baseName;
}

/**
 * Generate filename for all logs JSON export
 * @param {string} taskId - Task ID
 * @param {Date} date - Date for timestamp (defaults to current date)
 * @returns {string} Formatted filename (e.g., "task_abc12345_logs_2025-10-23T12-34-56.json")
 */
export function generateAllLogsFilename(taskId, date = new Date()) {
  const timestamp = formatTimestamp(date);
  const taskShortId = getTaskShortId(taskId);
  return `task_${taskShortId}_logs_${timestamp}.json`;
}

/**
 * Generate filename for single log file export
 * @param {string} taskId - Task ID
 * @param {string} logFilename - Original log filename
 * @param {Date} date - Date for timestamp (defaults to current date)
 * @param {string} extension - File extension (default: 'txt')
 * @returns {string} Formatted filename (e.g., "task_abc12345_snakemake_2025-10-23T12-34-56.txt")
 */
export function generateLogFilename(taskId, logFilename, date = new Date(), extension = 'txt') {
  const timestamp = formatTimestamp(date);
  const taskShortId = getTaskShortId(taskId);
  const baseName = getBaseName(logFilename);
  return `task_${taskShortId}_${baseName}_${timestamp}.${extension}`;
}

/**
 * Generate filename for execution manifest download
 * @param {string} taskId - Task ID
 * @param {Date} date - Date for timestamp (defaults to current date)
 * @returns {string} Formatted filename (e.g., "execution_manifest_abc12345_2025-10-23T12-34-56.json")
 */
export function generateManifestFilename(taskId, date = new Date()) {
  const timestamp = formatTimestamp(date);
  const taskShortId = getTaskShortId(taskId);
  return `execution_manifest_${taskShortId}_${timestamp}.json`;
}
