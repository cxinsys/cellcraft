import moment from "moment";

/**
 * Format bytes to human-readable format with customizable precision
 * @param {number} bytes - Byte size to format
 * @param {number} decimals - Number of decimal places (default: 2)
 * @returns {string} Formatted byte size (e.g., "1.5 MB", "2.34 GB")
 */
export function formatBytes(bytes, decimals = 2) {
  if (bytes === 0) return '0 Bytes';
  const k = 1024;
  const dm = decimals < 0 ? 0 : decimals;
  const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
  const i = Math.floor(Math.log(bytes) / Math.log(k));
  return parseFloat((bytes / Math.pow(k, i)).toFixed(dm)) + ' ' + sizes[i];
}

/**
 * Format file size in bytes to human-readable format
 * Wrapper for formatBytes with 2 decimal places
 * @param {number} bytes - File size in bytes
 * @returns {string} Formatted file size (e.g., "1.5 MB")
 */
export function formatFileSize(bytes) {
  return formatBytes(bytes, 2);
}

/**
 * Calculate time difference between two timestamps
 * @param {string|Date} startTime - Start time
 * @param {string|Date} endTime - End time
 * @returns {string} Formatted time difference (HH:MM:SS)
 */
export function getTimeDifference(startTime, endTime) {
  const start = new Date(startTime);
  const end = new Date(endTime);
  let diff = Math.abs(end - start); // Difference in milliseconds

  let hours = Math.floor(diff / 3600000);
  diff -= hours * 3600000;

  const minutes = Math.floor(diff / 60000);
  diff -= minutes * 60000;

  const seconds = Math.floor(diff / 1000);

  return `${hours.toString().padStart(2, "0")}:${minutes
    .toString()
    .padStart(2, "0")}:${seconds.toString().padStart(2, "0")}`;
}

/**
 * Calculate running time from start time to current time
 * @param {string|Date} startTime - Start time
 * @param {Date} currentTime - Current time
 * @returns {string} Formatted running time (HH:MM:SS)
 */
export function getRunningTime(startTime, currentTime) {
  const start = new Date(startTime);
  let diff = Math.abs(currentTime - start); // Difference in milliseconds

  let hours = Math.floor(diff / 3600000);
  diff -= hours * 3600000;

  const minutes = Math.floor(diff / 60000);
  diff -= minutes * 60000;

  const seconds = Math.floor(diff / 1000);

  return `${hours.toString().padStart(2, "0")}:${minutes
    .toString()
    .padStart(2, "0")}:${seconds.toString().padStart(2, "0")}`;
}

/**
 * Format date/time to standard format
 * @param {string|Date} dateTime - Date/time to format
 * @returns {string} Formatted date/time (YYYY.MM.DD-HH:mm) or "Not Yet Completed"
 */
export function formatDateTime(dateTime) {
  if (!dateTime || dateTime === null || dateTime === undefined) {
    return "Not Yet Completed";
  }
  const date = moment(dateTime).format("YYYY.MM.DD-HH:mm");
  if (date === "Invalid date") return "Not Yet Completed";
  return date;
}

/**
 * Format title with default value for null/empty
 * @param {string|null} title - Title to format
 * @returns {string} Formatted title or "Untitled"
 */
export function formatTitle(title) {
  if (title === null || title === "" || !title) return "Untitled";
  return title;
}

/**
 * Extract date from ISO string (cuts at 'T')
 * @param {string} isoString - ISO format date-time string (e.g., "2025-10-23T12:34:56")
 * @returns {string} Date portion only (e.g., "2025-10-23")
 */
export function cutDateFromISO(isoString) {
  if (!isoString || typeof isoString !== 'string') {
    return '';
  }
  return isoString.split("T")[0];
}

/**
 * Extract filename without extension
 * @param {string} fullName - Full filename with extension
 * @returns {string} Filename without extension
 */
export function extractFileName(fullName) {
  if (!fullName || typeof fullName !== 'string') {
    return '';
  }
  const parts = fullName.split(".");
  return parts.length > 1 ? parts.slice(0, -1).join(".") : fullName;
}

/**
 * Extract file extension from filename
 * @param {string} fullName - Full filename with extension
 * @returns {string} File extension (e.g., "txt", "h5ad")
 */
export function extractExtension(fullName) {
  if (!fullName || typeof fullName !== 'string') {
    return '';
  }
  const parts = fullName.split(".");
  return parts.length > 1 ? parts[parts.length - 1] : '';
}

/**
 * Format relative time from a past timestamp
 * @param {string|Date} timestamp - Past timestamp to format
 * @param {Date} currentTime - Current time (default: Date.now())
 * @returns {string} Relative time string (e.g., "Edited Recently", "Edited 5 hours ago")
 */
export function formatRelativeTime(timestamp, currentTime = null) {
  if (!timestamp) {
    return '';
  }

  const now = currentTime ? new Date(currentTime).getTime() : Date.now();
  const past = new Date(timestamp).getTime();

  if (isNaN(past)) {
    return '';
  }

  const time_diff = now - past;

  // Handle future timestamps (e.g., server clock ahead of client)
  if (time_diff < 0) {
    return "Edited Recently";
  }

  const hours = Math.floor(time_diff / (1000 * 60 * 60));

  if (hours === 0) {
    return "Edited Recently";
  } else if (hours === 1) {
    return "Edited 1 hour ago";
  } else {
    return `Edited ${hours} hours ago`;
  }
}
