/**
 * Create a Blob URL from data
 * @param {any} data - Data to create Blob from
 * @param {string} type - MIME type (default: 'application/octet-stream')
 * @returns {string} Blob URL
 */
export function createBlobURL(data, type = 'application/octet-stream') {
  const blob = new Blob([data], { type });
  return window.URL.createObjectURL(blob);
}

/**
 * Trigger browser download for given data
 * @param {any} data - Data to download (string, Blob, etc.)
 * @param {string} filename - Desired filename for download
 * @param {string} type - MIME type (optional)
 */
export function downloadFile(data, filename, type) {
  const blob = new Blob([data], type ? { type } : undefined);
  const url = window.URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;

  document.body.appendChild(link);
  link.click();

  // Clean up by removing the link and revoking the object URL
  document.body.removeChild(link);
  window.URL.revokeObjectURL(url);
}
