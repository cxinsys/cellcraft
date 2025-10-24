/**
 * Allowed file extensions for upload
 * @constant {string[]}
 */
export const ALLOWED_FILE_EXTENSIONS = ['h5ad', 'csv', 'txt'];

/**
 * Validate email format
 * @param {string} email - Email address to validate
 * @returns {boolean} True if valid email format
 */
function validateEmail(email) {
  var re =
    /^(([^<>()[\]\\.,;:\s@"]+(\.[^<>()[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/;
  return re.test(String(email).toLowerCase());
}

/**
 * Validate file extension against allowed types
 * @param {string} fileName - File name to validate
 * @param {string[]} allowedExtensions - Array of allowed extensions (default: ALLOWED_FILE_EXTENSIONS)
 * @returns {{isValid: boolean, extension: string, message: string}} Validation result
 */
export function validateFileExtension(fileName, allowedExtensions = ALLOWED_FILE_EXTENSIONS) {
  if (!fileName || typeof fileName !== 'string') {
    return {
      isValid: false,
      extension: '',
      message: 'Invalid file name'
    };
  }

  const extension = fileName.split('.').pop().toLowerCase();
  const isValid = allowedExtensions.includes(extension);

  return {
    isValid,
    extension,
    message: isValid
      ? ''
      : `Please upload ${allowedExtensions.join(', ')} file`
  };
}

export { validateEmail };
