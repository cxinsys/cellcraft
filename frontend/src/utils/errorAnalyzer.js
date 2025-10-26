/**
 * Error Analyzer
 *
 * Utilities for analyzing errors, extracting user-friendly messages,
 * and determining retry strategies.
 */

/**
 * Extract a user-friendly error message from an error object
 *
 * @param {Error|Object} error - Error object from API or network
 * @returns {string} User-friendly error message
 *
 * @example
 * const error = new Error('Network timeout');
 * const message = getErrorMessage(error);
 * // message: "Network error. Please check your internet connection."
 */
export function getErrorMessage(error) {
  if (!error) {
    return 'An unexpected error occurred';
  }

  // API error response with detailed information
  if (error.response) {
    // Extract detail from response data
    if (error.response.data?.detail) {
      return typeof error.response.data.detail === 'string'
        ? error.response.data.detail
        : error.response.data.detail.message || 'Unknown API error';
    }

    // HTTP status code specific messages
    const status = error.response.status;

    if (status === 404) {
      return 'Service not found. Please check if the server is running.';
    }

    if (status === 503) {
      return 'Service temporarily unavailable. Please try again later.';
    }

    if (status >= 500) {
      return 'Server error occurred. Please try again later.';
    }

    if (status === 401 || status === 403) {
      return 'Authentication failed. Please log in again.';
    }

    if (status === 400) {
      return 'Invalid request. Please check your input.';
    }

    return `Request failed with status ${status}`;
  }

  // Network error (request sent but no response received)
  if (error.request) {
    return 'Network error. Please check your internet connection.';
  }

  // JavaScript error or other error types
  return error.message || 'An unexpected error occurred';
}

/**
 * Determine if an error is retryable based on error type and status code
 *
 * @param {Error|Object} error - Error object from API or network
 * @returns {boolean} True if error is retryable
 *
 * @example
 * const error = { response: { status: 503 } };
 * if (canRetryError(error)) {
 *   // Retry the operation
 * }
 */
export function canRetryError(error) {
  if (!error) {
    return false;
  }

  // API response errors
  if (error.response) {
    const status = error.response.status;

    // Client errors (4xx) are generally not retryable
    // Except for specific cases like timeout (408) and rate limit (429)
    if (status >= 400 && status < 500) {
      return status === 408 || status === 429;
    }

    // Server errors (5xx) and timeout are retryable
    return status >= 500 || status === 408;
  }

  // Network errors (no response received) are retryable
  return !!error.request;
}

/**
 * Analyze an error and return comprehensive error information
 *
 * @param {Error|Object} error - Error object from API or network
 * @returns {Object} Error analysis result
 * @returns {string} result.message - User-friendly error message
 * @returns {boolean} result.canRetry - Whether the error is retryable
 * @returns {Date} result.timestamp - When the error occurred
 * @returns {string} [result.type] - Error type (api, network, client)
 * @returns {number} [result.statusCode] - HTTP status code if available
 *
 * @example
 * const error = { response: { status: 503 } };
 * const analysis = analyzeError(error);
 * // {
 * //   message: "Service temporarily unavailable...",
 * //   canRetry: true,
 * //   timestamp: Date,
 * //   type: "api",
 * //   statusCode: 503
 * // }
 */
export function analyzeError(error) {
  const analysis = {
    message: getErrorMessage(error),
    canRetry: canRetryError(error),
    timestamp: new Date()
  };

  if (!error) {
    analysis.type = 'unknown';
    return analysis;
  }

  // Determine error type
  if (error.response) {
    analysis.type = 'api';
    analysis.statusCode = error.response.status;
  } else if (error.request) {
    analysis.type = 'network';
  } else {
    analysis.type = 'client';
  }

  return analysis;
}

/**
 * Check if error is a timeout error
 *
 * @param {Error|Object} error - Error object
 * @returns {boolean} True if error is a timeout
 */
export function isTimeoutError(error) {
  if (!error) return false;

  return (
    error.code === 'ECONNABORTED' ||
    error.code === 'ETIMEDOUT' ||
    (error.response && error.response.status === 408) ||
    (error.message && error.message.toLowerCase().includes('timeout'))
  );
}

/**
 * Check if error is a network error
 *
 * @param {Error|Object} error - Error object
 * @returns {boolean} True if error is a network error
 */
export function isNetworkError(error) {
  if (!error) return false;

  return !!(error.request && !error.response);
}

/**
 * Check if error is an authentication error
 *
 * @param {Error|Object} error - Error object
 * @returns {boolean} True if error is authentication related
 */
export function isAuthError(error) {
  if (!error) return false;

  if (error.response) {
    const status = error.response.status;
    return status === 401 || status === 403;
  }

  return false;
}

/**
 * Get retry delay based on attempt number (exponential backoff)
 *
 * @param {number} attempt - Current attempt number (0-based)
 * @param {number} baseDelay - Base delay in milliseconds (default: 1000)
 * @param {number} maxDelay - Maximum delay in milliseconds (default: 30000)
 * @returns {number} Delay in milliseconds
 *
 * @example
 * const delay = getRetryDelay(2); // 3rd attempt
 * // delay: 4000 (1000 * 2^2)
 */
export function getRetryDelay(attempt, baseDelay = 1000, maxDelay = 30000) {
  const delay = baseDelay * Math.pow(2, attempt);
  return Math.min(delay, maxDelay);
}
