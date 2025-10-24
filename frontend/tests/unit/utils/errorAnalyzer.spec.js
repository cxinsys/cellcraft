import {
  getErrorMessage,
  canRetryError,
  analyzeError,
  isTimeoutError,
  isNetworkError,
  isAuthError,
  getRetryDelay
} from '@/utils/errorAnalyzer';

describe('errorAnalyzer', () => {
  describe('getErrorMessage', () => {
    it('should return default message for null/undefined error', () => {
      expect(getErrorMessage(null)).toBe('An unexpected error occurred');
      expect(getErrorMessage(undefined)).toBe('An unexpected error occurred');
    });

    it('should extract detail string from API response', () => {
      const error = {
        response: {
          data: {
            detail: 'Custom API error message'
          }
        }
      };

      expect(getErrorMessage(error)).toBe('Custom API error message');
    });

    it('should extract detail.message from API response object', () => {
      const error = {
        response: {
          data: {
            detail: {
              message: 'Nested error message'
            }
          }
        }
      };

      expect(getErrorMessage(error)).toBe('Nested error message');
    });

    it('should return appropriate message for 404 status', () => {
      const error = {
        response: {
          status: 404
        }
      };

      expect(getErrorMessage(error)).toBe('Service not found. Please check if the server is running.');
    });

    it('should return appropriate message for 503 status', () => {
      const error = {
        response: {
          status: 503
        }
      };

      expect(getErrorMessage(error)).toBe('Service temporarily unavailable. Please try again later.');
    });

    it('should return appropriate message for 500+ status', () => {
      const error = {
        response: {
          status: 500
        }
      };

      expect(getErrorMessage(error)).toBe('Server error occurred. Please try again later.');
    });

    it('should return appropriate message for 401 status', () => {
      const error = {
        response: {
          status: 401
        }
      };

      expect(getErrorMessage(error)).toBe('Authentication failed. Please log in again.');
    });

    it('should return appropriate message for 403 status', () => {
      const error = {
        response: {
          status: 403
        }
      };

      expect(getErrorMessage(error)).toBe('Authentication failed. Please log in again.');
    });

    it('should return appropriate message for 400 status', () => {
      const error = {
        response: {
          status: 400
        }
      };

      expect(getErrorMessage(error)).toBe('Invalid request. Please check your input.');
    });

    it('should return generic message for other status codes', () => {
      const error = {
        response: {
          status: 418
        }
      };

      expect(getErrorMessage(error)).toBe('Request failed with status 418');
    });

    it('should return network error message for request without response', () => {
      const error = {
        request: {}
      };

      expect(getErrorMessage(error)).toBe('Network error. Please check your internet connection.');
    });

    it('should return error message from Error object', () => {
      const error = new Error('Custom JavaScript error');

      expect(getErrorMessage(error)).toBe('Custom JavaScript error');
    });
  });

  describe('canRetryError', () => {
    it('should return false for null/undefined error', () => {
      expect(canRetryError(null)).toBe(false);
      expect(canRetryError(undefined)).toBe(false);
    });

    it('should return false for 4xx client errors (except 408 and 429)', () => {
      expect(canRetryError({ response: { status: 400 } })).toBe(false);
      expect(canRetryError({ response: { status: 401 } })).toBe(false);
      expect(canRetryError({ response: { status: 404 } })).toBe(false);
    });

    it('should return true for 408 timeout', () => {
      expect(canRetryError({ response: { status: 408 } })).toBe(true);
    });

    it('should return true for 429 rate limit', () => {
      expect(canRetryError({ response: { status: 429 } })).toBe(true);
    });

    it('should return true for 5xx server errors', () => {
      expect(canRetryError({ response: { status: 500 } })).toBe(true);
      expect(canRetryError({ response: { status: 503 } })).toBe(true);
    });

    it('should return true for network errors', () => {
      expect(canRetryError({ request: {} })).toBe(true);
    });
  });

  describe('analyzeError', () => {
    it('should analyze API error with status code', () => {
      const error = {
        response: {
          status: 500,
          data: { detail: 'Server error' }
        }
      };

      const analysis = analyzeError(error);

      expect(analysis.message).toBe('Server error');
      expect(analysis.canRetry).toBe(true);
      expect(analysis.type).toBe('api');
      expect(analysis.statusCode).toBe(500);
      expect(analysis.timestamp).toBeInstanceOf(Date);
    });

    it('should analyze network error', () => {
      const error = {
        request: {}
      };

      const analysis = analyzeError(error);

      expect(analysis.message).toBe('Network error. Please check your internet connection.');
      expect(analysis.canRetry).toBe(true);
      expect(analysis.type).toBe('network');
      expect(analysis.timestamp).toBeInstanceOf(Date);
    });

    it('should analyze client-side error', () => {
      const error = new Error('Client error');

      const analysis = analyzeError(error);

      expect(analysis.message).toBe('Client error');
      expect(analysis.canRetry).toBe(false);
      expect(analysis.type).toBe('client');
      expect(analysis.timestamp).toBeInstanceOf(Date);
    });

    it('should analyze null error', () => {
      const analysis = analyzeError(null);

      expect(analysis.message).toBe('An unexpected error occurred');
      expect(analysis.canRetry).toBe(false);
      expect(analysis.type).toBe('unknown');
      expect(analysis.timestamp).toBeInstanceOf(Date);
    });
  });

  describe('isTimeoutError', () => {
    it('should return false for null/undefined error', () => {
      expect(isTimeoutError(null)).toBe(false);
      expect(isTimeoutError(undefined)).toBe(false);
    });

    it('should detect ECONNABORTED code', () => {
      expect(isTimeoutError({ code: 'ECONNABORTED' })).toBe(true);
    });

    it('should detect ETIMEDOUT code', () => {
      expect(isTimeoutError({ code: 'ETIMEDOUT' })).toBe(true);
    });

    it('should detect 408 status code', () => {
      expect(isTimeoutError({ response: { status: 408 } })).toBe(true);
    });

    it('should detect timeout in error message', () => {
      expect(isTimeoutError({ message: 'Request timeout' })).toBe(true);
      expect(isTimeoutError({ message: 'Connection TIMEOUT' })).toBe(true);
    });

    it('should return false for non-timeout errors', () => {
      expect(isTimeoutError({ code: 'ENOTFOUND' })).toBe(false);
      expect(isTimeoutError({ response: { status: 500 } })).toBe(false);
    });
  });

  describe('isNetworkError', () => {
    it('should return false for null/undefined error', () => {
      expect(isNetworkError(null)).toBe(false);
      expect(isNetworkError(undefined)).toBe(false);
    });

    it('should detect network error (request without response)', () => {
      expect(isNetworkError({ request: {} })).toBe(true);
    });

    it('should return false when response exists', () => {
      expect(isNetworkError({ request: {}, response: {} })).toBe(false);
    });

    it('should return false for client errors', () => {
      expect(isNetworkError(new Error('Client error'))).toBe(false);
    });
  });

  describe('isAuthError', () => {
    it('should return false for null/undefined error', () => {
      expect(isAuthError(null)).toBe(false);
      expect(isAuthError(undefined)).toBe(false);
    });

    it('should detect 401 Unauthorized', () => {
      expect(isAuthError({ response: { status: 401 } })).toBe(true);
    });

    it('should detect 403 Forbidden', () => {
      expect(isAuthError({ response: { status: 403 } })).toBe(true);
    });

    it('should return false for other status codes', () => {
      expect(isAuthError({ response: { status: 404 } })).toBe(false);
      expect(isAuthError({ response: { status: 500 } })).toBe(false);
    });

    it('should return false for network errors', () => {
      expect(isAuthError({ request: {} })).toBe(false);
    });
  });

  describe('getRetryDelay', () => {
    it('should calculate exponential backoff', () => {
      expect(getRetryDelay(0)).toBe(1000);  // 1000 * 2^0 = 1000
      expect(getRetryDelay(1)).toBe(2000);  // 1000 * 2^1 = 2000
      expect(getRetryDelay(2)).toBe(4000);  // 1000 * 2^2 = 4000
      expect(getRetryDelay(3)).toBe(8000);  // 1000 * 2^3 = 8000
    });

    it('should respect max delay', () => {
      expect(getRetryDelay(10)).toBe(30000);  // Would be 1024000, capped at 30000
      expect(getRetryDelay(20)).toBe(30000);  // Would be huge, capped at 30000
    });

    it('should use custom base delay', () => {
      expect(getRetryDelay(0, 500)).toBe(500);   // 500 * 2^0 = 500
      expect(getRetryDelay(1, 500)).toBe(1000);  // 500 * 2^1 = 1000
      expect(getRetryDelay(2, 500)).toBe(2000);  // 500 * 2^2 = 2000
    });

    it('should use custom max delay', () => {
      expect(getRetryDelay(5, 1000, 5000)).toBe(5000);  // Would be 32000, capped at 5000
    });
  });
});
