import store from "@/store/index";

/**
 * Enhanced error response structure for better error handling
 */
class APIError extends Error {
  constructor(response, originalError) {
    const errorData = response?.data?.error || {};
    const message = errorData.message || response?.data?.detail || originalError.message || 'An error occurred';
    
    super(message);
    this.name = 'APIError';
    this.status = response?.status || 0;
    this.type = errorData.type || 'unknown';
    this.details = errorData.details || null;
    this.suggestedActions = errorData.suggested_actions || [];
    this.context = errorData.context || {};
    this.originalError = originalError;
  }

  /**
   * Get user-friendly error message
   */
  getUserMessage() {
    // Provide user-friendly messages for common error types
    switch (this.type) {
      case 'not_found':
        return this.message || 'The requested item was not found';
      case 'validation':
        return this.message || 'Invalid input provided';
      case 'plugin_error':
        return this.message || 'Plugin error occurred';
      case 'file_error':
        return this.message || 'File operation failed';
      case 'workflow_error':
        return this.message || 'Workflow error occurred';
      case 'task_error':
        return this.message || 'Task processing failed';
      case 'server_error':
        return this.message || 'Server error occurred';
      default:
        return this.message;
    }
  }

  /**
   * Get suggested actions for the user
   */
  getSuggestedActions() {
    if (this.suggestedActions && this.suggestedActions.length > 0) {
      return this.suggestedActions;
    }

    // Provide default suggestions based on status code
    switch (this.status) {
      case 400:
        return ['Check your input and try again', 'Ensure all required fields are provided'];
      case 401:
        return ['Please log in again', 'Check your authentication credentials'];
      case 403:
        return ['You may not have permission for this action', 'Contact administrator if needed'];
      case 404:
        return ['Check if the item exists', 'Refresh the page and try again'];
      case 429:
        return ['Please wait a moment and try again', 'Too many requests sent recently'];
      case 500:
        return ['Try again in a few moments', 'Contact support if the problem persists'];
      default:
        return ['Try again', 'Contact support if the problem continues'];
    }
  }
}

export function setInterceptors(instance) {
  // Add a request interceptor
  instance.interceptors.request.use(
    function (config) {
      // Add authorization header if token exists
      const token = store.state.token;
      if (token) {
        config.headers.Authorization = `Bearer ${token}`;
      }

      // Add request timestamp for debugging
      config.metadata = { startTime: Date.now() };
      
      return config;
    },
    function (error) {
      console.error('Request interceptor error:', error);
      return Promise.reject(new APIError(null, error));
    }
  );

  // Add a response interceptor
  instance.interceptors.response.use(
    function (response) {
      // Log response time for debugging
      if (response.config.metadata) {
        const duration = Date.now() - response.config.metadata.startTime;
        console.debug(`API call to ${response.config.url} completed in ${duration}ms`);
      }
      
      return response;
    },
    function (error) {
      console.error('API Error:', {
        url: error.config?.url,
        method: error.config?.method,
        status: error.response?.status,
        data: error.response?.data
      });

      // Handle network errors
      if (!error.response) {
        const networkError = new APIError(null, error);
        networkError.type = 'network_error';
        networkError.message = 'Network error - please check your internet connection';
        return Promise.reject(networkError);
      }

      // Handle authentication errors
      if (error.response.status === 401) {
        // Clear token and redirect to login if needed
        store.commit('clearToken');
        // You might want to redirect to login page here
        // router.push('/login');
      }

      // Create structured API error
      const apiError = new APIError(error.response, error);
      return Promise.reject(apiError);
    }
  );
  
  return instance;
}

/**
 * Global error handler utility
 */
export function handleAPIError(error, context = {}) {
  if (error instanceof APIError) {
    console.error('API Error Details:', {
      message: error.message,
      type: error.type,
      status: error.status,
      details: error.details,
      context: { ...error.context, ...context },
      suggestedActions: error.getSuggestedActions()
    });
    
    return {
      message: error.getUserMessage(),
      type: error.type,
      status: error.status,
      suggestedActions: error.getSuggestedActions(),
      details: error.details
    };
  } else {
    console.error('Unexpected error:', error);
    return {
      message: 'An unexpected error occurred',
      type: 'unknown',
      status: 0,
      suggestedActions: ['Try again', 'Contact support if the problem persists']
    };
  }
}

export { APIError };
