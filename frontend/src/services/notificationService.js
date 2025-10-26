/**
 * Notification Service
 *
 * Centralized service for displaying user notifications.
 * Currently uses browser alerts, but designed to be easily replaceable
 * with a toast/notification library in the future.
 */

/**
 * Notification Service Class
 *
 * Provides a consistent interface for displaying notifications to users.
 * Supports different notification types (success, error, warning, info)
 * and includes logging for debugging.
 */
export class NotificationService {
  /**
   * Create a new NotificationService instance
   *
   * @param {Object} [options] - Configuration options
   * @param {boolean} [options.enableLogging=true] - Enable console logging
   * @param {Function} [options.displayFunction] - Custom display function (for testing)
   */
  constructor(options = {}) {
    this.enableLogging = options.enableLogging !== false;
    this.displayFunction = options.displayFunction || this.defaultDisplay.bind(this);
  }

  /**
   * Default display function using browser alert
   *
   * @param {string} message - Message to display
   * @param {string} title - Notification title
   * @param {string} type - Notification type
   * @private
   */
  defaultDisplay(message, title, type) {
    const formattedMessage = title ? `${title}: ${message}` : message;
    window.alert(formattedMessage);
  }

  /**
   * Log notification to console
   *
   * @param {string} type - Notification type
   * @param {string} title - Notification title
   * @param {string} message - Notification message
   * @private
   */
  log(type, title, message) {
    if (!this.enableLogging) return;

    const logMessage = `${type.toUpperCase()} - ${title}: ${message}`;

    switch (type) {
      case 'success':
        console.log(logMessage);
        break;
      case 'error':
        console.error(logMessage);
        break;
      case 'warning':
        console.warn(logMessage);
        break;
      case 'info':
        console.info(logMessage);
        break;
      default:
        console.log(logMessage);
    }
  }

  /**
   * Display a success notification
   *
   * @param {string} message - Success message
   * @param {string} [title='Success'] - Notification title
   *
   * @example
   * notificationService.success('Plugin built successfully!', 'Build Complete');
   */
  success(message, title = 'Success') {
    this.log('success', title, message);
    this.displayFunction(message, title, 'success');
  }

  /**
   * Display an error notification
   *
   * @param {string} message - Error message
   * @param {string} [title='Error'] - Notification title
   *
   * @example
   * notificationService.error('Failed to load plugins', 'Load Error');
   */
  error(message, title = 'Error') {
    this.log('error', title, message);
    this.displayFunction(message, title, 'error');
  }

  /**
   * Display a warning notification
   *
   * @param {string} message - Warning message
   * @param {string} [title='Warning'] - Notification title
   *
   * @example
   * notificationService.warning('Some plugins failed to build', 'Build Warning');
   */
  warning(message, title = 'Warning') {
    this.log('warning', title, message);
    this.displayFunction(message, title, 'warning');
  }

  /**
   * Display an info notification
   *
   * @param {string} message - Info message
   * @param {string} [title='Info'] - Notification title
   *
   * @example
   * notificationService.info('Connection restored', 'Network Status');
   */
  info(message, title = 'Info') {
    this.log('info', title, message);
    this.displayFunction(message, title, 'info');
  }

  /**
   * Set a custom display function (useful for integrating toast libraries)
   *
   * @param {Function} displayFn - Custom display function (message, title, type) => void
   *
   * @example
   * // Integrate with a toast library
   * notificationService.setDisplayFunction((message, title, type) => {
   *   Toast.show({
   *     type: type,
   *     title: title,
   *     message: message
   *   });
   * });
   */
  setDisplayFunction(displayFn) {
    if (typeof displayFn === 'function') {
      this.displayFunction = displayFn;
    }
  }

  /**
   * Enable or disable console logging
   *
   * @param {boolean} enabled - Whether to enable logging
   */
  setLogging(enabled) {
    this.enableLogging = enabled;
  }
}

/**
 * Create a default notification service instance
 * This can be imported and used throughout the application
 */
export const notificationService = new NotificationService();
