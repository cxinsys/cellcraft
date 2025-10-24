/**
 * Dialog Service
 *
 * Abstraction layer for browser dialog APIs (confirm, alert).
 * Provides a testable interface for user confirmations and notifications.
 */

export class DialogService {
  /**
   * Show a confirmation dialog to the user
   *
   * @param {string} message - Message to display in the dialog
   * @returns {boolean} True if user confirmed, false if cancelled
   *
   * @example
   * const dialogService = new DialogService();
   * const confirmed = dialogService.confirm('Are you sure?');
   * if (confirmed) {
   *   // User clicked OK
   * }
   */
  confirm(message) {
    return window.confirm(message);
  }

  /**
   * Show an alert dialog to the user
   *
   * @param {string} message - Message to display in the alert
   *
   * @example
   * const dialogService = new DialogService();
   * dialogService.alert('Operation completed successfully');
   */
  alert(message) {
    window.alert(message);
  }

  /**
   * Show a formatted alert with title and message
   *
   * @param {string} message - Main message content
   * @param {string} [title='Alert'] - Alert title
   *
   * @example
   * dialogService.alertWithTitle('File saved successfully', 'Success');
   */
  alertWithTitle(message, title = 'Alert') {
    const formattedMessage = `${title}:\n${message}`;
    window.alert(formattedMessage);
  }
}

/**
 * Event-driven Dialog Service
 *
 * Alternative implementation that allows for event-based testing
 * and custom dialog implementations (e.g., Vue modals instead of browser dialogs)
 */
export class EventDialogService {
  constructor() {
    /**
     * Registered event listeners
     * @type {Object}
     */
    this.listeners = {};
  }

  /**
   * Show a confirmation dialog
   * If a custom handler is registered via onConfirm(), it will be used
   * Otherwise falls back to window.confirm()
   *
   * @param {string} message - Message to display
   * @returns {Promise<boolean>} Resolves to true if confirmed, false otherwise
   *
   * @example
   * const service = new EventDialogService();
   * const confirmed = await service.confirm('Delete this item?');
   */
  async confirm(message) {
    return new Promise((resolve) => {
      if (this.listeners.onConfirm) {
        // Use custom handler if registered
        this.listeners.onConfirm(message, resolve);
      } else {
        // Fall back to browser confirm
        resolve(window.confirm(message));
      }
    });
  }

  /**
   * Register a custom confirmation handler
   *
   * @param {Function} callback - Handler function (message, resolve) => void
   *
   * @example
   * service.onConfirm((message, resolve) => {
   *   // Show custom Vue modal
   *   this.showModal({
   *     message,
   *     onConfirm: () => resolve(true),
   *     onCancel: () => resolve(false)
   *   });
   * });
   */
  onConfirm(callback) {
    this.listeners.onConfirm = callback;
  }

  /**
   * Show an alert dialog
   * If a custom handler is registered via onAlert(), it will be used
   * Otherwise falls back to window.alert()
   *
   * @param {string} message - Message to display
   * @returns {Promise<void>} Resolves when alert is dismissed
   */
  async alert(message) {
    return new Promise((resolve) => {
      if (this.listeners.onAlert) {
        // Use custom handler if registered
        this.listeners.onAlert(message, resolve);
      } else {
        // Fall back to browser alert
        window.alert(message);
        resolve();
      }
    });
  }

  /**
   * Register a custom alert handler
   *
   * @param {Function} callback - Handler function (message, resolve) => void
   *
   * @example
   * service.onAlert((message, resolve) => {
   *   // Show custom Vue notification
   *   this.$notify({
   *     message,
   *     onClose: () => resolve()
   *   });
   * });
   */
  onAlert(callback) {
    this.listeners.onAlert = callback;
  }

  /**
   * Clear all registered listeners
   */
  clearListeners() {
    this.listeners = {};
  }
}
