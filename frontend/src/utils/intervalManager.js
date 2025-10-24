/**
 * Interval Manager
 *
 * Centralized management for setInterval timers with named keys.
 * Provides easy start/stop/cleanup for multiple concurrent intervals.
 */

export class IntervalManager {
  /**
   * Create a new IntervalManager instance
   */
  constructor() {
    /**
     * Map of interval keys to their interval IDs
     * @type {Map<string, number>}
     */
    this.intervals = new Map();
  }

  /**
   * Start a new interval with a unique key
   * If an interval with the same key already exists, it will be stopped first
   *
   * @param {string} key - Unique identifier for this interval
   * @param {Function} callback - Function to execute on each interval
   * @param {number} delay - Delay in milliseconds between executions
   * @returns {number} The interval ID
   *
   * @example
   * const manager = new IntervalManager();
   * manager.start('statusCheck', () => console.log('checking...'), 1000);
   */
  start(key, callback, delay) {
    // Stop existing interval with same key if it exists
    if (this.intervals.has(key)) {
      this.stop(key);
    }

    // Start new interval
    const intervalId = setInterval(callback, delay);
    this.intervals.set(key, intervalId);

    return intervalId;
  }

  /**
   * Stop a specific interval by its key
   *
   * @param {string} key - The unique key of the interval to stop
   * @returns {boolean} True if interval was found and stopped, false otherwise
   *
   * @example
   * manager.stop('statusCheck');
   */
  stop(key) {
    const intervalId = this.intervals.get(key);

    if (intervalId !== undefined) {
      clearInterval(intervalId);
      this.intervals.delete(key);
      return true;
    }

    return false;
  }

  /**
   * Stop all managed intervals and clear the registry
   *
   * @example
   * // Clean up all intervals on component destroy
   * beforeDestroy() {
   *   this.intervalManager.stopAll();
   * }
   */
  stopAll() {
    this.intervals.forEach((intervalId) => {
      clearInterval(intervalId);
    });
    this.intervals.clear();
  }

  /**
   * Check if an interval with the given key is currently active
   *
   * @param {string} key - The unique key to check
   * @returns {boolean} True if interval is active, false otherwise
   *
   * @example
   * if (manager.isActive('statusCheck')) {
   *   console.log('Status check is running');
   * }
   */
  isActive(key) {
    return this.intervals.has(key);
  }

  /**
   * Get the number of currently active intervals
   *
   * @returns {number} Count of active intervals
   */
  getActiveCount() {
    return this.intervals.size;
  }

  /**
   * Get all active interval keys
   *
   * @returns {Array<string>} Array of active interval keys
   */
  getActiveKeys() {
    return Array.from(this.intervals.keys());
  }
}
