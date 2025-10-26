import { IntervalManager } from '@/utils/intervalManager';

jest.useFakeTimers();

describe('IntervalManager', () => {
  let manager;

  beforeEach(() => {
    manager = new IntervalManager();
  });

  afterEach(() => {
    jest.clearAllTimers();
    manager.stopAll();
  });

  describe('start', () => {
    it('should start an interval with given key', () => {
      const callback = jest.fn();
      manager.start('test', callback, 1000);

      expect(callback).not.toHaveBeenCalled();

      jest.advanceTimersByTime(1000);
      expect(callback).toHaveBeenCalledTimes(1);

      jest.advanceTimersByTime(2000);
      expect(callback).toHaveBeenCalledTimes(3);
    });

    it('should return interval ID', () => {
      const callback = jest.fn();
      const intervalId = manager.start('test', callback, 1000);

      expect(typeof intervalId).toBe('number');
    });

    it('should replace existing interval with same key', () => {
      const callback1 = jest.fn();
      const callback2 = jest.fn();

      manager.start('test', callback1, 1000);
      manager.start('test', callback2, 1000);

      jest.advanceTimersByTime(1000);

      expect(callback1).not.toHaveBeenCalled();
      expect(callback2).toHaveBeenCalledTimes(1);
    });

    it('should handle multiple intervals with different keys', () => {
      const callback1 = jest.fn();
      const callback2 = jest.fn();

      manager.start('interval1', callback1, 1000);
      manager.start('interval2', callback2, 500);

      jest.advanceTimersByTime(1000);

      expect(callback1).toHaveBeenCalledTimes(1);
      expect(callback2).toHaveBeenCalledTimes(2); // Called at 500ms and 1000ms
    });
  });

  describe('stop', () => {
    it('should stop a specific interval by key', () => {
      const callback = jest.fn();
      manager.start('test', callback, 1000);

      jest.advanceTimersByTime(1000);
      expect(callback).toHaveBeenCalledTimes(1);

      manager.stop('test');
      jest.advanceTimersByTime(2000);

      // Should not be called again after stop
      expect(callback).toHaveBeenCalledTimes(1);
    });

    it('should return true when interval is found and stopped', () => {
      const callback = jest.fn();
      manager.start('test', callback, 1000);

      const result = manager.stop('test');

      expect(result).toBe(true);
    });

    it('should return false when interval key is not found', () => {
      const result = manager.stop('nonexistent');

      expect(result).toBe(false);
    });

    it('should not affect other intervals', () => {
      const callback1 = jest.fn();
      const callback2 = jest.fn();

      manager.start('interval1', callback1, 1000);
      manager.start('interval2', callback2, 1000);

      manager.stop('interval1');
      jest.advanceTimersByTime(1000);

      expect(callback1).not.toHaveBeenCalled();
      expect(callback2).toHaveBeenCalledTimes(1);
    });
  });

  describe('stopAll', () => {
    it('should stop all active intervals', () => {
      const callback1 = jest.fn();
      const callback2 = jest.fn();
      const callback3 = jest.fn();

      manager.start('interval1', callback1, 1000);
      manager.start('interval2', callback2, 1000);
      manager.start('interval3', callback3, 1000);

      manager.stopAll();
      jest.advanceTimersByTime(2000);

      expect(callback1).not.toHaveBeenCalled();
      expect(callback2).not.toHaveBeenCalled();
      expect(callback3).not.toHaveBeenCalled();
    });

    it('should clear the intervals map', () => {
      manager.start('interval1', jest.fn(), 1000);
      manager.start('interval2', jest.fn(), 1000);

      expect(manager.getActiveCount()).toBe(2);

      manager.stopAll();

      expect(manager.getActiveCount()).toBe(0);
    });
  });

  describe('isActive', () => {
    it('should return true for active interval', () => {
      manager.start('test', jest.fn(), 1000);

      expect(manager.isActive('test')).toBe(true);
    });

    it('should return false for inactive interval', () => {
      expect(manager.isActive('nonexistent')).toBe(false);
    });

    it('should return false after interval is stopped', () => {
      manager.start('test', jest.fn(), 1000);
      manager.stop('test');

      expect(manager.isActive('test')).toBe(false);
    });
  });

  describe('getActiveCount', () => {
    it('should return 0 for no active intervals', () => {
      expect(manager.getActiveCount()).toBe(0);
    });

    it('should return correct count of active intervals', () => {
      manager.start('interval1', jest.fn(), 1000);
      manager.start('interval2', jest.fn(), 1000);
      manager.start('interval3', jest.fn(), 1000);

      expect(manager.getActiveCount()).toBe(3);
    });

    it('should update count when intervals are stopped', () => {
      manager.start('interval1', jest.fn(), 1000);
      manager.start('interval2', jest.fn(), 1000);

      expect(manager.getActiveCount()).toBe(2);

      manager.stop('interval1');

      expect(manager.getActiveCount()).toBe(1);
    });
  });

  describe('getActiveKeys', () => {
    it('should return empty array for no active intervals', () => {
      expect(manager.getActiveKeys()).toEqual([]);
    });

    it('should return array of active interval keys', () => {
      manager.start('interval1', jest.fn(), 1000);
      manager.start('interval2', jest.fn(), 1000);

      const keys = manager.getActiveKeys();

      expect(keys).toHaveLength(2);
      expect(keys).toContain('interval1');
      expect(keys).toContain('interval2');
    });

    it('should not include stopped intervals', () => {
      manager.start('interval1', jest.fn(), 1000);
      manager.start('interval2', jest.fn(), 1000);
      manager.stop('interval1');

      const keys = manager.getActiveKeys();

      expect(keys).toEqual(['interval2']);
    });
  });
});
