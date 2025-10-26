import { BuildMonitor, BuildStatus } from '@/services/buildMonitorService';
import { IntervalManager } from '@/utils/intervalManager';

jest.useFakeTimers();

describe('BuildMonitor', () => {
  let intervalManager;
  let buildMonitor;
  let mockApiClient;

  beforeEach(() => {
    intervalManager = new IntervalManager();
    mockApiClient = {
      getBuildStatus: jest.fn()
    };
    buildMonitor = new BuildMonitor(intervalManager, mockApiClient);
  });

  afterEach(() => {
    jest.clearAllTimers();
    buildMonitor.stopAll();
  });

  describe('BuildStatus constants', () => {
    it('should have all required status constants', () => {
      expect(BuildStatus.PENDING).toBe('PENDING');
      expect(BuildStatus.RUNNING).toBe('RUNNING');
      expect(BuildStatus.SUCCESS).toBe('SUCCESS');
      expect(BuildStatus.FAILURE).toBe('FAILURE');
      expect(BuildStatus.REVOKED).toBe('REVOKED');
      expect(BuildStatus.ERROR).toBe('ERROR');
    });
  });

  describe('checkBuildStatus', () => {
    it('should return status on successful API call', async () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      const result = await buildMonitor.checkBuildStatus('task-123');

      expect(result).toEqual({
        status: BuildStatus.RUNNING,
        success: true
      });
      expect(mockApiClient.getBuildStatus).toHaveBeenCalledWith('task-123');
    });

    it('should return error status on API failure', async () => {
      const error = new Error('API error');
      mockApiClient.getBuildStatus.mockRejectedValue(error);

      const consoleSpy = jest.spyOn(console, 'error').mockImplementation(() => {});

      const result = await buildMonitor.checkBuildStatus('task-456');

      expect(result).toEqual({
        status: BuildStatus.ERROR,
        success: false,
        error
      });

      consoleSpy.mockRestore();
    });
  });

  describe('isTerminalStatus', () => {
    it.each([
      [BuildStatus.SUCCESS, true],
      [BuildStatus.FAILURE, true],
      [BuildStatus.REVOKED, true],
      [BuildStatus.ERROR, true],
      [BuildStatus.RUNNING, false],
      [BuildStatus.PENDING, false]
    ])('should return %s for status %s', (status, expected) => {
      expect(buildMonitor.isTerminalStatus(status)).toBe(expected);
    });
  });

  describe('startMonitoring', () => {
    it('should start monitoring and call onStatusUpdate', async () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      const onStatusUpdate = jest.fn();
      buildMonitor.startMonitoring('TENET', 'task-123', onStatusUpdate);

      // Initial check should be called immediately
      await Promise.resolve();
      expect(onStatusUpdate).toHaveBeenCalledWith({
        pluginName: 'TENET',
        taskId: 'task-123',
        status: BuildStatus.RUNNING,
        error: undefined
      });

      // Wait for interval
      jest.advanceTimersByTime(2000);
      await Promise.resolve();

      expect(onStatusUpdate).toHaveBeenCalledTimes(2);
    });

    it('should stop monitoring on terminal status', async () => {
      mockApiClient.getBuildStatus
        .mockResolvedValueOnce({ data: { state: BuildStatus.RUNNING } })
        .mockResolvedValueOnce({ data: { state: BuildStatus.SUCCESS } });

      const onStatusUpdate = jest.fn();
      buildMonitor.startMonitoring('TENET', 'task-123', onStatusUpdate);

      await Promise.resolve();
      expect(onStatusUpdate).toHaveBeenCalledTimes(1);

      jest.advanceTimersByTime(2000);
      await Promise.resolve();

      // Should have been called twice (initial + one interval)
      expect(onStatusUpdate).toHaveBeenCalledTimes(2);
      expect(onStatusUpdate).toHaveBeenLastCalledWith({
        pluginName: 'TENET',
        taskId: 'task-123',
        status: BuildStatus.SUCCESS,
        error: undefined
      });

      // No more calls after terminal status
      jest.advanceTimersByTime(4000);
      await Promise.resolve();
      expect(onStatusUpdate).toHaveBeenCalledTimes(2);
    });

    it('should use custom interval if provided', async () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      const onStatusUpdate = jest.fn();
      buildMonitor.startMonitoring('TENET', 'task-123', onStatusUpdate, { interval: 5000 });

      await Promise.resolve();
      expect(onStatusUpdate).toHaveBeenCalledTimes(1);

      // Should not be called yet after 2000ms
      jest.advanceTimersByTime(2000);
      await Promise.resolve();
      expect(onStatusUpdate).toHaveBeenCalledTimes(1);

      // Should be called after 5000ms
      jest.advanceTimersByTime(3000);
      await Promise.resolve();
      expect(onStatusUpdate).toHaveBeenCalledTimes(2);
    });

    it('should return monitor key', () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      const onStatusUpdate = jest.fn();
      const key = buildMonitor.startMonitoring('TENET', 'task-123', onStatusUpdate);

      expect(key).toBe('build-task-123');
    });

    it('should handle API errors during monitoring', async () => {
      const error = new Error('Network error');
      mockApiClient.getBuildStatus.mockRejectedValue(error);

      const consoleSpy = jest.spyOn(console, 'error').mockImplementation(() => {});

      const onStatusUpdate = jest.fn();
      buildMonitor.startMonitoring('TENET', 'task-123', onStatusUpdate);

      await Promise.resolve();

      expect(onStatusUpdate).toHaveBeenCalledWith({
        pluginName: 'TENET',
        taskId: 'task-123',
        status: BuildStatus.ERROR,
        error
      });

      // Should stop monitoring after error
      jest.advanceTimersByTime(4000);
      await Promise.resolve();
      expect(onStatusUpdate).toHaveBeenCalledTimes(1);

      consoleSpy.mockRestore();
    });
  });

  describe('stopMonitoring', () => {
    it('should stop monitoring for specific task', async () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      const onStatusUpdate = jest.fn();
      buildMonitor.startMonitoring('TENET', 'task-123', onStatusUpdate);

      await Promise.resolve();
      expect(onStatusUpdate).toHaveBeenCalledTimes(1);

      buildMonitor.stopMonitoring('task-123');

      jest.advanceTimersByTime(4000);
      await Promise.resolve();
      expect(onStatusUpdate).toHaveBeenCalledTimes(1);
    });

    it('should return true when task is found and stopped', async () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      buildMonitor.startMonitoring('TENET', 'task-123', jest.fn());

      const result = buildMonitor.stopMonitoring('task-123');

      expect(result).toBe(true);
    });

    it('should return false when task is not found', () => {
      const result = buildMonitor.stopMonitoring('nonexistent');

      expect(result).toBe(false);
    });
  });

  describe('stopAll', () => {
    it('should stop all active monitors', async () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      const callback1 = jest.fn();
      const callback2 = jest.fn();

      buildMonitor.startMonitoring('TENET', 'task-1', callback1);
      buildMonitor.startMonitoring('GENIE3', 'task-2', callback2);

      await Promise.resolve();

      buildMonitor.stopAll();

      jest.advanceTimersByTime(4000);
      await Promise.resolve();

      expect(callback1).toHaveBeenCalledTimes(1); // Only initial call
      expect(callback2).toHaveBeenCalledTimes(1); // Only initial call
    });
  });

  describe('getActiveMonitors', () => {
    it('should return empty array when no monitors', () => {
      expect(buildMonitor.getActiveMonitors()).toEqual([]);
    });

    it('should return array of active monitor info', async () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      buildMonitor.startMonitoring('TENET', 'task-1', jest.fn());
      buildMonitor.startMonitoring('GENIE3', 'task-2', jest.fn());

      const monitors = buildMonitor.getActiveMonitors();

      expect(monitors).toHaveLength(2);
      expect(monitors[0]).toMatchObject({
        taskId: 'task-1',
        pluginName: 'TENET',
        monitorKey: 'build-task-1'
      });
      expect(monitors[0].startTime).toBeDefined();
    });
  });

  describe('isMonitoring', () => {
    it('should return true for active monitor', async () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      buildMonitor.startMonitoring('TENET', 'task-123', jest.fn());

      expect(buildMonitor.isMonitoring('task-123')).toBe(true);
    });

    it('should return false for inactive monitor', () => {
      expect(buildMonitor.isMonitoring('nonexistent')).toBe(false);
    });

    it('should return false after monitor is stopped', async () => {
      mockApiClient.getBuildStatus.mockResolvedValue({
        data: { state: BuildStatus.RUNNING }
      });

      buildMonitor.startMonitoring('TENET', 'task-123', jest.fn());
      buildMonitor.stopMonitoring('task-123');

      expect(buildMonitor.isMonitoring('task-123')).toBe(false);
    });
  });
});
