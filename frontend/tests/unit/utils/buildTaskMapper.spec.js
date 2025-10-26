import { transformBuildTask, transformBuildTasks } from '@/utils/buildTaskMapper';

describe('buildTaskMapper', () => {
  describe('transformBuildTask', () => {
    it('should transform a single build task correctly', () => {
      const task = {
        task_id: '123',
        plugin_name: 'TENET',
        start_time: '2024-01-01T00:00:00',
        end_time: '2024-01-01T01:00:00',
        state: 'SUCCESS',
        error: null,
        info: { cpu: '100%' }
      };

      const durationCalculator = jest.fn(() => '1h 0m 0s');
      const result = transformBuildTask(task, durationCalculator);

      expect(result).toEqual({
        task_id: '123',
        plugin_name: 'TENET',
        start_time: '2024-01-01T00:00:00',
        end_time: '2024-01-01T01:00:00',
        running_time: '1h 0m 0s',
        status: 'SUCCESS',
        error: null,
        info: { cpu: '100%' }
      });

      expect(durationCalculator).toHaveBeenCalledWith(task);
    });

    it('should handle task with null values', () => {
      const task = {
        task_id: '456',
        plugin_name: 'GENIE3',
        start_time: null,
        end_time: null,
        state: 'PENDING',
        error: null,
        info: null
      };

      const durationCalculator = jest.fn(() => '-');
      const result = transformBuildTask(task, durationCalculator);

      expect(result).toEqual({
        task_id: '456',
        plugin_name: 'GENIE3',
        start_time: null,
        end_time: null,
        running_time: '-',
        status: 'PENDING',
        error: null,
        info: null
      });
    });

    it('should handle failed task with error message', () => {
      const task = {
        task_id: '789',
        plugin_name: 'FastTENET',
        start_time: '2024-01-01T00:00:00',
        end_time: '2024-01-01T00:05:00',
        state: 'FAILURE',
        error: 'Docker build failed',
        info: {}
      };

      const durationCalculator = jest.fn(() => '5m 0s');
      const result = transformBuildTask(task, durationCalculator);

      expect(result.status).toBe('FAILURE');
      expect(result.error).toBe('Docker build failed');
      expect(result.running_time).toBe('5m 0s');
    });
  });

  describe('transformBuildTasks', () => {
    it('should transform an array of build tasks', () => {
      const tasks = [
        {
          task_id: '1',
          plugin_name: 'TENET',
          start_time: '2024-01-01T00:00:00',
          end_time: '2024-01-01T01:00:00',
          state: 'SUCCESS',
          error: null,
          info: {}
        },
        {
          task_id: '2',
          plugin_name: 'GENIE3',
          start_time: '2024-01-01T00:00:00',
          end_time: null,
          state: 'RUNNING',
          error: null,
          info: {}
        }
      ];

      const durationCalculator = jest.fn()
        .mockReturnValueOnce('1h 0m 0s')
        .mockReturnValueOnce('30m 15s');

      const result = transformBuildTasks(tasks, durationCalculator);

      expect(result).toHaveLength(2);
      expect(result[0].task_id).toBe('1');
      expect(result[0].status).toBe('SUCCESS');
      expect(result[0].running_time).toBe('1h 0m 0s');

      expect(result[1].task_id).toBe('2');
      expect(result[1].status).toBe('RUNNING');
      expect(result[1].running_time).toBe('30m 15s');

      expect(durationCalculator).toHaveBeenCalledTimes(2);
    });

    it('should return empty array for null input', () => {
      const durationCalculator = jest.fn();
      const result = transformBuildTasks(null, durationCalculator);

      expect(result).toEqual([]);
      expect(durationCalculator).not.toHaveBeenCalled();
    });

    it('should return empty array for undefined input', () => {
      const durationCalculator = jest.fn();
      const result = transformBuildTasks(undefined, durationCalculator);

      expect(result).toEqual([]);
      expect(durationCalculator).not.toHaveBeenCalled();
    });

    it('should handle empty array', () => {
      const durationCalculator = jest.fn();
      const result = transformBuildTasks([], durationCalculator);

      expect(result).toEqual([]);
      expect(durationCalculator).not.toHaveBeenCalled();
    });

    it('should handle array with mixed task states', () => {
      const tasks = [
        { task_id: '1', plugin_name: 'P1', state: 'SUCCESS', start_time: '2024-01-01', end_time: '2024-01-02', error: null, info: {} },
        { task_id: '2', plugin_name: 'P2', state: 'FAILURE', start_time: '2024-01-01', end_time: '2024-01-01', error: 'Error', info: {} },
        { task_id: '3', plugin_name: 'P3', state: 'RUNNING', start_time: '2024-01-01', end_time: null, error: null, info: {} },
        { task_id: '4', plugin_name: 'P4', state: 'PENDING', start_time: null, end_time: null, error: null, info: {} }
      ];

      const durationCalculator = jest.fn((task) => {
        if (task.state === 'PENDING') return '-';
        if (task.state === 'RUNNING') return '10m';
        return '5m';
      });

      const result = transformBuildTasks(tasks, durationCalculator);

      expect(result).toHaveLength(4);
      expect(result[0].status).toBe('SUCCESS');
      expect(result[1].status).toBe('FAILURE');
      expect(result[2].status).toBe('RUNNING');
      expect(result[3].status).toBe('PENDING');
    });
  });
});
