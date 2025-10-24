import { describe, it, expect } from 'vitest';
import {
  TASK_STATUS,
  getStatusClass,
  isTaskCompleted,
  isTaskRunning
} from '@/utils/taskStatus';

describe('taskStatus.js', () => {
  describe('TASK_STATUS constants', () => {
    it('should have all required status constants', () => {
      expect(TASK_STATUS.SUCCESS).toBe('SUCCESS');
      expect(TASK_STATUS.FAILURE).toBe('FAILURE');
      expect(TASK_STATUS.REVOKED).toBe('REVOKED');
      expect(TASK_STATUS.RETRY).toBe('RETRY');
      expect(TASK_STATUS.RUNNING).toBe('RUNNING');
      expect(TASK_STATUS.PENDING).toBe('PENDING');
      expect(TASK_STATUS.INSTALLING).toBe('INSTALLING');
    });
  });

  describe('getStatusClass', () => {
    it('should return "status-success" for SUCCESS', () => {
      expect(getStatusClass('SUCCESS')).toBe('status-success');
    });

    it('should return "status-failure" for FAILURE', () => {
      expect(getStatusClass('FAILURE')).toBe('status-failure');
    });

    it('should return "status-failure" for REVOKED', () => {
      expect(getStatusClass('REVOKED')).toBe('status-failure');
    });

    it('should return "status-failure" for RETRY', () => {
      expect(getStatusClass('RETRY')).toBe('status-failure');
    });

    it('should return "status-running" for RUNNING', () => {
      expect(getStatusClass('RUNNING')).toBe('status-running');
    });

    it('should return "status-running" for PENDING', () => {
      expect(getStatusClass('PENDING')).toBe('status-running');
    });

    it('should return "status-running" for INSTALLING', () => {
      expect(getStatusClass('INSTALLING')).toBe('status-running');
    });

    it('should return empty string for unknown status', () => {
      expect(getStatusClass('UNKNOWN')).toBe('');
      expect(getStatusClass('INVALID')).toBe('');
      expect(getStatusClass('')).toBe('');
    });

    it('should handle null and undefined', () => {
      expect(getStatusClass(null)).toBe('');
      expect(getStatusClass(undefined)).toBe('');
    });
  });

  describe('isTaskCompleted', () => {
    it('should return true for SUCCESS', () => {
      expect(isTaskCompleted('SUCCESS')).toBe(true);
    });

    it('should return true for FAILURE', () => {
      expect(isTaskCompleted('FAILURE')).toBe(true);
    });

    it('should return true for REVOKED', () => {
      expect(isTaskCompleted('REVOKED')).toBe(true);
    });

    it('should return true for RETRY', () => {
      expect(isTaskCompleted('RETRY')).toBe(true);
    });

    it('should return false for RUNNING', () => {
      expect(isTaskCompleted('RUNNING')).toBe(false);
    });

    it('should return false for PENDING', () => {
      expect(isTaskCompleted('PENDING')).toBe(false);
    });

    it('should return false for INSTALLING', () => {
      expect(isTaskCompleted('INSTALLING')).toBe(false);
    });

    it('should return false for unknown status', () => {
      expect(isTaskCompleted('UNKNOWN')).toBe(false);
      expect(isTaskCompleted('INVALID')).toBe(false);
      expect(isTaskCompleted('')).toBe(false);
    });

    it('should handle null and undefined', () => {
      expect(isTaskCompleted(null)).toBe(false);
      expect(isTaskCompleted(undefined)).toBe(false);
    });
  });

  describe('isTaskRunning', () => {
    it('should return true for RUNNING', () => {
      expect(isTaskRunning('RUNNING')).toBe(true);
    });

    it('should return true for PENDING', () => {
      expect(isTaskRunning('PENDING')).toBe(true);
    });

    it('should return true for INSTALLING', () => {
      expect(isTaskRunning('INSTALLING')).toBe(true);
    });

    it('should return false for SUCCESS', () => {
      expect(isTaskRunning('SUCCESS')).toBe(false);
    });

    it('should return false for FAILURE', () => {
      expect(isTaskRunning('FAILURE')).toBe(false);
    });

    it('should return false for REVOKED', () => {
      expect(isTaskRunning('REVOKED')).toBe(false);
    });

    it('should return false for RETRY', () => {
      expect(isTaskRunning('RETRY')).toBe(false);
    });

    it('should return false for unknown status', () => {
      expect(isTaskRunning('UNKNOWN')).toBe(false);
      expect(isTaskRunning('INVALID')).toBe(false);
      expect(isTaskRunning('')).toBe(false);
    });

    it('should handle null and undefined', () => {
      expect(isTaskRunning(null)).toBe(false);
      expect(isTaskRunning(undefined)).toBe(false);
    });
  });

  describe('isTaskCompleted and isTaskRunning should be mutually exclusive', () => {
    const allStatuses = [
      'SUCCESS', 'FAILURE', 'REVOKED', 'RETRY',
      'RUNNING', 'PENDING', 'INSTALLING'
    ];

    it('should not have overlapping true values', () => {
      allStatuses.forEach(status => {
        const completed = isTaskCompleted(status);
        const running = isTaskRunning(status);

        // Either completed or running can be true, but not both
        expect(completed && running).toBe(false);

        // At least one should be true for known statuses
        expect(completed || running).toBe(true);
      });
    });
  });
});
