import {
  enrichPluginData,
  enrichPlugins,
  isBuildingStatus,
  isTerminalBuildStatus
} from '@/utils/pluginDataProcessor';

describe('pluginDataProcessor', () => {
  describe('enrichPluginData', () => {
    it('should enrich plugin with user association when user is included', () => {
      const plugin = {
        id: 1,
        name: 'TENET',
        description: 'GRN analysis tool',
        users: [
          { username: 'user1' },
          { username: 'user2' }
        ],
        latest_build: {
          task_id: 'task-123',
          status: 'RUNNING'
        },
        version: '1.0.0'
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.checked).toBe(true);
      expect(result.building).toBe(true);
      expect(result.buildTaskId).toBe('task-123');
      expect(result.buildStatus).toBe('RUNNING');
      expect(result.current_version).toBe('1.0.0');
      expect(result.imageExists).toBe(false);
    });

    it('should mark plugin as not checked when user is not associated', () => {
      const plugin = {
        id: 1,
        name: 'GENIE3',
        users: [{ username: 'user2' }],
        latest_build: null,
        version: '2.0.0'
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.checked).toBe(false);
      expect(result.building).toBe(false);
      expect(result.buildTaskId).toBe(null);
      expect(result.buildStatus).toBe(null);
    });

    it('should handle plugin with PENDING build status', () => {
      const plugin = {
        id: 1,
        name: 'FastTENET',
        users: [],
        latest_build: {
          task_id: 'task-456',
          status: 'PENDING'
        }
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.checked).toBe(false);
      expect(result.building).toBe(true); // PENDING is considered building
      expect(result.buildStatus).toBe('PENDING');
    });

    it('should handle plugin with SUCCESS build status', () => {
      const plugin = {
        id: 1,
        name: 'GRNBOOST2',
        users: [],
        latest_build: {
          task_id: 'task-789',
          status: 'SUCCESS'
        }
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.building).toBe(false); // SUCCESS is not building
      expect(result.buildStatus).toBe('SUCCESS');
    });

    it('should handle plugin without latest_build', () => {
      const plugin = {
        id: 1,
        name: 'LEAP',
        users: [],
        latest_build: null
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.building).toBe(false);
      expect(result.buildTaskId).toBe(null);
      expect(result.buildStatus).toBe(null);
    });

    it('should handle plugin without users array', () => {
      const plugin = {
        id: 1,
        name: 'Scribe',
        users: null,
        latest_build: null
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.checked).toBe(false);
    });

    it('should handle plugin with empty users array', () => {
      const plugin = {
        id: 1,
        name: 'Plugin',
        users: [],
        latest_build: null
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.checked).toBe(false);
    });

    it('should use plugin.version for current_version if available', () => {
      const plugin = {
        id: 1,
        name: 'Plugin',
        users: [],
        latest_build: null,
        version: '3.0.0'
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.current_version).toBe('3.0.0');
    });

    it('should use plugin.current_version if version is not available', () => {
      const plugin = {
        id: 1,
        name: 'Plugin',
        users: [],
        latest_build: null,
        current_version: '2.5.0'
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.current_version).toBe('2.5.0');
    });

    it('should default to "latest" if no version info available', () => {
      const plugin = {
        id: 1,
        name: 'Plugin',
        users: [],
        latest_build: null
      };

      const result = enrichPluginData(plugin, 'user1');

      expect(result.current_version).toBe('latest');
    });
  });

  describe('enrichPlugins', () => {
    it('should enrich multiple plugins', () => {
      const plugins = [
        {
          id: 1,
          name: 'TENET',
          users: [{ username: 'user1' }],
          latest_build: { task_id: 'task-1', status: 'SUCCESS' }
        },
        {
          id: 2,
          name: 'GENIE3',
          users: [],
          latest_build: { task_id: 'task-2', status: 'RUNNING' }
        }
      ];

      const result = enrichPlugins(plugins, 'user1');

      expect(result).toHaveLength(2);
      expect(result[0].checked).toBe(true);
      expect(result[0].building).toBe(false);
      expect(result[1].checked).toBe(false);
      expect(result[1].building).toBe(true);
    });

    it('should return empty array for null input', () => {
      const result = enrichPlugins(null, 'user1');
      expect(result).toEqual([]);
    });

    it('should return empty array for undefined input', () => {
      const result = enrichPlugins(undefined, 'user1');
      expect(result).toEqual([]);
    });

    it('should handle empty array', () => {
      const result = enrichPlugins([], 'user1');
      expect(result).toEqual([]);
    });
  });

  describe('isBuildingStatus', () => {
    it('should return true for RUNNING status', () => {
      expect(isBuildingStatus('RUNNING')).toBe(true);
    });

    it('should return true for PENDING status', () => {
      expect(isBuildingStatus('PENDING')).toBe(true);
    });

    it('should return false for SUCCESS status', () => {
      expect(isBuildingStatus('SUCCESS')).toBe(false);
    });

    it('should return false for FAILURE status', () => {
      expect(isBuildingStatus('FAILURE')).toBe(false);
    });

    it('should return false for null or undefined', () => {
      expect(isBuildingStatus(null)).toBe(false);
      expect(isBuildingStatus(undefined)).toBe(false);
    });
  });

  describe('isTerminalBuildStatus', () => {
    it('should return true for SUCCESS status', () => {
      expect(isTerminalBuildStatus('SUCCESS')).toBe(true);
    });

    it('should return true for FAILURE status', () => {
      expect(isTerminalBuildStatus('FAILURE')).toBe(true);
    });

    it('should return true for REVOKED status', () => {
      expect(isTerminalBuildStatus('REVOKED')).toBe(true);
    });

    it('should return true for ERROR status', () => {
      expect(isTerminalBuildStatus('ERROR')).toBe(true);
    });

    it('should return false for RUNNING status', () => {
      expect(isTerminalBuildStatus('RUNNING')).toBe(false);
    });

    it('should return false for PENDING status', () => {
      expect(isTerminalBuildStatus('PENDING')).toBe(false);
    });

    it('should return false for unknown status', () => {
      expect(isTerminalBuildStatus('UNKNOWN')).toBe(false);
    });

    it('should return false for null or undefined', () => {
      expect(isTerminalBuildStatus(null)).toBe(false);
      expect(isTerminalBuildStatus(undefined)).toBe(false);
    });
  });
});
