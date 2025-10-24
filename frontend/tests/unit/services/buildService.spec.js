import { BuildService } from '@/services/buildService';

describe('BuildService', () => {
  let service;
  let mockApiClient;

  beforeEach(() => {
    mockApiClient = {
      buildPluginDocker: jest.fn(),
      getBuildLogs: jest.fn(),
      checkPluginImage: jest.fn()
    };

    service = new BuildService(mockApiClient);
  });

  describe('constructor', () => {
    it('should create service with provided API client', () => {
      expect(service.apiClient).toBe(mockApiClient);
    });

    it('should create service with default API client', () => {
      const defaultService = new BuildService();

      expect(defaultService.apiClient).toBeDefined();
      expect(defaultService.apiClient.buildPluginDocker).toBeDefined();
    });
  });

  describe('buildPlugin', () => {
    it('should build plugin successfully', async () => {
      mockApiClient.buildPluginDocker.mockResolvedValue({
        data: { task_id: 'task-123' }
      });

      const result = await service.buildPlugin('TENET', false);

      expect(result).toEqual({
        success: true,
        plugin: 'TENET',
        taskId: 'task-123'
      });
      expect(mockApiClient.buildPluginDocker).toHaveBeenCalledWith('TENET', false);
    });

    it('should build plugin with GPU enabled', async () => {
      mockApiClient.buildPluginDocker.mockResolvedValue({
        data: { task_id: 'task-456' }
      });

      const result = await service.buildPlugin('GENIE3', true);

      expect(result).toEqual({
        success: true,
        plugin: 'GENIE3',
        taskId: 'task-456'
      });
      expect(mockApiClient.buildPluginDocker).toHaveBeenCalledWith('GENIE3', true);
    });

    it('should handle build failure', async () => {
      const error = new Error('Build failed');
      mockApiClient.buildPluginDocker.mockRejectedValue(error);

      const result = await service.buildPlugin('TENET', false);

      expect(result).toEqual({
        success: false,
        plugin: 'TENET',
        error: error
      });
    });

    it('should log error on failure', async () => {
      const consoleErrorSpy = jest.spyOn(console, 'error').mockImplementation(() => {});
      const error = new Error('Build failed');
      mockApiClient.buildPluginDocker.mockRejectedValue(error);

      await service.buildPlugin('TENET', false);

      expect(consoleErrorSpy).toHaveBeenCalledWith('Error building plugin TENET:', error);

      consoleErrorSpy.mockRestore();
    });
  });

  describe('buildMultiplePlugins', () => {
    it('should return empty array for empty input', async () => {
      const result = await service.buildMultiplePlugins([], false);

      expect(result).toEqual([]);
    });

    it('should return empty array for non-array input', async () => {
      const result = await service.buildMultiplePlugins(null, false);

      expect(result).toEqual([]);
    });

    it('should build multiple plugins successfully', async () => {
      mockApiClient.buildPluginDocker
        .mockResolvedValueOnce({ data: { task_id: 'task-1' } })
        .mockResolvedValueOnce({ data: { task_id: 'task-2' } });

      mockApiClient.checkPluginImage
        .mockResolvedValueOnce({ data: { image_exists: true } })
        .mockResolvedValueOnce({ data: { image_exists: false } });

      const plugins = [
        { name: 'TENET' },
        { name: 'GENIE3' }
      ];

      const results = await service.buildMultiplePlugins(plugins, false);

      expect(results).toHaveLength(2);
      expect(results[0]).toEqual({
        success: true,
        plugin: 'TENET',
        taskId: 'task-1',
        imageExists: true
      });
      expect(results[1]).toEqual({
        success: true,
        plugin: 'GENIE3',
        taskId: 'task-2',
        imageExists: false
      });
    });

    it('should handle partial failures', async () => {
      mockApiClient.buildPluginDocker
        .mockResolvedValueOnce({ data: { task_id: 'task-1' } })
        .mockRejectedValueOnce(new Error('Build failed'));

      mockApiClient.checkPluginImage
        .mockResolvedValueOnce({ data: { image_exists: true } });

      const plugins = [
        { name: 'TENET' },
        { name: 'GENIE3' }
      ];

      const results = await service.buildMultiplePlugins(plugins, false);

      expect(results).toHaveLength(2);
      expect(results[0].success).toBe(true);
      expect(results[1].success).toBe(false);
    });

    it('should call progress callback', async () => {
      mockApiClient.buildPluginDocker
        .mockResolvedValueOnce({ data: { task_id: 'task-1' } })
        .mockResolvedValueOnce({ data: { task_id: 'task-2' } });

      mockApiClient.checkPluginImage
        .mockResolvedValue({ data: { image_exists: false } });

      const plugins = [
        { name: 'TENET' },
        { name: 'GENIE3' }
      ];

      const onProgress = jest.fn();

      await service.buildMultiplePlugins(plugins, false, onProgress);

      expect(onProgress).toHaveBeenCalledWith(1, 2);
      expect(onProgress).toHaveBeenCalledWith(2, 2);
      expect(onProgress).toHaveBeenCalledTimes(2);
    });

    it('should call progress callback even on failure', async () => {
      mockApiClient.buildPluginDocker
        .mockResolvedValueOnce({ data: { task_id: 'task-1' } })
        .mockRejectedValueOnce(new Error('Build failed'));

      mockApiClient.checkPluginImage
        .mockResolvedValue({ data: { image_exists: false } });

      const plugins = [
        { name: 'TENET' },
        { name: 'GENIE3' }
      ];

      const onProgress = jest.fn();

      await service.buildMultiplePlugins(plugins, false, onProgress);

      expect(onProgress).toHaveBeenCalledTimes(2);
    });

    it('should handle image check failure gracefully', async () => {
      const consoleErrorSpy = jest.spyOn(console, 'error').mockImplementation(() => {});

      mockApiClient.buildPluginDocker
        .mockResolvedValueOnce({ data: { task_id: 'task-1' } });

      mockApiClient.checkPluginImage
        .mockRejectedValueOnce(new Error('Image check failed'));

      const plugins = [{ name: 'TENET' }];

      const results = await service.buildMultiplePlugins(plugins, false);

      expect(results[0]).toEqual({
        success: true,
        plugin: 'TENET',
        taskId: 'task-1',
        imageExists: false
      });

      expect(consoleErrorSpy).toHaveBeenCalled();

      consoleErrorSpy.mockRestore();
    });
  });

  describe('getBuildLogs', () => {
    it('should retrieve build logs successfully', async () => {
      mockApiClient.getBuildLogs.mockResolvedValue({
        data: { log_content: 'Build log content' }
      });

      const result = await service.getBuildLogs('TENET');

      expect(result).toEqual({
        logs: 'Build log content',
        success: true
      });
      expect(mockApiClient.getBuildLogs).toHaveBeenCalledWith('TENET');
    });

    it('should return default message for empty logs', async () => {
      mockApiClient.getBuildLogs.mockResolvedValue({
        data: { log_content: '' }
      });

      const result = await service.getBuildLogs('TENET');

      expect(result).toEqual({
        logs: 'No logs available',
        success: true
      });
    });

    it('should handle log retrieval failure', async () => {
      const error = new Error('Failed to get logs');
      mockApiClient.getBuildLogs.mockRejectedValue(error);

      const result = await service.getBuildLogs('TENET');

      expect(result).toEqual({
        logs: null,
        error: 'Failed to get logs',
        success: false
      });
    });

    it('should log error on failure', async () => {
      const consoleErrorSpy = jest.spyOn(console, 'error').mockImplementation(() => {});
      const error = new Error('Failed to get logs');
      mockApiClient.getBuildLogs.mockRejectedValue(error);

      await service.getBuildLogs('TENET');

      expect(consoleErrorSpy).toHaveBeenCalledWith('Error fetching build logs for plugin TENET:', error);

      consoleErrorSpy.mockRestore();
    });
  });

  describe('processBuildResults', () => {
    it('should process empty results', () => {
      const summary = service.processBuildResults([]);

      expect(summary).toEqual({
        total: 0,
        successful: 0,
        failed: 0,
        successfulPlugins: [],
        failedPlugins: []
      });
    });

    it('should process null input', () => {
      const summary = service.processBuildResults(null);

      expect(summary).toEqual({
        total: 0,
        successful: 0,
        failed: 0,
        successfulPlugins: [],
        failedPlugins: []
      });
    });

    it('should process all successful results', () => {
      const results = [
        { success: true, plugin: 'TENET' },
        { success: true, plugin: 'GENIE3' }
      ];

      const summary = service.processBuildResults(results);

      expect(summary).toEqual({
        total: 2,
        successful: 2,
        failed: 0,
        successfulPlugins: ['TENET', 'GENIE3'],
        failedPlugins: []
      });
    });

    it('should process all failed results', () => {
      const results = [
        { success: false, plugin: 'TENET' },
        { success: false, plugin: 'GENIE3' }
      ];

      const summary = service.processBuildResults(results);

      expect(summary).toEqual({
        total: 2,
        successful: 0,
        failed: 2,
        successfulPlugins: [],
        failedPlugins: ['TENET', 'GENIE3']
      });
    });

    it('should process mixed results', () => {
      const results = [
        { success: true, plugin: 'TENET' },
        { success: false, plugin: 'GENIE3' },
        { success: true, plugin: 'LEAP' }
      ];

      const summary = service.processBuildResults(results);

      expect(summary).toEqual({
        total: 3,
        successful: 2,
        failed: 1,
        successfulPlugins: ['TENET', 'LEAP'],
        failedPlugins: ['GENIE3']
      });
    });
  });

  describe('getBuildSummaryMessage', () => {
    it('should return info message for no builds', () => {
      const summary = { total: 0, successful: 0, failed: 0 };

      const message = service.getBuildSummaryMessage(summary);

      expect(message).toEqual({
        text: 'No plugins were built.',
        type: 'info'
      });
    });

    it('should return success message for all successful', () => {
      const summary = { total: 5, successful: 5, failed: 0 };

      const message = service.getBuildSummaryMessage(summary);

      expect(message).toEqual({
        text: 'All 5 plugin(s) were built successfully!',
        type: 'success'
      });
    });

    it('should return error message for all failed', () => {
      const summary = { total: 3, successful: 0, failed: 3 };

      const message = service.getBuildSummaryMessage(summary);

      expect(message).toEqual({
        text: 'All 3 plugin(s) failed to build.',
        type: 'error'
      });
    });

    it('should return warning message for mixed results', () => {
      const summary = { total: 5, successful: 3, failed: 2 };

      const message = service.getBuildSummaryMessage(summary);

      expect(message).toEqual({
        text: '3 plugin(s) built successfully, 2 plugin(s) failed.',
        type: 'warning'
      });
    });
  });

  describe('filterPluginsToBuild', () => {
    it('should return empty array for null input', () => {
      const result = service.filterPluginsToBuild(null);

      expect(result).toEqual([]);
    });

    it('should return empty array for non-array input', () => {
      const result = service.filterPluginsToBuild('not an array');

      expect(result).toEqual([]);
    });

    it('should filter out building plugins', () => {
      const plugins = [
        { name: 'TENET', building: false, imageExists: false },
        { name: 'GENIE3', building: true, imageExists: false },
        { name: 'LEAP', building: false, imageExists: false }
      ];

      const result = service.filterPluginsToBuild(plugins);

      expect(result).toHaveLength(2);
      expect(result.map(p => p.name)).toEqual(['TENET', 'LEAP']);
    });

    it('should filter out plugins with existing images', () => {
      const plugins = [
        { name: 'TENET', building: false, imageExists: true },
        { name: 'GENIE3', building: false, imageExists: false },
        { name: 'LEAP', building: false, imageExists: false }
      ];

      const result = service.filterPluginsToBuild(plugins);

      expect(result).toHaveLength(2);
      expect(result.map(p => p.name)).toEqual(['GENIE3', 'LEAP']);
    });

    it('should filter out both building and existing image plugins', () => {
      const plugins = [
        { name: 'TENET', building: false, imageExists: true },
        { name: 'GENIE3', building: true, imageExists: false },
        { name: 'LEAP', building: false, imageExists: false },
        { name: 'Scribe', building: false, imageExists: true }
      ];

      const result = service.filterPluginsToBuild(plugins);

      expect(result).toHaveLength(1);
      expect(result[0].name).toBe('LEAP');
    });

    it('should return all plugins if none are building or have images', () => {
      const plugins = [
        { name: 'TENET', building: false, imageExists: false },
        { name: 'GENIE3', building: false, imageExists: false }
      ];

      const result = service.filterPluginsToBuild(plugins);

      expect(result).toHaveLength(2);
    });

    it('should return empty array if all plugins are building or have images', () => {
      const plugins = [
        { name: 'TENET', building: true, imageExists: false },
        { name: 'GENIE3', building: false, imageExists: true }
      ];

      const result = service.filterPluginsToBuild(plugins);

      expect(result).toEqual([]);
    });
  });

  describe('delay', () => {
    it('should delay execution', async () => {
      jest.useFakeTimers();

      const promise = service.delay(1000);
      jest.advanceTimersByTime(1000);

      await promise;

      expect(true).toBe(true); // If we get here, delay worked

      jest.useRealTimers();
    });
  });

  describe('cancelBuildTask', () => {
    it('should cancel build task successfully', async () => {
      mockApiClient.cancelBuildTask.mockResolvedValue({ data: { message: 'Cancelled' } });

      const result = await service.cancelBuildTask('task-123');

      expect(result.success).toBe(true);
      expect(result.message).toBe('Build task has been cancelled.');
      expect(mockApiClient.cancelBuildTask).toHaveBeenCalledWith('task-123');
    });

    it('should handle cancellation failure', async () => {
      const error = new Error('Cancellation failed');
      mockApiClient.cancelBuildTask.mockRejectedValue(error);

      const result = await service.cancelBuildTask('task-123');

      expect(result.success).toBe(false);
      expect(result.error).toContain('Failed to cancel build task');
      expect(result.error).toContain('Cancellation failed');
    });

    it('should log error on failure', async () => {
      const consoleErrorSpy = jest.spyOn(console, 'error').mockImplementation(() => {});
      const error = new Error('Cancellation failed');
      mockApiClient.cancelBuildTask.mockRejectedValue(error);

      await service.cancelBuildTask('task-123');

      expect(consoleErrorSpy).toHaveBeenCalledWith('Failed to cancel build task:', error);

      consoleErrorSpy.mockRestore();
    });

    it('should handle API error with status code', async () => {
      const error = {
        response: {
          status: 404,
          data: { detail: 'Task not found' }
        }
      };
      mockApiClient.cancelBuildTask.mockRejectedValue(error);

      const result = await service.cancelBuildTask('task-123');

      expect(result.success).toBe(false);
      expect(result.error).toContain('Task not found');
    });
  });
});
