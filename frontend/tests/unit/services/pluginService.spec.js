import { PluginService } from '@/services/pluginService';
import * as api from '@/api/index';
import { vi } from 'vitest';

// Mock the API module
vi.mock('@/api/index');

describe('PluginService', () => {
  let service;

  beforeEach(() => {
    service = new PluginService();
    jest.clearAllMocks();
  });

  describe('getUserProfile', () => {
    it('should return user profile data', async () => {
      const mockProfile = {
        username: 'testuser',
        email: 'test@example.com',
        role: 'admin'
      };

      api.getUser.mockResolvedValue({ data: mockProfile });

      const result = await service.getUserProfile();

      expect(result).toEqual(mockProfile);
      expect(api.getUser).toHaveBeenCalledTimes(1);
    });

    it('should propagate API errors', async () => {
      const error = new Error('API Error');
      api.getUser.mockRejectedValue(error);

      await expect(service.getUserProfile()).rejects.toThrow('API Error');
    });
  });

  describe('getPluginsList', () => {
    it('should return array of plugins', async () => {
      const mockPlugins = [
        { id: 1, name: 'TENET', version: '1.0.0' },
        { id: 2, name: 'GENIE3', version: '2.0.0' }
      ];

      api.getPlugins.mockResolvedValue({
        data: { plugins: mockPlugins }
      });

      const result = await service.getPluginsList();

      expect(result).toEqual(mockPlugins);
      expect(api.getPlugins).toHaveBeenCalledTimes(1);
    });

    it('should handle empty plugin list', async () => {
      api.getPlugins.mockResolvedValue({
        data: { plugins: [] }
      });

      const result = await service.getPluginsList();

      expect(result).toEqual([]);
    });

    it('should propagate API errors', async () => {
      const error = new Error('Network Error');
      api.getPlugins.mockRejectedValue(error);

      await expect(service.getPluginsList()).rejects.toThrow('Network Error');
    });
  });

  describe('checkPluginImage', () => {
    it('should return true when image exists', async () => {
      api.checkPluginImage.mockResolvedValue({
        data: { image_exists: true }
      });

      const result = await service.checkPluginImage('TENET');

      expect(result).toBe(true);
      expect(api.checkPluginImage).toHaveBeenCalledWith('TENET');
    });

    it('should return false when image does not exist', async () => {
      api.checkPluginImage.mockResolvedValue({
        data: { image_exists: false }
      });

      const result = await service.checkPluginImage('GENIE3');

      expect(result).toBe(false);
      expect(api.checkPluginImage).toHaveBeenCalledWith('GENIE3');
    });

    it('should propagate API errors', async () => {
      const error = new Error('Image check failed');
      api.checkPluginImage.mockRejectedValue(error);

      await expect(service.checkPluginImage('FastTENET')).rejects.toThrow('Image check failed');
    });
  });

  describe('checkMultiplePluginImages', () => {
    it('should check images for multiple plugins', async () => {
      const plugins = [
        { id: 1, name: 'TENET' },
        { id: 2, name: 'GENIE3' }
      ];

      api.checkPluginImage
        .mockResolvedValueOnce({ data: { image_exists: true } })
        .mockResolvedValueOnce({ data: { image_exists: false } });

      const result = await service.checkMultiplePluginImages(plugins);

      expect(result).toEqual([
        { id: 1, name: 'TENET', imageExists: true },
        { id: 2, name: 'GENIE3', imageExists: false }
      ]);

      expect(api.checkPluginImage).toHaveBeenCalledTimes(2);
      expect(api.checkPluginImage).toHaveBeenCalledWith('TENET');
      expect(api.checkPluginImage).toHaveBeenCalledWith('GENIE3');
    });

    it('should handle API errors gracefully and set imageExists to false', async () => {
      const plugins = [
        { id: 1, name: 'TENET' },
        { id: 2, name: 'GENIE3' }
      ];

      const consoleSpy = jest.spyOn(console, 'error').mockImplementation(() => {});

      api.checkPluginImage
        .mockResolvedValueOnce({ data: { image_exists: true } })
        .mockRejectedValueOnce(new Error('Network error'));

      const result = await service.checkMultiplePluginImages(plugins);

      expect(result).toEqual([
        { id: 1, name: 'TENET', imageExists: true },
        { id: 2, name: 'GENIE3', imageExists: false }
      ]);

      expect(consoleSpy).toHaveBeenCalled();
      consoleSpy.mockRestore();
    });

    it('should handle empty plugin array', async () => {
      const result = await service.checkMultiplePluginImages([]);

      expect(result).toEqual([]);
      expect(api.checkPluginImage).not.toHaveBeenCalled();
    });

    it('should process all plugins in parallel', async () => {
      const plugins = [
        { name: 'plugin1' },
        { name: 'plugin2' },
        { name: 'plugin3' }
      ];

      api.checkPluginImage.mockResolvedValue({ data: { image_exists: true } });

      const startTime = Date.now();
      await service.checkMultiplePluginImages(plugins);
      const endTime = Date.now();

      // Should complete quickly since it's parallel (not 3x the time of single call)
      expect(endTime - startTime).toBeLessThan(1000);
      expect(api.checkPluginImage).toHaveBeenCalledTimes(3);
    });
  });

  describe('getCompletePluginData', () => {
    it('should fetch user profile and plugins in parallel', async () => {
      const mockProfile = { username: 'testuser' };
      const mockPlugins = [
        { id: 1, name: 'TENET' },
        { id: 2, name: 'GENIE3' }
      ];

      api.getUser.mockResolvedValue({ data: mockProfile });
      api.getPlugins.mockResolvedValue({ data: { plugins: mockPlugins } });

      const result = await service.getCompletePluginData('testuser');

      expect(result).toEqual({
        profile: mockProfile,
        plugins: mockPlugins
      });

      expect(api.getUser).toHaveBeenCalledTimes(1);
      expect(api.getPlugins).toHaveBeenCalledTimes(1);
    });

    it('should fail if user profile fetch fails', async () => {
      const error = new Error('User fetch failed');
      api.getUser.mockRejectedValue(error);
      api.getPlugins.mockResolvedValue({ data: { plugins: [] } });

      await expect(service.getCompletePluginData('testuser')).rejects.toThrow('User fetch failed');
    });

    it('should fail if plugins fetch fails', async () => {
      const error = new Error('Plugins fetch failed');
      api.getUser.mockResolvedValue({ data: { username: 'testuser' } });
      api.getPlugins.mockRejectedValue(error);

      await expect(service.getCompletePluginData('testuser')).rejects.toThrow('Plugins fetch failed');
    });

    it('should enrich plugins with user associations and build status', async () => {
      const mockProfile = { username: 'testuser' };
      const mockPlugins = [
        {
          id: 1,
          name: 'TENET',
          users: [{ username: 'testuser' }],
          latest_build: { task_id: 'task-1', status: 'RUNNING' }
        },
        {
          id: 2,
          name: 'GENIE3',
          users: [],
          latest_build: { task_id: 'task-2', status: 'SUCCESS' }
        }
      ];

      api.getUser.mockResolvedValue({ data: mockProfile });
      api.getPlugins.mockResolvedValue({ data: { plugins: mockPlugins } });

      const result = await service.getCompletePluginData('testuser');

      expect(result.plugins).toHaveLength(2);
      expect(result.plugins[0].checked).toBe(true); // User associated
      expect(result.plugins[0].building).toBe(true); // RUNNING status
      expect(result.plugins[1].checked).toBe(false); // Not associated
      expect(result.plugins[1].building).toBe(false); // SUCCESS status (not building)
    });
  });

  describe('associatePlugin', () => {
    it('should associate plugin successfully', async () => {
      const mockResponse = { message: 'Plugin associated' };
      api.associatePlugin.mockResolvedValue({ data: mockResponse });

      const result = await service.associatePlugin(123);

      expect(result).toEqual(mockResponse);
      expect(api.associatePlugin).toHaveBeenCalledWith(123);
    });

    it('should throw error on failure', async () => {
      const error = new Error('Association failed');
      api.associatePlugin.mockRejectedValue(error);

      await expect(service.associatePlugin(123)).rejects.toThrow('Association failed');
    });
  });

  describe('dissociatePlugin', () => {
    it('should dissociate plugin successfully', async () => {
      const mockResponse = { message: 'Plugin dissociated' };
      api.dissociatePlugin.mockResolvedValue({ data: mockResponse });

      const result = await service.dissociatePlugin(123);

      expect(result).toEqual(mockResponse);
      expect(api.dissociatePlugin).toHaveBeenCalledWith(123);
    });

    it('should throw error on failure', async () => {
      const error = new Error('Dissociation failed');
      api.dissociatePlugin.mockRejectedValue(error);

      await expect(service.dissociatePlugin(123)).rejects.toThrow('Dissociation failed');
    });
  });

  describe('getBuildTasks', () => {
    it('should get all build tasks when no plugin name provided', async () => {
      const mockTasks = [
        { task_id: 'task-1', plugin_name: 'TENET', status: 'RUNNING' },
        { task_id: 'task-2', plugin_name: 'GENIE3', status: 'SUCCESS' }
      ];
      api.getBuildTasks.mockResolvedValue({ data: { tasks: mockTasks } });

      const result = await service.getBuildTasks();

      expect(result).toEqual(mockTasks);
      expect(api.getBuildTasks).toHaveBeenCalledWith(null);
    });

    it('should get filtered build tasks for specific plugin', async () => {
      const mockTasks = [
        { task_id: 'task-1', plugin_name: 'TENET', status: 'RUNNING' }
      ];
      api.getBuildTasks.mockResolvedValue({ data: { tasks: mockTasks } });

      const result = await service.getBuildTasks('TENET');

      expect(result).toEqual(mockTasks);
      expect(api.getBuildTasks).toHaveBeenCalledWith('TENET');
    });

    it('should return empty array when no tasks in response', async () => {
      api.getBuildTasks.mockResolvedValue({ data: {} });

      const result = await service.getBuildTasks();

      expect(result).toEqual([]);
    });

    it('should throw error on failure', async () => {
      const error = new Error('Failed to fetch build tasks');
      api.getBuildTasks.mockRejectedValue(error);

      await expect(service.getBuildTasks()).rejects.toThrow('Failed to fetch build tasks');
    });
  });
});
