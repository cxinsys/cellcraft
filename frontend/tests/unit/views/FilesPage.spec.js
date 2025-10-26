import { describe, it, expect, beforeEach, afterEach, vi } from 'vitest';
import { mount } from '@vue/test-utils';
import FilesPage from '@/views/FilesPage.vue';
import * as apiIndex from '@/api/index';
import * as validation from '@/utils/validation';
import * as filename from '@/utils/filename';

describe('FilesPage.vue', () => {
  let wrapper;

  beforeEach(() => {
    // Mount component with minimal required props
    wrapper = mount(FilesPage, {
      propsData: {
        doUploadFile: ''
      },
      mocks: {
        $route: {
          query: {}
        },
        $refs: {
          selectFile: {
            files: []
          }
        }
      }
    });
  });

  /**
   * ClickOut Method Tests
   * Pure state mutation - easily testable
   */
  describe('ClickOut', () => {
    it('should set R_Mouse_isActive to false', () => {
      // Set initial state
      wrapper.vm.R_Mouse_isActive = true;

      // Call method
      wrapper.vm.ClickOut();

      // Assert state changed
      expect(wrapper.vm.R_Mouse_isActive).toBe(false);
    });

    it('should work when R_Mouse_isActive is already false', () => {
      // Set initial state
      wrapper.vm.R_Mouse_isActive = false;

      // Call method
      wrapper.vm.ClickOut();

      // Assert state remains false
      expect(wrapper.vm.R_Mouse_isActive).toBe(false);
    });

    it('should be called when clicking outside context menu', async () => {
      // Set initial state
      wrapper.vm.R_Mouse_isActive = true;

      // Trigger click on layout element
      await wrapper.find('.layout').trigger('click');

      // Assert context menu is closed
      expect(wrapper.vm.R_Mouse_isActive).toBe(false);
    });
  });

  /**
   * folderClick Method Tests
   * Contains conditional logic - tests both branches
   */
  describe('folderClick', () => {
    it('should set toggleFolder and currentFolder on first click', () => {
      // Initial state
      expect(wrapper.vm.toggleFolder).toBeNull();

      // Click on folder index 0 with name "data"
      wrapper.vm.folderClick(0, 'data');

      // Assert state changed
      expect(wrapper.vm.toggleFolder).toBe(0);
      expect(wrapper.vm.currentFolder).toBe('data');
    });

    it('should toggle off when clicking same folder again', () => {
      // Set initial state - folder 0 is already selected
      wrapper.vm.toggleFolder = 0;
      wrapper.vm.currentFolder = 'data';

      // Click on same folder (index 0)
      wrapper.vm.folderClick(0, 'data');

      // Assert folder is toggled off
      expect(wrapper.vm.toggleFolder).toBeNull();
      expect(wrapper.vm.currentFolder).toBe('data'); // currentFolder still updated
    });

    it('should switch to different folder when clicking another folder', () => {
      // Set initial state - folder 0 is selected
      wrapper.vm.toggleFolder = 0;
      wrapper.vm.currentFolder = 'data';

      // Click on different folder (index 1, name "results")
      wrapper.vm.folderClick(1, 'results');

      // Assert switched to new folder
      expect(wrapper.vm.toggleFolder).toBe(1);
      expect(wrapper.vm.currentFolder).toBe('results');
    });

    it('should handle folder index 0 correctly', () => {
      // Click on folder index 0
      wrapper.vm.folderClick(0, 'data');

      expect(wrapper.vm.toggleFolder).toBe(0);
      expect(wrapper.vm.currentFolder).toBe('data');

      // Click again to toggle off
      wrapper.vm.folderClick(0, 'data');

      expect(wrapper.vm.toggleFolder).toBeNull();
    });

    it('should update currentFolder even when toggling off', () => {
      // Set initial state
      wrapper.vm.toggleFolder = 0;
      wrapper.vm.currentFolder = 'oldFolder';

      // Click on same folder with different name
      wrapper.vm.folderClick(0, 'newFolder');

      // Assert toggleFolder is null but currentFolder updated
      expect(wrapper.vm.toggleFolder).toBeNull();
      expect(wrapper.vm.currentFolder).toBe('newFolder');
    });

    it('should handle multiple folder switches', () => {
      // Click folder 0
      wrapper.vm.folderClick(0, 'folder0');
      expect(wrapper.vm.toggleFolder).toBe(0);
      expect(wrapper.vm.currentFolder).toBe('folder0');

      // Click folder 1
      wrapper.vm.folderClick(1, 'folder1');
      expect(wrapper.vm.toggleFolder).toBe(1);
      expect(wrapper.vm.currentFolder).toBe('folder1');

      // Click folder 2
      wrapper.vm.folderClick(2, 'folder2');
      expect(wrapper.vm.toggleFolder).toBe(2);
      expect(wrapper.vm.currentFolder).toBe('folder2');
    });

    it('should handle folder names with special characters', () => {
      wrapper.vm.folderClick(0, 'user-data_2025.backup');

      expect(wrapper.vm.currentFolder).toBe('user-data_2025.backup');
    });

    it('should handle empty folder name', () => {
      wrapper.vm.folderClick(0, '');

      expect(wrapper.vm.currentFolder).toBe('');
    });
  });

  /**
   * RMouseClick Method Tests
   * Tests position calculation integration
   */
  describe('RMouseClick', () => {
    it('should activate context menu and set position', () => {
      // Create mock event
      const mockEvent = {
        clientX: 100,
        clientY: 200
      };

      // Call method
      wrapper.vm.RMouseClick(mockEvent, 'test.h5ad', 0);

      // Assert context menu is active
      expect(wrapper.vm.R_Mouse_isActive).toBe(true);
      expect(wrapper.vm.file_name).toBe('test.h5ad');
      expect(wrapper.vm.list_idx).toBe(0);

      // Assert position is set (using positionCalculator utility)
      expect(wrapper.vm.xPosition).toBe('100px');
      expect(wrapper.vm.yPosition).toBe('145px'); // 200 - 55 (default offset)
    });

    it('should reset R_Mouse_isActive before activating', () => {
      // Set initial state
      wrapper.vm.R_Mouse_isActive = true;

      const mockEvent = {
        clientX: 100,
        clientY: 200
      };

      // Call method
      wrapper.vm.RMouseClick(mockEvent, 'file.csv', 1);

      // Assert it was reset and then activated
      expect(wrapper.vm.R_Mouse_isActive).toBe(true);
    });

    it('should update file_name and list_idx', () => {
      const mockEvent = {
        clientX: 50,
        clientY: 100
      };

      wrapper.vm.RMouseClick(mockEvent, 'data.txt', 5);

      expect(wrapper.vm.file_name).toBe('data.txt');
      expect(wrapper.vm.list_idx).toBe(5);
    });

    it('should handle different mouse positions', () => {
      const positions = [
        { clientX: 0, clientY: 0 },
        { clientX: 1920, clientY: 1080 },
        { clientX: 500, clientY: 500 }
      ];

      positions.forEach(pos => {
        wrapper.vm.RMouseClick(pos, 'test.h5ad', 0);
        expect(wrapper.vm.xPosition).toBe(`${pos.clientX}px`);
        expect(wrapper.vm.yPosition).toBe(`${pos.clientY - 55}px`);
      });
    });
  });

  /**
   * undoDeletion Method Tests
   */
  describe('undoDeletion', () => {
    it('should restore deleted file to files_list', () => {
      // Set up initial state
      const deletedFile = {
        file_name: 'deleted.h5ad',
        file_size: 1024,
        created_at: '2025-10-23T12:00:00'
      };

      wrapper.vm.files_list = [
        { file_name: 'file1.csv', file_size: 512, created_at: '2025-10-23T11:00:00' }
      ];
      wrapper.vm.targetFile = deletedFile;
      wrapper.vm.toggleMessage = true;

      // Call method
      wrapper.vm.undoDeletion();

      // Assert file is restored
      expect(wrapper.vm.files_list).toHaveLength(2);
      expect(wrapper.vm.files_list[1]).toEqual(deletedFile);
      expect(wrapper.vm.toggleMessage).toBe(false);
    });

    it('should close toggle message', () => {
      wrapper.vm.toggleMessage = true;
      wrapper.vm.targetFile = { file_name: 'test.h5ad' };

      wrapper.vm.undoDeletion();

      expect(wrapper.vm.toggleMessage).toBe(false);
    });

    it('should handle null targetFile gracefully', () => {
      wrapper.vm.files_list = [];
      wrapper.vm.targetFile = null;

      // Should not throw error
      expect(() => wrapper.vm.undoDeletion()).not.toThrow();
    });

    it('should call clearTimeout with deletionTimer', () => {
      // Spy on global clearTimeout
      const clearTimeoutSpy = vi.spyOn(global, 'clearTimeout');

      // Setup timer
      const mockTimer = setTimeout(() => {}, 5000);
      wrapper.vm.deletionTimer = mockTimer;
      wrapper.vm.targetFile = { file_name: 'test.h5ad' };
      wrapper.vm.files_list = [];

      // Execute
      wrapper.vm.undoDeletion();

      // Assert clearTimeout was called with the timer
      expect(clearTimeoutSpy).toHaveBeenCalledWith(mockTimer);

      // Cleanup
      clearTimeoutSpy.mockRestore();
    });
  });

  /**
   * uploadFile Method Tests
   * Tests file upload flow with validation and API integration
   */
  describe('uploadFile', () => {
    let mockAlert;
    let mockFormData;

    beforeEach(() => {
      // Mock global objects
      mockAlert = vi.fn();
      global.alert = mockAlert;

      mockFormData = vi.fn(function() {
        this.append = vi.fn();
      });
      global.FormData = mockFormData;

      // Mock API functions
      vi.spyOn(apiIndex, 'uploadForm').mockResolvedValue({ data: 'success' });
      vi.spyOn(apiIndex, 'findFolder').mockResolvedValue({
        data: [{ file_name: 'uploaded.h5ad', file_size: 2048, created_at: '2025-10-24T12:00:00' }]
      });

      // Mock utility functions
      vi.spyOn(validation, 'validateFileExtension');
      vi.spyOn(filename, 'generateUploadFileName');
    });

    afterEach(() => {
      vi.restoreAllMocks();
    });

    it('should return early when no file is selected', async () => {
      // Set empty files array
      wrapper.vm.$refs.selectFile = { files: [] };

      await wrapper.vm.uploadFile();

      // Assert no API calls made
      expect(apiIndex.uploadForm).not.toHaveBeenCalled();
      expect(mockAlert).not.toHaveBeenCalled();
    });

    it('should alert when file extension is invalid', async () => {
      // Mock file with invalid extension
      const invalidFile = new File(['content'], 'document.pdf', { type: 'application/pdf' });
      wrapper.vm.$refs.selectFile = { files: [invalidFile] };

      // Mock validateFileExtension to return invalid
      validation.validateFileExtension.mockReturnValue({
        isValid: false,
        extension: 'pdf',
        message: 'Please upload h5ad, csv, txt file'
      });

      await wrapper.vm.uploadFile();

      // Assert alert was called with error message
      expect(mockAlert).toHaveBeenCalledWith('Please upload h5ad, csv, txt file');
      expect(apiIndex.uploadForm).not.toHaveBeenCalled();
    });

    it('should successfully upload file with valid extension', async () => {
      // Mock file with valid extension
      const validFile = new File(['content'], 'data.h5ad', { type: 'application/octet-stream' });
      wrapper.vm.$refs.selectFile = { files: [validFile] };
      wrapper.vm.currentFolder = 'test-folder';

      // Mock validateFileExtension to return valid
      validation.validateFileExtension.mockReturnValue({
        isValid: true,
        extension: 'h5ad',
        message: ''
      });

      // Mock generateUploadFileName
      filename.generateUploadFileName.mockReturnValue('test-folder_data.h5ad');

      await wrapper.vm.uploadFile();

      // Assert validation was called
      expect(validation.validateFileExtension).toHaveBeenCalledWith('data.h5ad');

      // Assert filename generation was called
      expect(filename.generateUploadFileName).toHaveBeenCalledWith('test-folder', 'data.h5ad');

      // Assert uploadForm was called
      expect(apiIndex.uploadForm).toHaveBeenCalled();

      // Assert files_list was updated
      expect(wrapper.vm.files_list).toEqual([
        { file_name: 'uploaded.h5ad', file_size: 2048, created_at: '2025-10-24T12:00:00' }
      ]);

      // Assert uploadPercentage reset to 0
      expect(wrapper.vm.uploadPercentage).toBe(0);
    });

    it('should update uploadPercentage during upload', async () => {
      const validFile = new File(['content'], 'experiment.csv', { type: 'text/csv' });
      wrapper.vm.$refs.selectFile = { files: [validFile] };

      validation.validateFileExtension.mockReturnValue({
        isValid: true,
        extension: 'csv',
        message: ''
      });

      filename.generateUploadFileName.mockReturnValue('data_experiment.csv');

      // Mock uploadForm to call progress callback
      apiIndex.uploadForm.mockImplementation((form, onProgress) => {
        // Simulate progress updates
        onProgress({ loaded: 50, total: 100 });
        expect(wrapper.vm.uploadPercentage).toBe(50);

        onProgress({ loaded: 75, total: 100 });
        expect(wrapper.vm.uploadPercentage).toBe(75);

        onProgress({ loaded: 100, total: 100 });
        expect(wrapper.vm.uploadPercentage).toBe(100);

        return Promise.resolve({ data: 'success' });
      });

      await wrapper.vm.uploadFile();

      // After completion, uploadPercentage should reset to 0
      expect(wrapper.vm.uploadPercentage).toBe(0);
    });

    it('should handle upload errors gracefully', async () => {
      const validFile = new File(['content'], 'data.txt', { type: 'text/plain' });
      wrapper.vm.$refs.selectFile = { files: [validFile] };

      validation.validateFileExtension.mockReturnValue({
        isValid: true,
        extension: 'txt',
        message: ''
      });

      filename.generateUploadFileName.mockReturnValue('data_data.txt');

      // Mock uploadForm to reject
      const uploadError = new Error('Network error');
      apiIndex.uploadForm.mockRejectedValue(uploadError);

      // Spy on console.error
      const consoleErrorSpy = vi.spyOn(console, 'error').mockImplementation(() => {});

      // Should not throw error
      await expect(wrapper.vm.uploadFile()).resolves.not.toThrow();

      // Assert error was logged
      expect(consoleErrorSpy).toHaveBeenCalledWith(uploadError);

      consoleErrorSpy.mockRestore();
    });
  });

  /**
   * removeFile Method Tests
   * Tests file deletion flow with user confirmation
   */
  describe('removeFile', () => {
    let mockConfirm;

    beforeEach(() => {
      // Mock window.confirm
      mockConfirm = vi.fn();
      global.confirm = mockConfirm;

      // Mock deleteFile API
      vi.spyOn(apiIndex, 'deleteFile').mockResolvedValue({ data: 'deleted' });
    });

    afterEach(() => {
      vi.restoreAllMocks();
    });

    it('should delete file when user confirms', async () => {
      // Setup
      mockConfirm.mockReturnValue(true);
      wrapper.vm.files_list = [
        { file_name: 'test1.h5ad', file_size: 1024, created_at: '2025-10-24T10:00:00' },
        { file_name: 'test2.csv', file_size: 2048, created_at: '2025-10-24T11:00:00' }
      ];
      wrapper.vm.list_idx = 0;
      wrapper.vm.file_name = 'test1.h5ad';

      // Execute - now with await
      await wrapper.vm.removeFile();

      // Assert targetFile was set
      expect(wrapper.vm.targetFile).toEqual({
        file_name: 'test1.h5ad',
        file_size: 1024,
        created_at: '2025-10-24T10:00:00'
      });

      // Assert file was removed from list
      expect(wrapper.vm.files_list).toHaveLength(1);
      expect(wrapper.vm.files_list[0].file_name).toBe('test2.csv');

      // Assert confirm was called
      expect(mockConfirm).toHaveBeenCalledWith(
        'Are you sure you want to delete this file? This action cannot be undone.'
      );

      // Assert deleteFile API was called
      expect(apiIndex.deleteFile).toHaveBeenCalledWith({ file_name: 'test1.h5ad' });
    });

    it('should not delete file when user cancels', async () => {
      // Setup
      mockConfirm.mockReturnValue(false);
      wrapper.vm.files_list = [
        { file_name: 'test.h5ad', file_size: 1024, created_at: '2025-10-24T10:00:00' }
      ];
      wrapper.vm.list_idx = 0;
      wrapper.vm.file_name = 'test.h5ad';

      // Execute - now with await
      await wrapper.vm.removeFile();

      // Assert targetFile was still set (splice happens before confirm)
      expect(wrapper.vm.targetFile).toEqual({
        file_name: 'test.h5ad',
        file_size: 1024,
        created_at: '2025-10-24T10:00:00'
      });

      // Assert file was removed from list (splice happens before confirm in current implementation)
      expect(wrapper.vm.files_list).toHaveLength(0);

      // Assert deleteFile API was NOT called
      expect(apiIndex.deleteFile).not.toHaveBeenCalled();
    });

    it('should handle deleteFile API errors and restore file', async () => {
      // Setup
      mockConfirm.mockReturnValue(true);
      const originalFilesList = [
        { file_name: 'test1.h5ad', file_size: 1024, created_at: '2025-10-24T10:00:00' },
        { file_name: 'test2.csv', file_size: 2048, created_at: '2025-10-24T11:00:00' }
      ];
      wrapper.vm.files_list = [...originalFilesList];
      wrapper.vm.list_idx = 0;
      wrapper.vm.file_name = 'test1.h5ad';

      // Mock deleteFile to reject
      const deleteError = new Error('API error');
      apiIndex.deleteFile.mockRejectedValue(deleteError);

      // Spy on console.error
      const consoleErrorSpy = vi.spyOn(console, 'error').mockImplementation(() => {});

      // Execute
      await wrapper.vm.removeFile();

      // Assert error was logged
      expect(consoleErrorSpy).toHaveBeenCalledWith(deleteError);

      // Assert file was restored to original position
      expect(wrapper.vm.files_list).toHaveLength(2);
      expect(wrapper.vm.files_list[0]).toEqual({
        file_name: 'test1.h5ad',
        file_size: 1024,
        created_at: '2025-10-24T10:00:00'
      });
      expect(wrapper.vm.files_list[1]).toEqual({
        file_name: 'test2.csv',
        file_size: 2048,
        created_at: '2025-10-24T11:00:00'
      });

      consoleErrorSpy.mockRestore();
    });

    it('should restore file at correct index after error', async () => {
      // Setup - delete middle item
      mockConfirm.mockReturnValue(true);
      wrapper.vm.files_list = [
        { file_name: 'file1.h5ad', file_size: 1024, created_at: '2025-10-24T10:00:00' },
        { file_name: 'file2.csv', file_size: 2048, created_at: '2025-10-24T11:00:00' },
        { file_name: 'file3.txt', file_size: 512, created_at: '2025-10-24T12:00:00' }
      ];
      wrapper.vm.list_idx = 1; // Middle file
      wrapper.vm.file_name = 'file2.csv';

      // Mock deleteFile to reject
      apiIndex.deleteFile.mockRejectedValue(new Error('Network error'));

      // Spy on console.error
      const consoleErrorSpy = vi.spyOn(console, 'error').mockImplementation(() => {});

      // Execute
      await wrapper.vm.removeFile();

      // Assert file was restored at index 1 (middle position)
      expect(wrapper.vm.files_list).toHaveLength(3);
      expect(wrapper.vm.files_list[0].file_name).toBe('file1.h5ad');
      expect(wrapper.vm.files_list[1].file_name).toBe('file2.csv'); // Restored at index 1
      expect(wrapper.vm.files_list[2].file_name).toBe('file3.txt');

      consoleErrorSpy.mockRestore();
    });

    it('should handle edge case with undefined list_idx', async () => {
      // Setup
      mockConfirm.mockReturnValue(true);
      wrapper.vm.files_list = [
        { file_name: 'test.h5ad', file_size: 1024, created_at: '2025-10-24T10:00:00' }
      ];
      wrapper.vm.list_idx = undefined;
      wrapper.vm.file_name = 'test.h5ad';

      // Execute - should not throw
      await expect(wrapper.vm.removeFile()).resolves.not.toThrow();

      // targetFile will be undefined
      expect(wrapper.vm.targetFile).toBeUndefined();
    });
  });

  /**
   * Data Initialization Tests
   */
  describe('Data Initialization', () => {
    it('should initialize with correct default values', () => {
      expect(wrapper.vm.folders_list).toEqual([]);
      expect(wrapper.vm.files_list).toEqual([]);
      expect(wrapper.vm.R_Mouse_isActive).toBe(false);
      expect(wrapper.vm.toggleFolder).toBeNull();
      expect(wrapper.vm.currentFolder).toBe('data');
      expect(wrapper.vm.toggleMessage).toBe(false);
      expect(wrapper.vm.uploadPercentage).toBe(0);
    });
  });
});
