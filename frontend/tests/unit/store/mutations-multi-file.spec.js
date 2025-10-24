/**
 * workflow/mutations.js 유닛 테스트 - 다중 파일 Mutations
 *
 * 테스트 범위:
 * - setWorkflowFiles
 * - setWorkflowSelectedFiles
 * - toggleWorkflowFileSelection
 * - removeWorkflowFiles
 */

import { describe, it, expect, beforeEach, vi } from 'vitest';
import mutations from '@/store/workflow/mutations';

describe('workflow/mutations - Multi File', () => {
  let mockState;

  beforeEach(() => {
    // Mock state 준비
    mockState = {
      workflow_info: {
        drawflow: {
          Home: {
            data: {
              // Node with multi-file format
              '1': {
                id: 1,
                name: 'InputFile',
                data: {
                  title: 'Multi File Node',
                  files: [
                    { name: 'file1.h5ad', selected: true, size: 1024 },
                    { name: 'file2.h5ad', selected: false, size: 2048 },
                    { name: 'file3.h5ad', selected: true, size: 3072 }
                  ]
                },
                inputs: {},
                outputs: {}
              },
              // Node with single file format (to be converted)
              '2': {
                id: 2,
                name: 'DataTable',
                data: {
                  title: 'Single File Node',
                  file: 'existing-file.h5ad'
                },
                inputs: {},
                outputs: {}
              },
              // Node without files
              '3': {
                id: 3,
                name: 'Visualization',
                data: {
                  title: 'Viz Node'
                },
                inputs: {},
                outputs: {}
              }
            }
          }
        }
      }
    };

    // console.error spy 설정
    vi.spyOn(console, 'error').mockImplementation(() => {});
  });

  describe('setWorkflowFiles', () => {
    it('should set multi-file array format', () => {
      const payload = {
        id: '3',
        files: [
          { name: 'new1.h5ad', selected: true, size: 500 },
          { name: 'new2.h5ad', selected: false, size: 600 }
        ]
      };

      mutations.setWorkflowFiles(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['3'].data;
      expect(nodeData.files).toEqual(payload.files);
      expect(nodeData.files).toHaveLength(2);
    });

    it('should clear old single file format', () => {
      const payload = {
        id: '2',
        files: [
          { name: 'converted.h5ad', selected: true, size: 1000 }
        ]
      };

      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBe('existing-file.h5ad');

      mutations.setWorkflowFiles(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['2'].data;
      expect(nodeData.files).toEqual(payload.files);
      expect(nodeData.file).toBeUndefined();
    });

    it('should replace existing multi-file array', () => {
      const payload = {
        id: '1',
        files: [
          { name: 'replaced1.h5ad', selected: false, size: 100 },
          { name: 'replaced2.h5ad', selected: false, size: 200 }
        ]
      };

      mutations.setWorkflowFiles(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.files).toEqual(payload.files);
      expect(nodeData.files).toHaveLength(2);
      expect(nodeData.files.some(f => f.name === 'file1.h5ad')).toBe(false);
    });

    it('should handle empty files array', () => {
      const payload = {
        id: '1',
        files: []
      };

      mutations.setWorkflowFiles(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.files).toEqual([]);
      expect(nodeData.files).toHaveLength(0);
    });

    it('should handle files without size property', () => {
      const payload = {
        id: '3',
        files: [
          { name: 'file1.h5ad', selected: true },
          { name: 'file2.h5ad', selected: false }
        ]
      };

      mutations.setWorkflowFiles(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['3'].data;
      expect(nodeData.files).toEqual(payload.files);
    });

    it('should handle large file arrays', () => {
      const largeFileArray = Array.from({ length: 100 }, (_, i) => ({
        name: `file${i}.h5ad`,
        selected: i % 2 === 0,
        size: i * 1000
      }));

      const payload = {
        id: '1',
        files: largeFileArray
      };

      mutations.setWorkflowFiles(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.files).toHaveLength(100);
      expect(nodeData.files[0].name).toBe('file0.h5ad');
      expect(nodeData.files[99].name).toBe('file99.h5ad');
    });

    it('should log error for invalid node ID', () => {
      const payload = {
        id: '999',
        files: [{ name: 'test.h5ad', selected: true }]
      };

      mutations.setWorkflowFiles(mockState, payload);

      expect(console.error).toHaveBeenCalledWith(
        'No object found with id: 999'
      );
    });

    it('should not crash for non-existent node', () => {
      const payload = {
        id: 'nonexistent',
        files: []
      };

      expect(() => {
        mutations.setWorkflowFiles(mockState, payload);
      }).not.toThrow();
    });
  });

  describe('setWorkflowSelectedFiles', () => {
    it('should update selected status for multiple files', () => {
      const payload = {
        id: '1',
        selectedFiles: ['file1.h5ad', 'file3.h5ad']
      };

      mutations.setWorkflowSelectedFiles(mockState, payload);

      const files = mockState.workflow_info.drawflow.Home.data['1'].data.files;
      expect(files[0].selected).toBe(true);  // file1.h5ad
      expect(files[1].selected).toBe(false); // file2.h5ad
      expect(files[2].selected).toBe(true);  // file3.h5ad
    });

    it('should deselect all when empty array provided', () => {
      const payload = {
        id: '1',
        selectedFiles: []
      };

      mutations.setWorkflowSelectedFiles(mockState, payload);

      const files = mockState.workflow_info.drawflow.Home.data['1'].data.files;
      expect(files.every(f => !f.selected)).toBe(true);
    });

    it('should select all files', () => {
      const payload = {
        id: '1',
        selectedFiles: ['file1.h5ad', 'file2.h5ad', 'file3.h5ad']
      };

      mutations.setWorkflowSelectedFiles(mockState, payload);

      const files = mockState.workflow_info.drawflow.Home.data['1'].data.files;
      expect(files.every(f => f.selected)).toBe(true);
    });

    it('should handle non-existent file names gracefully', () => {
      const payload = {
        id: '1',
        selectedFiles: ['file1.h5ad', 'nonexistent.h5ad', 'file3.h5ad']
      };

      mutations.setWorkflowSelectedFiles(mockState, payload);

      const files = mockState.workflow_info.drawflow.Home.data['1'].data.files;
      expect(files[0].selected).toBe(true);  // file1.h5ad
      expect(files[1].selected).toBe(false); // file2.h5ad
      expect(files[2].selected).toBe(true);  // file3.h5ad
    });

    it('should log error for non-multi-file nodes', () => {
      const payload = {
        id: '3',
        selectedFiles: ['test.h5ad']
      };

      mutations.setWorkflowSelectedFiles(mockState, payload);

      expect(console.error).toHaveBeenCalledWith(
        'No multi-file node found with id: 3'
      );
    });

    it('should handle node with single file format', () => {
      const payload = {
        id: '2',
        selectedFiles: ['existing-file.h5ad']
      };

      mutations.setWorkflowSelectedFiles(mockState, payload);

      expect(console.error).toHaveBeenCalledWith(
        'No multi-file node found with id: 2'
      );
    });

    it('should handle invalid node ID', () => {
      const payload = {
        id: '999',
        selectedFiles: ['test.h5ad']
      };

      mutations.setWorkflowSelectedFiles(mockState, payload);

      expect(console.error).toHaveBeenCalled();
    });
  });

  describe('toggleWorkflowFileSelection', () => {
    it('should toggle file selection state from true to false', () => {
      const payload = {
        id: '1',
        fileName: 'file1.h5ad'
      };

      const files = mockState.workflow_info.drawflow.Home.data['1'].data.files;
      expect(files[0].selected).toBe(true);

      mutations.toggleWorkflowFileSelection(mockState, payload);

      expect(files[0].selected).toBe(false);
    });

    it('should toggle file selection state from false to true', () => {
      const payload = {
        id: '1',
        fileName: 'file2.h5ad'
      };

      const files = mockState.workflow_info.drawflow.Home.data['1'].data.files;
      expect(files[1].selected).toBe(false);

      mutations.toggleWorkflowFileSelection(mockState, payload);

      expect(files[1].selected).toBe(true);
    });

    it('should toggle multiple times', () => {
      const payload = {
        id: '1',
        fileName: 'file1.h5ad'
      };

      const files = mockState.workflow_info.drawflow.Home.data['1'].data.files;
      const originalState = files[0].selected;

      mutations.toggleWorkflowFileSelection(mockState, payload);
      expect(files[0].selected).toBe(!originalState);

      mutations.toggleWorkflowFileSelection(mockState, payload);
      expect(files[0].selected).toBe(originalState);

      mutations.toggleWorkflowFileSelection(mockState, payload);
      expect(files[0].selected).toBe(!originalState);
    });

    it('should handle non-existent file name gracefully', () => {
      const payload = {
        id: '1',
        fileName: 'nonexistent.h5ad'
      };

      expect(() => {
        mutations.toggleWorkflowFileSelection(mockState, payload);
      }).not.toThrow();

      // Original files should remain unchanged
      const files = mockState.workflow_info.drawflow.Home.data['1'].data.files;
      expect(files[0].selected).toBe(true);
      expect(files[1].selected).toBe(false);
    });

    it('should log error for non-multi-file node', () => {
      const payload = {
        id: '3',
        fileName: 'test.h5ad'
      };

      mutations.toggleWorkflowFileSelection(mockState, payload);

      expect(console.error).toHaveBeenCalledWith(
        'No multi-file node found with id: 3'
      );
    });

    it('should log error for invalid node ID', () => {
      const payload = {
        id: '999',
        fileName: 'test.h5ad'
      };

      mutations.toggleWorkflowFileSelection(mockState, payload);

      expect(console.error).toHaveBeenCalled();
    });
  });

  describe('removeWorkflowFiles', () => {
    it('should clear multi-file array', () => {
      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.files).toHaveLength(3);

      mutations.removeWorkflowFiles(mockState, '1');

      expect(nodeData.files).toEqual([]);
      expect(nodeData.files).toHaveLength(0);
    });

    it('should clear single file format for complete cleanup', () => {
      const nodeData = mockState.workflow_info.drawflow.Home.data['2'].data;
      nodeData.files = [{ name: 'temp.h5ad', selected: true }];
      nodeData.file = 'backup.h5ad';

      mutations.removeWorkflowFiles(mockState, '2');

      expect(nodeData.files).toEqual([]);
      expect(nodeData.file).toBeNull();
    });

    it('should handle node without files', () => {
      expect(() => {
        mutations.removeWorkflowFiles(mockState, '3');
      }).not.toThrow();

      const nodeData = mockState.workflow_info.drawflow.Home.data['3'].data;
      expect(nodeData.files).toBeUndefined();
    });

    it('should handle node with already empty files array', () => {
      mockState.workflow_info.drawflow.Home.data['1'].data.files = [];

      mutations.removeWorkflowFiles(mockState, '1');

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.files).toEqual([]);
    });

    it('should handle invalid node ID gracefully', () => {
      expect(() => {
        mutations.removeWorkflowFiles(mockState, '999');
      }).not.toThrow();

      expect(console.error).toHaveBeenCalled();
    });

    it('should not affect other node properties', () => {
      const originalTitle = mockState.workflow_info.drawflow.Home.data['1'].data.title;

      mutations.removeWorkflowFiles(mockState, '1');

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.title).toBe(originalTitle);
    });
  });

  describe('Combined scenarios', () => {
    it('should handle setWorkflowFiles → setWorkflowSelectedFiles', () => {
      // First set files
      mutations.setWorkflowFiles(mockState, {
        id: '3',
        files: [
          { name: 'a.h5ad', selected: false, size: 100 },
          { name: 'b.h5ad', selected: false, size: 200 },
          { name: 'c.h5ad', selected: false, size: 300 }
        ]
      });

      // Then select specific files
      mutations.setWorkflowSelectedFiles(mockState, {
        id: '3',
        selectedFiles: ['a.h5ad', 'c.h5ad']
      });

      const files = mockState.workflow_info.drawflow.Home.data['3'].data.files;
      expect(files[0].selected).toBe(true);  // a.h5ad
      expect(files[1].selected).toBe(false); // b.h5ad
      expect(files[2].selected).toBe(true);  // c.h5ad
    });

    it('should handle setWorkflowFiles → toggleWorkflowFileSelection', () => {
      // Set files
      mutations.setWorkflowFiles(mockState, {
        id: '3',
        files: [
          { name: 'x.h5ad', selected: true },
          { name: 'y.h5ad', selected: false }
        ]
      });

      // Toggle selection
      mutations.toggleWorkflowFileSelection(mockState, {
        id: '3',
        fileName: 'x.h5ad'
      });

      mutations.toggleWorkflowFileSelection(mockState, {
        id: '3',
        fileName: 'y.h5ad'
      });

      const files = mockState.workflow_info.drawflow.Home.data['3'].data.files;
      expect(files[0].selected).toBe(false); // x.h5ad toggled
      expect(files[1].selected).toBe(true);  // y.h5ad toggled
    });

    it('should handle full workflow: set → select → toggle → remove', () => {
      const nodeId = '3';

      // 1. Set files
      mutations.setWorkflowFiles(mockState, {
        id: nodeId,
        files: [
          { name: 'file1.h5ad', selected: false, size: 100 },
          { name: 'file2.h5ad', selected: false, size: 200 }
        ]
      });

      let files = mockState.workflow_info.drawflow.Home.data[nodeId].data.files;
      expect(files).toHaveLength(2);

      // 2. Select files
      mutations.setWorkflowSelectedFiles(mockState, {
        id: nodeId,
        selectedFiles: ['file1.h5ad']
      });

      expect(files[0].selected).toBe(true);
      expect(files[1].selected).toBe(false);

      // 3. Toggle selection
      mutations.toggleWorkflowFileSelection(mockState, {
        id: nodeId,
        fileName: 'file2.h5ad'
      });

      expect(files[0].selected).toBe(true);
      expect(files[1].selected).toBe(true);

      // 4. Remove all
      mutations.removeWorkflowFiles(mockState, nodeId);

      // Re-fetch files reference after mutation (mutation assigns new array)
      files = mockState.workflow_info.drawflow.Home.data[nodeId].data.files;
      expect(files).toEqual([]);
    });

    it('should handle converting single file to multi file', () => {
      const nodeId = '2';
      const nodeData = mockState.workflow_info.drawflow.Home.data[nodeId].data;

      // Initially single file
      expect(nodeData.file).toBe('existing-file.h5ad');

      // Convert to multi-file
      mutations.setWorkflowFiles(mockState, {
        id: nodeId,
        files: [
          { name: 'converted1.h5ad', selected: true },
          { name: 'converted2.h5ad', selected: false }
        ]
      });

      expect(nodeData.file).toBeUndefined();
      expect(nodeData.files).toHaveLength(2);

      // Select files
      mutations.setWorkflowSelectedFiles(mockState, {
        id: nodeId,
        selectedFiles: ['converted1.h5ad', 'converted2.h5ad']
      });

      expect(nodeData.files.every(f => f.selected)).toBe(true);
    });
  });
});
