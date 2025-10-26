/**
 * workflow/getters.js 유닛 테스트 - File Query Getters
 *
 * 테스트 범위:
 * - getWorkflowNodeFileInfo
 * - getWorkflowNodeFilesInfo
 * - getSelectedWorkflowFiles
 * - getSelectedWorkflowFileNames
 * - getAlgorithmIdByTaskId
 */

import { describe, it, expect, beforeEach } from 'vitest';
import getters from '@/store/workflow/getters';

describe('workflow/getters - File Query', () => {
  let mockState;

  beforeEach(() => {
    // Mock state 준비
    mockState = {
      workflow_info: {
        drawflow: {
          Home: {
            data: {
              // Multi-file node with selections
              '1': {
                id: 1,
                name: 'InputFile',
                data: {
                  title: 'Multi-File Input',
                  files: [
                    { name: 'file1.h5ad', selected: true, size: 1024 },
                    { name: 'file2.h5ad', selected: false, size: 2048 },
                    { name: 'file3.h5ad', selected: true, size: 512 }
                  ]
                }
              },
              // Single-file node
              '2': {
                id: 2,
                name: 'DataTable',
                data: {
                  title: 'Single File',
                  file: 'single-data.h5ad'
                }
              },
              // Algorithm node with files object
              '3': {
                id: 3,
                name: 'Algorithm',
                class: 'Algorithm',
                data: {
                  title: 'Algorithm Node',
                  files: {
                    '1': 'file1.h5ad',
                    '2': 'file2.h5ad'
                  }
                }
              },
              // Node without files
              '4': {
                id: 4,
                name: 'Visualization',
                data: {
                  title: 'Viz Node'
                }
              },
              // Multi-file node with single selection
              '5': {
                id: 5,
                name: 'InputFile',
                data: {
                  files: [
                    { name: 'only-selected.h5ad', selected: true, size: 256 }
                  ]
                }
              },
              // Multi-file node with no selection
              '6': {
                id: 6,
                name: 'InputFile',
                data: {
                  files: [
                    { name: 'unselected1.h5ad', selected: false, size: 100 },
                    { name: 'unselected2.h5ad', selected: false, size: 200 }
                  ]
                }
              },
              // Empty multi-file array
              '7': {
                id: 7,
                name: 'InputFile',
                data: {
                  files: []
                }
              }
            }
          }
        }
      },
      taskAlgorithmMapping: {
        'task-123': '1',
        'task-456': '2',
        'task-789': '3'
      }
    };
  });

  describe('getWorkflowNodeFileInfo', () => {
    it('should return first selected file name for multi-file node', () => {
      const result = getters.getWorkflowNodeFileInfo(mockState);
      expect(result('1')).toBe('file1.h5ad');
    });

    it('should return single file name for single-file node', () => {
      const result = getters.getWorkflowNodeFileInfo(mockState);
      expect(result('2')).toBe('single-data.h5ad');
    });

    it('should return null for node without files', () => {
      const result = getters.getWorkflowNodeFileInfo(mockState);
      expect(result('4')).toBeNull();
    });

    it('should return null for multi-file node with no selection', () => {
      const result = getters.getWorkflowNodeFileInfo(mockState);
      expect(result('6')).toBeNull();
    });

    it('should return null for empty multi-file array', () => {
      const result = getters.getWorkflowNodeFileInfo(mockState);
      expect(result('7')).toBeNull();
    });

    it('should return only selected file name for single selection', () => {
      const result = getters.getWorkflowNodeFileInfo(mockState);
      expect(result('5')).toBe('only-selected.h5ad');
    });

    it('should return first selected file when multiple are selected', () => {
      const result = getters.getWorkflowNodeFileInfo(mockState);
      const fileName = result('1');
      expect(fileName).toBe('file1.h5ad');
    });

    it('should return null for non-existent node', () => {
      const result = getters.getWorkflowNodeFileInfo(mockState);
      expect(result('999')).toBeNull();
    });

    it('should handle node without data', () => {
      mockState.workflow_info.drawflow.Home.data['8'] = {
        id: 8,
        name: 'Test',
        data: null
      };
      const result = getters.getWorkflowNodeFileInfo(mockState);
      expect(result('8')).toBeNull();
    });
  });

  describe('getWorkflowNodeFilesInfo', () => {
    it('should return array of files for multi-file node', () => {
      const result = getters.getWorkflowNodeFilesInfo(mockState);
      const files = result('1');
      expect(files).toHaveLength(3);
      expect(files[0]).toEqual({ name: 'file1.h5ad', selected: true, size: 1024 });
    });

    it('should convert single file to array format', () => {
      const result = getters.getWorkflowNodeFilesInfo(mockState);
      const files = result('2');
      expect(files).toEqual([{
        name: 'single-data.h5ad',
        selected: true,
        size: 0
      }]);
    });

    it('should return algorithm files object', () => {
      const result = getters.getWorkflowNodeFilesInfo(mockState);
      const files = result('3');
      expect(files).toEqual({
        '1': 'file1.h5ad',
        '2': 'file2.h5ad'
      });
    });

    it('should return null for node without files', () => {
      const result = getters.getWorkflowNodeFilesInfo(mockState);
      expect(result('4')).toBeNull();
    });

    it('should return empty array for empty multi-file node', () => {
      const result = getters.getWorkflowNodeFilesInfo(mockState);
      const files = result('7');
      expect(files).toEqual([]);
    });

    it('should return null for non-existent node', () => {
      const result = getters.getWorkflowNodeFilesInfo(mockState);
      expect(result('999')).toBeNull();
    });

    it('should handle node without data', () => {
      mockState.workflow_info.drawflow.Home.data['9'] = {
        id: 9,
        name: 'Test',
        data: null
      };
      const result = getters.getWorkflowNodeFilesInfo(mockState);
      expect(result('9')).toBeNull();
    });

    it('should preserve file metadata in multi-file format', () => {
      const result = getters.getWorkflowNodeFilesInfo(mockState);
      const files = result('1');
      expect(files[0]).toHaveProperty('name');
      expect(files[0]).toHaveProperty('selected');
      expect(files[0]).toHaveProperty('size');
    });

    it('should return unmodified files array including unselected', () => {
      const result = getters.getWorkflowNodeFilesInfo(mockState);
      const files = result('6');
      expect(files).toHaveLength(2);
      expect(files.every(f => !f.selected)).toBe(true);
    });
  });

  describe('getSelectedWorkflowFiles', () => {
    it('should return only selected files from multi-file node', () => {
      const result = getters.getSelectedWorkflowFiles(mockState);
      const selected = result('1');
      expect(selected).toHaveLength(2);
      expect(selected[0].name).toBe('file1.h5ad');
      expect(selected[1].name).toBe('file3.h5ad');
    });

    it('should convert single file to array format with selected=true', () => {
      const result = getters.getSelectedWorkflowFiles(mockState);
      const selected = result('2');
      expect(selected).toEqual([{
        name: 'single-data.h5ad',
        selected: true,
        size: 0
      }]);
    });

    it('should return empty array for node without files', () => {
      const result = getters.getSelectedWorkflowFiles(mockState);
      expect(result('4')).toEqual([]);
    });

    it('should return empty array when no files are selected', () => {
      const result = getters.getSelectedWorkflowFiles(mockState);
      expect(result('6')).toEqual([]);
    });

    it('should return empty array for empty multi-file node', () => {
      const result = getters.getSelectedWorkflowFiles(mockState);
      expect(result('7')).toEqual([]);
    });

    it('should return single selected file in array', () => {
      const result = getters.getSelectedWorkflowFiles(mockState);
      const selected = result('5');
      expect(selected).toHaveLength(1);
      expect(selected[0].name).toBe('only-selected.h5ad');
    });

    it('should return empty array for non-existent node', () => {
      const result = getters.getSelectedWorkflowFiles(mockState);
      expect(result('999')).toEqual([]);
    });

    it('should handle node without data', () => {
      mockState.workflow_info.drawflow.Home.data['10'] = {
        id: 10,
        name: 'Test',
        data: null
      };
      const result = getters.getSelectedWorkflowFiles(mockState);
      expect(result('10')).toEqual([]);
    });

    it('should preserve file metadata for selected files', () => {
      const result = getters.getSelectedWorkflowFiles(mockState);
      const selected = result('1');
      expect(selected[0]).toHaveProperty('name');
      expect(selected[0]).toHaveProperty('selected');
      expect(selected[0]).toHaveProperty('size');
      expect(selected[0].selected).toBe(true);
    });

    it('should filter out unselected files', () => {
      const result = getters.getSelectedWorkflowFiles(mockState);
      const selected = result('1');
      const hasUnselected = selected.some(f => !f.selected);
      expect(hasUnselected).toBe(false);
    });
  });

  describe('getSelectedWorkflowFileNames', () => {
    it('should return array of selected file names', () => {
      const result = getters.getSelectedWorkflowFileNames(mockState);
      const names = result('1');
      expect(names).toEqual(['file1.h5ad', 'file3.h5ad']);
    });

    it('should return single file name in array for single-file node', () => {
      const result = getters.getSelectedWorkflowFileNames(mockState);
      const names = result('2');
      expect(names).toEqual(['single-data.h5ad']);
    });

    it('should return empty array for node without files', () => {
      const result = getters.getSelectedWorkflowFileNames(mockState);
      expect(result('4')).toEqual([]);
    });

    it('should return empty array when no files are selected', () => {
      const result = getters.getSelectedWorkflowFileNames(mockState);
      expect(result('6')).toEqual([]);
    });

    it('should return empty array for empty multi-file node', () => {
      const result = getters.getSelectedWorkflowFileNames(mockState);
      expect(result('7')).toEqual([]);
    });

    it('should return single name for single selection', () => {
      const result = getters.getSelectedWorkflowFileNames(mockState);
      const names = result('5');
      expect(names).toEqual(['only-selected.h5ad']);
    });

    it('should return empty array for non-existent node', () => {
      const result = getters.getSelectedWorkflowFileNames(mockState);
      expect(result('999')).toEqual([]);
    });

    it('should handle node without data', () => {
      mockState.workflow_info.drawflow.Home.data['11'] = {
        id: 11,
        name: 'Test',
        data: null
      };
      const result = getters.getSelectedWorkflowFileNames(mockState);
      expect(result('11')).toEqual([]);
    });

    it('should return only names, no metadata', () => {
      const result = getters.getSelectedWorkflowFileNames(mockState);
      const names = result('1');
      expect(names.every(item => typeof item === 'string')).toBe(true);
    });

    it('should maintain order of selected files', () => {
      const result = getters.getSelectedWorkflowFileNames(mockState);
      const names = result('1');
      expect(names[0]).toBe('file1.h5ad');
      expect(names[1]).toBe('file3.h5ad');
    });
  });

  describe('getAlgorithmIdByTaskId', () => {
    it('should return algorithm ID for valid task ID', () => {
      const result = getters.getAlgorithmIdByTaskId(mockState);
      expect(result('task-123')).toBe('1');
      expect(result('task-456')).toBe('2');
      expect(result('task-789')).toBe('3');
    });

    it('should return null for non-existent task ID', () => {
      const result = getters.getAlgorithmIdByTaskId(mockState);
      expect(result('task-999')).toBeNull();
    });

    it('should handle empty string task ID', () => {
      const result = getters.getAlgorithmIdByTaskId(mockState);
      expect(result('')).toBeNull();
    });

    it('should handle null task ID', () => {
      const result = getters.getAlgorithmIdByTaskId(mockState);
      expect(result(null)).toBeNull();
    });

    it('should handle undefined task ID', () => {
      const result = getters.getAlgorithmIdByTaskId(mockState);
      expect(result(undefined)).toBeNull();
    });

    it('should return correct mapping for multiple tasks', () => {
      const result = getters.getAlgorithmIdByTaskId(mockState);
      const mappings = [
        { taskId: 'task-123', expected: '1' },
        { taskId: 'task-456', expected: '2' },
        { taskId: 'task-789', expected: '3' }
      ];

      mappings.forEach(({ taskId, expected }) => {
        expect(result(taskId)).toBe(expected);
      });
    });

    it('should handle task ID with special characters', () => {
      mockState.taskAlgorithmMapping['task-abc-123-xyz'] = '10';
      const result = getters.getAlgorithmIdByTaskId(mockState);
      expect(result('task-abc-123-xyz')).toBe('10');
    });

    it('should handle numeric task IDs as strings', () => {
      mockState.taskAlgorithmMapping['12345'] = '99';
      const result = getters.getAlgorithmIdByTaskId(mockState);
      expect(result('12345')).toBe('99');
    });
  });

  describe('Combined scenarios', () => {
    it('should consistently identify files across different getters', () => {
      const fileInfo = getters.getWorkflowNodeFileInfo(mockState);
      const filesInfo = getters.getWorkflowNodeFilesInfo(mockState);
      const selectedFiles = getters.getSelectedWorkflowFiles(mockState);
      const selectedNames = getters.getSelectedWorkflowFileNames(mockState);

      // Node 1: multi-file
      expect(fileInfo('1')).toBe('file1.h5ad');
      expect(filesInfo('1')).toHaveLength(3);
      expect(selectedFiles('1')).toHaveLength(2);
      expect(selectedNames('1')).toEqual(['file1.h5ad', 'file3.h5ad']);

      // Node 2: single-file
      expect(fileInfo('2')).toBe('single-data.h5ad');
      expect(filesInfo('2')).toHaveLength(1);
      expect(selectedFiles('2')).toHaveLength(1);
      expect(selectedNames('2')).toEqual(['single-data.h5ad']);
    });

    it('should handle file format conversions correctly', () => {
      const filesInfo = getters.getWorkflowNodeFilesInfo(mockState);
      const selectedFiles = getters.getSelectedWorkflowFiles(mockState);

      // Single file converted to array format
      const singleAsArray = filesInfo('2');
      const selectedSingle = selectedFiles('2');

      expect(singleAsArray).toEqual(selectedSingle);
      expect(singleAsArray[0].selected).toBe(true);
    });

    it('should handle task-algorithm mapping integration', () => {
      const algorithmId = getters.getAlgorithmIdByTaskId(mockState);
      const nodeFiles = getters.getWorkflowNodeFilesInfo(mockState);

      // Get algorithm node ID from task
      const nodeId = algorithmId('task-789');
      expect(nodeId).toBe('3');

      // Get files from that algorithm node
      const files = nodeFiles(nodeId);
      expect(files).toEqual({
        '1': 'file1.h5ad',
        '2': 'file2.h5ad'
      });
    });

    it('should return empty/null consistently for nodes without files', () => {
      const fileInfo = getters.getWorkflowNodeFileInfo(mockState);
      const filesInfo = getters.getWorkflowNodeFilesInfo(mockState);
      const selectedFiles = getters.getSelectedWorkflowFiles(mockState);
      const selectedNames = getters.getSelectedWorkflowFileNames(mockState);

      expect(fileInfo('4')).toBeNull();
      expect(filesInfo('4')).toBeNull();
      expect(selectedFiles('4')).toEqual([]);
      expect(selectedNames('4')).toEqual([]);
    });
  });
});
