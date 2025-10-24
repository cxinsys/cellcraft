/**
 * workflow/getters.js 유닛 테스트 - Boolean & Format Getters
 *
 * 테스트 범위:
 * - isMultiFileNode
 * - hasWorkflowFiles
 * - hasSelectedWorkflowFiles
 * - getWorkflowFileFormat
 * - isAlgorithmNodeRunning
 */

import { describe, it, expect, beforeEach } from 'vitest';
import getters from '@/store/workflow/getters';

describe('workflow/getters - Boolean & Format', () => {
  let mockState;

  beforeEach(() => {
    // Mock state 준비
    mockState = {
      workflow_info: {
        drawflow: {
          Home: {
            data: {
              // Multi-file node (array format)
              '1': {
                id: 1,
                name: 'InputFile',
                data: {
                  files: [
                    { name: 'file1.h5ad', selected: true, size: 1024 },
                    { name: 'file2.h5ad', selected: false, size: 2048 }
                  ]
                }
              },
              // Single-file node (string format)
              '2': {
                id: 2,
                name: 'DataTable',
                data: {
                  file: 'single-file.h5ad'
                }
              },
              // Algorithm node with files object
              '3': {
                id: 3,
                name: 'Algorithm',
                data: {
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
                data: {}
              },
              // Multi-file node with no selection
              '5': {
                id: 5,
                name: 'InputFile',
                data: {
                  files: [
                    { name: 'file3.h5ad', selected: false, size: 512 },
                    { name: 'file4.h5ad', selected: false, size: 768 }
                  ]
                }
              },
              // Multi-file node with empty array
              '6': {
                id: 6,
                name: 'InputFile',
                data: {
                  files: []
                }
              }
            }
          }
        }
      },
      runningAlgorithmNodes: ['3', '7', '10']
    };
  });

  describe('isMultiFileNode', () => {
    it('should return true for multi-file nodes', () => {
      const result = getters.isMultiFileNode(mockState);
      expect(result('1')).toBe(true);
      expect(result('5')).toBe(true);
    });

    it('should return false for single-file nodes', () => {
      const result = getters.isMultiFileNode(mockState);
      expect(result('2')).toBe(false);
    });

    it('should return false for nodes without files', () => {
      const result = getters.isMultiFileNode(mockState);
      expect(result('4')).toBe(false);
    });

    it('should return false for algorithm files object', () => {
      const result = getters.isMultiFileNode(mockState);
      expect(result('3')).toBe(false);
    });

    it('should return true for empty multi-file array', () => {
      const result = getters.isMultiFileNode(mockState);
      expect(result('6')).toBe(true);
    });

    it('should return falsy for non-existent node', () => {
      const result = getters.isMultiFileNode(mockState);
      expect(result('999')).toBeFalsy();
    });
  });

  describe('hasWorkflowFiles', () => {
    it('should return true if node has multi-files', () => {
      const result = getters.hasWorkflowFiles(mockState);
      expect(result('1')).toBe(true);
      expect(result('5')).toBe(true);
    });

    it('should return true if node has single file', () => {
      const result = getters.hasWorkflowFiles(mockState);
      expect(result('2')).toBe(true);
    });

    it('should return false if no files', () => {
      const result = getters.hasWorkflowFiles(mockState);
      expect(result('4')).toBe(false);
    });

    it('should return false for empty multi-file array', () => {
      const result = getters.hasWorkflowFiles(mockState);
      expect(result('6')).toBe(false);
    });

    it('should return false for non-existent node', () => {
      const result = getters.hasWorkflowFiles(mockState);
      expect(result('999')).toBe(false);
    });

    it('should return false for node without data', () => {
      mockState.workflow_info.drawflow.Home.data['7'] = {
        id: 7,
        name: 'Test',
        data: null
      };
      const result = getters.hasWorkflowFiles(mockState);
      expect(result('7')).toBe(false);
    });
  });

  describe('hasSelectedWorkflowFiles', () => {
    it('should return true if any file is selected', () => {
      const result = getters.hasSelectedWorkflowFiles(mockState);
      expect(result('1')).toBe(true); // has selected file
    });

    it('should return false if no files selected', () => {
      const result = getters.hasSelectedWorkflowFiles(mockState);
      expect(result('5')).toBe(false); // no selected files
    });

    it('should return true for single file format', () => {
      const result = getters.hasSelectedWorkflowFiles(mockState);
      expect(result('2')).toBe(true); // single file is considered selected
    });

    it('should return false for nodes without files', () => {
      const result = getters.hasSelectedWorkflowFiles(mockState);
      expect(result('4')).toBe(false);
    });

    it('should return false for empty file array', () => {
      const result = getters.hasSelectedWorkflowFiles(mockState);
      expect(result('6')).toBe(false);
    });

    it('should return false for non-existent node', () => {
      const result = getters.hasSelectedWorkflowFiles(mockState);
      expect(result('999')).toBe(false);
    });

    it('should handle node without data', () => {
      mockState.workflow_info.drawflow.Home.data['8'] = {
        id: 8,
        name: 'Test',
        data: null
      };
      const result = getters.hasSelectedWorkflowFiles(mockState);
      expect(result('8')).toBe(false);
    });

    it('should handle all files selected in multi-file', () => {
      mockState.workflow_info.drawflow.Home.data['9'] = {
        id: 9,
        name: 'InputFile',
        data: {
          files: [
            { name: 'a.h5ad', selected: true },
            { name: 'b.h5ad', selected: true }
          ]
        }
      };
      const result = getters.hasSelectedWorkflowFiles(mockState);
      expect(result('9')).toBe(true);
    });
  });

  describe('getWorkflowFileFormat', () => {
    it('should return "multi" for array format', () => {
      const result = getters.getWorkflowFileFormat(mockState);
      expect(result('1')).toBe('multi');
      expect(result('5')).toBe('multi');
      expect(result('6')).toBe('multi');
    });

    it('should return "single" for string format', () => {
      const result = getters.getWorkflowFileFormat(mockState);
      expect(result('2')).toBe('single');
    });

    it('should return "algorithm" for object format', () => {
      const result = getters.getWorkflowFileFormat(mockState);
      expect(result('3')).toBe('algorithm');
    });

    it('should return "none" for no files', () => {
      const result = getters.getWorkflowFileFormat(mockState);
      expect(result('4')).toBe('none');
    });

    it('should return "none" for non-existent node', () => {
      const result = getters.getWorkflowFileFormat(mockState);
      expect(result('999')).toBe('none');
    });

    it('should return "none" for node without data', () => {
      mockState.workflow_info.drawflow.Home.data['10'] = {
        id: 10,
        name: 'Test',
        data: null
      };
      const result = getters.getWorkflowFileFormat(mockState);
      expect(result('10')).toBe('none');
    });

    it('should prioritize multi-file array over single file', () => {
      // 이론적으로는 두 형식이 동시에 존재할 수 없지만 검증
      mockState.workflow_info.drawflow.Home.data['11'] = {
        id: 11,
        name: 'Test',
        data: {
          files: [{ name: 'multi.h5ad', selected: true }],
          file: 'single.h5ad'
        }
      };
      const result = getters.getWorkflowFileFormat(mockState);
      expect(result('11')).toBe('multi');
    });

    it('should detect single file format when files is object (not array)', () => {
      mockState.workflow_info.drawflow.Home.data['12'] = {
        id: 12,
        name: 'Algorithm',
        data: {
          files: { '1': 'file.h5ad' },  // Object format (algorithm), not array
          file: 'single-file.h5ad'
        }
      };
      const result = getters.getWorkflowFileFormat(mockState);
      // files object is not an array, so single file format is detected
      expect(result('12')).toBe('single');
    });
  });

  describe('isAlgorithmNodeRunning', () => {
    it('should return true if node ID is in running list', () => {
      const result = getters.isAlgorithmNodeRunning(mockState);
      expect(result('3')).toBe(true);
      expect(result('7')).toBe(true);
      expect(result('10')).toBe(true);
    });

    it('should return false otherwise', () => {
      const result = getters.isAlgorithmNodeRunning(mockState);
      expect(result('1')).toBe(false);
      expect(result('2')).toBe(false);
      expect(result('999')).toBe(false);
    });

    it('should handle string ID conversion', () => {
      const result = getters.isAlgorithmNodeRunning(mockState);
      expect(result('3')).toBe(true);
      expect(result(3)).toBe(true); // number conversion
    });

    it('should return false for number when only string exists', () => {
      mockState.runningAlgorithmNodes = ['100'];
      const result = getters.isAlgorithmNodeRunning(mockState);
      expect(result('100')).toBe(true);
      expect(result(100)).toBe(true); // String() conversion
    });

    it('should handle empty running list', () => {
      mockState.runningAlgorithmNodes = [];
      const result = getters.isAlgorithmNodeRunning(mockState);
      expect(result('3')).toBe(false);
      expect(result('999')).toBe(false);
    });

    it('should handle null/undefined node ID', () => {
      const result = getters.isAlgorithmNodeRunning(mockState);
      expect(result(null)).toBe(false);
      expect(result(undefined)).toBe(false);
    });

    it('should handle special string node IDs', () => {
      mockState.runningAlgorithmNodes = ['node-abc-123', 'test_node_456'];
      const result = getters.isAlgorithmNodeRunning(mockState);
      expect(result('node-abc-123')).toBe(true);
      expect(result('test_node_456')).toBe(true);
      expect(result('other')).toBe(false);
    });
  });

  describe('Combined scenarios', () => {
    it('should handle workflow with mixed file formats', () => {
      const hasFiles = getters.hasWorkflowFiles(mockState);
      const format = getters.getWorkflowFileFormat(mockState);
      const hasSelected = getters.hasSelectedWorkflowFiles(mockState);

      // Node 1: multi-file with selection
      expect(hasFiles('1')).toBe(true);
      expect(format('1')).toBe('multi');
      expect(hasSelected('1')).toBe(true);

      // Node 2: single file
      expect(hasFiles('2')).toBe(true);
      expect(format('2')).toBe('single');
      expect(hasSelected('2')).toBe(true);

      // Node 4: no files
      expect(hasFiles('4')).toBe(false);
      expect(format('4')).toBe('none');
      expect(hasSelected('4')).toBe(false);
    });

    it('should correctly identify running vs non-running algorithm nodes', () => {
      const isRunning = getters.isAlgorithmNodeRunning(mockState);
      const runningNodes = mockState.runningAlgorithmNodes;

      runningNodes.forEach(nodeId => {
        expect(isRunning(nodeId)).toBe(true);
      });

      ['1', '2', '4', '5', '6'].forEach(nodeId => {
        expect(isRunning(nodeId)).toBe(false);
      });
    });
  });
});
