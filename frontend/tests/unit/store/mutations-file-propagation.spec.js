/**
 * workflow/mutations.js 유닛 테스트 - File Propagation Mutations
 *
 * 테스트 범위:
 * - shareWorkflowFiles (파일 전파)
 * - removeWorkflowFile (파일 삭제 전파)
 * - removeWorkflowFiles (전체 삭제 전파)
 */

import { describe, it, expect, beforeEach, vi } from 'vitest';
import mutations from '@/store/workflow/mutations';

describe('workflow/mutations - File Propagation', () => {
  let mockState;

  beforeEach(() => {
    // Mock state with connected workflow graph
    mockState = {
      workflow_info: {
        drawflow: {
          Home: {
            data: {
              // InputFile → DataTable → Algorithm → Visualization
              '1': {
                id: 1,
                name: 'InputFile',
                data: {
                  title: 'Input',
                  file: 'input.h5ad'
                },
                inputs: {},
                outputs: {
                  'output_1': {
                    connections: [{ node: '2', output: 'input_1' }]
                  }
                }
              },
              '2': {
                id: 2,
                name: 'DataTable',
                data: {
                  title: 'Table'
                },
                inputs: {
                  'input_1': {
                    connections: [{ node: '1', input: 'output_1' }]
                  }
                },
                outputs: {
                  'output_1': {
                    connections: [{ node: '3', output: 'input_1' }]
                  }
                }
              },
              '3': {
                id: 3,
                name: 'Algorithm',
                data: {
                  title: 'TENET',
                  files: {}
                },
                inputs: {
                  'input_1': {
                    connections: [{ node: '2', input: 'output_1' }]
                  }
                },
                outputs: {
                  'output_1': {
                    connections: [{ node: '4', output: 'input_1' }]
                  }
                }
              },
              '4': {
                id: 4,
                name: 'Visualization',
                data: {
                  title: 'Viz'
                },
                inputs: {
                  'input_1': {
                    connections: [{ node: '3', input: 'output_1' }]
                  }
                },
                outputs: {}
              },
              // Isolated node with no connections
              '5': {
                id: 5,
                name: 'InputFile',
                data: {
                  title: 'Isolated',
                  file: 'isolated.h5ad'
                },
                inputs: {},
                outputs: {}
              },
              // Multi-file node
              '6': {
                id: 6,
                name: 'InputFile',
                data: {
                  title: 'Multi',
                  files: [
                    { name: 'file1.h5ad', selected: true, size: 100 },
                    { name: 'file2.h5ad', selected: true, size: 200 },
                    { name: 'file3.h5ad', selected: false, size: 300 }
                  ]
                },
                inputs: {},
                outputs: {
                  'output_1': {
                    connections: [{ node: '7', output: 'input_1' }]
                  }
                }
              },
              '7': {
                id: 7,
                name: 'DataTable',
                data: {
                  title: 'Table for Multi'
                },
                inputs: {
                  'input_1': {
                    connections: [{ node: '6', input: 'output_1' }]
                  }
                },
                outputs: {}
              },
              // Branching workflow: 8 → 9 and 8 → 10
              '8': {
                id: 8,
                name: 'InputFile',
                data: {
                  title: 'Branch Source',
                  file: 'branch.h5ad'
                },
                inputs: {},
                outputs: {
                  'output_1': {
                    connections: [
                      { node: '9', output: 'input_1' },
                      { node: '10', output: 'input_1' }
                    ]
                  }
                }
              },
              '9': {
                id: 9,
                name: 'DataTable',
                data: { title: 'Branch 1' },
                inputs: {
                  'input_1': {
                    connections: [{ node: '8', input: 'output_1' }]
                  }
                },
                outputs: {}
              },
              '10': {
                id: 10,
                name: 'DataTable',
                data: { title: 'Branch 2' },
                inputs: {
                  'input_1': {
                    connections: [{ node: '8', input: 'output_1' }]
                  }
                },
                outputs: {}
              }
            }
          }
        }
      }
    };

    // console spy
    vi.spyOn(console, 'error').mockImplementation(() => {});
  });

  describe('shareWorkflowFiles', () => {
    it('should propagate single file through connected nodes', () => {
      mutations.shareWorkflowFiles(mockState, '1');

      // Check propagation: 1 → 2 → 3 (Algorithm)
      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBe('input.h5ad');
      expect(mockState.workflow_info.drawflow.Home.data['3'].data.files['1'])
        .toBe('input.h5ad');
    });

    it('should stop propagation at Algorithm node', () => {
      mutations.shareWorkflowFiles(mockState, '1');

      // Algorithm node should receive the file in files object
      const algorithmFiles = mockState.workflow_info.drawflow.Home.data['3'].data.files;
      expect(algorithmFiles['1']).toBe('input.h5ad');

      // Visualization node after Algorithm should not receive the file
      expect(mockState.workflow_info.drawflow.Home.data['4'].data.file)
        .toBeUndefined();
    });

    it('should propagate multi-file through connected nodes', () => {
      mutations.shareWorkflowFiles(mockState, '6');

      // Selected files should be propagated
      const table = mockState.workflow_info.drawflow.Home.data['7'];
      expect(table.data.files).toBeDefined();
      expect(table.data.files).toHaveLength(2); // Only selected files
      expect(table.data.files[0].name).toBe('file1.h5ad');
      expect(table.data.files[1].name).toBe('file2.h5ad');
    });

    it('should handle branching workflow', () => {
      mutations.shareWorkflowFiles(mockState, '8');

      // Both branches should receive the file
      expect(mockState.workflow_info.drawflow.Home.data['9'].data.file)
        .toBe('branch.h5ad');
      expect(mockState.workflow_info.drawflow.Home.data['10'].data.file)
        .toBe('branch.h5ad');
    });

    it('should not propagate for isolated node', () => {
      mutations.shareWorkflowFiles(mockState, '5');

      // Isolated node should not trigger any errors
      expect(console.error).not.toHaveBeenCalled();
    });

    it('should handle node without output connections', () => {
      mutations.shareWorkflowFiles(mockState, '4');

      // Should not throw and should not propagate anything
      expect(() => mutations.shareWorkflowFiles(mockState, '4')).not.toThrow();
    });

    it('should log error for non-existent node', () => {
      mutations.shareWorkflowFiles(mockState, '999');

      expect(console.error).toHaveBeenCalledWith(
        'No node found with id: 999'
      );
    });

    it('should log error when no files are selected in multi-file node', () => {
      mockState.workflow_info.drawflow.Home.data['6'].data.files.forEach(f => {
        f.selected = false;
      });

      mutations.shareWorkflowFiles(mockState, '6');

      expect(console.error).toHaveBeenCalledWith(
        'No file(s) selected in node with id: 6'
      );
    });

    it('should log error when node has no file', () => {
      delete mockState.workflow_info.drawflow.Home.data['1'].data.file;

      mutations.shareWorkflowFiles(mockState, '1');

      expect(console.error).toHaveBeenCalledWith(
        'No file(s) selected in node with id: 1'
      );
    });

    it('should convert multi-file to multi-file format in target', () => {
      mutations.shareWorkflowFiles(mockState, '6');

      const table = mockState.workflow_info.drawflow.Home.data['7'];
      expect(Array.isArray(table.data.files)).toBe(true);
      expect(table.data.files[0]).toHaveProperty('name');
      expect(table.data.files[0]).toHaveProperty('selected');
      expect(table.data.files[0]).toHaveProperty('size');
    });

    it('should clear old single file when converting to multi-file', () => {
      mockState.workflow_info.drawflow.Home.data['7'].data.file = 'old-file.h5ad';

      mutations.shareWorkflowFiles(mockState, '6');

      expect(mockState.workflow_info.drawflow.Home.data['7'].data.file)
        .toBeUndefined();
    });

    it('should propagate single selected file as single file', () => {
      // Node 6 has only one selected file
      mockState.workflow_info.drawflow.Home.data['6'].data.files[1].selected = false;

      mutations.shareWorkflowFiles(mockState, '6');

      const table = mockState.workflow_info.drawflow.Home.data['7'];
      expect(table.data.file).toBe('file1.h5ad');
    });

    it('should store multi-file array in Algorithm node files mapping', () => {
      // Use existing chain: Node 1 → Node 2 → Node 3 (Algorithm)
      // Modify Node 1 to have multi-file format
      mockState.workflow_info.drawflow.Home.data['1'].data = {
        title: 'Multi Input',
        files: [
          { name: 'file1.h5ad', selected: true, size: 100 },
          { name: 'file2.h5ad', selected: true, size: 200 }
        ]
      };

      mutations.shareWorkflowFiles(mockState, '1');

      // Check Algorithm node receives multi-file array in files mapping
      const algorithmFiles = mockState.workflow_info.drawflow.Home.data['3'].data.files;
      expect(algorithmFiles['1']).toBeDefined();
      expect(Array.isArray(algorithmFiles['1'])).toBe(true);
      expect(algorithmFiles['1']).toEqual(['file1.h5ad', 'file2.h5ad']);
    });

    it('should handle multi-hop propagation', () => {
      // Create chain: 1 → 2 → 11 → 12
      mockState.workflow_info.drawflow.Home.data['11'] = {
        id: 11,
        name: 'DataTable',
        data: { title: 'Middle' },
        inputs: {
          'input_1': {
            connections: [{ node: '2', input: 'output_1' }]
          }
        },
        outputs: {
          'output_1': {
            connections: [{ node: '12', output: 'input_1' }]
          }
        }
      };
      mockState.workflow_info.drawflow.Home.data['12'] = {
        id: 12,
        name: 'DataTable',
        data: { title: 'End' },
        inputs: {
          'input_1': {
            connections: [{ node: '11', input: 'output_1' }]
          }
        },
        outputs: {}
      };
      mockState.workflow_info.drawflow.Home.data['2'].outputs.output_1.connections = [
        { node: '11', output: 'input_1' }
      ];

      mutations.shareWorkflowFiles(mockState, '1');

      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBe('input.h5ad');
      expect(mockState.workflow_info.drawflow.Home.data['11'].data.file)
        .toBe('input.h5ad');
      expect(mockState.workflow_info.drawflow.Home.data['12'].data.file)
        .toBe('input.h5ad');
    });
  });

  describe('removeWorkflowFile', () => {
    beforeEach(() => {
      // Set files first
      mutations.shareWorkflowFiles(mockState, '1');
    });

    it('should remove file and propagate removal', () => {
      mutations.removeWorkflowFile(mockState, '1');

      // File should be removed from DataTable
      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBeNull();

      // Algorithm files object should be updated
      expect(mockState.workflow_info.drawflow.Home.data['3'].data.files['1'])
        .toBeUndefined();
    });

    it('should handle branching removal', () => {
      mutations.shareWorkflowFiles(mockState, '8');
      mutations.removeWorkflowFile(mockState, '8');

      // Both branches should have null files
      expect(mockState.workflow_info.drawflow.Home.data['9'].data.file)
        .toBeNull();
      expect(mockState.workflow_info.drawflow.Home.data['10'].data.file)
        .toBeNull();
    });

    it('should log error for node without file', () => {
      delete mockState.workflow_info.drawflow.Home.data['1'].data.file;
      mockState.workflow_info.drawflow.Home.data['1'].data.files = [];

      mutations.removeWorkflowFile(mockState, '1');

      expect(console.error).toHaveBeenCalledWith(
        'No file(s) found in node with id: 1'
      );
    });

    it('should log error for non-existent node', () => {
      mutations.removeWorkflowFile(mockState, '999');

      expect(console.error).toHaveBeenCalledWith(
        'No node found with id: 999'
      );
    });

    it('should handle node with no output connections', () => {
      expect(() => mutations.removeWorkflowFile(mockState, '4')).not.toThrow();
    });

    it('should clear multi-file arrays during removal', () => {
      mutations.shareWorkflowFiles(mockState, '6');
      mutations.removeWorkflowFile(mockState, '6');

      const table = mockState.workflow_info.drawflow.Home.data['7'];
      expect(table.data.files).toEqual([]);
    });
  });

  describe('removeWorkflowFiles', () => {
    beforeEach(() => {
      // Set files first
      mutations.shareWorkflowFiles(mockState, '6');
    });

    it('should clear multi-file array', () => {
      mutations.removeWorkflowFiles(mockState, '6');

      const files = mockState.workflow_info.drawflow.Home.data['6'].data.files;
      expect(files).toEqual([]);
    });

    it('should clear single file format', () => {
      mockState.workflow_info.drawflow.Home.data['1'].data.file = 'test.h5ad';

      mutations.removeWorkflowFiles(mockState, '1');

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
        .toBeNull();
    });

    it('should propagate removal through connected nodes', () => {
      mutations.removeWorkflowFiles(mockState, '6');

      const table = mockState.workflow_info.drawflow.Home.data['7'];
      expect(table.data.files).toEqual([]);
      expect(table.data.file).toBeNull();
    });

    it('should log error for non-existent node', () => {
      mutations.removeWorkflowFiles(mockState, '999');

      expect(console.error).toHaveBeenCalled();
    });

    it('should handle already empty files array', () => {
      mockState.workflow_info.drawflow.Home.data['6'].data.files = [];

      expect(() => mutations.removeWorkflowFiles(mockState, '6')).not.toThrow();
    });

    it('should remove source node mapping from Algorithm files object', () => {
      // Setup: Node 1 → Node 2 → Node 3 (Algorithm)
      mutations.shareWorkflowFiles(mockState, '1');

      // Verify Algorithm has the mapping
      let algorithmFiles = mockState.workflow_info.drawflow.Home.data['3'].data.files;
      expect(algorithmFiles['1']).toBe('input.h5ad');

      // Remove files from Node 1
      mutations.removeWorkflowFiles(mockState, '1');

      // Check Algorithm.files['1'] is removed
      algorithmFiles = mockState.workflow_info.drawflow.Home.data['3'].data.files;
      expect(algorithmFiles['1']).toBeUndefined();
      expect(algorithmFiles).toEqual({}); // Should be empty object
    });
  });

  describe('Combined propagation scenarios', () => {
    it('should handle share → remove → share cycle', () => {
      // Share
      mutations.shareWorkflowFiles(mockState, '1');
      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBe('input.h5ad');

      // Remove
      mutations.removeWorkflowFile(mockState, '1');
      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBeNull();

      // Share again
      mutations.shareWorkflowFiles(mockState, '1');
      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBe('input.h5ad');
    });

    it('should handle multiple sources to same Algorithm node', () => {
      // Connect node 8 to algorithm 3
      mockState.workflow_info.drawflow.Home.data['8'].outputs.output_2 = {
        connections: [{ node: '3', output: 'input_2' }]
      };

      mutations.shareWorkflowFiles(mockState, '1');
      mutations.shareWorkflowFiles(mockState, '8');

      const algorithmFiles = mockState.workflow_info.drawflow.Home.data['3'].data.files;
      expect(algorithmFiles['1']).toBe('input.h5ad');
      expect(algorithmFiles['8']).toBe('branch.h5ad');
    });

    it('should handle complex multi-branch propagation', () => {
      mutations.shareWorkflowFiles(mockState, '8');

      expect(mockState.workflow_info.drawflow.Home.data['9'].data.file)
        .toBe('branch.h5ad');
      expect(mockState.workflow_info.drawflow.Home.data['10'].data.file)
        .toBe('branch.h5ad');

      mutations.removeWorkflowFile(mockState, '8');

      expect(mockState.workflow_info.drawflow.Home.data['9'].data.file)
        .toBeNull();
      expect(mockState.workflow_info.drawflow.Home.data['10'].data.file)
        .toBeNull();
    });
  });
});
