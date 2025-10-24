/**
 * workflow/utils.js 유닛 테스트
 *
 * 테스트 범위:
 * - isAlgorithmNode
 * - getNodeFromState
 * - getNodeOutputConnections
 * - traverseGraphBFS
 * - propagateFileToConnectedNodes
 * - removeFileFromConnectedNodes
 */

import { describe, it, expect, beforeEach, vi } from 'vitest';
import {
  isAlgorithmNode,
  getNodeFromState,
  getNodeOutputConnections,
  traverseGraphBFS,
  propagateFileToConnectedNodes,
  removeFileFromConnectedNodes
} from '@/store/workflow/utils';

describe('workflow/utils', () => {
  let mockState;

  beforeEach(() => {
    // Console spy 설정
    vi.spyOn(console, 'error').mockImplementation(() => {});

    // Mock state 준비 - 연결된 워크플로우 그래프
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
              // Branching: 5 → 6, 7
              '5': {
                id: 5,
                name: 'InputFile',
                data: {
                  file: 'branch.h5ad'
                },
                inputs: {},
                outputs: {
                  'output_1': {
                    connections: [
                      { node: '6', output: 'input_1' },
                      { node: '7', output: 'input_1' }
                    ]
                  }
                }
              },
              '6': {
                id: 6,
                name: 'DataTable',
                data: {},
                inputs: {
                  'input_1': {
                    connections: [{ node: '5', input: 'output_1' }]
                  }
                },
                outputs: {}
              },
              '7': {
                id: 7,
                name: 'DataTable',
                data: {},
                inputs: {
                  'input_1': {
                    connections: [{ node: '5', input: 'output_1' }]
                  }
                },
                outputs: {}
              }
            }
          }
        }
      }
    };
  });

  describe('isAlgorithmNode', () => {
    it('should return true for Algorithm node', () => {
      const node = { name: 'Algorithm', data: {} };
      expect(isAlgorithmNode(node)).toBe(true);
    });

    it('should return false for non-Algorithm node', () => {
      const node = { name: 'DataTable', data: {} };
      expect(isAlgorithmNode(node)).toBe(false);
    });

    it('should return false for null node', () => {
      expect(isAlgorithmNode(null)).toBe(false);
    });

    it('should return false for undefined node', () => {
      expect(isAlgorithmNode(undefined)).toBe(false);
    });
  });

  describe('getNodeFromState', () => {
    it('should return node for valid ID', () => {
      const node = getNodeFromState(mockState, '1');
      expect(node).toBeDefined();
      expect(node.id).toBe(1);
      expect(node.name).toBe('InputFile');
    });

    it('should return null for non-existent node ID', () => {
      const node = getNodeFromState(mockState, '999');
      expect(node).toBeNull();
      expect(console.error).toHaveBeenCalledWith('No node found with id: 999');
    });

    it('should return null for invalid state structure', () => {
      const invalidState = {};
      const node = getNodeFromState(invalidState, '1');
      expect(node).toBeNull();
      expect(console.error).toHaveBeenCalledWith('Invalid state structure');
    });

    it('should return null for null state', () => {
      const node = getNodeFromState(null, '1');
      expect(node).toBeNull();
    });

    it('should handle node with all properties', () => {
      const node = getNodeFromState(mockState, '2');
      expect(node.name).toBe('DataTable');
      expect(node.inputs).toBeDefined();
      expect(node.outputs).toBeDefined();
    });
  });

  describe('getNodeOutputConnections', () => {
    it('should return empty array for node without outputs', () => {
      const node = { id: 1, name: 'Test', outputs: {} };
      const connections = getNodeOutputConnections(node);
      expect(connections).toEqual([]);
    });

    it('should return single connection', () => {
      const node = mockState.workflow_info.drawflow.Home.data['1'];
      const connections = getNodeOutputConnections(node);
      expect(connections).toHaveLength(1);
      expect(connections[0]).toEqual({
        nodeId: '2',
        outputKey: 'output_1'
      });
    });

    it('should return multiple connections from single output', () => {
      const node = mockState.workflow_info.drawflow.Home.data['5'];
      const connections = getNodeOutputConnections(node);
      expect(connections).toHaveLength(2);
      expect(connections[0].nodeId).toBe('6');
      expect(connections[1].nodeId).toBe('7');
    });

    it('should return empty array for null node', () => {
      const connections = getNodeOutputConnections(null);
      expect(connections).toEqual([]);
    });

    it('should return empty array for node without outputs property', () => {
      const node = { id: 1, name: 'Test' };
      const connections = getNodeOutputConnections(node);
      expect(connections).toEqual([]);
    });
  });

  describe('traverseGraphBFS', () => {
    it('should traverse simple linear graph', () => {
      const visited = [];
      traverseGraphBFS(mockState, '1', (current, target) => {
        visited.push(target.id);
        return 'continue';
      });

      expect(visited).toEqual([2, 3]); // Stops at Algorithm
    });

    it('should stop at Algorithm node by default', () => {
      const visited = [];
      traverseGraphBFS(mockState, '1', (current, target) => {
        visited.push(target.id);
        return 'continue';
      });

      expect(visited).not.toContain(4); // Visualization not reached
    });

    it('should call onAlgorithmNode callback', () => {
      const algorithmNodes = [];
      traverseGraphBFS(
        mockState,
        '1',
        () => 'continue',
        {
          onAlgorithmNode: (current, target) => {
            algorithmNodes.push(target.id);
          }
        }
      );

      expect(algorithmNodes).toEqual([3]);
    });

    it('should handle branching graph', () => {
      const visited = [];
      traverseGraphBFS(mockState, '5', (current, target) => {
        visited.push(target.id);
        return 'continue';
      });

      expect(visited).toHaveLength(2);
      expect(visited).toContain(6);
      expect(visited).toContain(7);
    });

    it('should stop traversal when callback returns "stop"', () => {
      const visited = [];
      traverseGraphBFS(mockState, '5', (current, target) => {
        visited.push(target.id);
        return 'stop'; // Stop immediately
      });

      expect(visited).toHaveLength(1); // Only first node visited
    });

    it('should skip branch when callback returns "skip"', () => {
      const visited = [];
      traverseGraphBFS(mockState, '5', (current, target) => {
        if (target.id === 6) {
          return 'skip';
        }
        visited.push(target.id);
        return 'continue';
      });

      expect(visited).toEqual([7]); // Node 6 skipped
    });

    it('should handle circular references with visited set', () => {
      // Create circular reference: 8 → 9 → 10 → 8
      mockState.workflow_info.drawflow.Home.data['8'] = {
        id: 8,
        name: 'DataTable',
        data: {},
        inputs: {
          'input_1': {
            connections: [{ node: '10', input: 'output_1' }]
          }
        },
        outputs: {
          'output_1': {
            connections: [{ node: '9', output: 'input_1' }]
          }
        }
      };
      mockState.workflow_info.drawflow.Home.data['9'] = {
        id: 9,
        name: 'DataTable',
        data: {},
        inputs: {
          'input_1': {
            connections: [{ node: '8', input: 'output_1' }]
          }
        },
        outputs: {
          'output_1': {
            connections: [{ node: '10', output: 'input_1' }]
          }
        }
      };
      mockState.workflow_info.drawflow.Home.data['10'] = {
        id: 10,
        name: 'DataTable',
        data: {},
        inputs: {
          'input_1': {
            connections: [{ node: '9', input: 'output_1' }]
          }
        },
        outputs: {
          'output_1': {
            connections: [{ node: '8', output: 'input_1' }]
          }
        }
      };

      const visited = [];
      traverseGraphBFS(mockState, '8', (current, target) => {
        visited.push(target.id);
        return 'continue';
      });

      // Should not infinite loop
      expect(visited.length).toBeLessThanOrEqual(3);
    });

    it('should handle node without output connections', () => {
      const visited = [];
      traverseGraphBFS(mockState, '4', (current, target) => {
        visited.push(target.id);
        return 'continue';
      });

      expect(visited).toEqual([]); // No outputs to traverse
    });

    it('should handle non-existent start node', () => {
      const visited = [];
      traverseGraphBFS(mockState, '999', (current, target) => {
        visited.push(target.id);
        return 'continue';
      });

      expect(visited).toEqual([]);
      expect(console.error).toHaveBeenCalledWith('No node found with id: 999');
    });
  });

  describe('propagateFileToConnectedNodes', () => {
    it('should propagate single file to connected nodes', () => {
      propagateFileToConnectedNodes(mockState, '1', 'test.h5ad', false);

      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file).toBe('test.h5ad');
    });

    it('should propagate single file to Algorithm node files object', () => {
      propagateFileToConnectedNodes(mockState, '1', 'test.h5ad', false);

      const algorithmNode = mockState.workflow_info.drawflow.Home.data['3'];
      expect(algorithmNode.data.files['1']).toBe('test.h5ad');
    });

    it('should propagate multi-file as array format', () => {
      const files = ['file1.h5ad', 'file2.h5ad'];
      propagateFileToConnectedNodes(mockState, '1', files, true);

      const targetNode = mockState.workflow_info.drawflow.Home.data['2'];
      expect(Array.isArray(targetNode.data.files)).toBe(true);
      expect(targetNode.data.files).toHaveLength(2);
      expect(targetNode.data.files[0]).toEqual({
        name: 'file1.h5ad',
        selected: true,
        size: 0
      });
    });

    it('should clear single file when propagating multi-file', () => {
      mockState.workflow_info.drawflow.Home.data['2'].data.file = 'old.h5ad';

      const files = ['file1.h5ad', 'file2.h5ad'];
      propagateFileToConnectedNodes(mockState, '1', files, true);

      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file).toBeUndefined();
    });

    it('should propagate to branching nodes', () => {
      propagateFileToConnectedNodes(mockState, '5', 'branch.h5ad', false);

      expect(mockState.workflow_info.drawflow.Home.data['6'].data.file).toBe('branch.h5ad');
      expect(mockState.workflow_info.drawflow.Home.data['7'].data.file).toBe('branch.h5ad');
    });

    it('should handle non-existent start node', () => {
      propagateFileToConnectedNodes(mockState, '999', 'test.h5ad', false);

      expect(console.error).toHaveBeenCalled();
    });
  });

  describe('removeFileFromConnectedNodes', () => {
    beforeEach(() => {
      // Set up files first
      mockState.workflow_info.drawflow.Home.data['2'].data.file = 'test.h5ad';
      mockState.workflow_info.drawflow.Home.data['3'].data.files = { '1': 'test.h5ad' };
    });

    it('should remove single file from connected nodes', () => {
      removeFileFromConnectedNodes(mockState, '1');

      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file).toBeNull();
    });

    it('should remove multi-file arrays from connected nodes', () => {
      mockState.workflow_info.drawflow.Home.data['2'].data.files = [
        { name: 'file1.h5ad', selected: true, size: 100 }
      ];

      removeFileFromConnectedNodes(mockState, '1');

      expect(mockState.workflow_info.drawflow.Home.data['2'].data.files).toEqual([]);
    });

    it('should remove from Algorithm node files object', () => {
      removeFileFromConnectedNodes(mockState, '1');

      const algorithmNode = mockState.workflow_info.drawflow.Home.data['3'];
      expect(algorithmNode.data.files['1']).toBeUndefined();
    });

    it('should remove from branching nodes', () => {
      mockState.workflow_info.drawflow.Home.data['6'].data.file = 'branch.h5ad';
      mockState.workflow_info.drawflow.Home.data['7'].data.file = 'branch.h5ad';

      removeFileFromConnectedNodes(mockState, '5');

      expect(mockState.workflow_info.drawflow.Home.data['6'].data.file).toBeNull();
      expect(mockState.workflow_info.drawflow.Home.data['7'].data.file).toBeNull();
    });

    it('should handle non-existent start node', () => {
      removeFileFromConnectedNodes(mockState, '999');

      expect(console.error).toHaveBeenCalled();
    });
  });

  describe('Integration scenarios', () => {
    it('should propagate then remove file correctly', () => {
      // Propagate
      propagateFileToConnectedNodes(mockState, '1', 'test.h5ad', false);
      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file).toBe('test.h5ad');

      // Remove
      removeFileFromConnectedNodes(mockState, '1');
      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file).toBeNull();
    });

    it('should handle multiple file propagations to same Algorithm node', () => {
      propagateFileToConnectedNodes(mockState, '1', 'file1.h5ad', false);

      // Add another source node
      mockState.workflow_info.drawflow.Home.data['11'] = {
        id: 11,
        name: 'InputFile',
        data: { file: 'file2.h5ad' },
        inputs: {},
        outputs: {
          'output_1': {
            connections: [{ node: '3', output: 'input_2' }]
          }
        }
      };

      propagateFileToConnectedNodes(mockState, '11', 'file2.h5ad', false);

      const algorithmNode = mockState.workflow_info.drawflow.Home.data['3'];
      expect(algorithmNode.data.files['1']).toBe('file1.h5ad');
      expect(algorithmNode.data.files['11']).toBe('file2.h5ad');
    });
  });
});
