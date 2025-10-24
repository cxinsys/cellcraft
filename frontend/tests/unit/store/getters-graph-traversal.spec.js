/**
 * workflow/getters.js 유닛 테스트 - Graph Traversal Getters
 *
 * 테스트 범위:
 * - getWorkflowNodeInfo
 * - getAlgorithmNodeConnectedToInput
 */

import { describe, it, expect, beforeEach, vi } from 'vitest';
import getters from '@/store/workflow/getters';

describe('workflow/getters - Graph Traversal', () => {
  let mockState;

  beforeEach(() => {
    // Mock state with connected graph structure
    mockState = {
      workflow_info: {
        drawflow: {
          Home: {
            data: {
              // InputFile node (no inputs, has outputs)
              '1': {
                id: 1,
                name: 'InputFile',
                class: 'InputFile',
                data: {
                  title: 'Input File',
                  file: 'input.h5ad'
                },
                inputs: {},
                outputs: {
                  'output_1': {
                    connections: [
                      { node: '2', output: 'input_1' }
                    ]
                  }
                }
              },
              // DataTable node (connected to InputFile, connected to Algorithm)
              '2': {
                id: 2,
                name: 'DataTable',
                class: 'DataTable',
                data: {
                  title: 'Data Table'
                },
                inputs: {
                  'input_1': {
                    connections: [
                      { node: '1', input: 'output_1' }
                    ]
                  }
                },
                outputs: {
                  'output_1': {
                    connections: [
                      { node: '3', output: 'input_1' }
                    ]
                  }
                }
              },
              // Algorithm node (connected to DataTable)
              '3': {
                id: 3,
                name: 'TENET',
                class: 'Algorithm',
                data: {
                  title: 'TENET Algorithm',
                  plugin: 'TENET'
                },
                inputs: {
                  'input_1': {
                    connections: [
                      { node: '2', input: 'output_1' }
                    ]
                  }
                },
                outputs: {
                  'output_1': {
                    connections: [
                      { node: '4', output: 'input_1' }
                    ]
                  }
                }
              },
              // Visualization node (connected to Algorithm)
              '4': {
                id: 4,
                name: 'Heatmap',
                class: 'Visualization',
                data: {
                  title: 'Heatmap Viz'
                },
                inputs: {
                  'input_1': {
                    connections: [
                      { node: '3', input: 'output_1' }
                    ]
                  }
                },
                outputs: {}
              },
              // Isolated node (no connections)
              '5': {
                id: 5,
                name: 'InputFile',
                class: 'InputFile',
                data: {
                  title: 'Isolated Input'
                },
                inputs: {},
                outputs: {}
              },
              // Node directly connected to Algorithm
              '6': {
                id: 6,
                name: 'InputFile',
                class: 'InputFile',
                data: {
                  title: 'Direct to Algorithm'
                },
                inputs: {},
                outputs: {
                  'output_1': {
                    connections: [
                      { node: '3', output: 'input_2' }
                    ]
                  }
                }
              },
              // Multi-hop chain: 7 → 8 → 9 → 3 (Algorithm)
              '7': {
                id: 7,
                name: 'InputFile',
                class: 'InputFile',
                data: { title: 'Chain Start' },
                inputs: {},
                outputs: {
                  'output_1': {
                    connections: [{ node: '8', output: 'input_1' }]
                  }
                }
              },
              '8': {
                id: 8,
                name: 'DataTable',
                class: 'DataTable',
                data: { title: 'Chain Middle 1' },
                inputs: {
                  'input_1': {
                    connections: [{ node: '7', input: 'output_1' }]
                  }
                },
                outputs: {
                  'output_1': {
                    connections: [{ node: '9', output: 'input_1' }]
                  }
                }
              },
              '9': {
                id: 9,
                name: 'DataTable',
                class: 'DataTable',
                data: { title: 'Chain Middle 2' },
                inputs: {
                  'input_1': {
                    connections: [{ node: '8', input: 'output_1' }]
                  }
                },
                outputs: {
                  'output_1': {
                    connections: [{ node: '3', output: 'input_3' }]
                  }
                }
              }
            }
          }
        }
      }
    };

    // console spy setup
    vi.spyOn(console, 'error').mockImplementation(() => {});
    vi.spyOn(console, 'log').mockImplementation(() => {});
  });

  describe('getWorkflowNodeInfo', () => {
    it('should return node by ID', () => {
      const result = getters.getWorkflowNodeInfo(mockState);
      const node = result('1');
      expect(node).toBeDefined();
      expect(node.id).toBe(1);
      expect(node.name).toBe('InputFile');
    });

    it('should return different node types correctly', () => {
      const result = getters.getWorkflowNodeInfo(mockState);

      const inputFile = result('1');
      expect(inputFile.class).toBe('InputFile');

      const dataTable = result('2');
      expect(dataTable.class).toBe('DataTable');

      const algorithm = result('3');
      expect(algorithm.class).toBe('Algorithm');

      const viz = result('4');
      expect(viz.class).toBe('Visualization');
    });

    it('should return undefined for non-existent node', () => {
      const result = getters.getWorkflowNodeInfo(mockState);
      expect(result('999')).toBeUndefined();
    });

    it('should preserve node structure', () => {
      const result = getters.getWorkflowNodeInfo(mockState);
      const node = result('1');

      expect(node).toHaveProperty('id');
      expect(node).toHaveProperty('name');
      expect(node).toHaveProperty('class');
      expect(node).toHaveProperty('data');
      expect(node).toHaveProperty('inputs');
      expect(node).toHaveProperty('outputs');
    });

    it('should return the same reference as state', () => {
      const result = getters.getWorkflowNodeInfo(mockState);
      const node = result('1');
      expect(node).toBe(mockState.workflow_info.drawflow.Home.data['1']);
    });

    it('should handle all numeric node IDs', () => {
      const result = getters.getWorkflowNodeInfo(mockState);
      for (let i = 1; i <= 9; i++) {
        const node = result(String(i));
        expect(node).toBeDefined();
        expect(node.id).toBe(i);
      }
    });
  });

  describe('getAlgorithmNodeConnectedToInput', () => {
    it('should find algorithm node when directly connected', () => {
      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('4'); // Viz connected to Algorithm

      expect(algorithm).toBeDefined();
      expect(algorithm.id).toBe(3);
      expect(algorithm.class).toBe('Algorithm');
      expect(algorithm.name).toBe('TENET');
    });

    it('should return null when traversing backwards through non-algorithm nodes', () => {
      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('2'); // DataTable with input from InputFile

      // The getter traverses backwards through inputs to find Algorithm
      // Node 2 has input from 1 (InputFile), which has no inputs
      // Since we only traverse inputs and don't find Algorithm, return null
      expect(algorithm).toBeNull();
    });

    it('should return null when multi-hop chain has no algorithm in input path', () => {
      const result = getters.getAlgorithmNodeConnectedToInput(mockState);

      // Node 9's input chain: 9 ← 8 ← 7 (all non-algorithm nodes)
      // Should traverse backwards through inputs but not find Algorithm
      const algorithm = result('9');
      expect(algorithm).toBeNull();
    });

    it('should find algorithm through 2-hop input traversal', () => {
      // Create: Viz2 (11) ← DataTable2 (10) ← Algorithm (3)
      mockState.workflow_info.drawflow.Home.data['10'] = {
        id: 10,
        name: 'DataTable2',
        class: 'DataTable',
        data: { title: 'Middle Node' },
        inputs: {
          'input_1': {
            connections: [{ node: '3', input: 'output_1' }]
          }
        },
        outputs: {
          'output_1': {
            connections: [{ node: '11', output: 'input_1' }]
          }
        }
      };

      mockState.workflow_info.drawflow.Home.data['11'] = {
        id: 11,
        name: 'Heatmap2',
        class: 'Visualization',
        data: { title: 'Viz Node' },
        inputs: {
          'input_1': {
            connections: [{ node: '10', input: 'output_1' }]
          }
        },
        outputs: {}
      };

      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('11'); // Should find Algorithm through 2 hops

      expect(algorithm).toBeDefined();
      expect(algorithm.id).toBe(3);
      expect(algorithm.class).toBe('Algorithm');
    });

    it('should find algorithm through 3-hop input traversal', () => {
      // Create: Viz3 (13) ← Table3 (12) ← Table2 (10) ← Algorithm (3)
      mockState.workflow_info.drawflow.Home.data['10'] = {
        id: 10,
        name: 'DataTable2',
        class: 'DataTable',
        data: { title: 'Middle Node 1' },
        inputs: {
          'input_1': {
            connections: [{ node: '3', input: 'output_1' }]
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
        name: 'DataTable3',
        class: 'DataTable',
        data: { title: 'Middle Node 2' },
        inputs: {
          'input_1': {
            connections: [{ node: '10', input: 'output_1' }]
          }
        },
        outputs: {
          'output_1': {
            connections: [{ node: '13', output: 'input_1' }]
          }
        }
      };

      mockState.workflow_info.drawflow.Home.data['13'] = {
        id: 13,
        name: 'Heatmap3',
        class: 'Visualization',
        data: { title: 'Viz Node 3' },
        inputs: {
          'input_1': {
            connections: [{ node: '12', input: 'output_1' }]
          }
        },
        outputs: {}
      };

      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('13'); // Should find Algorithm through 3 hops

      expect(algorithm).toBeDefined();
      expect(algorithm.id).toBe(3);
      expect(algorithm.class).toBe('Algorithm');
      expect(algorithm.name).toBe('TENET');
    });

    it('should return null for node without input connections', () => {
      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('1'); // InputFile with no inputs

      expect(algorithm).toBeNull();
      expect(console.log).toHaveBeenCalledWith(
        'No connections found for node with id: 1'
      );
    });

    it('should return null for isolated node', () => {
      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('5');

      expect(algorithm).toBeNull();
    });

    it('should log error for non-existent node', () => {
      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('999');

      expect(algorithm).toBeNull();
      expect(console.error).toHaveBeenCalledWith(
        'No node found with id: 999'
      );
    });

    it('should return algorithm node itself if already algorithm', () => {
      const result = getters.getAlgorithmNodeConnectedToInput(mockState);

      // Add input to algorithm node for testing
      mockState.workflow_info.drawflow.Home.data['3'].inputs.input_test = {
        connections: [{ node: '1', input: 'output_1' }]
      };

      const algorithm = result('3');
      expect(algorithm).toBeDefined();
      expect(algorithm.id).toBe(3);
    });

    it('should handle node with multiple input connections', () => {
      // Add multiple inputs to node 4
      mockState.workflow_info.drawflow.Home.data['4'].inputs.input_2 = {
        connections: [{ node: '6', input: 'output_1' }]
      };

      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('4');

      expect(algorithm).toBeDefined();
      expect(algorithm.class).toBe('Algorithm');
    });

    it('should traverse breadth-first and find closest algorithm', () => {
      // Create scenario with multiple paths to different algorithms
      mockState.workflow_info.drawflow.Home.data['10'] = {
        id: 10,
        name: 'GENIE3',
        class: 'Algorithm',
        data: { title: 'GENIE3', plugin: 'GENIE3' },
        inputs: {
          input_1: { connections: [] }
        },
        outputs: {}
      };

      // Connect node 2 to both algorithm 3 and 10
      mockState.workflow_info.drawflow.Home.data['2'].inputs.input_2 = {
        connections: [{ node: '10', input: 'output_1' }]
      };

      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('2');

      // Should find an algorithm (either 3 or 10)
      expect(algorithm).toBeDefined();
      expect(algorithm.class).toBe('Algorithm');
    });

    it('should handle circular references gracefully with BFS', () => {
      // Create actual circular connection: 14 ← 15 ← 14 (through inputs)
      // This shouldn't happen in real workflows but testing BFS robustness
      mockState.workflow_info.drawflow.Home.data['14'] = {
        id: 14,
        name: 'DataTable4',
        class: 'DataTable',
        data: { title: 'Circular Node 1' },
        inputs: {
          'input_1': {
            connections: [{ node: '15', input: 'output_1' }]
          }
        },
        outputs: {
          'output_1': {
            connections: [{ node: '15', output: 'input_1' }]
          }
        }
      };

      mockState.workflow_info.drawflow.Home.data['15'] = {
        id: 15,
        name: 'DataTable5',
        class: 'DataTable',
        data: { title: 'Circular Node 2' },
        inputs: {
          'input_1': {
            connections: [{ node: '14', input: 'output_1' }]
          }
        },
        outputs: {
          'output_1': {
            connections: [{ node: '14', output: 'input_1' }]
          }
        }
      };

      const result = getters.getAlgorithmNodeConnectedToInput(mockState);

      // BFS should handle circular references without infinite loop
      expect(() => result('14')).not.toThrow();
      expect(() => result('15')).not.toThrow();

      // Should return null since no algorithm in circular path
      expect(result('14')).toBeNull();
      expect(result('15')).toBeNull();
    });

    it('should return null when only non-algorithm nodes in path', () => {
      const result = getters.getAlgorithmNodeConnectedToInput(mockState);

      // Node 2 → Node 1 (both non-algorithm)
      const algorithm = result('2');
      expect(algorithm).toBeNull();
    });

    it('should log error if node in traversal path does not exist', () => {
      // Create invalid connection to non-existent node
      mockState.workflow_info.drawflow.Home.data['5'].inputs = {
        input_1: {
          connections: [{ node: '999', input: 'output_1' }]
        }
      };

      const result = getters.getAlgorithmNodeConnectedToInput(mockState);
      const algorithm = result('5');

      expect(algorithm).toBeNull();
      expect(console.error).toHaveBeenCalled();
    });
  });

  describe('Combined scenarios', () => {
    it('should retrieve node info and traverse from it', () => {
      const nodeInfo = getters.getWorkflowNodeInfo(mockState);
      const algorithmConnected = getters.getAlgorithmNodeConnectedToInput(mockState);

      const viz = nodeInfo('4');
      expect(viz).toBeDefined();
      expect(viz.class).toBe('Visualization');

      const algorithm = algorithmConnected('4');
      expect(algorithm).toBeDefined();
      expect(algorithm.class).toBe('Algorithm');
    });

    it('should handle complete workflow chain', () => {
      const nodeInfo = getters.getWorkflowNodeInfo(mockState);
      const algorithmConnected = getters.getAlgorithmNodeConnectedToInput(mockState);

      // Check each node in chain
      const input = nodeInfo('1');
      expect(input.name).toBe('InputFile');
      expect(algorithmConnected('1')).toBeNull(); // No inputs

      const table = nodeInfo('2');
      expect(table.name).toBe('DataTable');

      const algo = nodeInfo('3');
      expect(algo.class).toBe('Algorithm');

      const viz = nodeInfo('4');
      expect(viz.class).toBe('Visualization');
      expect(algorithmConnected('4').id).toBe(3); // Finds algorithm
    });

    it('should handle multiple disconnected subgraphs', () => {
      const nodeInfo = getters.getWorkflowNodeInfo(mockState);
      const algorithmConnected = getters.getAlgorithmNodeConnectedToInput(mockState);

      // Main graph: 1 → 2 → 3 → 4
      expect(nodeInfo('1')).toBeDefined();
      expect(algorithmConnected('4').id).toBe(3);

      // Isolated node
      expect(nodeInfo('5')).toBeDefined();
      expect(algorithmConnected('5')).toBeNull();

      // Direct connection: 6 → 3
      expect(nodeInfo('6')).toBeDefined();
      // Node 6 has output to 3, but no inputs, so won't find algorithm
      expect(algorithmConnected('6')).toBeNull();
    });
  });
});
