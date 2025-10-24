import { describe, it, expect, beforeEach } from 'vitest';
import { HierarchicalLayout } from '@/utils/HierarchicalLayout';

describe('HierarchicalLayout.js', () => {
  let simpleNodes, simpleEdges;
  let branchingNodes, branchingEdges;
  let cyclicNodes, cyclicEdges;

  beforeEach(() => {
    // Simple linear DAG: A -> B -> C
    simpleNodes = [
      { id: 'A', label: 'Node A' },
      { id: 'B', label: 'Node B' },
      { id: 'C', label: 'Node C' }
    ];
    simpleEdges = [
      { source: 'A', target: 'B' },
      { source: 'B', target: 'C' }
    ];

    // Branching DAG: A -> B, A -> C, B -> D, C -> D
    branchingNodes = [
      { id: 'A', label: 'Root' },
      { id: 'B', label: 'Branch 1' },
      { id: 'C', label: 'Branch 2' },
      { id: 'D', label: 'Merge' }
    ];
    branchingEdges = [
      { source: 'A', target: 'B' },
      { source: 'A', target: 'C' },
      { source: 'B', target: 'D' },
      { source: 'C', target: 'D' }
    ];

    // Cyclic graph: A -> B -> C -> A
    cyclicNodes = [
      { id: 'A', label: 'Node A' },
      { id: 'B', label: 'Node B' },
      { id: 'C', label: 'Node C' }
    ];
    cyclicEdges = [
      { source: 'A', target: 'B' },
      { source: 'B', target: 'C' },
      { source: 'C', target: 'A' }
    ];
  });

  describe('constructor', () => {
    it('should initialize with default options', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);

      expect(layout.nodes).toEqual(simpleNodes);
      expect(layout.edges).toEqual(simpleEdges);
      expect(layout.options.direction).toBe('TB');
      expect(layout.options.rankSep).toBe(180);
      expect(layout.options.nodeSep).toBe(250);
      expect(layout.levels).toEqual({});
      expect(layout.positions).toEqual({});
    });

    it('should accept custom options', () => {
      const customOptions = {
        direction: 'LR',
        rankSep: 200,
        nodeSep: 300
      };
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges, customOptions);

      expect(layout.options.direction).toBe('LR');
      expect(layout.options.rankSep).toBe(200);
      expect(layout.options.nodeSep).toBe(300);
    });
  });

  describe('assignLevels', () => {
    it('should assign levels for simple linear DAG', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      const levels = layout.assignLevels();

      expect(levels['A']).toBe(0);
      expect(levels['B']).toBe(1);
      expect(levels['C']).toBe(2);
    });

    it('should handle branching DAG correctly', () => {
      const layout = new HierarchicalLayout(branchingNodes, branchingEdges);
      const levels = layout.assignLevels();

      expect(levels['A']).toBe(0);
      expect(levels['B']).toBe(1);
      expect(levels['C']).toBe(1);
      expect(levels['D']).toBe(2);
    });

    it('should handle cyclic graphs gracefully', () => {
      const layout = new HierarchicalLayout(cyclicNodes, cyclicEdges);
      const levels = layout.assignLevels();

      // All nodes should be assigned level 0 due to cycle
      expect(levels['A']).toBe(0);
      expect(levels['B']).toBeDefined();
      expect(levels['C']).toBeDefined();
    });

    it('should handle isolated nodes', () => {
      const nodes = [
        { id: 'A', label: 'Connected' },
        { id: 'B', label: 'Isolated' }
      ];
      const edges = [];

      const layout = new HierarchicalLayout(nodes, edges);
      const levels = layout.assignLevels();

      expect(levels['A']).toBe(0);
      expect(levels['B']).toBe(0);
    });
  });

  describe('calculatePositions', () => {
    it('should calculate positions for TB (top-bottom) direction', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      const positions = layout.calculatePositions('TB');

      expect(positions['A'].y).toBe(0);
      expect(positions['B'].y).toBe(180); // rankSep
      expect(positions['C'].y).toBe(360); // 2 * rankSep
      expect(positions['A'].level).toBe(0);
      expect(positions['B'].level).toBe(1);
      expect(positions['C'].level).toBe(2);
    });

    it('should calculate positions for LR (left-right) direction', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      const positions = layout.calculatePositions('LR');

      expect(positions['A'].x).toBe(0);
      expect(positions['B'].x).toBe(180);
      expect(positions['C'].x).toBe(360);
    });

    it('should handle branching with proper spacing', () => {
      const layout = new HierarchicalLayout(branchingNodes, branchingEdges);
      const positions = layout.calculatePositions('TB');

      // B and C should be at same level but different x positions
      expect(positions['B'].level).toBe(1);
      expect(positions['C'].level).toBe(1);
      expect(positions['B'].x).not.toBe(positions['C'].x);
      expect(positions['D'].level).toBe(2);
    });
  });

  describe('calculateDynamicSpacing', () => {
    it('should reduce spacing for many nodes', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      const spacing = layout.calculateDynamicSpacing(10);

      expect(spacing).toBeLessThan(layout.options.nodeSep);
      expect(spacing).toBeGreaterThanOrEqual(150); // minSpacing
    });

    it('should increase spacing for few nodes', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      const spacing = layout.calculateDynamicSpacing(2);

      expect(spacing).toBeGreaterThan(layout.options.nodeSep);
      expect(spacing).toBeLessThanOrEqual(400); // maxSpacing
    });

    it('should use base spacing for moderate node count', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      const spacing = layout.calculateDynamicSpacing(4);

      expect(spacing).toBe(layout.options.nodeSep);
    });
  });

  describe('minimizeCrossings and reorderLevel', () => {
    it('should minimize edge crossings through reordering', () => {
      const layout = new HierarchicalLayout(branchingNodes, branchingEdges);
      layout.assignLevels();

      // Setup level nodes
      layout.levelNodes = {
        0: ['A'],
        1: ['C', 'B'], // Intentionally wrong order
        2: ['D']
      };

      layout.minimizeCrossings();

      // After optimization, B and C should be in optimal order
      expect(layout.levelNodes[1]).toBeDefined();
      expect(layout.levelNodes[1].length).toBe(2);
    });

    it('should not reorder when already optimal', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      layout.assignLevels();
      layout.levelNodes = {
        0: ['A'],
        1: ['B'],
        2: ['C']
      };

      const originalOrder = [...layout.levelNodes[1]];
      layout.minimizeCrossings();

      expect(layout.levelNodes[1]).toEqual(originalOrder);
    });
  });

  describe('arraysEqual', () => {
    it('should return true for equal arrays', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);

      expect(layout.arraysEqual([1, 2, 3], [1, 2, 3])).toBe(true);
      expect(layout.arraysEqual(['A', 'B'], ['A', 'B'])).toBe(true);
    });

    it('should return false for different arrays', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);

      expect(layout.arraysEqual([1, 2, 3], [3, 2, 1])).toBe(false);
      expect(layout.arraysEqual([1, 2], [1, 2, 3])).toBe(false);
    });
  });

  describe('changeDirection', () => {
    it('should change direction and recalculate positions', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      layout.calculatePositions('TB');

      const tbY = layout.positions['B'].y;
      const newPositions = layout.changeDirection('LR');

      expect(layout.options.direction).toBe('LR');
      expect(newPositions['B'].x).toBe(tbY); // Y becomes X in LR
    });
  });

  describe('scaleToFit', () => {
    it('should scale layout to fit container', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      layout.calculatePositions('TB');

      const result = layout.scaleToFit(800, 600, 50);

      expect(result.positions).toBeDefined();
      expect(result.scale).toBeDefined();
      expect(result.bounds).toBeDefined();
      expect(result.scale).toBeLessThanOrEqual(1);
    });

    it('should handle empty positions gracefully', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);

      const result = layout.scaleToFit(800, 600);

      expect(result.positions).toBeDefined();
      expect(Object.keys(result.positions).length).toBeGreaterThan(0);
    });

    it('should center layout in container', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      layout.calculatePositions('TB');

      const result = layout.scaleToFit(800, 600, 50);
      const xValues = Object.values(result.positions).map(pos => pos.x);
      const yValues = Object.values(result.positions).map(pos => pos.y);

      // Scaled layout should be centered (approximately)
      const avgX = xValues.reduce((a, b) => a + b, 0) / xValues.length;
      const avgY = yValues.reduce((a, b) => a + b, 0) / yValues.length;

      expect(Math.abs(avgX)).toBeLessThan(100);
      expect(Math.abs(avgY)).toBeLessThan(100);
    });
  });

  describe('findCriticalPath', () => {
    it('should find longest path in linear DAG', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      layout.assignLevels();

      const path = layout.findCriticalPath();

      expect(path).toEqual(['A', 'B', 'C']);
    });

    it('should find critical path in branching DAG', () => {
      const layout = new HierarchicalLayout(branchingNodes, branchingEdges);
      layout.assignLevels();

      const path = layout.findCriticalPath();

      expect(path.length).toBe(3); // A -> (B or C) -> D
      expect(path[0]).toBe('A');
      expect(path[2]).toBe('D');
    });
  });

  describe('topologicalSort', () => {
    it('should return topologically sorted nodes', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      const sorted = layout.topologicalSort();

      expect(sorted).toEqual(['A', 'B', 'C']);
    });

    it('should handle branching DAG', () => {
      const layout = new HierarchicalLayout(branchingNodes, branchingEdges);
      const sorted = layout.topologicalSort();

      expect(sorted[0]).toBe('A');
      expect(sorted[3]).toBe('D');
      // B and C can be in any order
      expect(sorted.includes('B')).toBe(true);
      expect(sorted.includes('C')).toBe(true);
    });
  });

  describe('getLayoutInfo', () => {
    it('should return complete layout information', () => {
      const layout = new HierarchicalLayout(simpleNodes, simpleEdges);
      layout.calculatePositions();

      const info = layout.getLayoutInfo();

      expect(info.nodes).toEqual(simpleNodes);
      expect(info.edges).toEqual(simpleEdges);
      expect(info.positions).toBeDefined();
      expect(info.levels).toBeDefined();
      expect(info.levelNodes).toBeDefined();
      expect(info.options).toBeDefined();
      expect(info.criticalPath).toBeDefined();
    });
  });

  describe('handleSpecialCases', () => {
    it('should handle GRN nodes reordering', () => {
      const nodes = [
        { id: 'node1', label: 'Start' },
        { id: 'GRN_NumLinks', label: 'GRN NumLinks' },
        { id: 'GRN_FDR', label: 'GRN FDR' }
      ];
      const edges = [
        { source: 'node1', target: 'GRN_NumLinks' },
        { source: 'node1', target: 'GRN_FDR' }
      ];

      const layout = new HierarchicalLayout(nodes, edges);
      layout.assignLevels();
      layout.levelNodes = {
        0: ['node1'],
        2: ['GRN_NumLinks', 'GRN_FDR']
      };

      layout.handleSpecialCases();

      // GRN_FDR should be before GRN_NumLinks
      const level2 = layout.levelNodes[2];
      expect(level2.indexOf('GRN_FDR')).toBeLessThan(level2.indexOf('GRN_NumLinks'));
    });
  });

  describe('edge cases', () => {
    it('should handle empty graph', () => {
      const layout = new HierarchicalLayout([], []);
      const levels = layout.assignLevels();
      const positions = layout.calculatePositions();

      expect(Object.keys(levels).length).toBe(0);
      expect(Object.keys(positions).length).toBe(0);
    });

    it('should handle single node without edges', () => {
      const nodes = [{ id: 'A', label: 'Alone' }];
      const layout = new HierarchicalLayout(nodes, []);

      const levels = layout.assignLevels();
      const positions = layout.calculatePositions();

      expect(levels['A']).toBe(0);
      expect(positions['A']).toBeDefined();
      expect(positions['A'].x).toBe(0);
      expect(positions['A'].y).toBe(0);
    });

    it('should handle disconnected components', () => {
      const nodes = [
        { id: 'A', label: 'Graph 1 Start' },
        { id: 'B', label: 'Graph 1 End' },
        { id: 'C', label: 'Graph 2 Start' },
        { id: 'D', label: 'Graph 2 End' }
      ];
      const edges = [
        { source: 'A', target: 'B' },
        { source: 'C', target: 'D' }
      ];

      const layout = new HierarchicalLayout(nodes, edges);
      const levels = layout.assignLevels();

      expect(levels['A']).toBe(0);
      expect(levels['B']).toBe(1);
      expect(levels['C']).toBe(0);
      expect(levels['D']).toBe(1);
    });
  });
});
