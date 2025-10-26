/**
 * workflow/getters.js 유닛 테스트 - 기본 Getters
 *
 * 테스트 범위:
 * - getTitle, getThumbnail
 * - getWorkflowInfo
 * - getWorkflowVisualizationNodeInfo
 * - getRunningAlgorithmNodes
 * - getTaskAlgorithmMapping
 */

import { describe, it, expect, beforeEach } from 'vitest';
import getters from '@/store/workflow/getters';

describe('workflow/getters - Basic', () => {
  let mockState;

  beforeEach(() => {
    // Mock state 준비
    mockState = {
      title: 'Test Workflow',
      thumbnail: 'data:image/png;base64,iVBORw0KGgoAAAANSU...',
      workflow_info: {
        id: 1,
        name: 'My Workflow',
        drawflow: {
          Home: {
            data: {
              '1': {
                id: 1,
                name: 'InputFile',
                class: 'InputFile',
                data: { file: 'test.h5ad' },
                inputs: {},
                outputs: {}
              },
              '2': {
                id: 2,
                name: 'Algorithm',
                class: 'Algorithm',
                data: { plugin: 'TENET' },
                inputs: {},
                outputs: {}
              },
              '3': {
                id: 3,
                name: 'Visualization',
                class: 'Visualization',
                data: { type: 'heatmap' },
                inputs: {},
                outputs: {}
              },
              '4': {
                id: 4,
                name: 'Visualization',
                class: 'Visualization',
                data: { type: 'scatterplot' },
                inputs: {},
                outputs: {}
              }
            }
          }
        }
      },
      runningAlgorithmNodes: ['2', '5'],
      taskAlgorithmMapping: {
        'task-123': '2',
        'task-456': '5'
      }
    };
  });

  describe('getTitle', () => {
    it('should return current title from state', () => {
      const result = getters.getTitle(mockState);
      expect(result).toBe('Test Workflow');
    });

    it('should return null when title is not set', () => {
      mockState.title = null;
      const result = getters.getTitle(mockState);
      expect(result).toBeNull();
    });

    it('should return empty string when title is empty', () => {
      mockState.title = '';
      const result = getters.getTitle(mockState);
      expect(result).toBe('');
    });
  });

  describe('getThumbnail', () => {
    it('should return current thumbnail from state', () => {
      const result = getters.getThumbnail(mockState);
      expect(result).toBe('data:image/png;base64,iVBORw0KGgoAAAANSU...');
    });

    it('should return null when thumbnail is not set', () => {
      mockState.thumbnail = null;
      const result = getters.getThumbnail(mockState);
      expect(result).toBeNull();
    });

    it('should handle different thumbnail formats', () => {
      mockState.thumbnail = 'https://example.com/thumbnail.png';
      const result = getters.getThumbnail(mockState);
      expect(result).toBe('https://example.com/thumbnail.png');
    });
  });

  describe('getWorkflowInfo', () => {
    it('should return entire workflow_info object', () => {
      const result = getters.getWorkflowInfo(mockState);
      expect(result).toEqual(mockState.workflow_info);
      expect(result.id).toBe(1);
      expect(result.name).toBe('My Workflow');
    });

    it('should return null when workflow_info is not set', () => {
      mockState.workflow_info = null;
      const result = getters.getWorkflowInfo(mockState);
      expect(result).toBeNull();
    });

    it('should preserve drawflow structure', () => {
      const result = getters.getWorkflowInfo(mockState);
      expect(result.drawflow).toBeDefined();
      expect(result.drawflow.Home).toBeDefined();
      expect(result.drawflow.Home.data).toBeDefined();
    });

    it('should return the same reference to workflow_info', () => {
      const result = getters.getWorkflowInfo(mockState);
      expect(result).toBe(mockState.workflow_info);
    });
  });

  describe('getWorkflowVisualizationNodeInfo', () => {
    it('should filter nodes with class "Visualization"', () => {
      const result = getters.getWorkflowVisualizationNodeInfo(mockState);
      expect(result).toHaveLength(2);
      expect(result.every(node => node.class === 'Visualization')).toBe(true);
    });

    it('should return all visualization nodes with correct data', () => {
      const result = getters.getWorkflowVisualizationNodeInfo(mockState);

      const heatmapNode = result.find(node => node.data.type === 'heatmap');
      expect(heatmapNode).toBeDefined();
      expect(heatmapNode.id).toBe(3);

      const scatterplotNode = result.find(node => node.data.type === 'scatterplot');
      expect(scatterplotNode).toBeDefined();
      expect(scatterplotNode.id).toBe(4);
    });

    it('should return empty array when no visualization nodes exist', () => {
      // Remove visualization nodes
      delete mockState.workflow_info.drawflow.Home.data['3'];
      delete mockState.workflow_info.drawflow.Home.data['4'];

      const result = getters.getWorkflowVisualizationNodeInfo(mockState);
      expect(result).toEqual([]);
    });

    it('should not include Algorithm or InputFile nodes', () => {
      const result = getters.getWorkflowVisualizationNodeInfo(mockState);

      const hasNonVisualization = result.some(
        node => node.class !== 'Visualization'
      );
      expect(hasNonVisualization).toBe(false);
    });

    it('should handle empty drawflow data', () => {
      mockState.workflow_info.drawflow.Home.data = {};
      const result = getters.getWorkflowVisualizationNodeInfo(mockState);
      expect(result).toEqual([]);
    });
  });

  describe('getRunningAlgorithmNodes', () => {
    it('should return array of running algorithm node IDs', () => {
      const result = getters.getRunningAlgorithmNodes(mockState);
      expect(result).toEqual(['2', '5']);
    });

    it('should return empty array when no running nodes', () => {
      mockState.runningAlgorithmNodes = [];
      const result = getters.getRunningAlgorithmNodes(mockState);
      expect(result).toEqual([]);
    });

    it('should preserve array reference', () => {
      const result = getters.getRunningAlgorithmNodes(mockState);
      expect(result).toBe(mockState.runningAlgorithmNodes);
    });

    it('should handle single running node', () => {
      mockState.runningAlgorithmNodes = ['10'];
      const result = getters.getRunningAlgorithmNodes(mockState);
      expect(result).toEqual(['10']);
      expect(result).toHaveLength(1);
    });

    it('should handle multiple running nodes', () => {
      mockState.runningAlgorithmNodes = ['1', '2', '3', '4', '5'];
      const result = getters.getRunningAlgorithmNodes(mockState);
      expect(result).toHaveLength(5);
    });
  });

  describe('getTaskAlgorithmMapping', () => {
    it('should return task-to-algorithm mapping object', () => {
      const result = getters.getTaskAlgorithmMapping(mockState);
      expect(result).toEqual({
        'task-123': '2',
        'task-456': '5'
      });
    });

    it('should return empty object when no mappings exist', () => {
      mockState.taskAlgorithmMapping = {};
      const result = getters.getTaskAlgorithmMapping(mockState);
      expect(result).toEqual({});
    });

    it('should preserve object reference', () => {
      const result = getters.getTaskAlgorithmMapping(mockState);
      expect(result).toBe(mockState.taskAlgorithmMapping);
    });

    it('should handle single mapping', () => {
      mockState.taskAlgorithmMapping = { 'task-999': '10' };
      const result = getters.getTaskAlgorithmMapping(mockState);
      expect(result).toEqual({ 'task-999': '10' });
      expect(Object.keys(result)).toHaveLength(1);
    });

    it('should handle multiple mappings', () => {
      mockState.taskAlgorithmMapping = {
        'task-1': 'node-1',
        'task-2': 'node-2',
        'task-3': 'node-3'
      };
      const result = getters.getTaskAlgorithmMapping(mockState);
      expect(Object.keys(result)).toHaveLength(3);
    });

    it('should maintain key-value relationships', () => {
      const result = getters.getTaskAlgorithmMapping(mockState);
      expect(result['task-123']).toBe('2');
      expect(result['task-456']).toBe('5');
    });
  });
});
