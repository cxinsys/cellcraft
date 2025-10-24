/**
 * workflow/mutations.js 유닛 테스트 - 단순 Mutations
 *
 * 테스트 범위:
 * - setTitle, clearTitle
 * - setThumbnail, clearThumbnail
 * - setWorkflow, clearWorkflow
 * - updateWorkflowNodeTitle
 * - Running Algorithm Nodes 관리
 * - Task Algorithm Mapping 관리
 */

import { describe, it, expect, beforeEach, vi } from 'vitest';
import mutations from '@/store/workflow/mutations';

describe('workflow/mutations - Simple', () => {
  let mockState;

  beforeEach(() => {
    // Mock state 준비
    mockState = {
      title: 'Initial Title',
      thumbnail: null,
      workflow_info: {
        id: 1,
        name: 'Test Workflow',
        drawflow: {
          Home: {
            data: {
              '1': {
                id: 1,
                name: 'InputFile',
                data: { title: 'Original Title', file: 'test.h5ad' }
              },
              '2': {
                id: 2,
                name: 'Algorithm',
                data: { title: 'Algorithm Node' }
              }
            }
          }
        }
      },
      runningAlgorithmNodes: [],
      taskAlgorithmMapping: {}
    };

    // console.error spy 설정
    vi.spyOn(console, 'error').mockImplementation(() => {});
  });

  describe('setTitle', () => {
    it('should update state.title', () => {
      mutations.setTitle(mockState, 'New Workflow Title');
      expect(mockState.title).toBe('New Workflow Title');
    });

    it('should handle empty string', () => {
      mutations.setTitle(mockState, '');
      expect(mockState.title).toBe('');
    });

    it('should handle special characters', () => {
      const specialTitle = 'Test <>&"\' Workflow';
      mutations.setTitle(mockState, specialTitle);
      expect(mockState.title).toBe(specialTitle);
    });

    it('should handle very long titles', () => {
      const longTitle = 'A'.repeat(1000);
      mutations.setTitle(mockState, longTitle);
      expect(mockState.title).toBe(longTitle);
      expect(mockState.title).toHaveLength(1000);
    });
  });

  describe('clearTitle', () => {
    it('should reset title to "Untitled"', () => {
      mockState.title = 'Custom Title';
      mutations.clearTitle(mockState);
      expect(mockState.title).toBe('Untitled');
    });

    it('should reset from empty string to "Untitled"', () => {
      mockState.title = '';
      mutations.clearTitle(mockState);
      expect(mockState.title).toBe('Untitled');
    });

    it('should reset from null to "Untitled"', () => {
      mockState.title = null;
      mutations.clearTitle(mockState);
      expect(mockState.title).toBe('Untitled');
    });
  });

  describe('setThumbnail', () => {
    it('should update state.thumbnail', () => {
      const thumbnailData = 'data:image/png;base64,iVBORw0KGgo...';
      mutations.setThumbnail(mockState, thumbnailData);
      expect(mockState.thumbnail).toBe(thumbnailData);
    });

    it('should handle URL format', () => {
      const thumbnailUrl = 'https://example.com/thumbnail.png';
      mutations.setThumbnail(mockState, thumbnailUrl);
      expect(mockState.thumbnail).toBe(thumbnailUrl);
    });

    it('should overwrite existing thumbnail', () => {
      mockState.thumbnail = 'old-thumbnail';
      mutations.setThumbnail(mockState, 'new-thumbnail');
      expect(mockState.thumbnail).toBe('new-thumbnail');
    });

    it('should handle null value', () => {
      mutations.setThumbnail(mockState, null);
      expect(mockState.thumbnail).toBeNull();
    });
  });

  describe('clearThumbnail', () => {
    it('should reset thumbnail to null', () => {
      mockState.thumbnail = 'data:image/png;base64,...';
      mutations.clearThumbnail(mockState);
      expect(mockState.thumbnail).toBeNull();
    });

    it('should clear from URL format', () => {
      mockState.thumbnail = 'https://example.com/thumbnail.png';
      mutations.clearThumbnail(mockState);
      expect(mockState.thumbnail).toBeNull();
    });
  });

  describe('setWorkflow', () => {
    it('should update state.workflow_info', () => {
      const newWorkflow = {
        id: 2,
        name: 'New Workflow',
        drawflow: { Home: { data: {} } }
      };
      mutations.setWorkflow(mockState, newWorkflow);
      expect(mockState.workflow_info).toEqual(newWorkflow);
    });

    it('should replace entire workflow object', () => {
      const newWorkflow = { id: 5, name: 'Replaced' };
      mutations.setWorkflow(mockState, newWorkflow);
      expect(mockState.workflow_info.id).toBe(5);
      expect(mockState.workflow_info.name).toBe('Replaced');
    });

    it('should handle complex workflow structure', () => {
      const complexWorkflow = {
        id: 10,
        name: 'Complex',
        drawflow: {
          Home: {
            data: {
              '1': { id: 1, name: 'Node1' },
              '2': { id: 2, name: 'Node2' }
            }
          }
        },
        metadata: { created: '2024-01-01', author: 'test' }
      };
      mutations.setWorkflow(mockState, complexWorkflow);
      expect(mockState.workflow_info).toEqual(complexWorkflow);
      expect(mockState.workflow_info.metadata).toBeDefined();
    });
  });

  describe('clearWorkflow', () => {
    it('should reset workflow_info to null', () => {
      mutations.clearWorkflow(mockState);
      expect(mockState.workflow_info).toBeNull();
    });

    it('should clear complex workflow', () => {
      mockState.workflow_info = {
        id: 999,
        name: 'Complex',
        drawflow: { /* complex structure */ }
      };
      mutations.clearWorkflow(mockState);
      expect(mockState.workflow_info).toBeNull();
    });
  });

  describe('updateWorkflowNodeTitle', () => {
    it('should update node title by ID', () => {
      mutations.updateWorkflowNodeTitle(mockState, {
        nodeId: '1',
        newTitle: 'Updated Input File'
      });
      expect(mockState.workflow_info.drawflow.Home.data['1'].data.title)
        .toBe('Updated Input File');
    });

    it('should update different node', () => {
      mutations.updateWorkflowNodeTitle(mockState, {
        nodeId: '2',
        newTitle: 'Updated Algorithm'
      });
      expect(mockState.workflow_info.drawflow.Home.data['2'].data.title)
        .toBe('Updated Algorithm');
    });

    it('should log error for invalid node ID', () => {
      mutations.updateWorkflowNodeTitle(mockState, {
        nodeId: '999',
        newTitle: 'Invalid'
      });
      expect(console.error).toHaveBeenCalledWith(
        'No object found with id: 999'
      );
    });

    it('should not crash when updating non-existent node', () => {
      expect(() => {
        mutations.updateWorkflowNodeTitle(mockState, {
          nodeId: 'nonexistent',
          newTitle: 'Test'
        });
      }).not.toThrow();
    });

    it('should handle empty string as new title', () => {
      mutations.updateWorkflowNodeTitle(mockState, {
        nodeId: '1',
        newTitle: ''
      });
      expect(mockState.workflow_info.drawflow.Home.data['1'].data.title)
        .toBe('');
    });
  });

  describe('setRunningAlgorithmNodes', () => {
    it('should replace running algorithm nodes array', () => {
      mutations.setRunningAlgorithmNodes(mockState, ['1', '2', '3']);
      expect(mockState.runningAlgorithmNodes).toEqual(['1', '2', '3']);
    });

    it('should create new array (not mutate reference)', () => {
      const newNodes = ['5', '6'];
      mutations.setRunningAlgorithmNodes(mockState, newNodes);
      expect(mockState.runningAlgorithmNodes).toEqual(['5', '6']);
      expect(mockState.runningAlgorithmNodes).not.toBe(newNodes);
    });

    it('should handle empty array', () => {
      mockState.runningAlgorithmNodes = ['1', '2'];
      mutations.setRunningAlgorithmNodes(mockState, []);
      expect(mockState.runningAlgorithmNodes).toEqual([]);
    });

    it('should handle single node', () => {
      mutations.setRunningAlgorithmNodes(mockState, ['10']);
      expect(mockState.runningAlgorithmNodes).toEqual(['10']);
    });
  });

  describe('addRunningAlgorithmNode', () => {
    it('should add node ID if not present', () => {
      mutations.addRunningAlgorithmNode(mockState, '5');
      expect(mockState.runningAlgorithmNodes).toContain('5');
    });

    it('should not duplicate existing node ID', () => {
      mockState.runningAlgorithmNodes = ['1', '2'];
      mutations.addRunningAlgorithmNode(mockState, '1');
      expect(mockState.runningAlgorithmNodes).toEqual(['1', '2']);
      expect(mockState.runningAlgorithmNodes).toHaveLength(2);
    });

    it('should add to empty array', () => {
      mutations.addRunningAlgorithmNode(mockState, 'first');
      expect(mockState.runningAlgorithmNodes).toEqual(['first']);
    });

    it('should add multiple different nodes', () => {
      mutations.addRunningAlgorithmNode(mockState, '1');
      mutations.addRunningAlgorithmNode(mockState, '2');
      mutations.addRunningAlgorithmNode(mockState, '3');
      expect(mockState.runningAlgorithmNodes).toEqual(['1', '2', '3']);
    });

    it('should maintain order when adding', () => {
      mutations.addRunningAlgorithmNode(mockState, 'a');
      mutations.addRunningAlgorithmNode(mockState, 'b');
      mutations.addRunningAlgorithmNode(mockState, 'c');
      expect(mockState.runningAlgorithmNodes[0]).toBe('a');
      expect(mockState.runningAlgorithmNodes[1]).toBe('b');
      expect(mockState.runningAlgorithmNodes[2]).toBe('c');
    });
  });

  describe('removeRunningAlgorithmNode', () => {
    it('should remove node ID from array', () => {
      mockState.runningAlgorithmNodes = ['1', '2', '3'];
      mutations.removeRunningAlgorithmNode(mockState, '2');
      expect(mockState.runningAlgorithmNodes).toEqual(['1', '3']);
    });

    it('should handle non-existent node ID gracefully', () => {
      mockState.runningAlgorithmNodes = ['1', '2'];
      mutations.removeRunningAlgorithmNode(mockState, '999');
      expect(mockState.runningAlgorithmNodes).toEqual(['1', '2']);
    });

    it('should handle empty array', () => {
      mutations.removeRunningAlgorithmNode(mockState, '1');
      expect(mockState.runningAlgorithmNodes).toEqual([]);
    });

    it('should remove first node', () => {
      mockState.runningAlgorithmNodes = ['1', '2', '3'];
      mutations.removeRunningAlgorithmNode(mockState, '1');
      expect(mockState.runningAlgorithmNodes).toEqual(['2', '3']);
    });

    it('should remove last node', () => {
      mockState.runningAlgorithmNodes = ['1', '2', '3'];
      mutations.removeRunningAlgorithmNode(mockState, '3');
      expect(mockState.runningAlgorithmNodes).toEqual(['1', '2']);
    });

    it('should handle single node removal', () => {
      mockState.runningAlgorithmNodes = ['only'];
      mutations.removeRunningAlgorithmNode(mockState, 'only');
      expect(mockState.runningAlgorithmNodes).toEqual([]);
    });
  });

  describe('clearRunningAlgorithmNodes', () => {
    it('should empty the array', () => {
      mockState.runningAlgorithmNodes = ['1', '2', '3'];
      mutations.clearRunningAlgorithmNodes(mockState);
      expect(mockState.runningAlgorithmNodes).toEqual([]);
    });

    it('should handle already empty array', () => {
      mutations.clearRunningAlgorithmNodes(mockState);
      expect(mockState.runningAlgorithmNodes).toEqual([]);
    });

    it('should create new empty array', () => {
      mockState.runningAlgorithmNodes = ['1', '2'];
      mutations.clearRunningAlgorithmNodes(mockState);
      expect(mockState.runningAlgorithmNodes).toHaveLength(0);
    });
  });

  describe('setTaskAlgorithmMapping', () => {
    it('should set task-algorithm mapping object', () => {
      const mapping = { 'task-1': 'node-1', 'task-2': 'node-2' };
      mutations.setTaskAlgorithmMapping(mockState, mapping);
      expect(mockState.taskAlgorithmMapping).toEqual(mapping);
    });

    it('should create new object (not mutate reference)', () => {
      const mapping = { 'task-1': 'node-1' };
      mutations.setTaskAlgorithmMapping(mockState, mapping);
      expect(mockState.taskAlgorithmMapping).toEqual(mapping);
      expect(mockState.taskAlgorithmMapping).not.toBe(mapping);
    });

    it('should handle empty object', () => {
      mockState.taskAlgorithmMapping = { 'old': 'mapping' };
      mutations.setTaskAlgorithmMapping(mockState, {});
      expect(mockState.taskAlgorithmMapping).toEqual({});
    });

    it('should replace existing mapping', () => {
      mockState.taskAlgorithmMapping = { 'task-1': 'node-1' };
      mutations.setTaskAlgorithmMapping(mockState, { 'task-2': 'node-2' });
      expect(mockState.taskAlgorithmMapping).toEqual({ 'task-2': 'node-2' });
      expect(mockState.taskAlgorithmMapping['task-1']).toBeUndefined();
    });

    it('should handle complex mapping structure', () => {
      const complexMapping = {
        'task-123-abc': 'node-456-def',
        'task-789-xyz': 'node-012-uvw'
      };
      mutations.setTaskAlgorithmMapping(mockState, complexMapping);
      expect(mockState.taskAlgorithmMapping).toEqual(complexMapping);
    });
  });

  describe('clearTaskAlgorithmMapping', () => {
    it('should clear task-algorithm mapping', () => {
      mockState.taskAlgorithmMapping = {
        'task-1': 'node-1',
        'task-2': 'node-2'
      };
      mutations.clearTaskAlgorithmMapping(mockState);
      expect(mockState.taskAlgorithmMapping).toEqual({});
    });

    it('should handle already empty mapping', () => {
      mutations.clearTaskAlgorithmMapping(mockState);
      expect(mockState.taskAlgorithmMapping).toEqual({});
    });

    it('should create new empty object', () => {
      mockState.taskAlgorithmMapping = { 'task': 'node' };
      mutations.clearTaskAlgorithmMapping(mockState);
      expect(Object.keys(mockState.taskAlgorithmMapping)).toHaveLength(0);
    });
  });
});
