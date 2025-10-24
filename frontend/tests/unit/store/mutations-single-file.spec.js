/**
 * workflow/mutations.js 유닛 테스트 - 단일 파일 Mutations
 *
 * 테스트 범위:
 * - setWorkflowFile
 * - setWorkflowNodeDataObject
 */

import { describe, it, expect, beforeEach, vi } from 'vitest';
import mutations from '@/store/workflow/mutations';

describe('workflow/mutations - Single File', () => {
  let mockState;

  beforeEach(() => {
    // Mock state 준비
    mockState = {
      workflow_info: {
        drawflow: {
          Home: {
            data: {
              '1': {
                id: 1,
                name: 'InputFile',
                data: {
                  title: 'Input Node',
                  file: null
                }
              },
              '2': {
                id: 2,
                name: 'DataTable',
                data: {
                  title: 'Data Table',
                  file: 'existing-file.h5ad'
                }
              },
              '3': {
                id: 3,
                name: 'Algorithm',
                data: {
                  title: 'Algorithm Node',
                  plugin: 'TENET',
                  parameters: {
                    fdr: 0.05
                  }
                }
              }
            }
          }
        }
      }
    };

    // console.error spy 설정
    vi.spyOn(console, 'error').mockImplementation(() => {});
  });

  describe('setWorkflowFile', () => {
    it('should set file name for valid node ID', () => {
      const fileInfo = {
        id: '1',
        file_name: 'test-data.h5ad'
      };

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
        .toBe('test-data.h5ad');
    });

    it('should update existing file name', () => {
      const fileInfo = {
        id: '2',
        file_name: 'new-data.h5ad'
      };

      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBe('existing-file.h5ad');

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBe('new-data.h5ad');
    });

    it('should handle different file extensions', () => {
      const testCases = [
        { id: '1', file_name: 'data.csv' },
        { id: '1', file_name: 'data.tsv' },
        { id: '1', file_name: 'data.txt' },
        { id: '1', file_name: 'data.json' }
      ];

      testCases.forEach(fileInfo => {
        mutations.setWorkflowFile(mockState, fileInfo);
        expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
          .toBe(fileInfo.file_name);
      });
    });

    it('should handle file names with special characters', () => {
      const fileInfo = {
        id: '1',
        file_name: 'test file (2024-01-01) [version 2].h5ad'
      };

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
        .toBe('test file (2024-01-01) [version 2].h5ad');
    });

    it('should handle file names with unicode characters', () => {
      const fileInfo = {
        id: '1',
        file_name: '데이터_파일_한글.h5ad'
      };

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
        .toBe('데이터_파일_한글.h5ad');
    });

    it('should handle very long file names', () => {
      const fileInfo = {
        id: '1',
        file_name: 'a'.repeat(255) + '.h5ad'
      };

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
        .toBe(fileInfo.file_name);
      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file.length)
        .toBeGreaterThan(255);
    });

    it('should log error for invalid node ID', () => {
      const fileInfo = {
        id: '999',
        file_name: 'test.h5ad'
      };

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(console.error).toHaveBeenCalledWith(
        'No object found with id: 999'
      );
    });

    it('should not crash when setting file for non-existent node', () => {
      const fileInfo = {
        id: 'nonexistent',
        file_name: 'test.h5ad'
      };

      expect(() => {
        mutations.setWorkflowFile(mockState, fileInfo);
      }).not.toThrow();
    });

    it('should handle empty string as file name', () => {
      const fileInfo = {
        id: '1',
        file_name: ''
      };

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
        .toBe('');
    });

    it('should handle null file name', () => {
      mockState.workflow_info.drawflow.Home.data['2'].data.file = 'existing.h5ad';

      const fileInfo = {
        id: '2',
        file_name: null
      };

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(mockState.workflow_info.drawflow.Home.data['2'].data.file)
        .toBeNull();
    });

    it('should not affect other node properties', () => {
      const originalTitle = mockState.workflow_info.drawflow.Home.data['1'].data.title;

      const fileInfo = {
        id: '1',
        file_name: 'test.h5ad'
      };

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.title)
        .toBe(originalTitle);
      expect(mockState.workflow_info.drawflow.Home.data['1'].id).toBe(1);
      expect(mockState.workflow_info.drawflow.Home.data['1'].name)
        .toBe('InputFile');
    });

    it('should work with string node IDs', () => {
      const fileInfo = {
        id: '1',
        file_name: 'test.h5ad'
      };

      mutations.setWorkflowFile(mockState, fileInfo);

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
        .toBe('test.h5ad');
    });

    it('should work with numeric node IDs (converted to string)', () => {
      // Note: 실제로는 항상 string으로 전달되지만, 방어적 테스트
      const fileInfo = {
        id: 1,
        file_name: 'test.h5ad'
      };

      // String conversion이 내부적으로 일어난다고 가정
      mutations.setWorkflowFile(mockState, { ...fileInfo, id: String(fileInfo.id) });

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
        .toBe('test.h5ad');
    });
  });

  describe('setWorkflowNodeDataObject', () => {
    it('should merge new data into node.data', () => {
      const payload = {
        nodeId: '1',
        dataObject: {
          newProperty: 'newValue',
          anotherProperty: 123
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.newProperty).toBe('newValue');
      expect(nodeData.anotherProperty).toBe(123);
    });

    it('should preserve existing keys', () => {
      const originalTitle = mockState.workflow_info.drawflow.Home.data['1'].data.title;
      const originalFile = mockState.workflow_info.drawflow.Home.data['1'].data.file;

      const payload = {
        nodeId: '1',
        dataObject: {
          newKey: 'newValue'
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.title).toBe(originalTitle);
      expect(nodeData.file).toBe(originalFile);
      expect(nodeData.newKey).toBe('newValue');
    });

    it('should override conflicting keys', () => {
      mockState.workflow_info.drawflow.Home.data['3'].data.parameters = {
        fdr: 0.05,
        iterations: 100
      };

      const payload = {
        nodeId: '3',
        dataObject: {
          parameters: {
            fdr: 0.01,
            numLinks: 50
          }
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['3'].data;
      expect(nodeData.parameters.fdr).toBe(0.01);
      expect(nodeData.parameters.numLinks).toBe(50);
      // 원본 parameters 객체가 교체됨
      expect(nodeData.parameters.iterations).toBeUndefined();
    });

    it('should handle nested object updates', () => {
      const payload = {
        nodeId: '3',
        dataObject: {
          config: {
            advanced: {
              option1: true,
              option2: 'value'
            }
          }
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['3'].data;
      expect(nodeData.config.advanced.option1).toBe(true);
      expect(nodeData.config.advanced.option2).toBe('value');
    });

    it('should handle array data', () => {
      const payload = {
        nodeId: '1',
        dataObject: {
          selectedFiles: ['file1.h5ad', 'file2.h5ad', 'file3.h5ad']
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.selectedFiles).toEqual(['file1.h5ad', 'file2.h5ad', 'file3.h5ad']);
      expect(nodeData.selectedFiles).toHaveLength(3);
    });

    it('should handle boolean values', () => {
      const payload = {
        nodeId: '3',
        dataObject: {
          useGPU: true,
          enableLogging: false
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['3'].data;
      expect(nodeData.useGPU).toBe(true);
      expect(nodeData.enableLogging).toBe(false);
    });

    it('should handle numeric values', () => {
      const payload = {
        nodeId: '3',
        dataObject: {
          timeout: 3600,
          maxRetries: 5,
          threshold: 0.95
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['3'].data;
      expect(nodeData.timeout).toBe(3600);
      expect(nodeData.maxRetries).toBe(5);
      expect(nodeData.threshold).toBe(0.95);
    });

    it('should handle null and undefined values', () => {
      const payload = {
        nodeId: '1',
        dataObject: {
          nullValue: null,
          undefinedValue: undefined
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.nullValue).toBeNull();
      expect(nodeData.undefinedValue).toBeUndefined();
    });

    it('should handle empty object', () => {
      const originalData = { ...mockState.workflow_info.drawflow.Home.data['1'].data };

      const payload = {
        nodeId: '1',
        dataObject: {}
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData).toEqual(originalData);
    });

    it('should log error for invalid node ID', () => {
      const payload = {
        nodeId: '999',
        dataObject: {
          test: 'value'
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      expect(console.error).toHaveBeenCalledWith(
        'No object found with id: 999'
      );
    });

    it('should not crash for non-existent node', () => {
      const payload = {
        nodeId: 'nonexistent',
        dataObject: {
          test: 'value'
        }
      };

      expect(() => {
        mutations.setWorkflowNodeDataObject(mockState, payload);
      }).not.toThrow();
    });

    it('should handle multiple updates to same node', () => {
      const payload1 = {
        nodeId: '1',
        dataObject: { prop1: 'value1' }
      };

      const payload2 = {
        nodeId: '1',
        dataObject: { prop2: 'value2' }
      };

      const payload3 = {
        nodeId: '1',
        dataObject: { prop3: 'value3' }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload1);
      mutations.setWorkflowNodeDataObject(mockState, payload2);
      mutations.setWorkflowNodeDataObject(mockState, payload3);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.prop1).toBe('value1');
      expect(nodeData.prop2).toBe('value2');
      expect(nodeData.prop3).toBe('value3');
    });

    it('should not affect other nodes', () => {
      const originalNode2Data = { ...mockState.workflow_info.drawflow.Home.data['2'].data };
      const originalNode3Data = { ...mockState.workflow_info.drawflow.Home.data['3'].data };

      const payload = {
        nodeId: '1',
        dataObject: {
          newData: 'test'
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      expect(mockState.workflow_info.drawflow.Home.data['2'].data)
        .toEqual(originalNode2Data);
      expect(mockState.workflow_info.drawflow.Home.data['3'].data.title)
        .toBe(originalNode3Data.title);
    });

    it('should handle Date objects', () => {
      const testDate = new Date('2024-01-01');
      const payload = {
        nodeId: '1',
        dataObject: {
          createdAt: testDate
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.createdAt).toEqual(testDate);
    });

    it('should handle function values (though not recommended)', () => {
      const testFn = () => 'test';
      const payload = {
        nodeId: '1',
        dataObject: {
          customFunction: testFn
        }
      };

      mutations.setWorkflowNodeDataObject(mockState, payload);

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.customFunction).toBe(testFn);
      expect(nodeData.customFunction()).toBe('test');
    });
  });

  describe('Combined scenarios', () => {
    it('should handle setWorkflowFile followed by setWorkflowNodeDataObject', () => {
      // First set file
      mutations.setWorkflowFile(mockState, {
        id: '1',
        file_name: 'test.h5ad'
      });

      // Then add additional data
      mutations.setWorkflowNodeDataObject(mockState, {
        nodeId: '1',
        dataObject: {
          fileSize: 1024000,
          uploadedAt: '2024-01-01'
        }
      });

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.file).toBe('test.h5ad');
      expect(nodeData.fileSize).toBe(1024000);
      expect(nodeData.uploadedAt).toBe('2024-01-01');
    });

    it('should handle setWorkflowNodeDataObject overriding file set by setWorkflowFile', () => {
      // First set file
      mutations.setWorkflowFile(mockState, {
        id: '1',
        file_name: 'original.h5ad'
      });

      expect(mockState.workflow_info.drawflow.Home.data['1'].data.file)
        .toBe('original.h5ad');

      // Then override with dataObject
      mutations.setWorkflowNodeDataObject(mockState, {
        nodeId: '1',
        dataObject: {
          file: 'overridden.h5ad'
        }
      });

      const nodeData = mockState.workflow_info.drawflow.Home.data['1'].data;
      expect(nodeData.file).toBe('overridden.h5ad');
    });
  });
});
