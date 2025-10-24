/**
 * workflow/actions.js 유닛 테스트
 *
 * 테스트 범위:
 * - compileNodes 액션의 API 호출 검증
 * - 올바른 데이터 전달 확인
 * - 성공/실패 시나리오 처리
 */

import { describe, it, expect, vi, beforeEach } from 'vitest';
import actions from '@/store/workflow/actions';
import * as api from '@/api/index';

// API 모듈 전체를 모킹
vi.mock('@/api/index', () => ({
  exportData: vi.fn()
}));

describe('workflow/actions', () => {
  describe('compileNodes', () => {
    let mockContext;
    let mockLinkedNodes;

    beforeEach(() => {
      // 각 테스트마다 mock 초기화
      vi.clearAllMocks();

      // Mock 데이터 준비
      mockLinkedNodes = {
        nodes: [
          { id: 1, type: 'InputFile', data: { file: 'test.h5ad' } },
          { id: 2, type: 'Algorithm', data: { plugin: 'TENET' } }
        ],
        connections: [
          { from: 1, to: 2 }
        ]
      };

      // Vuex context mock
      mockContext = {
        state: {
          linked_nodes: mockLinkedNodes
        },
        commit: vi.fn(),
        dispatch: vi.fn(),
        getters: {}
      };
    });

    it('should call exportData with stringified linked_nodes', async () => {
      // Arrange
      const expectedResponse = { data: { success: true, workflow_id: 123 } };
      api.exportData.mockResolvedValue(expectedResponse);

      // Act
      await actions.compileNodes(mockContext);

      // Assert
      expect(api.exportData).toHaveBeenCalledTimes(1);
      expect(api.exportData).toHaveBeenCalledWith(
        JSON.stringify(mockLinkedNodes)
      );
    });

    it('should return the response from exportData', async () => {
      // Arrange
      const expectedResponse = {
        data: {
          success: true,
          workflow_id: 456,
          message: 'Workflow compiled successfully'
        }
      };
      api.exportData.mockResolvedValue(expectedResponse);

      // Act
      const result = await actions.compileNodes(mockContext);

      // Assert
      expect(result).toEqual(expectedResponse);
    });

    it('should handle API errors properly', async () => {
      // Arrange
      const mockError = new Error('API Error: Failed to compile workflow');
      api.exportData.mockRejectedValue(mockError);

      // Act & Assert
      await expect(actions.compileNodes(mockContext)).rejects.toThrow(
        'API Error: Failed to compile workflow'
      );
    });

    it('should handle empty linked_nodes', async () => {
      // Arrange
      mockContext.state.linked_nodes = { nodes: [], connections: [] };
      const expectedResponse = { data: { success: true, workflow_id: null } };
      api.exportData.mockResolvedValue(expectedResponse);

      // Act
      await actions.compileNodes(mockContext);

      // Assert
      expect(api.exportData).toHaveBeenCalledWith(
        JSON.stringify({ nodes: [], connections: [] })
      );
    });

    it('should handle null linked_nodes gracefully', async () => {
      // Arrange
      mockContext.state.linked_nodes = null;
      const expectedResponse = { data: { success: false } };
      api.exportData.mockResolvedValue(expectedResponse);

      // Act
      await actions.compileNodes(mockContext);

      // Assert
      expect(api.exportData).toHaveBeenCalledWith(
        JSON.stringify(null)
      );
    });

    it('should preserve async behavior and await exportData', async () => {
      // Arrange
      let callbackExecuted = false;
      api.exportData.mockImplementation(() => {
        return new Promise((resolve) => {
          setTimeout(() => {
            callbackExecuted = true;
            resolve({ data: { success: true } });
          }, 10);
        });
      });

      // Act
      await actions.compileNodes(mockContext);

      // Assert
      expect(callbackExecuted).toBe(true);
    });

    it('should handle network timeout errors', async () => {
      // Arrange
      const timeoutError = new Error('Network timeout');
      timeoutError.code = 'ETIMEDOUT';
      api.exportData.mockRejectedValue(timeoutError);

      // Act & Assert
      await expect(actions.compileNodes(mockContext)).rejects.toMatchObject({
        message: 'Network timeout',
        code: 'ETIMEDOUT'
      });
    });

    it('should handle 500 server errors', async () => {
      // Arrange
      const serverError = new Error('Internal Server Error');
      serverError.response = { status: 500 };
      api.exportData.mockRejectedValue(serverError);

      // Act & Assert
      await expect(actions.compileNodes(mockContext)).rejects.toMatchObject({
        message: 'Internal Server Error',
        response: { status: 500 }
      });
    });
  });
});
