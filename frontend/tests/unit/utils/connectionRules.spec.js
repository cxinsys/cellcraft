/**
 * connectionRules.js 유닛 테스트
 *
 * 테스트 범위:
 * - connectionRules 객체 구조 검증
 * - 각 노드 타입별 허용 연결 검증
 * - canConnect 함수 (있는 경우) 동작 검증
 */

import { describe, it, expect } from 'vitest';
import { connectionRules, canConnect } from '@/utils/connectionRules';

describe('connectionRules', () => {
  describe('기본 구조 검증', () => {
    it('connectionRules 객체가 정의되어 있어야 함', () => {
      expect(connectionRules).toBeDefined();
      expect(typeof connectionRules).toBe('object');
    });

    it('모든 노드 타입이 배열로 정의되어 있어야 함', () => {
      Object.values(connectionRules).forEach(rule => {
        expect(Array.isArray(rule)).toBe(true);
      });
    });
  });

  describe('InputFile 노드 연결 규칙', () => {
    it('InputFile은 DataTable, ScatterPlot, Algorithm에 연결 가능해야 함', () => {
      const allowedConnections = connectionRules.InputFile;

      expect(allowedConnections).toContain('DataTable');
      expect(allowedConnections).toContain('ScatterPlot');
      expect(allowedConnections).toContain('Algorithm');
    });

    it('InputFile은 Visualization에 직접 연결할 수 없어야 함', () => {
      const allowedConnections = connectionRules.InputFile;

      expect(allowedConnections).not.toContain('Visualization');
    });
  });

  describe('DataTable 노드 연결 규칙', () => {
    it('DataTable은 ScatterPlot, Algorithm에 연결 가능해야 함', () => {
      const allowedConnections = connectionRules.DataTable;

      expect(allowedConnections).toContain('ScatterPlot');
      expect(allowedConnections).toContain('Algorithm');
    });
  });

  describe('Algorithm 노드 연결 규칙', () => {
    it('Algorithm은 ResultFiles, Visualization에 연결 가능해야 함', () => {
      const allowedConnections = connectionRules.Algorithm;

      expect(allowedConnections).toContain('ResultFiles');
      expect(allowedConnections).toContain('Visualization');
    });
  });

  describe('ResultFiles 노드 연결 규칙', () => {
    it('ResultFiles는 Visualization에만 연결 가능해야 함', () => {
      const allowedConnections = connectionRules.ResultFiles;

      expect(allowedConnections).toContain('Visualization');
      expect(allowedConnections).toHaveLength(1);
    });
  });

  describe('Visualization 노드 연결 규칙', () => {
    it('Visualization은 아무 노드에도 연결할 수 없어야 함 (종단 노드)', () => {
      const allowedConnections = connectionRules.Visualization;

      expect(allowedConnections).toHaveLength(0);
    });
  });

  describe('canConnect 함수 테스트 (함수가 있는 경우)', () => {
    it('유효한 연결을 true로 반환해야 함', () => {
      if (typeof canConnect === 'function') {
        expect(canConnect('InputFile', 'Algorithm')).toBe(true);
        expect(canConnect('Algorithm', 'ResultFiles')).toBe(true);
        expect(canConnect('ResultFiles', 'Visualization')).toBe(true);
      }
    });

    it('무효한 연결을 false로 반환해야 함', () => {
      if (typeof canConnect === 'function') {
        expect(canConnect('InputFile', 'Visualization')).toBe(false);
        expect(canConnect('Visualization', 'Algorithm')).toBe(false);
        expect(canConnect('ResultFiles', 'InputFile')).toBe(false);
      }
    });

    it('존재하지 않는 노드 타입은 false로 반환해야 함', () => {
      if (typeof canConnect === 'function') {
        expect(canConnect('NonExistent', 'Algorithm')).toBe(false);
        expect(canConnect('InputFile', 'NonExistent')).toBe(false);
      }
    });
  });

  describe('GRN 워크플로우 시나리오', () => {
    it('일반적인 GRN 분석 플로우가 허용되어야 함', () => {
      // InputFile → Algorithm → ResultFiles → Visualization
      const flow = [
        { from: 'InputFile', to: 'Algorithm' },
        { from: 'Algorithm', to: 'ResultFiles' },
        { from: 'ResultFiles', to: 'Visualization' }
      ];

      flow.forEach(({ from, to }) => {
        const allowedConnections = connectionRules[from];
        expect(allowedConnections).toContain(to);
      });
    });

    it('데이터 시각화를 거친 플로우가 허용되어야 함', () => {
      // InputFile → DataTable → ScatterPlot → Algorithm
      const flow = [
        { from: 'InputFile', to: 'DataTable' },
        { from: 'DataTable', to: 'ScatterPlot' }
      ];

      flow.forEach(({ from, to }) => {
        const allowedConnections = connectionRules[from];
        expect(allowedConnections).toContain(to);
      });
    });
  });
});
