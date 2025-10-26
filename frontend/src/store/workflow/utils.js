/**
 * Workflow Graph Utilities
 *
 * 공통 그래프 순회 및 노드 조작 유틸리티
 * mutations.js의 중복 코드를 제거하기 위한 헬퍼 함수들
 */

import { NODE_TYPES, NODE_PROPERTIES } from './constants';

/**
 * 노드가 Algorithm 타입인지 확인
 * @param {Object} node - Workflow node
 * @returns {boolean}
 */
export function isAlgorithmNode(node) {
  if (!node) {
    return false;
  }
  return node[NODE_PROPERTIES.NAME] === NODE_TYPES.ALGORITHM;
}

/**
 * State에서 노드를 안전하게 가져오기
 * @param {Object} state - Vuex state
 * @param {string} nodeId - Node ID
 * @returns {Object|null} - Node 또는 null
 */
export function getNodeFromState(state, nodeId) {
  if (!state || !state.workflow_info || !state.workflow_info.drawflow || !state.workflow_info.drawflow.Home) {
    console.error('Invalid state structure');
    return null;
  }

  const node = state.workflow_info.drawflow.Home.data[nodeId];
  if (!node) {
    console.error(`No node found with id: ${nodeId}`);
    return null;
  }

  return node;
}

/**
 * 노드의 출력 연결을 가져오기
 * @param {Object} node - Workflow node
 * @returns {Array} - Array of {nodeId, outputKey}
 */
export function getNodeOutputConnections(node) {
  if (!node || !node.outputs) {
    return [];
  }

  const connections = [];
  Object.keys(node.outputs).forEach(outputKey => {
    if (node.outputs[outputKey].connections) {
      node.outputs[outputKey].connections.forEach(connection => {
        connections.push({
          nodeId: connection.node,
          outputKey
        });
      });
    }
  });

  return connections;
}

/**
 * BFS 그래프 순회 (출력 방향)
 *
 * @param {Object} state - Vuex state
 * @param {string} startNodeId - 시작 노드 ID
 * @param {Function} onNodeVisit - 각 노드 방문 시 호출되는 콜백
 *   콜백 시그니처: (currentNode, targetNode, state) => 'continue' | 'stop' | 'skip'
 *   - 'continue': 다음 노드로 계속 진행
 *   - 'stop': 전체 순회 중단
 *   - 'skip': 현재 브랜치만 건너뛰기 (다른 브랜치는 계속)
 * @param {Object} options - 순회 옵션
 * @param {boolean} options.stopAtAlgorithm - Algorithm 노드에서 중단 (기본: true)
 * @param {Function} options.onAlgorithmNode - Algorithm 노드 만났을 때 호출되는 콜백
 */
export function traverseGraphBFS(state, startNodeId, onNodeVisit, options = {}) {
  const {
    stopAtAlgorithm = true,
    onAlgorithmNode = null
  } = options;

  const startNode = getNodeFromState(state, startNodeId);
  if (!startNode) {
    return;
  }

  // 출력 연결이 없으면 종료
  const hasOutputs = Object.keys(startNode.outputs || {}).some(
    outputKey => startNode.outputs[outputKey].connections.length > 0
  );
  if (!hasOutputs) {
    return;
  }

  // BFS 큐 초기화
  let currentNodes = [startNodeId];
  const visited = new Set(); // 순환 참조 방지

  while (currentNodes.length > 0) {
    const nextNodes = [];

    for (const currentNodeId of currentNodes) {
      // 이미 방문한 노드는 건너뛰기
      if (visited.has(currentNodeId)) {
        continue;
      }
      visited.add(currentNodeId);

      const currentNode = getNodeFromState(state, currentNodeId);
      if (!currentNode) {
        continue;
      }

      // Algorithm 노드 체크
      if (isAlgorithmNode(currentNode)) {
        if (stopAtAlgorithm) {
          return;
        }
        continue;
      }

      // 현재 노드의 출력 연결 순회
      const connections = getNodeOutputConnections(currentNode);

      for (const { nodeId } of connections) {
        const targetNode = getNodeFromState(state, nodeId);
        if (!targetNode) {
          continue;
        }

        // Algorithm 노드 특수 처리
        const isAlgorithm = isAlgorithmNode(targetNode);
        if (isAlgorithm && onAlgorithmNode) {
          onAlgorithmNode(currentNode, targetNode, state);
        }

        // 콜백 실행 (Algorithm 노드도 방문)
        const result = onNodeVisit(currentNode, targetNode, state);

        if (result === 'stop') {
          return;
        }

        if (result === 'skip') {
          continue;
        }

        // Algorithm 노드는 nextNodes에 추가하지 않음 (전파 중단)
        if (isAlgorithm && stopAtAlgorithm) {
          continue;
        }

        // 다음 레벨에 추가
        if (!visited.has(nodeId)) {
          nextNodes.push(nodeId);
        }
      }
    }

    currentNodes = nextNodes;
  }
}

/**
 * 파일을 연결된 노드들에 전파
 *
 * @param {Object} state - Vuex state
 * @param {string} sourceNodeId - 소스 노드 ID
 * @param {string|Array} filesToShare - 전파할 파일 (단일 파일명 또는 파일명 배열)
 * @param {boolean} isMultiFile - 멀티파일 여부
 */
export function propagateFileToConnectedNodes(state, sourceNodeId, filesToShare, isMultiFile = false) {
  traverseGraphBFS(
    state,
    sourceNodeId,
    (currentNode, targetNode) => {
      // Algorithm 노드는 onAlgorithmNode 콜백에서 이미 처리되므로 건너뛰기
      if (isAlgorithmNode(targetNode)) {
        return 'continue';
      }

      // 일반 노드에 파일 전파
      if (isMultiFile && Array.isArray(filesToShare)) {
        // 멀티파일: 배열 형식으로 변환
        targetNode.data.files = filesToShare.map(name => ({
          name,
          selected: true,
          size: 0
        }));
        // 기존 단일 파일 제거
        if (targetNode.data.file) {
          delete targetNode.data.file;
        }
      } else {
        // 단일 파일: 하위 호환성 유지
        targetNode.data.file = filesToShare;
      }

      return 'continue';
    },
    {
      stopAtAlgorithm: true,
      onAlgorithmNode: (currentNode, targetNode) => {
        // Algorithm 노드: files 객체에 소스 노드 ID를 키로 저장
        if (!targetNode.data.files) {
          targetNode.data.files = {};
        }
        targetNode.data.files[sourceNodeId] = filesToShare;
      }
    }
  );
}

/**
 * 연결된 노드들에서 파일 제거
 *
 * @param {Object} state - Vuex state
 * @param {string} sourceNodeId - 소스 노드 ID
 */
export function removeFileFromConnectedNodes(state, sourceNodeId) {
  traverseGraphBFS(
    state,
    sourceNodeId,
    (currentNode, targetNode) => {
      // Algorithm 노드는 onAlgorithmNode 콜백에서 이미 처리되므로 건너뛰기
      if (isAlgorithmNode(targetNode)) {
        return 'continue';
      }

      // 단일 파일 및 멀티파일 모두 제거
      targetNode.data.file = null;
      if (targetNode.data.files && Array.isArray(targetNode.data.files)) {
        targetNode.data.files = [];
      }

      return 'continue';
    },
    {
      stopAtAlgorithm: true,
      onAlgorithmNode: (currentNode, targetNode) => {
        // Algorithm 노드: files 객체에서 소스 노드 ID 삭제
        if (targetNode.data.files) {
          delete targetNode.data.files[sourceNodeId];
        }
      }
    }
  );
}
