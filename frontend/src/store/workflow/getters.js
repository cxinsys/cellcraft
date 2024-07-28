export default {
  getTitle(state) {
    return state.title;
  },
  getThumbnail(state) {
    return state.thumbnail;
  },
  getWorkflowInfo(state) {
    return state.workflow_info;
  },
  getWorkflowNodeInfo(state) {
    return (id) => {
      return state.workflow_info.drawflow.Home.data[id];
    };
  },
  getWorkflowNodeFileInfo(state) {
    return (id) => {
        const node = state.workflow_info.drawflow.Home.data[id];
        if (node && node.data && node.data.file) {
            return node.data.file;
        } else {
            return null;
        }
    };
  },
  getWorkflowNodeFilesInfo(state) {
    return (id) => {
        const node = state.workflow_info.drawflow.Home.data[id];
        if (node && node.data && node.data.files) {
            return node.data.files;
        } else {
            return null;
        }
    };
  },
  getNodeInfo(state) {
    return (id) => {
      return state.nodes.find((node) => node.id === id);
    };
  },
  getCurrentNode(state) {
    return state.current_node;
  },
  getCurrentNodeInfo(state) {
    return state.nodes.find((node) => node.id === state.current_node);
  },
  getNodes(state) {
    return state.nodes;
  },
  getLinkedNodes(state) {
    return state.linked_nodes;
  },
  getCurrentLinkedNodes(state) {
    return state.linked_nodes.filter((node) =>
      node.connection.includes(state.current_node)
    );
  },
  getCurrentFile(state) {
    return state.linked_nodes.find((node) =>
      node.connection.includes(state.current_node)
    );
  },
  getNodesAlgorithmOptions(state) {
    return (nodeIndex) => {
      return state.linked_nodes[nodeIndex].algorithmOptions;
    };
  },
  // 현재 노드의 인덱스를 찾는 getter 함수
  getCurrentNodeIndex(state) {
    // 현재 노드의 ID를 찾습니다.
    const currentNodeId = state.current_node;
    // linked_nodes 배열을 순회하며 현재 노드의 인덱스를 찾습니다.
    const nodeIndex = state.linked_nodes.findIndex((node) =>
      node.connection.find((nodeId) => nodeId === currentNodeId)
    );
    // 인덱스를 반환합니다. 노드를 찾지 못한 경우 -1을 반환합니다.
    return nodeIndex;
  },
};
