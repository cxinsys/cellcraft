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
  getWorkflowVisualizationNodeInfo(state) {
    const nodes = state.workflow_info.drawflow.Home.data;
    return Object.values(nodes).filter(node => node.class === "Visualization");
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
  // id를 기반으로 node의 inputs들을 타고 들어가서 순회하면서 연결된 algorithm node를 찾는다.
  // algorithm node를 찾으면 해당 node의 id를 반환한다.
  getAlgorithmNodeConnectedToInput(state) {
    return (id) => {
      const node = state.workflow_info.drawflow.Home.data[id];
      if (!node) {
        console.error(`No node found with id: ${id}`);
        return null;
      }

      if (!Object.keys(node.inputs).some(inputKey => node.inputs[inputKey].connections.length > 0)) {
        console.log(`No connections found for node with id: ${id}`);
        return null;
      }

      let currentNodes = [id];
      while (currentNodes.length > 0) {
        const nextNodes = [];
        for (const currentNodeId of currentNodes) {
          const currentNode = state.workflow_info.drawflow.Home.data[currentNodeId];

          if (!currentNode) {
            console.error(`No node found with id: ${currentNodeId}`);
            continue;
          }

          // Check if the current node is of type "Algorithm"
          if (currentNode.class === "Algorithm") {
            return state.workflow_info.drawflow.Home.data[currentNodeId];
          }

          Object.keys(currentNode.inputs).forEach(inputKey => {
            currentNode.inputs[inputKey].connections.forEach(connection => {
              nextNodes.push(connection.node);
            });
          });
        }
        currentNodes = nextNodes;
      }

      return null;
    };
  }, 
};
