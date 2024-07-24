function isAlgorithmOptionsEmpty(options) {
  const defaultOptions = {
    algorithm: "",
    optionName: "",
    commonOptions: {
      annotationColumn: "",
      pseudotimeColumn: "",
      clusterColumn: [],
    },
    tenetOptions: null,
    optionFilePath: null,
  };

  return JSON.stringify(options) === JSON.stringify(defaultOptions);
}
export default {
  setTitle(state, title) {
    state.title = title;
  },
  setThumbnail(state, thumbnail) {
    state.thumbnail = thumbnail;
  },
  clearTitle(state) {
    state.title = "Untitled";
  },
  clearThumbnail(state) {
    state.thumbnail = null;
  },
  setWorkflow(state, workflow_info) {
    state.workflow_info = workflow_info;
  },
  clearWorkflow(state) {
    state.workflow_info = null;
  },
  setWorkflowFile(state, file_info) {
    if (state.workflow_info.drawflow.Home.data[file_info.id]) {
      console.log("setWorkflowFile this node data : " + state.workflow_info.drawflow.Home.data[file_info.id].data);
      console.log("file : " + file_info.file_name);
      state.workflow_info.drawflow.Home.data[file_info.id].data.file = file_info.file_name;
    } else {
      console.error(`No object found with id: ${file_info.id}`);
    }
  },
  shareWorkflowFile(state, id) {
    const node = state.workflow_info.drawflow.Home.data[id];
    console.log("implement shareWorkflowFile this node : " + node);

    if (!node) {
        console.error(`No node found with id: ${id}`);
        return;
    }

    if (node.name === 'Algorithm') {
        console.log(`Node with id: ${id} is of type 'Algorithm'. Function execution stopped.`);
        return;
    }

    const file_name = node.data.file;
    if (!file_name) {
        console.error(`No file found in node with id: ${id}`);
        return;
    }

    if (!Object.keys(node.outputs).some(outputKey => node.outputs[outputKey].connections.length > 0)) {
        console.log(`No connections found for node with id: ${id}`);
        return;
    }

    let files = {};

    // Iterate over outputs to find connections and add file data to connected nodes
    Object.keys(node.outputs).forEach(outputKey => {
        node.outputs[outputKey].connections.forEach(connection => {
            const targetNode = state.workflow_info.drawflow.Home.data[connection.node];
            
            if (targetNode) {
                if (!targetNode.data.files) {
                    targetNode.data.files = {};
                }
                targetNode.data.files[id] = file_name;

                // If target node has its own file, add it to the files object
                if (targetNode.data.file) {
                    files[connection.node] = targetNode.data.file;
                }

                // Recursively call shareWorkflowFile on the target node
                // this.shareWorkflowFile(state, connection.node);
            }
        });
    });

    // After updating connected nodes, add the files object to the current node's data
    if (Object.keys(files).length > 0) {
        node.data.files = files;
    }
  },
  updateWorkflowNodeTitle(state, { nodeId, newTitle }) {
    if (state.workflow_info.drawflow.Home.data[nodeId]) {
      state.workflow_info.drawflow.Home.data[nodeId].data.title = newTitle;
    } else {
      console.error(`No object found with id: ${nodeId}`);
    }
  },
  createNode(state, node) {
    state.nodes.push(node);
  },
  deleteNode(state, node) {
    state.nodes.forEach((obj, idx) => {
      if (obj.id === node.id) {
        state.nodes.splice(idx, 1);
      }
    });
  },
  changeFile(state, file) {
    state.nodes.forEach((obj) => {
      if (obj.id === state.current_node) {
        obj.file = file;
      }
    });
  },
  setNodes(state, nodes) {
    state.nodes = nodes;
  },
  updateNodeTitle(state, { nodeId, newTitle }) {
    const node = state.nodes.find((node) => node.id === nodeId);
    if (node) {
      node.title = newTitle;
    }
  },
  clearNodes(state) {
    state.nodes = [];
  },
  createConnection(state, connection_info) {
    // 경우 : 1. 연결된 노드가 없을 때
    if (state.linked_nodes.length === 0) {
      state.linked_nodes.push(connection_info);
    } else {
      for (let i = 0; i < state.linked_nodes.length; i++) {
        if (
          state.linked_nodes[i].connection.at(-1) ===
          connection_info.connection[0]
        ) {
          state.linked_nodes[i].connection.push(connection_info.connection[1]);
          state.linked_nodes[i].lastNode = connection_info.lastNode;
          return 0;
        } else if (
          state.linked_nodes[i].connection.at(0) ===
          connection_info.connection[1]
        ) {
          state.linked_nodes[i].connection.splice(
            0,
            0,
            connection_info.connection[0]
          );
          return 0;
        } else if (
          state.linked_nodes[i].connection.indexOf(
            connection_info.connection[0]
          ) != -1
        ) {
          const based_connection = state.linked_nodes[i].connection.slice(
            0,
            state.linked_nodes[i].connection.indexOf(
              connection_info.connection[0]
            )
          );
          state.linked_nodes.push({
            connection: [...based_connection, ...connection_info.connection],
            file: "",
            lastNode: connection_info.lastNode,
            algorithmOptions: {
              algorithm: "TENET",
              optionName: "Untitled",
              commonOptions: {
                annotationColumn: "",
                pseudotimeColumn: "",
                clusterColumn: [],
              },
              tenetOptions: null,
              optionFilePath: null,
            },
          });
          return 0;
        }
      }
      state.linked_nodes.push(connection_info);
    }
  },
  // 코드 설명
  // 1. state.linked_nodes를 순회하면서
  // 2. connection[0]이 있는지 확인
  // 3. 있으면 connection[1]을 찾아서
  // 4. connection[1]을 제거
  // 5. connection[1]이 없으면 connection[0]을 제거
  // 6. connection[0]이 없으면 아무것도 안함
  deleteConnection(state, connection) {
    state.linked_nodes.forEach((ele, idx) => {
      ele.connection.forEach((item, idx) => {
        if (
          item === connection[0] &&
          ele.connection[idx + 1] === connection[1]
        ) {
          if (ele.connection.length == 2) {
            ele.connection.splice(idx, 2);
          } else {
            ele.connection.splice(idx + 1, 1);
          }
        }
      });
      if (ele.connection.length === 0) {
        state.linked_nodes.splice(idx, 1);
      }
    });
  },
  setSelectedIndices(state, selectedList) {
    state.linked_nodes.forEach((obj) => {
      if (obj.connection.includes(state.current_node)) {
        obj.group = selectedList[0];
        obj.selectedIndices = selectedList[1];
      }
    });
  },
  shareConnectionFile(state) {
    // state.linked_nodes.forEach((ele) => {
    //   ele.forEach((item) => {
    //   })
    // })
    /** 노드들을 검사해 */
    state.nodes.forEach((node) => {
      // 노드의 이름이 File이고 file에 무언가 데이터가 있다면
      if (node.file != "" && node.name === "File") {
        //해당 노드와 연결된 노드들에 file_name을 부여
        state.linked_nodes.forEach((obj) => {
          if (obj.connection.includes(node.id)) {
            obj.file = node.file;
          }
        });
      }
    });
  },
  setLinkedNodes(state, linked_nodes) {
    state.linked_nodes = linked_nodes;
    // linked_nodes를 순회해서 각 아이템 안에 algorithmOptions가 없으면 추가합니다.
    state.linked_nodes.forEach((node) => {
      if (!node.algorithmOptions) {
        node.algorithmOptions = {
          algorithm: "TENET",
          optionName: "Untitled",
          commonOptions: {
            annotationColumn: "",
            pseudotimeColumn: "",
            clusterColumn: [],
          },
          tenetOptions: null,
          optionFilePath: null,
        };
      }
    });
  },
  clearLinkedNodes(state) {
    state.linked_nodes = [];
  },
  changeNode(state, id) {
    state.current_node = id;
  },
  clearCurrentNode(state) {
    state.current_node = 0;
  },
  setLinkedNodeAlgorithm(state, { nodeIndex, algorithm }) {
    if (state.linked_nodes[nodeIndex]) {
      state.linked_nodes[nodeIndex].algorithmOptions.algorithm = algorithm;
    }
  },
  setLinkedNodeCommonOptions(state, { nodeIndex, commonOptions }) {
    if (state.linked_nodes[nodeIndex]) {
      state.linked_nodes[nodeIndex].algorithmOptions.commonOptions =
        commonOptions;
    }
  },
  setLinkedNodeTenetOptions(state, { nodeIndex, tenetOptions }) {
    if (state.linked_nodes[nodeIndex]) {
      state.linked_nodes[nodeIndex].algorithmOptions.tenetOptions =
        tenetOptions;
    }
  },
  setLinkedNodeOptionFilePath(state, { nodeIndex, optionFilePath }) {
    if (state.linked_nodes[nodeIndex]) {
      state.linked_nodes[nodeIndex].algorithmOptions.optionFilePath =
        optionFilePath;
    }
  },
  // linked_nodes를 순회하면서 만약 중복된 connection을 가진 linked_nodes 2개 이상이 있다면 하나만 남기고 다 제거
  removeDuplicateLinkedNodes(state) {
    const newLinkedNodes = [];
    state.linked_nodes.forEach((node) => {
      const isDuplicated = newLinkedNodes.some((newNode) => {
        return newNode.connection.join() === node.connection.join();
      });
      if (!isDuplicated) {
        newLinkedNodes.push(node);
      }
    });
    state.linked_nodes = newLinkedNodes;
  },
  // 연결된 노드들 돌면서 connection에 포함된 노드ID들이 실제로 존재하는 노드의 ID인지 체크
  // 존재하지 않는다면 해당 노드를 제거
  removeInvalidLinkedNodes(state) {
    const newLinkedNodes = [];
    state.linked_nodes.forEach((node) => {
      const isInvalid = node.connection.some((nodeId) => {
        return !state.nodes.some((node) => node.id === nodeId);
      });
      if (!isInvalid) {
        newLinkedNodes.push(node);
      }
    });
    state.linked_nodes = newLinkedNodes;
  },
  // 연결된 노드가 2개 이상 있을 때만 실행
  // 연결된 노드를 돌면서 connection에 포함된 노드ID를 통해 algorithm node가 있는지 체크
  // algorithm node가 있다면 다른 연결된 노드에도 해당 algorithm node가 있는지 체크
  // 있다면 algortihmOptions가 채워져 있는 쪽 내용을 복사해 비어있는 쪽 채워넣기
  // 연결된 노드가 1개 이하라면 아무것도 하지 않음
  fillAlgorithmOptions(state) {
    if (state.linked_nodes.length > 1) {
      state.linked_nodes.forEach((node) => {
        // algorithmOptions가 비어있는지 확인
        if (!isAlgorithmOptionsEmpty(node.algorithmOptions)) {
          node.connection.forEach((nodeId) => {
            const algorithmNode = state.nodes.find(
              (node) => node.id === nodeId
            );
            if (algorithmNode && algorithmNode.name === "Algorithm") {
              const otherLinkedNode = state.linked_nodes.find(
                (linkedNode) =>
                  linkedNode.connection.includes(algorithmNode.id) &&
                  linkedNode.connection.length > 1 &&
                  linkedNode.connection.join() !== node.connection.join()
              );
              if (
                otherLinkedNode &&
                isAlgorithmOptionsEmpty(otherLinkedNode.algorithmOptions)
              ) {
                otherLinkedNode.algorithmOptions = node.algorithmOptions;
              }
            }
          });
        }
      });
    }
  },
};
