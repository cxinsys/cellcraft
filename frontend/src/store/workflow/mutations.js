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

    const file_name = node.data.file;
    if (!file_name) {
        console.error(`No file found in node with id: ${id}`);
        return;
    }

    if (!Object.keys(node.outputs).some(outputKey => node.outputs[outputKey].connections.length > 0)) {
        console.log(`No connections found for node with id: ${id}`);
        return;
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
          if (currentNode.name === 'Algorithm') {
              console.log(`Node with id: ${currentNodeId} is of type 'Algorithm'. Stopping.`);
              return;
          }
  
          // Iterate over the outputs to find connections
          Object.keys(currentNode.outputs).forEach(outputKey => {
              currentNode.outputs[outputKey].connections.forEach(connection => {
                  const targetNode = state.workflow_info.drawflow.Home.data[connection.node];
  
                  if (targetNode) {
                      if (targetNode.name === 'Algorithm') {
                          console.log(`Node with id: ${targetNode.id} is of type 'Algorithm'. Stopping.`);
                          if (!targetNode.data.files) {
                              targetNode.data.files = {};
                          }
                          targetNode.data.files[id] = file_name;
                          return;
                      }

                      targetNode.data.file = file_name;
  
                      // Add the connected node to the next nodes to process
                      nextNodes.push(connection.node);
                  }
              });
          });
      }
  
      currentNodes = nextNodes;
    }
  },
  removeWorkflowFile(state, id) {
    const node = state.workflow_info.drawflow.Home.data[id];
    if (!node) {
        console.error(`No node found with id: ${id}`);
        return;
    }

    if (!node.data.file) {
        console.error(`No file found in node with id: ${id}`);
        return;
    }

    if (!Object.keys(node.outputs).some(outputKey => node.outputs[outputKey].connections.length > 0)) {
        console.log(`No connections found for node with id: ${id}`);
        return;
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
            if (currentNode.name === 'Algorithm') {
                console.log(`Node with id: ${currentNodeId} is of type 'Algorithm'. Stopping.`);
                return;
            }

            // Iterate over the outputs to find connections
            Object.keys(currentNode.outputs).forEach(outputKey => {
                currentNode.outputs[outputKey].connections.forEach(connection => {
                    const targetNode = state.workflow_info.drawflow.Home.data[connection.node];

                    if (targetNode) {
                        if (targetNode.name === 'Algorithm') {
                            console.log(`Node with id: ${targetNode.id} is of type 'Algorithm'. Stopping.`);
                            if (targetNode.data.files) {
                                delete targetNode.data.files[id];
                            }
                            return;
                        }

                        targetNode.data.file = null;

                        // Add the connected node to the next nodes to process
                        nextNodes.push(connection.node);
                    }
                });
            });
        }

        currentNodes = nextNodes;
    }
  },
  updateWorkflowNodeTitle(state, { nodeId, newTitle }) {
    if (state.workflow_info.drawflow.Home.data[nodeId]) {
      state.workflow_info.drawflow.Home.data[nodeId].data.title = newTitle;
    } else {
      console.error(`No object found with id: ${nodeId}`);
    }
  },
  setWorkflowNodeDataObject(state, { nodeId, dataObject }) {
    if (state.workflow_info.drawflow.Home.data[nodeId]) {
      // data에 dataObject의 key-value를 추가한다. (기존 key는 유지하고 key : value가 없으면 추가)
      state.workflow_info.drawflow.Home.data[nodeId].data = {
        ...state.workflow_info.drawflow.Home.data[nodeId].data,
        ...dataObject,
      };
    } else {
      console.error(`No object found with id: ${nodeId}`);
    }
  },
};
