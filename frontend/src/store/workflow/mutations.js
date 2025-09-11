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
      state.workflow_info.drawflow.Home.data[file_info.id].data.file = file_info.file_name;
    } else {
      console.error(`No object found with id: ${file_info.id}`);
    }
  },
  shareWorkflowFile(state, id) {
    const node = state.workflow_info.drawflow.Home.data[id];

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
              return;
          }
  
          // Iterate over the outputs to find connections
          Object.keys(currentNode.outputs).forEach(outputKey => {
              currentNode.outputs[outputKey].connections.forEach(connection => {
                  const targetNode = state.workflow_info.drawflow.Home.data[connection.node];
  
                  if (targetNode) {
                      if (targetNode.name === 'Algorithm') {
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

    // Check both single and multi-file formats
    const hasSingleFile = node.data.file;
    const hasMultiFiles = node.data.files && Array.isArray(node.data.files) && node.data.files.length > 0;
    
    if (!hasSingleFile && !hasMultiFiles) {
        console.error(`No file(s) found in node with id: ${id}`);
        return;
    }

    if (!Object.keys(node.outputs).some(outputKey => node.outputs[outputKey].connections.length > 0)) {
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

                        // Clear both single and multi-file formats
                        targetNode.data.file = null;
                        if (targetNode.data.files && Array.isArray(targetNode.data.files)) {
                          targetNode.data.files = [];
                        }

                        // Add the connected node to the next nodes to process
                        nextNodes.push(connection.node);
                    }
                });
            });
        }

        currentNodes = nextNodes;
    }
  },
  
  removeWorkflowFiles(state, id) {
    const node = state.workflow_info.drawflow.Home.data[id];
    if (!node) {
        console.error(`No node found with id: ${id}`);
        return;
    }

    // Clear multi-file format
    if (node.data.files && Array.isArray(node.data.files)) {
      node.data.files = [];
    }
    
    // Clear single file format for complete cleanup
    if (node.data.file) {
      node.data.file = null;
    }
    
    // Also propagate removal using existing logic
    // This ensures connected nodes are also cleared
    if (Object.keys(node.outputs).some(outputKey => node.outputs[outputKey].connections.length > 0)) {
      // Use existing removal logic by temporarily restoring a dummy file
      const originalFile = node.data.file;
      node.data.file = 'temp';
      
      // Call existing removal logic
      let currentNodes = [id];
      while (currentNodes.length > 0) {
          const nextNodes = [];
          for (const currentNodeId of currentNodes) {
              const currentNode = state.workflow_info.drawflow.Home.data[currentNodeId];

              if (!currentNode || currentNode.name === 'Algorithm') {
                  continue;
              }

              Object.keys(currentNode.outputs).forEach(outputKey => {
                  currentNode.outputs[outputKey].connections.forEach(connection => {
                      const targetNode = state.workflow_info.drawflow.Home.data[connection.node];

                      if (targetNode) {
                          if (targetNode.name === 'Algorithm') {
                              if (targetNode.data.files) {
                                  delete targetNode.data.files[id];
                              }
                              return;
                          }

                          targetNode.data.file = null;
                          if (targetNode.data.files && Array.isArray(targetNode.data.files)) {
                            targetNode.data.files = [];
                          }
                          nextNodes.push(connection.node);
                      }
                  });
              });
          }
          currentNodes = nextNodes;
      }
      
      // Clean up the temporary file
      node.data.file = originalFile;
    }
  },
  updateWorkflowNodeTitle(state, { nodeId, newTitle }) {
    if (state.workflow_info.drawflow.Home.data[nodeId]) {
      state.workflow_info.drawflow.Home.data[nodeId].data.title = newTitle;
    } else {
      console.error(`No object found with id: ${nodeId}`);
    }
  },
  
  // NEW: Multi-file support mutations (backward compatible)
  setWorkflowFiles(state, { id, files }) {
    const node = state.workflow_info.drawflow.Home.data[id];
    if (node) {
      // Set multi-file format: array of { name, size, selected }
      node.data.files = files;
      // Clear old single file format to avoid confusion
      if (node.data.file) {
        delete node.data.file;
      }
    } else {
      console.error(`No object found with id: ${id}`);
    }
  },
  
  setWorkflowSelectedFiles(state, { id, selectedFiles }) {
    const node = state.workflow_info.drawflow.Home.data[id];
    if (!node || !node.data.files || !Array.isArray(node.data.files)) {
      console.error(`No multi-file node found with id: ${id}`);
      return;
    }
    
    // Update selected status for all files
    node.data.files.forEach(file => {
      file.selected = selectedFiles.includes(file.name);
    });
  },
  
  toggleWorkflowFileSelection(state, { id, fileName }) {
    const node = state.workflow_info.drawflow.Home.data[id];
    if (!node || !node.data.files || !Array.isArray(node.data.files)) {
      console.error(`No multi-file node found with id: ${id}`);
      return;
    }
    
    const file = node.data.files.find(f => f.name === fileName);
    if (file) {
      file.selected = !file.selected;
    }
  },
  
  shareWorkflowFiles(state, id) {
    const node = state.workflow_info.drawflow.Home.data[id];

    if (!node) {
        console.error(`No node found with id: ${id}`);
        return;
    }

    // Determine what files to share - handle both single and multi-file formats
    let filesToShare = null;
    let isMultiFile = false;
    
    if (node.data.files && Array.isArray(node.data.files)) {
      // Multi-file format: get selected files
      const selectedFiles = node.data.files.filter(f => f.selected).map(f => f.name);
      if (selectedFiles.length > 0) {
        filesToShare = selectedFiles.length === 1 ? selectedFiles[0] : selectedFiles;
        isMultiFile = selectedFiles.length > 1;
      }
    } else if (node.data.file) {
      // Single file format (backward compatibility)
      filesToShare = node.data.file;
    }
    
    if (!filesToShare) {
        console.error(`No file(s) selected in node with id: ${id}`);
        return;
    }

    if (!Object.keys(node.outputs).some(outputKey => node.outputs[outputKey].connections.length > 0)) {
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
              return;
          }
  
          // Iterate over the outputs to find connections
          Object.keys(currentNode.outputs).forEach(outputKey => {
              currentNode.outputs[outputKey].connections.forEach(connection => {
                  const targetNode = state.workflow_info.drawflow.Home.data[connection.node];
  
                  if (targetNode) {
                      if (targetNode.name === 'Algorithm') {
                          if (!targetNode.data.files) {
                              targetNode.data.files = {};
                          }
                          // Store files from this source node
                          targetNode.data.files[id] = filesToShare;
                          return;
                      }

                      // For non-algorithm nodes, propagate in compatible format
                      if (isMultiFile) {
                        // If source has multiple files, convert to multi-file format
                        targetNode.data.files = filesToShare.map(name => ({ 
                          name, 
                          selected: true, 
                          size: 0 
                        }));
                        // Clear old single file
                        if (targetNode.data.file) {
                          delete targetNode.data.file;
                        }
                      } else {
                        // Single file: maintain backward compatibility
                        targetNode.data.file = filesToShare;
                      }
  
                      // Add the connected node to the next nodes to process
                      nextNodes.push(connection.node);
                  }
              });
          });
      }
  
      currentNodes = nextNodes;
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
  
  // NEW: Running algorithm nodes management
  setRunningAlgorithmNodes(state, nodeIds) {
    state.runningAlgorithmNodes = [...nodeIds];
  },
  
  addRunningAlgorithmNode(state, nodeId) {
    if (!state.runningAlgorithmNodes.includes(nodeId)) {
      state.runningAlgorithmNodes.push(nodeId);
    }
  },
  
  removeRunningAlgorithmNode(state, nodeId) {
    const index = state.runningAlgorithmNodes.indexOf(nodeId);
    if (index > -1) {
      state.runningAlgorithmNodes.splice(index, 1);
    }
  },
  
  clearRunningAlgorithmNodes(state) {
    state.runningAlgorithmNodes = [];
  },
  
  setTaskAlgorithmMapping(state, mapping) {
    state.taskAlgorithmMapping = { ...mapping };
  },
  
  clearTaskAlgorithmMapping(state) {
    state.taskAlgorithmMapping = {};
  },
};
