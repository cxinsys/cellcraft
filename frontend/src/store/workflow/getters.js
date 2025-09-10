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
        if (!node || !node.data) {
            return null;
        }
        
        // NEW: Check for multi-file format first
        if (node.data.files && Array.isArray(node.data.files)) {
            const selectedFiles = node.data.files.filter(f => f.selected);
            if (selectedFiles.length === 1) {
                // Single selected file - return as string for backward compatibility
                return selectedFiles[0].name;
            } else if (selectedFiles.length > 1) {
                // Multiple selected files - return first selected for backward compatibility
                return selectedFiles[0].name;
            } else {
                // No files selected
                return null;
            }
        }
        
        // EXISTING: Single file fallback (unchanged)
        if (node.data.file) {
            return node.data.file;
        }
        
        return null;
    };
  },
  getWorkflowNodeFilesInfo(state) {
    return (id) => {
        const node = state.workflow_info.drawflow.Home.data[id];
        if (!node || !node.data) {
            return null;
        }
        
        // NEW: Return array format for multi-file nodes
        if (node.data.files && Array.isArray(node.data.files)) {
            return node.data.files;
        }
        
        // COMPATIBILITY: Convert single file to array format
        if (node.data.file) {
            return [{ 
                name: node.data.file, 
                selected: true, 
                size: 0 
            }];
        }
        
        // EXISTING: Handle algorithm files (files object format)
        if (node.data.files && typeof node.data.files === 'object' && !Array.isArray(node.data.files)) {
            return node.data.files;
        }
        
        return null;
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
  
  // NEW: Multi-file support getters
  isMultiFileNode(state) {
    return (id) => {
      const node = state.workflow_info.drawflow.Home.data[id];
      return node && node.data && Array.isArray(node.data.files);
    };
  },
  
  getSelectedWorkflowFiles(state) {
    return (id) => {
      const node = state.workflow_info.drawflow.Home.data[id];
      if (!node || !node.data) {
        return [];
      }
      
      // Multi-file format
      if (node.data.files && Array.isArray(node.data.files)) {
        return node.data.files.filter(f => f.selected);
      }
      
      // Single file format - convert to array
      if (node.data.file) {
        return [{ name: node.data.file, selected: true, size: 0 }];
      }
      
      return [];
    };
  },
  
  getSelectedWorkflowFileNames(state) {
    return (id) => {
      const node = state.workflow_info.drawflow.Home.data[id];
      if (!node || !node.data) {
        return [];
      }
      
      // Multi-file format
      if (node.data.files && Array.isArray(node.data.files)) {
        return node.data.files.filter(f => f.selected).map(f => f.name);
      }
      
      // Single file format
      if (node.data.file) {
        return [node.data.file];
      }
      
      return [];
    };
  },
  
  hasWorkflowFiles(state) {
    return (id) => {
      const node = state.workflow_info.drawflow.Home.data[id];
      if (!node || !node.data) {
        return false;
      }
      
      // Check multi-file format
      if (node.data.files && Array.isArray(node.data.files)) {
        return node.data.files.length > 0;
      }
      
      // Check single file format
      return !!node.data.file;
    };
  },
  
  hasSelectedWorkflowFiles(state) {
    return (id) => {
      const node = state.workflow_info.drawflow.Home.data[id];
      if (!node || !node.data) {
        return false;
      }
      
      // Check multi-file format
      if (node.data.files && Array.isArray(node.data.files)) {
        return node.data.files.some(f => f.selected);
      }
      
      // Check single file format
      return !!node.data.file;
    };
  },
  
  getWorkflowFileFormat(state) {
    return (id) => {
      const node = state.workflow_info.drawflow.Home.data[id];
      if (!node || !node.data) {
        return 'none';
      }
      
      if (node.data.files && Array.isArray(node.data.files)) {
        return 'multi';
      }
      
      if (node.data.file) {
        return 'single';
      }
      
      if (node.data.files && typeof node.data.files === 'object') {
        return 'algorithm';
      }
      
      return 'none';
    };
  },
  
  // NEW: Running algorithm nodes getters
  getRunningAlgorithmNodes(state) {
    return state.runningAlgorithmNodes;
  },
  
  isAlgorithmNodeRunning(state) {
    return (nodeId) => {
      return state.runningAlgorithmNodes.includes(String(nodeId));
    };
  },
  
  getTaskAlgorithmMapping(state) {
    return state.taskAlgorithmMapping;
  },
  
  getAlgorithmIdByTaskId(state) {
    return (taskId) => {
      return state.taskAlgorithmMapping[taskId] || null;
    };
  },
};
