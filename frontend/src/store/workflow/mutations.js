import {
  getNodeFromState,
  propagateFileToConnectedNodes,
  removeFileFromConnectedNodes
} from './utils';
import { DEFAULT_VALUES } from './constants';

export default {
  setTitle(state, title) {
    state.title = title;
  },
  setThumbnail(state, thumbnail) {
    state.thumbnail = thumbnail;
  },
  clearTitle(state) {
    state.title = DEFAULT_VALUES.TITLE;
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
      state.workflow_info.drawflow.Home.data[file_info.id].data.fileSource = file_info.source || "user";
    } else {
      console.error(`No object found with id: ${file_info.id}`);
    }
  },
  shareWorkflowFile(state, id) {
    const node = getNodeFromState(state, id);
    if (!node) {
      return;
    }

    const file_name = node.data.file;
    if (!file_name) {
      console.error(`No file found in node with id: ${id}`);
      return;
    }

    propagateFileToConnectedNodes(state, id, file_name, false);
  },
  removeWorkflowFile(state, id) {
    const node = getNodeFromState(state, id);
    if (!node) {
      return;
    }

    // Check both single and multi-file formats
    const hasSingleFile = node.data.file;
    const hasMultiFiles = node.data.files && Array.isArray(node.data.files) && node.data.files.length > 0;

    if (!hasSingleFile && !hasMultiFiles) {
      console.error(`No file(s) found in node with id: ${id}`);
      return;
    }

    removeFileFromConnectedNodes(state, id);
  },
  
  removeWorkflowFiles(state, id) {
    const node = getNodeFromState(state, id);
    if (!node) {
      return;
    }

    // Clear multi-file format
    if (node.data.files && Array.isArray(node.data.files)) {
      node.data.files = [];
    }

    // Clear single file format
    if (node.data.file) {
      node.data.file = null;
    }

    // Propagate removal to connected nodes
    removeFileFromConnectedNodes(state, id);
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
    const node = getNodeFromState(state, id);
    if (!node) {
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

    propagateFileToConnectedNodes(state, id, filesToShare, isMultiFile);
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
