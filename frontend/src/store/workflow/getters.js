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
};
