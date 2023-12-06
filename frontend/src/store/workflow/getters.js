export default {
  getTitle(state) {
    return state.title;
  },
  getThumbnail(state) {
    return state.thumbnail;
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
  getAlgorithmOptions(state) {
    return state.algorithmOptions;
  },
};
