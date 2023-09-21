export default {
  getTitle(state) {
    return state.title;
  },
  getNodeInfo(state) {
    return (id) => {
      return state.nodes.find((node) => node.id === id);
    };
  },
  getCurrentNode(state) {
    return state.current_node;
  },
  getNodes(state) {
    return state.nodes;
  },
  getLinkedNodes(state) {
    return state.linked_nodes;
  },
  getCurrentFile(state) {
    return state.linked_nodes.find((node) =>
      node.connection.includes(state.current_node)
    );
  },
};
