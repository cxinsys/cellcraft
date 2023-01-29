export default {
  getNodeInfo(state) {
    return (id) => {
      return state.nodes.find((node) => node.id === id);
    };
  },
  getCurrentNode(state) {
    return state.current_node;
  },
  getCurrentFile(state) {
    return state.linked_nodes.find((node) =>
      node.connection.includes(state.current_node)
    );
  },
};
