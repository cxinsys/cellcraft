export default {
  createNode(state, node) {
    state.nodes.push(node);
    console.log("yes");
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
  createConnection(state, connection_info) {
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
          });
          return 0;
        }
      }
      state.linked_nodes.push(connection_info);
    }
  },
  deleteConnection(state, connection) {
    state.linked_nodes.forEach((ele, idx) => {
      ele.forEach((item, idx) => {
        if (item === connection[0] && ele[idx + 1] === connection[1]) {
          if (ele.length == 2) {
            ele.splice(idx, 2);
          } else {
            ele.splice(idx + 1, 1);
          }
        }
      });
      if (ele.length === 0) {
        state.linked_nodes.splice(idx, 1);
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
  changeNode(state, id) {
    state.current_node = id;
  },
};
