import { exportData } from "@/api/index";
export default {
  async compileNodes(context) {
    return await exportData(JSON.stringify(context.state.linked_nodes));
  },
};
