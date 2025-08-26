import state from "@/store/workflow/state";
import getters from "@/store/workflow/getters";
import actions from "@/store/workflow/actions";
import mutations from "@/store/workflow/mutations";
export default {
  state: {
    ...state,
  },
  getters: {
    ...getters,
  },
  actions: {
    ...actions,
  },
  mutations: {
    ...mutations,
  },
};
