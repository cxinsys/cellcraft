import Vue from "vue";
import Vuex from "vuex";
import { getAuthFromCookie, saveAuthToCookie } from "@/utils/cookies";
import { loginUser } from "@/api/index";

Vue.use(Vuex);

export default new Vuex.Store({
  state: {
    token: getAuthFromCookie() || "",
    userInfo: "",
  },
  getters: {
    isLogin(state) {
      return state.token !== "";
    },
  },
  mutations: {
    setToken(state, token) {
      state.token = token;
    },
    clearToken(state) {
      state.token = "";
    },
    setUserInfo(state, userInfo) {
      state.userInfo = userInfo;
    },
  },
  actions: {
    async LOGIN({ commit }, userData) {
      const response = await loginUser(userData);
      console.log(response);
      commit("setToken", response.data.access_token);
      saveAuthToCookie(response.data.access_token);
    },
  },
});
