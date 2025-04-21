import Vue from "vue";
import Vuex from "vuex";
import { getAuthFromCookie, saveAuthToCookie, getUserInfoFromCookie, saveUserInfoToCookie } from "@/utils/cookies";
import { loginUser } from "@/api/index";
import createPersistedState from "vuex-persistedstate";
import workflow from "./workflow/workflow";

Vue.use(Vuex);

export default new Vuex.Store({
  modules: {
    workflow,
  },
  plugins: [
    createPersistedState({
      paths: ["workflow", "userInfo"],
    }),
  ],
  state: {
    token: getAuthFromCookie() || "",
    userInfo: getUserInfoFromCookie() || { is_superuser: false },
  },
  getters: {
    isLogin(state) {
      return state.token !== "";
    },
    isSuperUser(state) {
      return state.userInfo && state.userInfo.is_superuser === true;
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
      saveUserInfoToCookie(userInfo);
    },
    clearUserInfo(state) {
      state.userInfo = { is_superuser: false };
    },
  },
  actions: {
    async LOGIN({ commit }, userData) {
      const response = await loginUser(userData);
      commit("setToken", response.data.access_token);
      commit("setUserInfo", response.data.user_info);
      saveAuthToCookie(response.data.access_token);
    },
    LOGOUT({ commit }) {
      commit("clearToken");
      commit("clearUserInfo");
    },
  },
});
