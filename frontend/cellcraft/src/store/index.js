import Vue from 'vue'
import Vuex from 'vuex'

Vue.use(Vuex)

export default new Vuex.Store({
  state: {
    token: '',
    userInfo: ''
  },
  getters: {
    isLogin (state) {
      return state.token !== ''
    }
  },
  mutations: {
    setToken (state, token) {
      state.token = token
    },
    clearToken (state, token) {
      state.token = ''
    },
    setUserInfo (state, userInfo) {
      state.userInfo = userInfo
    }
  }
})
