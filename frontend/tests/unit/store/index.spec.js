/**
 * Root Store 유닛 테스트 - Authentication & Authorization
 *
 * 테스트 범위:
 * - 초기 상태 (쿠키 로딩)
 * - Getters: isLogin, isSuperUser
 * - Mutations: setToken, clearToken, setUserInfo, clearUserInfo
 * - Actions: LOGIN, LOGOUT
 */

import { describe, it, expect, beforeEach, vi } from 'vitest';
import Vuex from 'vuex';
import Vue from 'vue';

// Mock external dependencies
vi.mock('@/utils/cookies', () => ({
  getAuthFromCookie: vi.fn(),
  saveAuthToCookie: vi.fn(),
  getUserInfoFromCookie: vi.fn(),
  saveUserInfoToCookie: vi.fn()
}));

vi.mock('@/api/index', () => ({
  loginUser: vi.fn()
}));

// Import after mocking
import {
  getAuthFromCookie,
  saveAuthToCookie,
  getUserInfoFromCookie,
  saveUserInfoToCookie
} from '@/utils/cookies';
import { loginUser } from '@/api/index';

Vue.use(Vuex);

describe('Root Store - Authentication', () => {
  let store;

  beforeEach(() => {
    // Reset all mocks
    vi.clearAllMocks();

    // Default mock implementations
    getAuthFromCookie.mockReturnValue(null);
    getUserInfoFromCookie.mockReturnValue(null);

    // Create fresh store instance for each test
    store = new Vuex.Store({
      state: {
        token: getAuthFromCookie() || "",
        userInfo: getUserInfoFromCookie() || { is_superuser: false }
      },
      getters: {
        isLogin(state) {
          return state.token !== "";
        },
        isSuperUser(state) {
          return state.userInfo && state.userInfo.is_superuser === true;
        }
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
        }
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
        }
      }
    });
  });

  describe('Initial State', () => {
    it('should initialize with empty token when no cookie exists', () => {
      expect(store.state.token).toBe("");
    });

    it('should initialize with token from cookie when available', () => {
      getAuthFromCookie.mockReturnValue('test-token-123');

      store = new Vuex.Store({
        state: {
          token: getAuthFromCookie() || "",
          userInfo: getUserInfoFromCookie() || { is_superuser: false }
        },
        getters: {},
        mutations: {},
        actions: {}
      });

      expect(store.state.token).toBe('test-token-123');
    });

    it('should initialize with default userInfo when no cookie exists', () => {
      expect(store.state.userInfo).toEqual({ is_superuser: false });
    });

    it('should initialize with userInfo from cookie when available', () => {
      const mockUserInfo = {
        id: 1,
        username: 'testuser',
        is_superuser: true
      };
      getUserInfoFromCookie.mockReturnValue(mockUserInfo);

      store = new Vuex.Store({
        state: {
          token: getAuthFromCookie() || "",
          userInfo: getUserInfoFromCookie() || { is_superuser: false }
        },
        getters: {},
        mutations: {},
        actions: {}
      });

      expect(store.state.userInfo).toEqual(mockUserInfo);
    });
  });

  describe('Getters', () => {
    describe('isLogin', () => {
      it('should return false when token is empty string', () => {
        expect(store.getters.isLogin).toBe(false);
      });

      it('should return true when token exists', () => {
        store.state.token = 'valid-token';
        expect(store.getters.isLogin).toBe(true);
      });

      it('should return true for any non-empty token', () => {
        store.state.token = 'x';
        expect(store.getters.isLogin).toBe(true);
      });

      it('should be reactive to token changes', () => {
        expect(store.getters.isLogin).toBe(false);

        store.state.token = 'new-token';
        expect(store.getters.isLogin).toBe(true);

        store.state.token = '';
        expect(store.getters.isLogin).toBe(false);
      });
    });

    describe('isSuperUser', () => {
      it('should return false when userInfo.is_superuser is false', () => {
        expect(store.getters.isSuperUser).toBe(false);
      });

      it('should return true when userInfo.is_superuser is true', () => {
        store.state.userInfo = { is_superuser: true };
        expect(store.getters.isSuperUser).toBe(true);
      });

      it('should return false when userInfo is null', () => {
        store.state.userInfo = null;
        expect(store.getters.isSuperUser).toBeFalsy();
      });

      it('should return false when userInfo is undefined', () => {
        store.state.userInfo = undefined;
        expect(store.getters.isSuperUser).toBeFalsy();
      });

      it('should return false when is_superuser property is missing', () => {
        store.state.userInfo = { username: 'test' };
        expect(store.getters.isSuperUser).toBe(false);
      });

      it('should be reactive to userInfo changes', () => {
        expect(store.getters.isSuperUser).toBe(false);

        store.state.userInfo = { is_superuser: true };
        expect(store.getters.isSuperUser).toBe(true);

        store.state.userInfo = { is_superuser: false };
        expect(store.getters.isSuperUser).toBe(false);
      });
    });
  });

  describe('Mutations', () => {
    describe('setToken', () => {
      it('should set token in state', () => {
        store.commit('setToken', 'new-token-123');
        expect(store.state.token).toBe('new-token-123');
      });

      it('should overwrite existing token', () => {
        store.state.token = 'old-token';
        store.commit('setToken', 'new-token');
        expect(store.state.token).toBe('new-token');
      });

      it('should handle empty string token', () => {
        store.state.token = 'existing-token';
        store.commit('setToken', '');
        expect(store.state.token).toBe('');
      });

      it('should not call any cookie functions', () => {
        store.commit('setToken', 'token');
        expect(saveAuthToCookie).not.toHaveBeenCalled();
        expect(saveUserInfoToCookie).not.toHaveBeenCalled();
      });
    });

    describe('clearToken', () => {
      it('should clear token to empty string', () => {
        store.state.token = 'existing-token';
        store.commit('clearToken');
        expect(store.state.token).toBe('');
      });

      it('should work when token is already empty', () => {
        store.state.token = '';
        store.commit('clearToken');
        expect(store.state.token).toBe('');
      });
    });

    describe('setUserInfo', () => {
      it('should set userInfo in state', () => {
        const userInfo = {
          id: 1,
          username: 'testuser',
          is_superuser: false
        };

        store.commit('setUserInfo', userInfo);
        expect(store.state.userInfo).toEqual(userInfo);
      });

      it('should call saveUserInfoToCookie with userInfo', () => {
        const userInfo = {
          id: 2,
          username: 'admin',
          is_superuser: true
        };

        store.commit('setUserInfo', userInfo);
        expect(saveUserInfoToCookie).toHaveBeenCalledWith(userInfo);
        expect(saveUserInfoToCookie).toHaveBeenCalledTimes(1);
      });

      it('should overwrite existing userInfo', () => {
        store.state.userInfo = { username: 'old' };

        const newUserInfo = { username: 'new', is_superuser: true };
        store.commit('setUserInfo', newUserInfo);

        expect(store.state.userInfo).toEqual(newUserInfo);
      });

      it('should handle partial userInfo objects', () => {
        const partialInfo = { is_superuser: true };
        store.commit('setUserInfo', partialInfo);
        expect(store.state.userInfo).toEqual(partialInfo);
      });
    });

    describe('clearUserInfo', () => {
      it('should reset userInfo to default state', () => {
        store.state.userInfo = {
          id: 1,
          username: 'testuser',
          is_superuser: true
        };

        store.commit('clearUserInfo');
        expect(store.state.userInfo).toEqual({ is_superuser: false });
      });

      it('should work when userInfo is already at default', () => {
        store.state.userInfo = { is_superuser: false };
        store.commit('clearUserInfo');
        expect(store.state.userInfo).toEqual({ is_superuser: false });
      });

      it('should not call saveUserInfoToCookie', () => {
        store.commit('clearUserInfo');
        expect(saveUserInfoToCookie).not.toHaveBeenCalled();
      });
    });
  });

  describe('Actions', () => {
    describe('LOGIN', () => {
      it('should successfully login with valid credentials', async () => {
        const userData = {
          username: 'testuser',
          password: 'testpass'
        };

        const mockResponse = {
          data: {
            access_token: 'jwt-token-123',
            user_info: {
              id: 1,
              username: 'testuser',
              is_superuser: false
            }
          }
        };

        loginUser.mockResolvedValue(mockResponse);

        await store.dispatch('LOGIN', userData);

        expect(loginUser).toHaveBeenCalledWith(userData);
        expect(store.state.token).toBe('jwt-token-123');
        expect(store.state.userInfo).toEqual(mockResponse.data.user_info);
        expect(saveAuthToCookie).toHaveBeenCalledWith('jwt-token-123');
      });

      it('should handle superuser login', async () => {
        const mockResponse = {
          data: {
            access_token: 'admin-token',
            user_info: {
              id: 1,
              username: 'admin',
              is_superuser: true
            }
          }
        };

        loginUser.mockResolvedValue(mockResponse);

        await store.dispatch('LOGIN', {});

        expect(store.state.userInfo.is_superuser).toBe(true);
        expect(store.getters.isSuperUser).toBe(true);
      });

      it('should update state and call cookie functions', async () => {
        const mockResponse = {
          data: {
            access_token: 'token',
            user_info: { username: 'user', is_superuser: false }
          }
        };

        loginUser.mockResolvedValue(mockResponse);

        await store.dispatch('LOGIN', {});

        // Verify state was updated
        expect(store.state.token).toBe('token');
        expect(store.state.userInfo).toEqual({ username: 'user', is_superuser: false });

        // Verify cookie function was called
        expect(saveAuthToCookie).toHaveBeenCalledWith('token');
        expect(saveUserInfoToCookie).toHaveBeenCalledWith({ username: 'user', is_superuser: false });
      });

      it('should propagate API errors', async () => {
        const apiError = new Error('Invalid credentials');
        loginUser.mockRejectedValue(apiError);

        await expect(
          store.dispatch('LOGIN', { username: 'wrong', password: 'wrong' })
        ).rejects.toThrow('Invalid credentials');
      });

      it('should not modify state on login failure', async () => {
        const originalToken = store.state.token;
        const originalUserInfo = store.state.userInfo;

        loginUser.mockRejectedValue(new Error('Login failed'));

        try {
          await store.dispatch('LOGIN', {});
        } catch (e) {
          // Expected to fail
        }

        expect(store.state.token).toBe(originalToken);
        expect(store.state.userInfo).toEqual(originalUserInfo);
      });
    });

    describe('LOGOUT', () => {
      it('should clear token and userInfo', () => {
        // Set up logged-in state
        store.state.token = 'active-token';
        store.state.userInfo = {
          id: 1,
          username: 'testuser',
          is_superuser: true
        };

        store.dispatch('LOGOUT');

        expect(store.state.token).toBe('');
        expect(store.state.userInfo).toEqual({ is_superuser: false });
      });

      it('should clear both token and userInfo on logout', () => {
        // Set up logged-in state
        store.state.token = 'active-token';
        store.state.userInfo = { id: 1, username: 'test', is_superuser: true };

        store.dispatch('LOGOUT');

        // Verify both were cleared
        expect(store.state.token).toBe('');
        expect(store.state.userInfo).toEqual({ is_superuser: false });
      });

      it('should work when already logged out', () => {
        store.state.token = '';
        store.state.userInfo = { is_superuser: false };

        store.dispatch('LOGOUT');

        expect(store.state.token).toBe('');
        expect(store.state.userInfo).toEqual({ is_superuser: false });
      });

      it('should update isLogin getter to false', () => {
        store.state.token = 'active-token';
        expect(store.getters.isLogin).toBe(true);

        store.dispatch('LOGOUT');

        expect(store.getters.isLogin).toBe(false);
      });

      it('should update isSuperUser getter to false', () => {
        store.state.userInfo = { is_superuser: true };
        expect(store.getters.isSuperUser).toBe(true);

        store.dispatch('LOGOUT');

        expect(store.getters.isSuperUser).toBe(false);
      });
    });
  });

  describe('Integration Scenarios', () => {
    it('should handle complete login → logout cycle', async () => {
      // Initial state - not logged in
      expect(store.getters.isLogin).toBe(false);
      expect(store.getters.isSuperUser).toBe(false);

      // Login
      const mockResponse = {
        data: {
          access_token: 'session-token',
          user_info: {
            id: 1,
            username: 'testuser',
            is_superuser: true
          }
        }
      };
      loginUser.mockResolvedValue(mockResponse);

      await store.dispatch('LOGIN', { username: 'test', password: 'pass' });

      // Verify logged in state
      expect(store.getters.isLogin).toBe(true);
      expect(store.getters.isSuperUser).toBe(true);
      expect(store.state.token).toBe('session-token');

      // Logout
      store.dispatch('LOGOUT');

      // Verify logged out state
      expect(store.getters.isLogin).toBe(false);
      expect(store.getters.isSuperUser).toBe(false);
      expect(store.state.token).toBe('');
    });

    it('should handle multiple login attempts', async () => {
      // First login
      loginUser.mockResolvedValueOnce({
        data: {
          access_token: 'token1',
          user_info: { username: 'user1', is_superuser: false }
        }
      });

      await store.dispatch('LOGIN', { username: 'user1' });
      expect(store.state.token).toBe('token1');
      expect(store.state.userInfo.username).toBe('user1');

      // Second login (different user)
      loginUser.mockResolvedValueOnce({
        data: {
          access_token: 'token2',
          user_info: { username: 'user2', is_superuser: true }
        }
      });

      await store.dispatch('LOGIN', { username: 'user2' });
      expect(store.state.token).toBe('token2');
      expect(store.state.userInfo.username).toBe('user2');
      expect(store.getters.isSuperUser).toBe(true);
    });

    it('should maintain state consistency after failed login followed by successful login', async () => {
      // Failed login
      loginUser.mockRejectedValueOnce(new Error('Failed'));

      try {
        await store.dispatch('LOGIN', { username: 'wrong' });
      } catch (e) {
        // Expected
      }

      expect(store.getters.isLogin).toBe(false);

      // Successful login
      loginUser.mockResolvedValueOnce({
        data: {
          access_token: 'success-token',
          user_info: { username: 'correct' }
        }
      });

      await store.dispatch('LOGIN', { username: 'correct' });

      expect(store.getters.isLogin).toBe(true);
      expect(store.state.token).toBe('success-token');
    });
  });
});
