/**
 * Vitest 전역 설정 파일
 * 모든 테스트 실행 전에 자동으로 로드됩니다.
 */

import { config } from '@vue/test-utils';
import { vi } from 'vitest';

/**
 * Vue Test Utils 전역 설정 (Vue 2.x + @vue/test-utils 1.x)
 */
config.mocks = {
  // 라우터 mock
  $router: {
    push: vi.fn(),
    replace: vi.fn(),
    go: vi.fn(),
    back: vi.fn()
  },
  $route: {
    path: '/',
    params: {},
    query: {},
    name: 'home'
  },
  // 스토어 mock (필요 시 개별 테스트에서 override)
  $store: {
    state: {},
    getters: {},
    commit: vi.fn(),
    dispatch: vi.fn()
  }
};

/**
 * 전역 스텁 컴포넌트
 */
config.stubs = {
  // 외부 라이브러리 컴포넌트 스텁
  'router-link': true,
  'router-view': true
};

/**
 * 전역 테스트 유틸리티
 */

// 로컬 스토리지 mock
const localStorageMock = {
  getItem: vi.fn(),
  setItem: vi.fn(),
  removeItem: vi.fn(),
  clear: vi.fn(),
};
global.localStorage = localStorageMock;

// 세션 스토리지 mock
const sessionStorageMock = {
  getItem: vi.fn(),
  setItem: vi.fn(),
  removeItem: vi.fn(),
  clear: vi.fn(),
};
global.sessionStorage = sessionStorageMock;

/**
 * 각 테스트 후 정리
 */
afterEach(() => {
  // Mock 초기화
  vi.clearAllMocks();

  // 스토리지 초기화
  localStorage.clear();
  sessionStorage.clear();
});

/**
 * 테스트 타임아웃 설정
 */
vi.setConfig({ testTimeout: 10000 });

/**
 * 콘솔 에러/경고 필터링 (노이즈 감소)
 */
const originalError = console.error;
const originalWarn = console.warn;

beforeAll(() => {
  console.error = (...args) => {
    // Vue 관련 알려진 경고 무시
    if (
      typeof args[0] === 'string' &&
      (args[0].includes('[Vue warn]') ||
       args[0].includes('Not implemented: HTMLFormElement.prototype.submit'))
    ) {
      return;
    }
    originalError.call(console, ...args);
  };

  console.warn = (...args) => {
    // 테스트 환경 관련 경고 무시
    if (
      typeof args[0] === 'string' &&
      args[0].includes('Avoided redundant navigation')
    ) {
      return;
    }
    originalWarn.call(console, ...args);
  };
});

afterAll(() => {
  console.error = originalError;
  console.warn = originalWarn;
});

/**
 * 공통 테스트 헬퍼 export
 */
export const createMockStore = (overrides = {}) => ({
  state: {
    workflow_info: {
      id: null,
      drawflow: { Home: { data: {} } }
    },
    token: '',
    userInfo: { is_superuser: false },
    ...overrides.state
  },
  getters: {
    isAuthenticated: () => false,
    ...overrides.getters
  },
  mutations: {
    SET_TOKEN: vi.fn(),
    SET_USER_INFO: vi.fn(),
    ...overrides.mutations
  },
  actions: {
    login: vi.fn(),
    logout: vi.fn(),
    ...overrides.actions
  },
  commit: vi.fn(),
  dispatch: vi.fn()
});

export const createMockRouter = (overrides = {}) => ({
  push: vi.fn(),
  replace: vi.fn(),
  go: vi.fn(),
  back: vi.fn(),
  currentRoute: {
    path: '/',
    params: {},
    query: {},
    name: 'home',
    ...overrides.currentRoute
  },
  ...overrides
});

/**
 * API 응답 mock 헬퍼
 */
export const mockApiSuccess = (data) => Promise.resolve({ data });
export const mockApiError = (error) => Promise.reject({ response: { data: error } });

/**
 * 비동기 테스트 헬퍼
 */
export const flushPromises = () => new Promise(resolve => setImmediate(resolve));
