import { vi, describe, it, expect, beforeEach } from 'vitest';

vi.mock('@/api/index', () => ({
  changePassword: vi.fn(),
  loginUser: vi.fn(),
}));

import store from '@/store/index';
import { passwordChangeGuard, CHANGE_PASSWORD_PATH } from '@/router/index';

describe('password change store state', () => {
  it('sets and clears the must_change_password flag', () => {
    store.commit('setUserInfo', { is_superuser: false, email: 'a@b.com', is_active: true, must_change_password: false });
    expect(store.getters.mustChangePassword).toBe(false);

    store.commit('setPasswordChangeRequired');
    expect(store.getters.mustChangePassword).toBe(true);
    expect(store.state.userInfo.email).toBe('a@b.com');

    store.commit('setPasswordChanged');
    expect(store.getters.mustChangePassword).toBe(false);
  });
});

describe('passwordChangeGuard', () => {
  beforeEach(() => {
    store.commit('setToken', 'token');
    store.commit('setPasswordChangeRequired');
  });

  it('redirects to the change page while a change is required', () => {
    const next = vi.fn();
    passwordChangeGuard({ path: '/projects' }, {}, next);
    expect(next).toHaveBeenCalledWith(CHANGE_PASSWORD_PATH);
  });

  it.each([CHANGE_PASSWORD_PATH, '/login'])('allows %s', (path) => {
    const next = vi.fn();
    passwordChangeGuard({ path }, {}, next);
    expect(next).toHaveBeenCalledWith();
  });

  it('does not affect users without the flag', () => {
    store.commit('setPasswordChanged');
    const next = vi.fn();
    passwordChangeGuard({ path: '/projects' }, {}, next);
    expect(next).toHaveBeenCalledWith();
  });
});
