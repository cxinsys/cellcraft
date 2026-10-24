import { shallowMount } from '@vue/test-utils';
import { vi, describe, it, expect, beforeEach } from 'vitest';

vi.mock('@/api/index', () => ({
  changePassword: vi.fn(),
  loginUser: vi.fn(),
}));

import { changePassword } from '@/api/index';
import { APIError } from '@/api/common/interceptors';
import ChangePasswordForm, { MIN_PASSWORD_LENGTH } from '@/components/ChangePasswordForm.vue';

function mountForm(getters = {}) {
  const $store = {
    getters: { mustChangePassword: true, isSuperUser: false, ...getters },
    commit: vi.fn(),
  };
  const $router = { push: vi.fn() };
  const wrapper = shallowMount(ChangePasswordForm, { mocks: { $store, $router } });
  return { wrapper, $store, $router };
}

describe('ChangePasswordForm', () => {
  beforeEach(() => {
    changePassword.mockReset();
  });

  it('shows the required notice when a change is required', () => {
    const { wrapper } = mountForm();
    expect(wrapper.find('.login__notice').exists()).toBe(true);
  });

  it.each([
    ['short', 'short', 'current-pass', `at least ${MIN_PASSWORD_LENGTH}`],
    ['new-password-1', 'new-password-2', 'current-pass', 'do not match'],
    ['current-pass', 'current-pass', 'current-pass', 'different'],
  ])('rejects invalid input locally (%s)', async (newPassword, confirmPassword, currentPassword, message) => {
    const { wrapper } = mountForm();
    await wrapper.setData({ currentPassword, newPassword, confirmPassword });
    await wrapper.vm.submitForm();

    expect(changePassword).not.toHaveBeenCalled();
    expect(wrapper.vm.errorMessage).toContain(message);
  });

  it('changes the password, clears the flag and leaves the page', async () => {
    changePassword.mockResolvedValue({ data: {} });
    const { wrapper, $store, $router } = mountForm();
    await wrapper.setData({ currentPassword: 'temp-pass', newPassword: 'new-password-1', confirmPassword: 'new-password-1' });
    await wrapper.vm.submitForm();

    expect(changePassword).toHaveBeenCalledWith('temp-pass', 'new-password-1');
    expect($store.commit).toHaveBeenCalledWith('setPasswordChanged');
    expect($router.push).toHaveBeenCalledWith('/projects');
  });

  it('shows the server error message', async () => {
    const response = { status: 400, data: { detail: 'Incorrect current password' } };
    changePassword.mockRejectedValue(new APIError(response, new Error('Request failed with status code 400')));
    const { wrapper, $store } = mountForm();
    await wrapper.setData({ currentPassword: 'wrong', newPassword: 'new-password-1', confirmPassword: 'new-password-1' });
    await wrapper.vm.submitForm();

    expect(wrapper.vm.errorMessage).toBe('Incorrect current password');
    expect($store.commit).not.toHaveBeenCalled();
  });
});


describe('APIError message', () => {
  it('uses a string detail', () => {
    expect(new APIError({ status: 400, data: { detail: 'Bad input' } }, new Error('x')).message).toBe('Bad input');
  });

  it('uses the message of an object detail', () => {
    const response = { status: 403, data: { detail: { code: 'password_change_required', message: 'Password change required' } } };
    expect(new APIError(response, new Error('x')).message).toBe('Password change required');
  });
});
