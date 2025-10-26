import { describe, it, expect, beforeEach, vi } from 'vitest';
import {
  saveAuthToCookie,
  getAuthFromCookie,
  saveUserInfoToCookie,
  getUserInfoFromCookie,
  deleteCookie
} from '@/utils/cookies';

describe('cookies.js', () => {
  beforeEach(() => {
    // Clear cookies before each test
    Object.defineProperty(document, 'cookie', {
      writable: true,
      value: ''
    });
  });

  describe('saveAuthToCookie and getAuthFromCookie', () => {
    it('should save and retrieve auth token', () => {
      const token = 'test-auth-token-12345';
      saveAuthToCookie(token);

      const retrieved = getAuthFromCookie();
      expect(retrieved).toBe(token);
    });

    it('should return empty string when no auth cookie exists', () => {
      const retrieved = getAuthFromCookie();
      expect(retrieved).toBe('');
    });

    it('should handle auth token with special characters', () => {
      const token = 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJzdWIiOiIxMjM0NTY3ODkwIiwibmFtZSI6IkpvaG4gRG9lIiwiaWF0IjoxNTE2MjM5MDIyfQ';
      saveAuthToCookie(token);

      const retrieved = getAuthFromCookie();
      expect(retrieved).toBe(token);
    });
  });

  describe('saveUserInfoToCookie and getUserInfoFromCookie', () => {
    it('should save and retrieve user info object', () => {
      const userInfo = {
        username: 'testuser',
        email: 'test@example.com',
        role: 'admin'
      };

      saveUserInfoToCookie(userInfo);
      const retrieved = getUserInfoFromCookie();

      expect(retrieved).toEqual(userInfo);
    });

    it('should handle user info with nested objects', () => {
      const userInfo = {
        username: 'testuser',
        profile: {
          firstName: 'Test',
          lastName: 'User',
          preferences: { theme: 'dark' }
        }
      };

      saveUserInfoToCookie(userInfo);
      const retrieved = getUserInfoFromCookie();

      expect(retrieved).toEqual(userInfo);
    });

    it('should return null when no user info cookie exists', () => {
      const retrieved = getUserInfoFromCookie();
      expect(retrieved).toBeNull();
    });

    it('should handle malformed JSON gracefully', () => {
      // Manually set malformed cookie
      Object.defineProperty(document, 'cookie', {
        writable: true,
        value: 'scap_user={invalid-json'
      });

      const retrieved = getUserInfoFromCookie();
      expect(retrieved).toBeNull();
    });
  });

  describe('deleteCookie', () => {
    it('should delete auth cookie', () => {
      saveAuthToCookie('test-token');
      expect(getAuthFromCookie()).toBe('test-token');

      deleteCookie();

      // After deletion, cookie should be expired (empty string)
      const retrieved = getAuthFromCookie();
      expect(retrieved).toBe('');
    });

    it('should delete user info cookie', () => {
      const userInfo = { username: 'testuser' };
      saveUserInfoToCookie(userInfo);
      expect(getUserInfoFromCookie()).toEqual(userInfo);

      deleteCookie();

      // After deletion, cookie should be null
      const retrieved = getUserInfoFromCookie();
      expect(retrieved).toBeNull();
    });

    it('should delete both cookies simultaneously', () => {
      saveAuthToCookie('test-token');
      saveUserInfoToCookie({ username: 'testuser' });

      deleteCookie();

      expect(getAuthFromCookie()).toBe('');
      expect(getUserInfoFromCookie()).toBeNull();
    });
  });

  describe('edge cases', () => {
    it('should handle empty string auth token', () => {
      saveAuthToCookie('');
      const retrieved = getAuthFromCookie();
      expect(retrieved).toBe('');
    });

    it('should handle empty object user info', () => {
      saveUserInfoToCookie({});
      const retrieved = getUserInfoFromCookie();
      expect(retrieved).toEqual({});
    });

    it('should handle null user info gracefully', () => {
      saveUserInfoToCookie(null);
      const retrieved = getUserInfoFromCookie();
      // null stringifies to 'null' then parses back
      expect(retrieved).toBe(null);
    });
  });
});
