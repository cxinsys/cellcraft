import { describe, it, expect } from 'vitest';
import { validateEmail, validateFileExtension, ALLOWED_FILE_EXTENSIONS } from '@/utils/validation';

describe('validation.js', () => {
  describe('validateEmail - valid emails', () => {
    it('should validate simple email addresses', () => {
      expect(validateEmail('user@example.com')).toBe(true);
      expect(validateEmail('test.user@example.com')).toBe(true);
    });

    it('should validate emails with subdomains', () => {
      expect(validateEmail('user@mail.example.com')).toBe(true);
      expect(validateEmail('admin@dev.company.co.uk')).toBe(true);
    });

    it('should validate emails with numbers and hyphens', () => {
      expect(validateEmail('user123@example.com')).toBe(true);
      expect(validateEmail('test-user@my-company.com')).toBe(true);
    });

    it('should validate emails with plus addressing', () => {
      expect(validateEmail('user+tag@example.com')).toBe(true);
      expect(validateEmail('test+123@domain.org')).toBe(true);
    });
  });

  describe('validateEmail - invalid emails', () => {
    it('should reject emails without @ symbol', () => {
      expect(validateEmail('userexample.com')).toBe(false);
      expect(validateEmail('user.example.com')).toBe(false);
    });

    it('should reject emails without domain', () => {
      expect(validateEmail('user@')).toBe(false);
      expect(validateEmail('@example.com')).toBe(false);
    });

    it('should reject emails with special characters in wrong places', () => {
      expect(validateEmail('user@exam ple.com')).toBe(false);
      expect(validateEmail('us er@example.com')).toBe(false);
      expect(validateEmail('user@@example.com')).toBe(false);
    });
  });

  describe('validateEmail - edge cases', () => {
    it('should handle null and undefined gracefully', () => {
      expect(validateEmail(null)).toBe(false);
      expect(validateEmail(undefined)).toBe(false);
    });

    it('should handle empty strings', () => {
      expect(validateEmail('')).toBe(false);
      expect(validateEmail('   ')).toBe(false);
    });

    it('should be case insensitive', () => {
      expect(validateEmail('User@Example.COM')).toBe(true);
      expect(validateEmail('TEST@DOMAIN.ORG')).toBe(true);
      expect(validateEmail('MixedCase@Test.Co.UK')).toBe(true);
    });
  });

  describe('validateFileExtension - valid file extensions', () => {
    it('should validate h5ad files', () => {
      const result = validateFileExtension('data.h5ad');
      expect(result.isValid).toBe(true);
      expect(result.extension).toBe('h5ad');
      expect(result.message).toBe('');
    });

    it('should validate csv files', () => {
      const result = validateFileExtension('experiment.csv');
      expect(result.isValid).toBe(true);
      expect(result.extension).toBe('csv');
      expect(result.message).toBe('');
    });

    it('should validate txt files', () => {
      const result = validateFileExtension('notes.txt');
      expect(result.isValid).toBe(true);
      expect(result.extension).toBe('txt');
      expect(result.message).toBe('');
    });

    it('should be case insensitive for extensions', () => {
      expect(validateFileExtension('data.H5AD').isValid).toBe(true);
      expect(validateFileExtension('data.CSV').isValid).toBe(true);
      expect(validateFileExtension('data.TXT').isValid).toBe(true);
      expect(validateFileExtension('data.H5ad').extension).toBe('h5ad');
    });

    it('should handle filenames with multiple dots', () => {
      const result = validateFileExtension('my.data.file.h5ad');
      expect(result.isValid).toBe(true);
      expect(result.extension).toBe('h5ad');
    });
  });

  describe('validateFileExtension - invalid file extensions', () => {
    it('should reject unsupported file types', () => {
      const result = validateFileExtension('document.pdf');
      expect(result.isValid).toBe(false);
      expect(result.extension).toBe('pdf');
      expect(result.message).toBe('Please upload h5ad, csv, txt file');
    });

    it('should reject executable files', () => {
      const result = validateFileExtension('script.exe');
      expect(result.isValid).toBe(false);
      expect(result.message).toBe('Please upload h5ad, csv, txt file');
    });

    it('should reject image files', () => {
      expect(validateFileExtension('image.png').isValid).toBe(false);
      expect(validateFileExtension('photo.jpg').isValid).toBe(false);
    });
  });

  describe('validateFileExtension - edge cases', () => {
    it('should handle null and undefined gracefully', () => {
      expect(validateFileExtension(null).isValid).toBe(false);
      expect(validateFileExtension(null).message).toBe('Invalid file name');
      expect(validateFileExtension(undefined).isValid).toBe(false);
    });

    it('should handle empty strings', () => {
      const result = validateFileExtension('');
      expect(result.isValid).toBe(false);
      expect(result.message).toBe('Invalid file name');
    });

    it('should handle filenames without extensions', () => {
      const result = validateFileExtension('filename');
      expect(result.isValid).toBe(false);
      expect(result.extension).toBe('filename');
    });

    it('should allow custom allowed extensions', () => {
      const customExtensions = ['json', 'yaml'];
      const result = validateFileExtension('config.json', customExtensions);
      expect(result.isValid).toBe(true);
      expect(result.extension).toBe('json');
    });

    it('should reject files with custom allowed extensions', () => {
      const customExtensions = ['json', 'yaml'];
      const result = validateFileExtension('data.h5ad', customExtensions);
      expect(result.isValid).toBe(false);
      expect(result.message).toBe('Please upload json, yaml file');
    });
  });

  describe('ALLOWED_FILE_EXTENSIONS constant', () => {
    it('should export the allowed extensions array', () => {
      expect(ALLOWED_FILE_EXTENSIONS).toEqual(['h5ad', 'csv', 'txt']);
    });
  });
});
