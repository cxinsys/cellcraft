import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import {
  formatBytes,
  formatFileSize,
  getTimeDifference,
  getRunningTime,
  formatDateTime,
  formatTitle,
  cutDateFromISO,
  extractFileName,
  extractExtension,
  formatRelativeTime
} from '@/utils/formatters';

describe('formatters.js', () => {
  /**
   * formatBytes Tests
   */
  describe('formatBytes', () => {
    it('should return "0 Bytes" for 0 bytes', () => {
      expect(formatBytes(0)).toBe('0 Bytes');
      expect(formatBytes(0, 0)).toBe('0 Bytes');
      expect(formatBytes(0, 3)).toBe('0 Bytes');
    });

    it('should format bytes correctly with default decimals', () => {
      expect(formatBytes(500)).toBe('500 Bytes');
      expect(formatBytes(1023)).toBe('1023 Bytes');
    });

    it('should format KB correctly', () => {
      expect(formatBytes(1024)).toBe('1 KB');
      expect(formatBytes(2048)).toBe('2 KB');
      expect(formatBytes(1536)).toBe('1.5 KB');
    });

    it('should format MB correctly', () => {
      expect(formatBytes(1048576)).toBe('1 MB');
      expect(formatBytes(2097152)).toBe('2 MB');
      expect(formatBytes(1572864)).toBe('1.5 MB');
    });

    it('should format GB correctly', () => {
      expect(formatBytes(1073741824)).toBe('1 GB');
      expect(formatBytes(2147483648)).toBe('2 GB');
      expect(formatBytes(1610612736)).toBe('1.5 GB');
    });

    it('should format TB correctly', () => {
      expect(formatBytes(1099511627776)).toBe('1 TB');
      expect(formatBytes(2199023255552)).toBe('2 TB');
      expect(formatBytes(1649267441664)).toBe('1.5 TB');
    });

    it('should handle custom decimals parameter', () => {
      // 0 decimals
      expect(formatBytes(1536, 0)).toBe('2 KB');
      expect(formatBytes(1572864, 0)).toBe('2 MB');

      // 1 decimal
      expect(formatBytes(1536, 1)).toBe('1.5 KB');
      expect(formatBytes(1572864, 1)).toBe('1.5 MB');

      // 3 decimals
      expect(formatBytes(1536, 3)).toBe('1.5 KB');
      expect(formatBytes(512486, 3)).toBe('500.475 KB');
    });

    it('should handle negative decimals as 0', () => {
      expect(formatBytes(1536, -1)).toBe('2 KB');
      expect(formatBytes(1536, -5)).toBe('2 KB');
    });

    it('should handle edge case: 1023 bytes to 1024 bytes boundary', () => {
      expect(formatBytes(1023)).toBe('1023 Bytes');
      expect(formatBytes(1024)).toBe('1 KB');
    });

    it('should format large TB values correctly', () => {
      // 10 TB
      expect(formatBytes(10995116277760)).toBe('10 TB');
      // 100 TB
      expect(formatBytes(109951162777600)).toBe('100 TB');
    });

    it('should maintain precision with different decimals', () => {
      const bytes = 1234567890; // ~1.15 GB

      expect(formatBytes(bytes, 0)).toBe('1 GB');
      expect(formatBytes(bytes, 1)).toBe('1.1 GB');
      expect(formatBytes(bytes, 2)).toBe('1.15 GB');
      expect(formatBytes(bytes, 3)).toBe('1.15 GB');
    });
  });

  /**
   * formatFileSize Tests (wrapper for formatBytes)
   */
  describe('formatFileSize', () => {
    it('should return "0 Bytes" for 0 bytes', () => {
      expect(formatFileSize(0)).toBe('0 Bytes');
    });

    it('should format bytes correctly', () => {
      expect(formatFileSize(500)).toBe('500 Bytes');
      expect(formatFileSize(1023)).toBe('1023 Bytes');
    });

    it('should format KB correctly', () => {
      expect(formatFileSize(1024)).toBe('1 KB');
      expect(formatFileSize(2048)).toBe('2 KB');
      expect(formatFileSize(1536)).toBe('1.5 KB');
    });

    it('should format MB correctly', () => {
      expect(formatFileSize(1048576)).toBe('1 MB');
      expect(formatFileSize(2097152)).toBe('2 MB');
      expect(formatFileSize(1572864)).toBe('1.5 MB');
    });

    it('should format GB correctly', () => {
      expect(formatFileSize(1073741824)).toBe('1 GB');
      expect(formatFileSize(2147483648)).toBe('2 GB');
    });

    it('should use 2 decimals by default (wrapper behavior)', () => {
      // formatFileSize should call formatBytes with decimals=2
      expect(formatFileSize(1234567890)).toBe(formatBytes(1234567890, 2));
    });

    it('should support TB through formatBytes', () => {
      // Even though old formatFileSize didn't have TB, the new wrapper supports it
      expect(formatFileSize(1099511627776)).toBe('1 TB');
    });
  });

  describe('getTimeDifference', () => {
    it('should calculate time difference correctly for same time', () => {
      const time = new Date('2024-01-01T12:00:00');
      expect(getTimeDifference(time, time)).toBe('00:00:00');
    });

    it('should calculate time difference for seconds', () => {
      const start = new Date('2024-01-01T12:00:00');
      const end = new Date('2024-01-01T12:00:30');
      expect(getTimeDifference(start, end)).toBe('00:00:30');
    });

    it('should calculate time difference for minutes', () => {
      const start = new Date('2024-01-01T12:00:00');
      const end = new Date('2024-01-01T12:05:00');
      expect(getTimeDifference(start, end)).toBe('00:05:00');
    });

    it('should calculate time difference for hours', () => {
      const start = new Date('2024-01-01T12:00:00');
      const end = new Date('2024-01-01T14:30:45');
      expect(getTimeDifference(start, end)).toBe('02:30:45');
    });

    it('should handle string input dates', () => {
      const start = '2024-01-01T12:00:00';
      const end = '2024-01-01T13:15:30';
      expect(getTimeDifference(start, end)).toBe('01:15:30');
    });

    it('should use absolute difference (order-independent)', () => {
      const start = new Date('2024-01-01T14:00:00');
      const end = new Date('2024-01-01T12:00:00');
      expect(getTimeDifference(start, end)).toBe('02:00:00');
      expect(getTimeDifference(end, start)).toBe('02:00:00');
    });
  });

  describe('getRunningTime', () => {
    it('should calculate running time from start to current', () => {
      const start = new Date('2024-01-01T12:00:00');
      const current = new Date('2024-01-01T12:05:30');
      expect(getRunningTime(start, current)).toBe('00:05:30');
    });

    it('should handle string start time', () => {
      const start = '2024-01-01T12:00:00';
      const current = new Date('2024-01-01T13:00:00');
      expect(getRunningTime(start, current)).toBe('01:00:00');
    });

    it('should calculate hours correctly', () => {
      const start = new Date('2024-01-01T10:00:00');
      const current = new Date('2024-01-01T15:30:45');
      expect(getRunningTime(start, current)).toBe('05:30:45');
    });

    it('should use absolute difference', () => {
      const start = new Date('2024-01-01T14:00:00');
      const current = new Date('2024-01-01T12:00:00');
      expect(getRunningTime(start, current)).toBe('02:00:00');
    });
  });

  describe('formatDateTime', () => {
    it('should format valid date correctly', () => {
      const dateTime = '2024-01-15T14:30:00';
      expect(formatDateTime(dateTime)).toBe('2024.01.15-14:30');
    });

    it('should return "Not Yet Completed" for invalid date', () => {
      expect(formatDateTime('invalid-date')).toBe('Not Yet Completed');
      expect(formatDateTime(null)).toBe('Not Yet Completed');
      expect(formatDateTime(undefined)).toBe('Not Yet Completed');
    });

    it('should handle Date objects', () => {
      const date = new Date('2024-03-20T09:15:00');
      expect(formatDateTime(date)).toBe('2024.03.20-09:15');
    });

    it('should format different times correctly', () => {
      expect(formatDateTime('2024-12-31T23:59:00')).toBe('2024.12.31-23:59');
      expect(formatDateTime('2024-01-01T00:00:00')).toBe('2024.01.01-00:00');
    });
  });

  describe('formatTitle', () => {
    it('should return "Untitled" for null', () => {
      expect(formatTitle(null)).toBe('Untitled');
    });

    it('should return "Untitled" for empty string', () => {
      expect(formatTitle('')).toBe('Untitled');
    });

    it('should return "Untitled" for undefined', () => {
      expect(formatTitle(undefined)).toBe('Untitled');
    });

    it('should return the title as-is for valid strings', () => {
      expect(formatTitle('My Workflow')).toBe('My Workflow');
      expect(formatTitle('Test 123')).toBe('Test 123');
    });

    it('should handle whitespace-only strings', () => {
      // Current implementation treats whitespace as valid
      expect(formatTitle('   ')).toBe('   ');
    });
  });

  /**
   * cutDateFromISO Tests
   */
  describe('cutDateFromISO', () => {
    it('should extract date from ISO string', () => {
      expect(cutDateFromISO('2025-10-23T12:34:56')).toBe('2025-10-23');
      expect(cutDateFromISO('2024-01-15T00:00:00')).toBe('2024-01-15');
      expect(cutDateFromISO('2024-12-31T23:59:59')).toBe('2024-12-31');
    });

    it('should handle ISO strings with timezone', () => {
      expect(cutDateFromISO('2025-10-23T12:34:56Z')).toBe('2025-10-23');
      expect(cutDateFromISO('2025-10-23T12:34:56+09:00')).toBe('2025-10-23');
      expect(cutDateFromISO('2025-10-23T12:34:56-05:00')).toBe('2025-10-23');
    });

    it('should handle ISO strings with milliseconds', () => {
      expect(cutDateFromISO('2025-10-23T12:34:56.789')).toBe('2025-10-23');
      expect(cutDateFromISO('2025-10-23T12:34:56.123Z')).toBe('2025-10-23');
    });

    it('should handle edge cases gracefully', () => {
      expect(cutDateFromISO(null)).toBe('');
      expect(cutDateFromISO(undefined)).toBe('');
      expect(cutDateFromISO('')).toBe('');
    });

    it('should handle strings without T separator', () => {
      expect(cutDateFromISO('2025-10-23')).toBe('2025-10-23');
      expect(cutDateFromISO('invalid-date')).toBe('invalid-date');
    });

    it('should handle non-string inputs', () => {
      expect(cutDateFromISO(123)).toBe('');
      expect(cutDateFromISO({})).toBe('');
      expect(cutDateFromISO([])).toBe('');
    });
  });

  /**
   * extractFileName Tests
   */
  describe('extractFileName', () => {
    it('should extract filename without extension', () => {
      expect(extractFileName('data.h5ad')).toBe('data');
      expect(extractFileName('experiment.csv')).toBe('experiment');
      expect(extractFileName('notes.txt')).toBe('notes');
    });

    it('should handle filenames with multiple dots', () => {
      expect(extractFileName('my.data.file.h5ad')).toBe('my.data.file');
      expect(extractFileName('backup.2025.10.23.csv')).toBe('backup.2025.10.23');
      expect(extractFileName('project.v1.0.txt')).toBe('project.v1.0');
    });

    it('should return original filename if no extension', () => {
      expect(extractFileName('filename')).toBe('filename');
      expect(extractFileName('README')).toBe('README');
    });

    it('should handle edge cases gracefully', () => {
      expect(extractFileName(null)).toBe('');
      expect(extractFileName(undefined)).toBe('');
      expect(extractFileName('')).toBe('');
    });

    it('should handle files starting with dot', () => {
      expect(extractFileName('.gitignore')).toBe('');
      expect(extractFileName('.env.production')).toBe('.env');
    });

    it('should handle non-string inputs', () => {
      expect(extractFileName(123)).toBe('');
      expect(extractFileName({})).toBe('');
      expect(extractFileName([])).toBe('');
    });

    it('should preserve case', () => {
      expect(extractFileName('MyFile.H5AD')).toBe('MyFile');
      expect(extractFileName('DATA.CSV')).toBe('DATA');
    });
  });

  /**
   * extractExtension Tests
   */
  describe('extractExtension', () => {
    it('should extract file extension', () => {
      expect(extractExtension('data.h5ad')).toBe('h5ad');
      expect(extractExtension('experiment.csv')).toBe('csv');
      expect(extractExtension('notes.txt')).toBe('txt');
    });

    it('should handle filenames with multiple dots', () => {
      expect(extractExtension('my.data.file.h5ad')).toBe('h5ad');
      expect(extractExtension('backup.2025.10.23.csv')).toBe('csv');
    });

    it('should return empty string for filenames without extension', () => {
      expect(extractExtension('filename')).toBe('');
      expect(extractExtension('README')).toBe('');
    });

    it('should handle edge cases gracefully', () => {
      expect(extractExtension(null)).toBe('');
      expect(extractExtension(undefined)).toBe('');
      expect(extractExtension('')).toBe('');
    });

    it('should handle files starting with dot', () => {
      expect(extractExtension('.gitignore')).toBe('gitignore');
      expect(extractExtension('.env.production')).toBe('production');
    });

    it('should preserve case', () => {
      expect(extractExtension('file.H5AD')).toBe('H5AD');
      expect(extractExtension('data.CSV')).toBe('CSV');
      expect(extractExtension('notes.TXT')).toBe('TXT');
    });

    it('should handle non-string inputs', () => {
      expect(extractExtension(123)).toBe('');
      expect(extractExtension({})).toBe('');
      expect(extractExtension([])).toBe('');
    });

    it('should work consistently with extractFileName', () => {
      const filename = 'my.project.data.h5ad';
      const name = extractFileName(filename);
      const ext = extractExtension(filename);
      expect(`${name}.${ext}`).toBe(filename);
    });
  });

  /**
   * formatRelativeTime Tests
   */
  describe('formatRelativeTime', () => {
    it('should return "Edited Recently" for timestamps within the last hour', () => {
      const now = new Date('2025-10-24T12:00:00');
      const recent = new Date('2025-10-24T11:30:00'); // 30 minutes ago

      expect(formatRelativeTime(recent, now)).toBe('Edited Recently');
    });

    it('should return "Edited 1 hour ago" for exactly 1 hour', () => {
      const now = new Date('2025-10-24T12:00:00');
      const oneHourAgo = new Date('2025-10-24T11:00:00');

      expect(formatRelativeTime(oneHourAgo, now)).toBe('Edited 1 hour ago');
    });

    it('should return "Edited X hours ago" for multiple hours', () => {
      const now = new Date('2025-10-24T12:00:00');
      const fiveHoursAgo = new Date('2025-10-24T07:00:00');
      const tenHoursAgo = new Date('2025-10-24T02:00:00');

      expect(formatRelativeTime(fiveHoursAgo, now)).toBe('Edited 5 hours ago');
      expect(formatRelativeTime(tenHoursAgo, now)).toBe('Edited 10 hours ago');
    });

    it('should handle edge cases gracefully', () => {
      expect(formatRelativeTime(null)).toBe('');
      expect(formatRelativeTime(undefined)).toBe('');
      expect(formatRelativeTime('')).toBe('');
    });

    it('should return empty string for invalid timestamps', () => {
      const now = new Date('2025-10-24T12:00:00');

      expect(formatRelativeTime('invalid-date', now)).toBe('');
      expect(formatRelativeTime('not a date', now)).toBe('');
    });

    it('should use current time when currentTime is not provided', () => {
      const recent = new Date(Date.now() - 30 * 60 * 1000); // 30 minutes ago
      const result = formatRelativeTime(recent);

      expect(result).toBe('Edited Recently');
    });

    it('should handle ISO string timestamps', () => {
      const now = new Date('2025-10-24T12:00:00');
      const isoString = '2025-10-24T08:00:00';

      expect(formatRelativeTime(isoString, now)).toBe('Edited 4 hours ago');
    });

    it('should handle Date objects', () => {
      const now = new Date('2025-10-24T12:00:00');
      const pastDate = new Date('2025-10-24T06:00:00');

      expect(formatRelativeTime(pastDate, now)).toBe('Edited 6 hours ago');
    });

    it('should handle future timestamps gracefully', () => {
      const now = new Date('2025-10-24T10:00:00').getTime();
      const future = new Date('2025-10-24T15:00:00');  // 5 hours in the future

      const result = formatRelativeTime(future, now);

      expect(result).toBe('Edited Recently');
    });

    it('should not show negative hours for future dates', () => {
      const now = new Date('2025-10-24T10:00:00').getTime();
      const futureByOneHour = new Date('2025-10-24T11:00:00');
      const futureByTenHours = new Date('2025-10-24T20:00:00');

      const result1 = formatRelativeTime(futureByOneHour, now);
      const result2 = formatRelativeTime(futureByTenHours, now);

      // Should not contain negative numbers
      expect(result1).not.toContain('-');
      expect(result2).not.toContain('-');

      // Should return "Edited Recently" for future timestamps
      expect(result1).toBe('Edited Recently');
      expect(result2).toBe('Edited Recently');
    });
  });
});
