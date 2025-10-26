import { describe, it, expect, beforeEach } from 'vitest';
import {
  formatTimestamp,
  getTaskShortId,
  getBaseName,
  generateAllLogsFilename,
  generateLogFilename,
  generateManifestFilename,
  extractFileExtension,
  generateUploadFileName
} from '@/utils/filename';

describe('filename.js', () => {
  let fixedDate;

  beforeEach(() => {
    // Create a fixed date for consistent testing
    fixedDate = new Date('2025-10-23T12:34:56.789Z');
  });

  /**
   * formatTimestamp Tests
   */
  describe('formatTimestamp', () => {
    it('should format date to ISO string without colons and dots', () => {
      const result = formatTimestamp(fixedDate);

      expect(result).toBe('2025-10-23T12-34-56');
      expect(result).toMatch(/^\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}$/);
    });

    it('should use current date when no date provided', () => {
      const result = formatTimestamp();

      expect(result).toMatch(/^\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}$/);
      expect(typeof result).toBe('string');
    });

    it('should handle different dates correctly', () => {
      const date1 = new Date('2025-01-01T00:00:00.000Z');
      const date2 = new Date('2025-12-31T23:59:59.999Z');

      expect(formatTimestamp(date1)).toBe('2025-01-01T00-00-00');
      expect(formatTimestamp(date2)).toBe('2025-12-31T23-59-59');
    });

    it('should not include milliseconds', () => {
      const result = formatTimestamp(fixedDate);

      expect(result).not.toContain('.');
      expect(result).not.toContain('789');
    });
  });

  /**
   * getTaskShortId Tests
   */
  describe('getTaskShortId', () => {
    it('should extract first 8 characters by default', () => {
      const taskId = 'abc123def456ghi789';
      const result = getTaskShortId(taskId);

      expect(result).toBe('abc123de');
      expect(result.length).toBe(8);
    });

    it('should extract specified number of characters', () => {
      const taskId = 'abc123def456ghi789';

      expect(getTaskShortId(taskId, 4)).toBe('abc1');
      expect(getTaskShortId(taskId, 12)).toBe('abc123def456');
      expect(getTaskShortId(taskId, 20)).toBe(taskId); // Full ID when length > ID length
    });

    it('should handle short task IDs', () => {
      expect(getTaskShortId('abc', 8)).toBe('abc');
      expect(getTaskShortId('ab', 8)).toBe('ab');
    });

    it('should handle empty string', () => {
      expect(getTaskShortId('', 8)).toBe('');
    });

    it('should handle null and undefined gracefully', () => {
      expect(getTaskShortId(null, 8)).toBe('');
      expect(getTaskShortId(undefined, 8)).toBe('');
    });

    it('should handle non-string input', () => {
      expect(getTaskShortId(123, 8)).toBe('');
      expect(getTaskShortId({}, 8)).toBe('');
    });
  });

  /**
   * getBaseName Tests
   */
  describe('getBaseName', () => {
    it('should remove file extension', () => {
      expect(getBaseName('snakemake.log')).toBe('snakemake');
      expect(getBaseName('error.txt')).toBe('error');
      expect(getBaseName('data.json')).toBe('data');
    });

    it('should handle filenames with multiple dots', () => {
      expect(getBaseName('my.file.name.log')).toBe('my.file.name');
      expect(getBaseName('archive.tar.gz')).toBe('archive.tar');
    });

    it('should return filename if no extension', () => {
      expect(getBaseName('README')).toBe('README');
      expect(getBaseName('logfile')).toBe('logfile');
    });

    it('should handle empty string', () => {
      expect(getBaseName('')).toBe('');
    });

    it('should handle null and undefined gracefully', () => {
      expect(getBaseName(null)).toBe('');
      expect(getBaseName(undefined)).toBe('');
    });

    it('should handle filenames starting with dot', () => {
      expect(getBaseName('.gitignore')).toBe('.gitignore');
      expect(getBaseName('.env.local')).toBe('.env');
    });
  });

  /**
   * generateAllLogsFilename Tests
   */
  describe('generateAllLogsFilename', () => {
    it('should generate correct filename for all logs', () => {
      const taskId = 'abc123def456ghi789';
      const result = generateAllLogsFilename(taskId, fixedDate);

      expect(result).toBe('task_abc123de_logs_2025-10-23T12-34-56.json');
    });

    it('should use current date when not provided', () => {
      const taskId = 'abc123def456';
      const result = generateAllLogsFilename(taskId);

      expect(result).toMatch(/^task_abc123de_logs_\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.json$/);
    });

    it('should handle short task IDs', () => {
      const result = generateAllLogsFilename('abc', fixedDate);

      expect(result).toBe('task_abc_logs_2025-10-23T12-34-56.json');
    });

    it('should format with consistent structure', () => {
      const taskId = 'test1234567890';
      const result = generateAllLogsFilename(taskId, fixedDate);

      expect(result).toContain('task_');
      expect(result).toContain('_logs_');
      expect(result).toMatch(/\.json$/);
    });
  });

  /**
   * generateLogFilename Tests
   */
  describe('generateLogFilename', () => {
    it('should generate correct filename for single log with default extension', () => {
      const taskId = 'abc123def456ghi789';
      const logFilename = 'snakemake.log';
      const result = generateLogFilename(taskId, logFilename, fixedDate);

      expect(result).toBe('task_abc123de_snakemake_2025-10-23T12-34-56.txt');
    });

    it('should use specified extension', () => {
      const taskId = 'abc123def456';
      const logFilename = 'error.log';
      const result = generateLogFilename(taskId, logFilename, fixedDate, 'log');

      expect(result).toBe('task_abc123de_error_2025-10-23T12-34-56.log');
    });

    it('should use current date when not provided', () => {
      const taskId = 'abc123def456';
      const logFilename = 'output.log';
      const result = generateLogFilename(taskId, logFilename);

      expect(result).toMatch(/^task_abc123de_output_\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.txt$/);
    });

    it('should handle log filenames with multiple extensions', () => {
      const taskId = 'abc123def456';
      const logFilename = 'archive.tar.gz';
      const result = generateLogFilename(taskId, logFilename, fixedDate);

      expect(result).toBe('task_abc123de_archive.tar_2025-10-23T12-34-56.txt');
    });

    it('should handle log filename without extension', () => {
      const taskId = 'abc123def456';
      const logFilename = 'stdout';
      const result = generateLogFilename(taskId, logFilename, fixedDate);

      expect(result).toBe('task_abc123de_stdout_2025-10-23T12-34-56.txt');
    });

    it('should format with consistent structure', () => {
      const taskId = 'test1234567890';
      const logFilename = 'log.txt';
      const result = generateLogFilename(taskId, logFilename, fixedDate);

      expect(result).toContain('task_');
      expect(result).toContain('_log_');
      expect(result).toMatch(/_\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\./);
    });
  });

  /**
   * generateManifestFilename Tests
   */
  describe('generateManifestFilename', () => {
    it('should generate correct filename for execution manifest', () => {
      const taskId = 'abc123def456ghi789';
      const result = generateManifestFilename(taskId, fixedDate);

      expect(result).toBe('execution_manifest_abc123de_2025-10-23T12-34-56.json');
    });

    it('should use current date when not provided', () => {
      const taskId = 'abc123def456';
      const result = generateManifestFilename(taskId);

      expect(result).toMatch(/^execution_manifest_abc123de_\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.json$/);
    });

    it('should handle short task IDs', () => {
      const result = generateManifestFilename('xyz', fixedDate);

      expect(result).toBe('execution_manifest_xyz_2025-10-23T12-34-56.json');
    });

    it('should format with consistent structure', () => {
      const taskId = 'test1234567890';
      const result = generateManifestFilename(taskId, fixedDate);

      expect(result).toContain('execution_manifest_');
      expect(result).toMatch(/\.json$/);
      expect(result).toMatch(/_\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.json$/);
    });
  });

  /**
   * Integration Tests - Test combined functionality
   */
  describe('Integration Tests', () => {
    it('should generate unique filenames for same task at different times', () => {
      const taskId = 'abc123def456';
      const date1 = new Date('2025-10-23T12:00:00.000Z');
      const date2 = new Date('2025-10-23T13:00:00.000Z');

      const filename1 = generateAllLogsFilename(taskId, date1);
      const filename2 = generateAllLogsFilename(taskId, date2);

      expect(filename1).not.toBe(filename2);
      expect(filename1).toContain('12-00-00');
      expect(filename2).toContain('13-00-00');
    });

    it('should generate consistent filenames for same inputs', () => {
      const taskId = 'abc123def456';
      const logFilename = 'test.log';

      const result1 = generateLogFilename(taskId, logFilename, fixedDate);
      const result2 = generateLogFilename(taskId, logFilename, fixedDate);

      expect(result1).toBe(result2);
    });

    it('should work with real-world task ID patterns', () => {
      // UUID-like task ID
      const taskId = '550e8400-e29b-41d4-a716-446655440000';

      const allLogs = generateAllLogsFilename(taskId, fixedDate);
      const singleLog = generateLogFilename(taskId, 'snakemake.log', fixedDate);
      const manifest = generateManifestFilename(taskId, fixedDate);

      expect(allLogs).toContain('550e8400');
      expect(singleLog).toContain('550e8400');
      expect(manifest).toContain('550e8400');
    });
  });

  /**
   * extractFileExtension Tests
   */
  describe('extractFileExtension', () => {
    it('should extract extension from filename', () => {
      expect(extractFileExtension('data.h5ad')).toBe('h5ad');
      expect(extractFileExtension('experiment.csv')).toBe('csv');
      expect(extractFileExtension('notes.txt')).toBe('txt');
    });

    it('should return lowercase extension', () => {
      expect(extractFileExtension('DATA.H5AD')).toBe('h5ad');
      expect(extractFileExtension('File.CSV')).toBe('csv');
      expect(extractFileExtension('Document.TXT')).toBe('txt');
    });

    it('should handle filenames with multiple dots', () => {
      expect(extractFileExtension('my.data.file.h5ad')).toBe('h5ad');
      expect(extractFileExtension('backup.2025.10.23.csv')).toBe('csv');
    });

    it('should return empty string for filenames without extension', () => {
      expect(extractFileExtension('filename')).toBe('');
      expect(extractFileExtension('README')).toBe('');
    });

    it('should handle edge cases gracefully', () => {
      expect(extractFileExtension(null)).toBe('');
      expect(extractFileExtension(undefined)).toBe('');
      expect(extractFileExtension('')).toBe('');
    });

    it('should handle files starting with dot', () => {
      expect(extractFileExtension('.gitignore')).toBe('gitignore');
      expect(extractFileExtension('.env.production')).toBe('production');
    });
  });

  /**
   * generateUploadFileName Tests
   */
  describe('generateUploadFileName', () => {
    it('should generate filename with folder prefix', () => {
      expect(generateUploadFileName('data', 'experiment.h5ad')).toBe('data_experiment.h5ad');
      expect(generateUploadFileName('results', 'output.csv')).toBe('results_output.csv');
      expect(generateUploadFileName('logs', 'error.txt')).toBe('logs_error.txt');
    });

    it('should handle different folder names', () => {
      expect(generateUploadFileName('user_data', 'file.h5ad')).toBe('user_data_file.h5ad');
      expect(generateUploadFileName('temp', 'upload.csv')).toBe('temp_upload.csv');
    });

    it('should preserve original filename extension', () => {
      const result = generateUploadFileName('data', 'my.complex.file.h5ad');
      expect(result).toBe('data_my.complex.file.h5ad');
      expect(result).toContain('.h5ad');
    });

    it('should handle edge cases gracefully', () => {
      expect(generateUploadFileName(null, 'file.h5ad')).toBe('file.h5ad');
      expect(generateUploadFileName(undefined, 'file.csv')).toBe('file.csv');
      expect(generateUploadFileName('', 'file.txt')).toBe('file.txt');
      expect(generateUploadFileName('folder', null)).toBe('');
      expect(generateUploadFileName('folder', undefined)).toBe('');
      expect(generateUploadFileName('folder', '')).toBe('');
    });

    it('should handle special characters in folder names', () => {
      expect(generateUploadFileName('data-2025', 'file.h5ad')).toBe('data-2025_file.h5ad');
      expect(generateUploadFileName('user.backup', 'file.csv')).toBe('user.backup_file.csv');
    });

    it('should not add double prefix if already present', () => {
      const result = generateUploadFileName('data', 'data_file.h5ad');
      expect(result).toBe('data_data_file.h5ad');
    });
  });
});
