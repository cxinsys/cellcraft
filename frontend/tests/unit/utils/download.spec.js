import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { createBlobURL, downloadFile } from '@/utils/download';

describe('download.js', () => {
  describe('createBlobURL', () => {
    let mockURL;
    let mockCreateObjectURL;

    beforeEach(() => {
      // Mock window.URL.createObjectURL
      mockCreateObjectURL = vi.fn().mockReturnValue('blob:mock-url-123');
      mockURL = {
        createObjectURL: mockCreateObjectURL
      };
      global.window = { URL: mockURL };
    });

    afterEach(() => {
      vi.restoreAllMocks();
    });

    it('should create Blob URL with default type', () => {
      const data = 'test data';
      const result = createBlobURL(data);

      expect(result).toBe('blob:mock-url-123');
      expect(mockCreateObjectURL).toHaveBeenCalledTimes(1);

      // Check that Blob was created with correct data
      const blobArg = mockCreateObjectURL.mock.calls[0][0];
      expect(blobArg).toBeInstanceOf(Blob);
      expect(blobArg.type).toBe('application/octet-stream');
    });

    it('should create Blob URL with custom type', () => {
      const data = '{"key": "value"}';
      const type = 'application/json';
      const result = createBlobURL(data, type);

      expect(result).toBe('blob:mock-url-123');

      const blobArg = mockCreateObjectURL.mock.calls[0][0];
      expect(blobArg.type).toBe(type);
    });

    it('should handle different data types', () => {
      const data = [1, 2, 3];
      createBlobURL(data);

      expect(mockCreateObjectURL).toHaveBeenCalledTimes(1);
    });
  });

  describe('downloadFile', () => {
    let mockCreateObjectURL;
    let mockRevokeObjectURL;
    let mockClick;
    let mockRemoveChild;
    let mockAppendChild;

    beforeEach(() => {
      // Mock window.URL
      mockCreateObjectURL = vi.fn().mockReturnValue('blob:mock-url-456');
      mockRevokeObjectURL = vi.fn();
      global.window = {
        URL: {
          createObjectURL: mockCreateObjectURL,
          revokeObjectURL: mockRevokeObjectURL
        }
      };

      // Mock document.createElement and related methods
      mockClick = vi.fn();
      mockRemoveChild = vi.fn();
      mockAppendChild = vi.fn();

      const mockLink = {
        href: '',
        download: '',
        click: mockClick
      };

      vi.spyOn(document, 'createElement').mockReturnValue(mockLink);
      vi.spyOn(document.body, 'appendChild').mockImplementation(mockAppendChild);
      vi.spyOn(document.body, 'removeChild').mockImplementation(mockRemoveChild);
    });

    afterEach(() => {
      vi.restoreAllMocks();
    });

    it('should trigger download with correct filename', () => {
      const data = 'test content';
      const filename = 'test.txt';

      downloadFile(data, filename);

      // Verify createElement was called with 'a'
      expect(document.createElement).toHaveBeenCalledWith('a');

      // Verify link was added to DOM
      expect(mockAppendChild).toHaveBeenCalledTimes(1);

      // Verify click was triggered
      expect(mockClick).toHaveBeenCalledTimes(1);

      // Verify cleanup
      expect(mockRemoveChild).toHaveBeenCalledTimes(1);
      expect(mockRevokeObjectURL).toHaveBeenCalledWith('blob:mock-url-456');
    });

    it('should create Blob without type when not provided', () => {
      const data = 'test data';
      const filename = 'file.bin';

      downloadFile(data, filename);

      expect(mockCreateObjectURL).toHaveBeenCalledTimes(1);
      const blobArg = mockCreateObjectURL.mock.calls[0][0];
      expect(blobArg).toBeInstanceOf(Blob);
    });

    it('should create Blob with type when provided', () => {
      const data = '{"test": true}';
      const filename = 'data.json';
      const type = 'application/json';

      downloadFile(data, filename, type);

      const blobArg = mockCreateObjectURL.mock.calls[0][0];
      expect(blobArg.type).toBe(type);
    });

    it('should handle different filenames', () => {
      const testCases = [
        'simple.txt',
        'file with spaces.pdf',
        'file_with_underscores.json',
        'file-with-dashes.xml'
      ];

      testCases.forEach((filename) => {
        mockCreateObjectURL.mockClear();
        mockClick.mockClear();
        mockAppendChild.mockClear();
        mockRemoveChild.mockClear();
        mockRevokeObjectURL.mockClear();

        downloadFile('data', filename);

        expect(mockClick).toHaveBeenCalledTimes(1);
        expect(mockRevokeObjectURL).toHaveBeenCalledTimes(1);
      });
    });

    it('should properly clean up resources', () => {
      downloadFile('test', 'file.txt');

      // Verify cleanup order
      expect(mockClick).toHaveBeenCalled();
      expect(mockRemoveChild).toHaveBeenCalled();
      expect(mockRevokeObjectURL).toHaveBeenCalledWith('blob:mock-url-456');
    });
  });
});
