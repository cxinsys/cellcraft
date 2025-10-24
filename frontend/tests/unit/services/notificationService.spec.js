import { NotificationService } from '@/services/notificationService';

describe('NotificationService', () => {
  let service;
  let mockDisplayFunction;
  let consoleLogSpy;
  let consoleErrorSpy;
  let consoleWarnSpy;
  let consoleInfoSpy;

  beforeEach(() => {
    mockDisplayFunction = jest.fn();
    consoleLogSpy = jest.spyOn(console, 'log').mockImplementation(() => {});
    consoleErrorSpy = jest.spyOn(console, 'error').mockImplementation(() => {});
    consoleWarnSpy = jest.spyOn(console, 'warn').mockImplementation(() => {});
    consoleInfoSpy = jest.spyOn(console, 'info').mockImplementation(() => {});
  });

  afterEach(() => {
    consoleLogSpy.mockRestore();
    consoleErrorSpy.mockRestore();
    consoleWarnSpy.mockRestore();
    consoleInfoSpy.mockRestore();
  });

  describe('constructor', () => {
    it('should create service with default options', () => {
      service = new NotificationService();

      expect(service.enableLogging).toBe(true);
      expect(service.displayFunction).toBeDefined();
    });

    it('should create service with custom display function', () => {
      service = new NotificationService({
        displayFunction: mockDisplayFunction
      });

      expect(service.displayFunction).toBe(mockDisplayFunction);
    });

    it('should create service with logging disabled', () => {
      service = new NotificationService({
        enableLogging: false
      });

      expect(service.enableLogging).toBe(false);
    });
  });

  describe('defaultDisplay', () => {
    it('should call window.alert with formatted message', () => {
      const alertSpy = jest.spyOn(window, 'alert').mockImplementation(() => {});
      service = new NotificationService();

      service.defaultDisplay('Test message', 'Test Title', 'info');

      expect(alertSpy).toHaveBeenCalledWith('Test Title: Test message');

      alertSpy.mockRestore();
    });

    it('should call window.alert with message only if no title', () => {
      const alertSpy = jest.spyOn(window, 'alert').mockImplementation(() => {});
      service = new NotificationService();

      service.defaultDisplay('Test message', '', 'info');

      expect(alertSpy).toHaveBeenCalledWith('Test message');

      alertSpy.mockRestore();
    });
  });

  describe('log', () => {
    beforeEach(() => {
      service = new NotificationService();
    });

    it('should log success messages to console.log', () => {
      service.log('success', 'Success Title', 'Success message');

      expect(consoleLogSpy).toHaveBeenCalledWith('SUCCESS - Success Title: Success message');
    });

    it('should log error messages to console.error', () => {
      service.log('error', 'Error Title', 'Error message');

      expect(consoleErrorSpy).toHaveBeenCalledWith('ERROR - Error Title: Error message');
    });

    it('should log warning messages to console.warn', () => {
      service.log('warning', 'Warning Title', 'Warning message');

      expect(consoleWarnSpy).toHaveBeenCalledWith('WARNING - Warning Title: Warning message');
    });

    it('should log info messages to console.info', () => {
      service.log('info', 'Info Title', 'Info message');

      expect(consoleInfoSpy).toHaveBeenCalledWith('INFO - Info Title: Info message');
    });

    it('should log unknown types to console.log', () => {
      service.log('unknown', 'Unknown Title', 'Unknown message');

      expect(consoleLogSpy).toHaveBeenCalledWith('UNKNOWN - Unknown Title: Unknown message');
    });

    it('should not log when logging is disabled', () => {
      service = new NotificationService({ enableLogging: false });

      service.log('success', 'Title', 'Message');

      expect(consoleLogSpy).not.toHaveBeenCalled();
    });
  });

  describe('success', () => {
    beforeEach(() => {
      service = new NotificationService({
        displayFunction: mockDisplayFunction
      });
    });

    it('should call display function with success type', () => {
      service.success('Operation successful', 'Success');

      expect(mockDisplayFunction).toHaveBeenCalledWith(
        'Operation successful',
        'Success',
        'success'
      );
    });

    it('should use default title if not provided', () => {
      service.success('Operation successful');

      expect(mockDisplayFunction).toHaveBeenCalledWith(
        'Operation successful',
        'Success',
        'success'
      );
    });

    it('should log the message', () => {
      service.success('Operation successful', 'Success');

      expect(consoleLogSpy).toHaveBeenCalledWith('SUCCESS - Success: Operation successful');
    });
  });

  describe('error', () => {
    beforeEach(() => {
      service = new NotificationService({
        displayFunction: mockDisplayFunction
      });
    });

    it('should call display function with error type', () => {
      service.error('Operation failed', 'Error');

      expect(mockDisplayFunction).toHaveBeenCalledWith(
        'Operation failed',
        'Error',
        'error'
      );
    });

    it('should use default title if not provided', () => {
      service.error('Operation failed');

      expect(mockDisplayFunction).toHaveBeenCalledWith(
        'Operation failed',
        'Error',
        'error'
      );
    });

    it('should log the message', () => {
      service.error('Operation failed', 'Error');

      expect(consoleErrorSpy).toHaveBeenCalledWith('ERROR - Error: Operation failed');
    });
  });

  describe('warning', () => {
    beforeEach(() => {
      service = new NotificationService({
        displayFunction: mockDisplayFunction
      });
    });

    it('should call display function with warning type', () => {
      service.warning('Warning message', 'Warning');

      expect(mockDisplayFunction).toHaveBeenCalledWith(
        'Warning message',
        'Warning',
        'warning'
      );
    });

    it('should use default title if not provided', () => {
      service.warning('Warning message');

      expect(mockDisplayFunction).toHaveBeenCalledWith(
        'Warning message',
        'Warning',
        'warning'
      );
    });

    it('should log the message', () => {
      service.warning('Warning message', 'Warning');

      expect(consoleWarnSpy).toHaveBeenCalledWith('WARNING - Warning: Warning message');
    });
  });

  describe('info', () => {
    beforeEach(() => {
      service = new NotificationService({
        displayFunction: mockDisplayFunction
      });
    });

    it('should call display function with info type', () => {
      service.info('Info message', 'Information');

      expect(mockDisplayFunction).toHaveBeenCalledWith(
        'Info message',
        'Information',
        'info'
      );
    });

    it('should use default title if not provided', () => {
      service.info('Info message');

      expect(mockDisplayFunction).toHaveBeenCalledWith(
        'Info message',
        'Info',
        'info'
      );
    });

    it('should log the message', () => {
      service.info('Info message', 'Information');

      expect(consoleInfoSpy).toHaveBeenCalledWith('INFO - Information: Info message');
    });
  });

  describe('setDisplayFunction', () => {
    beforeEach(() => {
      service = new NotificationService();
    });

    it('should update display function', () => {
      const newDisplayFn = jest.fn();

      service.setDisplayFunction(newDisplayFn);

      expect(service.displayFunction).toBe(newDisplayFn);
    });

    it('should ignore non-function values', () => {
      const originalFn = service.displayFunction;

      service.setDisplayFunction('not a function');

      expect(service.displayFunction).toBe(originalFn);
    });

    it('should use new display function for notifications', () => {
      const newDisplayFn = jest.fn();
      service.setDisplayFunction(newDisplayFn);

      service.success('Test message');

      expect(newDisplayFn).toHaveBeenCalledWith('Test message', 'Success', 'success');
    });
  });

  describe('setLogging', () => {
    beforeEach(() => {
      service = new NotificationService();
    });

    it('should enable logging', () => {
      service.setLogging(true);

      expect(service.enableLogging).toBe(true);
    });

    it('should disable logging', () => {
      service.setLogging(false);

      expect(service.enableLogging).toBe(false);
    });

    it('should prevent logging when disabled', () => {
      service.setLogging(false);

      service.success('Test message');

      expect(consoleLogSpy).not.toHaveBeenCalled();
    });
  });

  describe('integration', () => {
    it('should work with custom toast library integration', () => {
      const mockToast = {
        show: jest.fn()
      };

      service = new NotificationService({
        displayFunction: (message, title, type) => {
          mockToast.show({
            type: type,
            title: title,
            message: message
          });
        }
      });

      service.success('Plugin built successfully', 'Build Complete');

      expect(mockToast.show).toHaveBeenCalledWith({
        type: 'success',
        title: 'Build Complete',
        message: 'Plugin built successfully'
      });
    });
  });
});
