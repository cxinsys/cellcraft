import { DialogService, EventDialogService } from '@/services/dialogService';

describe('DialogService', () => {
  let service;
  let confirmSpy;
  let alertSpy;

  beforeEach(() => {
    service = new DialogService();
    confirmSpy = jest.spyOn(window, 'confirm').mockReturnValue(true);
    alertSpy = jest.spyOn(window, 'alert').mockImplementation(() => {});
  });

  afterEach(() => {
    confirmSpy.mockRestore();
    alertSpy.mockRestore();
  });

  describe('confirm', () => {
    it('should call window.confirm with message', () => {
      const message = 'Are you sure?';
      service.confirm(message);

      expect(confirmSpy).toHaveBeenCalledWith(message);
    });

    it('should return true when user confirms', () => {
      confirmSpy.mockReturnValue(true);

      const result = service.confirm('Test message');

      expect(result).toBe(true);
    });

    it('should return false when user cancels', () => {
      confirmSpy.mockReturnValue(false);

      const result = service.confirm('Test message');

      expect(result).toBe(false);
    });
  });

  describe('alert', () => {
    it('should call window.alert with message', () => {
      const message = 'Success!';
      service.alert(message);

      expect(alertSpy).toHaveBeenCalledWith(message);
    });
  });

  describe('alertWithTitle', () => {
    it('should format message with title', () => {
      service.alertWithTitle('File saved', 'Success');

      expect(alertSpy).toHaveBeenCalledWith('Success:\nFile saved');
    });

    it('should use default title if not provided', () => {
      service.alertWithTitle('Some message');

      expect(alertSpy).toHaveBeenCalledWith('Alert:\nSome message');
    });
  });
});

describe('EventDialogService', () => {
  let service;
  let confirmSpy;
  let alertSpy;

  beforeEach(() => {
    service = new EventDialogService();
    confirmSpy = jest.spyOn(window, 'confirm').mockReturnValue(true);
    alertSpy = jest.spyOn(window, 'alert').mockImplementation(() => {});
  });

  afterEach(() => {
    confirmSpy.mockRestore();
    alertSpy.mockRestore();
  });

  describe('confirm', () => {
    it('should use window.confirm by default', async () => {
      confirmSpy.mockReturnValue(true);

      const result = await service.confirm('Are you sure?');

      expect(result).toBe(true);
      expect(confirmSpy).toHaveBeenCalledWith('Are you sure?');
    });

    it('should use custom handler if registered', async () => {
      const customHandler = jest.fn((message, resolve) => resolve(true));
      service.onConfirm(customHandler);

      const result = await service.confirm('Custom confirm?');

      expect(result).toBe(true);
      expect(customHandler).toHaveBeenCalledWith('Custom confirm?', expect.any(Function));
      expect(confirmSpy).not.toHaveBeenCalled();
    });

    it('should handle custom handler returning false', async () => {
      const customHandler = jest.fn((message, resolve) => resolve(false));
      service.onConfirm(customHandler);

      const result = await service.confirm('Custom confirm?');

      expect(result).toBe(false);
    });

    it('should return promise that resolves with user choice', async () => {
      confirmSpy.mockReturnValue(false);

      const result = await service.confirm('Test?');

      expect(result).toBe(false);
    });
  });

  describe('alert', () => {
    it('should use window.alert by default', async () => {
      await service.alert('Test message');

      expect(alertSpy).toHaveBeenCalledWith('Test message');
    });

    it('should use custom handler if registered', async () => {
      const customHandler = jest.fn((message, resolve) => resolve());
      service.onAlert(customHandler);

      await service.alert('Custom alert');

      expect(customHandler).toHaveBeenCalledWith('Custom alert', expect.any(Function));
      expect(alertSpy).not.toHaveBeenCalled();
    });

    it('should return promise that resolves when dismissed', async () => {
      const promise = service.alert('Test');

      expect(promise).toBeInstanceOf(Promise);
      await expect(promise).resolves.toBeUndefined();
    });
  });

  describe('onConfirm', () => {
    it('should register custom confirm handler', () => {
      const handler = jest.fn();
      service.onConfirm(handler);

      expect(service.listeners.onConfirm).toBe(handler);
    });

    it('should replace existing handler', () => {
      const handler1 = jest.fn();
      const handler2 = jest.fn();

      service.onConfirm(handler1);
      service.onConfirm(handler2);

      expect(service.listeners.onConfirm).toBe(handler2);
    });
  });

  describe('onAlert', () => {
    it('should register custom alert handler', () => {
      const handler = jest.fn();
      service.onAlert(handler);

      expect(service.listeners.onAlert).toBe(handler);
    });

    it('should replace existing handler', () => {
      const handler1 = jest.fn();
      const handler2 = jest.fn();

      service.onAlert(handler1);
      service.onAlert(handler2);

      expect(service.listeners.onAlert).toBe(handler2);
    });
  });

  describe('clearListeners', () => {
    it('should clear all registered listeners', () => {
      service.onConfirm(jest.fn());
      service.onAlert(jest.fn());

      service.clearListeners();

      expect(service.listeners).toEqual({});
    });

    it('should fall back to window methods after clearing', async () => {
      const customHandler = jest.fn((message, resolve) => resolve(true));
      service.onConfirm(customHandler);

      service.clearListeners();

      confirmSpy.mockReturnValue(true);
      const result = await service.confirm('Test?');

      expect(result).toBe(true);
      expect(confirmSpy).toHaveBeenCalled();
      expect(customHandler).not.toHaveBeenCalled();
    });
  });
});
