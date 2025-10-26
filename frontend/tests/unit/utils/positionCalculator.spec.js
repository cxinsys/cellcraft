import { describe, it, expect } from 'vitest';
import { calculateContextMenuPosition, CONTEXT_MENU_OFFSET_Y } from '@/utils/positionCalculator';

describe('positionCalculator.js', () => {
  /**
   * calculateContextMenuPosition Tests
   */
  describe('calculateContextMenuPosition', () => {
    it('should calculate position with default offset', () => {
      const result = calculateContextMenuPosition(100, 200);

      expect(result.x).toBe('100px');
      expect(result.y).toBe('145px'); // 200 - 55 = 145
    });

    it('should calculate position with custom offset', () => {
      const result = calculateContextMenuPosition(100, 200, 30);

      expect(result.x).toBe('100px');
      expect(result.y).toBe('170px'); // 200 - 30 = 170
    });

    it('should handle zero coordinates', () => {
      const result = calculateContextMenuPosition(0, 0);

      expect(result.x).toBe('0px');
      expect(result.y).toBe('-55px'); // 0 - 55 = -55
    });

    it('should handle negative Y offset resulting in negative position', () => {
      const result = calculateContextMenuPosition(100, 20);

      expect(result.x).toBe('100px');
      expect(result.y).toBe('-35px'); // 20 - 55 = -35
    });

    it('should handle large coordinates', () => {
      const result = calculateContextMenuPosition(1920, 1080);

      expect(result.x).toBe('1920px');
      expect(result.y).toBe('1025px'); // 1080 - 55 = 1025
    });

    it('should handle zero offset', () => {
      const result = calculateContextMenuPosition(100, 200, 0);

      expect(result.x).toBe('100px');
      expect(result.y).toBe('200px'); // 200 - 0 = 200
    });

    it('should handle negative offset (moves menu down)', () => {
      const result = calculateContextMenuPosition(100, 200, -10);

      expect(result.x).toBe('100px');
      expect(result.y).toBe('210px'); // 200 - (-10) = 210
    });

    it('should handle floating point coordinates', () => {
      const result = calculateContextMenuPosition(100.5, 200.7);

      expect(result.x).toBe('100.5px');
      expect(result.y).toBe('145.7px'); // 200.7 - 55 = 145.7
    });

    it('should handle edge case with very small numbers', () => {
      const result = calculateContextMenuPosition(0.1, 0.1);

      expect(result.x).toBe('0.1px');
      expect(result.y).toBe('-54.9px'); // 0.1 - 55 = -54.9
    });
  });

  describe('calculateContextMenuPosition - edge cases', () => {
    it('should return default position for null coordinates', () => {
      const result = calculateContextMenuPosition(null, null);

      expect(result.x).toBe('0px');
      expect(result.y).toBe('0px');
    });

    it('should return default position for undefined coordinates', () => {
      const result = calculateContextMenuPosition(undefined, undefined);

      expect(result.x).toBe('0px');
      expect(result.y).toBe('0px');
    });

    it('should return default position for non-number clientX', () => {
      const result = calculateContextMenuPosition('100', 200);

      expect(result.x).toBe('0px');
      expect(result.y).toBe('0px');
    });

    it('should return default position for non-number clientY', () => {
      const result = calculateContextMenuPosition(100, '200');

      expect(result.x).toBe('0px');
      expect(result.y).toBe('0px');
    });

    it('should return default position for NaN coordinates', () => {
      const result = calculateContextMenuPosition(NaN, NaN);

      expect(result.x).toBe('0px');
      expect(result.y).toBe('0px');
    });

    it('should return default position for object coordinates', () => {
      const result = calculateContextMenuPosition({}, {});

      expect(result.x).toBe('0px');
      expect(result.y).toBe('0px');
    });

    it('should return default position for array coordinates', () => {
      const result = calculateContextMenuPosition([100], [200]);

      expect(result.x).toBe('0px');
      expect(result.y).toBe('0px');
    });
  });

  describe('calculateContextMenuPosition - real-world scenarios', () => {
    it('should position menu near top of screen', () => {
      // Simulating click near top of screen
      const result = calculateContextMenuPosition(500, 60);

      expect(result.x).toBe('500px');
      expect(result.y).toBe('5px'); // Menu positioned slightly below click
    });

    it('should position menu in center of typical screen', () => {
      // Simulating click in center of 1920x1080 screen
      const result = calculateContextMenuPosition(960, 540);

      expect(result.x).toBe('960px');
      expect(result.y).toBe('485px');
    });

    it('should position menu near bottom of screen', () => {
      // Simulating click near bottom of screen
      const result = calculateContextMenuPosition(800, 1050);

      expect(result.x).toBe('800px');
      expect(result.y).toBe('995px');
    });

    it('should position menu at exact edge of screen', () => {
      // Simulating click at exact screen edge
      const result = calculateContextMenuPosition(0, 768);

      expect(result.x).toBe('0px');
      expect(result.y).toBe('713px');
    });
  });

  describe('CONTEXT_MENU_OFFSET_Y constant', () => {
    it('should export the default offset value', () => {
      expect(CONTEXT_MENU_OFFSET_Y).toBe(55);
    });

    it('should be a number', () => {
      expect(typeof CONTEXT_MENU_OFFSET_Y).toBe('number');
    });
  });
});
