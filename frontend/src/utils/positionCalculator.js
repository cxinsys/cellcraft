/**
 * Position calculation utilities for UI elements
 * Provides functions for calculating positions of context menus and tooltips
 */

/**
 * Default offset for context menu Y position
 * @constant {number}
 */
export const CONTEXT_MENU_OFFSET_Y = 55;

/**
 * Calculate position for context menu
 * @param {number} clientX - Mouse click X coordinate
 * @param {number} clientY - Mouse click Y coordinate
 * @param {number} offsetY - Y-axis offset (default: CONTEXT_MENU_OFFSET_Y)
 * @returns {{x: string, y: string}} Position object with pixel values
 */
export function calculateContextMenuPosition(clientX, clientY, offsetY = CONTEXT_MENU_OFFSET_Y) {
  if (typeof clientX !== 'number' || typeof clientY !== 'number' ||
      Number.isNaN(clientX) || Number.isNaN(clientY)) {
    return { x: '0px', y: '0px' };
  }

  return {
    x: `${clientX}px`,
    y: `${clientY - offsetY}px`
  };
}
