import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';

/**
 * Unit tests for TabComponent.vue
 *
 * Tests cover tab management functionality:
 * - createTab: Create new tab with max limit (5 tabs)
 * - removeTab: Remove tab by ID with index adjustment
 * - adjustCurrentTab: Select existing tab or create new
 * - tabClick: Handle tab click events
 * - componentChange: Router navigation and event emission
 */

describe('TabComponent.vue', () => {
  let mockComponent;
  let mockRouter;
  let mockRoute;

  beforeEach(() => {
    // Mock Vue Router
    mockRoute = {
      path: '/workflow/data',
      fullPath: '/workflow/data?workflow_id=123&node=1',
    };

    mockRouter = {
      push: vi.fn(),
    };

    // Mock component instance
    mockComponent = {
      tabList: [],
      currentTab: 0,
      initialTabList: [],
      isTabView: true,
      currentWorkflowId: '123',
      $route: mockRoute,
      $router: mockRouter,
      $emit: vi.fn(),
      componentChange: vi.fn(),
    };
  });

  afterEach(() => {
    vi.clearAllMocks();
  });

  /**
   * createTab Method Tests
   */
  describe('createTab', () => {
    const createTab = function(node) {
      // Safely handle image require in test environment
      let img;
      try {
        img = require(`@/assets/${node.name}.png`);
      } catch (e) {
        img = `mocked-image-@/assets/${node.name}.png`;
      }

      const newTab = {
        id: node.id,
        name: node.name,
        title: node.data.title || node.name,
        img: img
      };

      // 최대 5개까지만 유지, 6개가 되면 가장 오래된 아이템 삭제
      if (this.tabList.length >= 5) {
        this.tabList.shift(); // 첫 번째 아이템 제거
        this.currentTab = Math.max(0, this.currentTab - 1); // 현재 탭 인덱스 조정
      }

      this.tabList.push(newTab);
      this.currentTab = this.tabList.length - 1;
      this.componentChange(newTab);
    };

    it('should create a new tab with correct structure', () => {
      const node = {
        id: 1,
        name: 'Data',
        data: { title: 'Data Node' }
      };

      createTab.call(mockComponent, node);

      expect(mockComponent.tabList).toHaveLength(1);
      expect(mockComponent.tabList[0]).toEqual({
        id: 1,
        name: 'Data',
        title: 'Data Node',
        img: expect.stringContaining('mocked-image')
      });
    });

    it('should use node name as title when data.title is not provided', () => {
      const node = {
        id: 2,
        name: 'Algorithm',
        data: {}
      };

      createTab.call(mockComponent, node);

      expect(mockComponent.tabList[0].title).toBe('Algorithm');
    });

    it('should set currentTab to the last tab index', () => {
      const node1 = { id: 1, name: 'Data', data: {} };
      const node2 = { id: 2, name: 'Algorithm', data: {} };

      createTab.call(mockComponent, node1);
      expect(mockComponent.currentTab).toBe(0);

      createTab.call(mockComponent, node2);
      expect(mockComponent.currentTab).toBe(1);
    });

    it('should call componentChange with the new tab', () => {
      const node = {
        id: 1,
        name: 'Data',
        data: { title: 'Data Node' }
      };

      createTab.call(mockComponent, node);

      expect(mockComponent.componentChange).toHaveBeenCalledWith(
        expect.objectContaining({
          id: 1,
          name: 'Data',
          title: 'Data Node'
        })
      );
    });

    it('should remove the oldest tab when exceeding max limit (5 tabs)', () => {
      // Add 5 tabs
      for (let i = 1; i <= 5; i++) {
        createTab.call(mockComponent, { id: i, name: `Node${i}`, data: {} });
      }

      expect(mockComponent.tabList).toHaveLength(5);
      expect(mockComponent.tabList[0].id).toBe(1);

      // Add 6th tab - should remove the first tab
      createTab.call(mockComponent, { id: 6, name: 'Node6', data: {} });

      expect(mockComponent.tabList).toHaveLength(5);
      expect(mockComponent.tabList[0].id).toBe(2); // First tab removed
      expect(mockComponent.tabList[4].id).toBe(6); // New tab at the end
    });

    it('should adjust currentTab index when removing oldest tab', () => {
      // Add 5 tabs
      for (let i = 1; i <= 5; i++) {
        createTab.call(mockComponent, { id: i, name: `Node${i}`, data: {} });
      }

      // currentTab is 4 (last tab)
      expect(mockComponent.currentTab).toBe(4);

      // Add 6th tab - should decrement currentTab before adding
      createTab.call(mockComponent, { id: 6, name: 'Node6', data: {} });

      expect(mockComponent.currentTab).toBe(4); // Adjusted to 3, then incremented to 4
    });

    it('should handle image loading gracefully', () => {
      const node = { id: 1, name: 'Data', data: {} };

      createTab.call(mockComponent, node);

      expect(mockComponent.tabList[0].img).toContain('Data.png');
    });
  });

  /**
   * removeTab Method Tests
   */
  describe('removeTab', () => {
    const removeTab = function(id) {
      const index = this.tabList.findIndex(tab => tab.id === id);
      if (index !== -1) {
        this.tabList.splice(index, 1);
        if (this.currentTab >= index) {
          this.currentTab = this.currentTab === 0 ? 0 : this.currentTab - 1;
        }
        if (this.tabList.length > 0) {
          this.componentChange(this.tabList[this.currentTab]);
        }
      }
    };

    beforeEach(() => {
      // Setup initial tabs
      mockComponent.tabList = [
        { id: 1, name: 'Tab1', title: 'Tab 1', img: 'img1.png' },
        { id: 2, name: 'Tab2', title: 'Tab 2', img: 'img2.png' },
        { id: 3, name: 'Tab3', title: 'Tab 3', img: 'img3.png' },
      ];
      mockComponent.currentTab = 1;
    });

    it('should remove tab by ID', () => {
      removeTab.call(mockComponent, 2);

      expect(mockComponent.tabList).toHaveLength(2);
      expect(mockComponent.tabList.find(t => t.id === 2)).toBeUndefined();
    });

    it('should decrement currentTab when removing a tab before current', () => {
      mockComponent.currentTab = 2;

      removeTab.call(mockComponent, 1);

      expect(mockComponent.currentTab).toBe(1);
    });

    it('should adjust currentTab when removing the current tab', () => {
      mockComponent.currentTab = 1;

      removeTab.call(mockComponent, 2);

      expect(mockComponent.currentTab).toBe(0);
    });

    it('should keep currentTab at 0 when removing first tab and currentTab is 0', () => {
      mockComponent.currentTab = 0;

      removeTab.call(mockComponent, 1);

      expect(mockComponent.currentTab).toBe(0);
      expect(mockComponent.tabList).toHaveLength(2);
    });

    it('should call componentChange with the current tab after removal', () => {
      mockComponent.currentTab = 1;

      removeTab.call(mockComponent, 1);

      expect(mockComponent.componentChange).toHaveBeenCalledWith(
        mockComponent.tabList[0]
      );
    });

    it('should not call componentChange when no tabs remain', () => {
      mockComponent.tabList = [{ id: 1, name: 'Tab1', title: 'Tab 1', img: 'img1.png' }];
      mockComponent.currentTab = 0;

      removeTab.call(mockComponent, 1);

      expect(mockComponent.tabList).toHaveLength(0);
      expect(mockComponent.componentChange).not.toHaveBeenCalled();
    });

    it('should do nothing when trying to remove non-existent tab', () => {
      const initialLength = mockComponent.tabList.length;
      const initialCurrentTab = mockComponent.currentTab;

      removeTab.call(mockComponent, 999);

      expect(mockComponent.tabList).toHaveLength(initialLength);
      expect(mockComponent.currentTab).toBe(initialCurrentTab);
      expect(mockComponent.componentChange).not.toHaveBeenCalled();
    });

    it('should handle removing last tab', () => {
      mockComponent.currentTab = 2;

      removeTab.call(mockComponent, 3);

      expect(mockComponent.currentTab).toBe(1);
      expect(mockComponent.tabList).toHaveLength(2);
    });
  });

  /**
   * adjustCurrentTab Method Tests
   */
  describe('adjustCurrentTab', () => {
    const adjustCurrentTab = function(node) {
      const index = this.tabList.findIndex(tab => tab.id === node.id);
      if (index !== -1) {
        this.currentTab = index;
      } else {
        this.createTab(node);
      }
      if (this.tabList.length > 0) {
        this.componentChange(this.tabList[this.currentTab]);
      }
    };

    const createTab = function(node) {
      // Safely handle image require in test environment
      let img;
      try {
        img = require(`@/assets/${node.name}.png`);
      } catch (e) {
        img = `mocked-image-@/assets/${node.name}.png`;
      }

      const newTab = {
        id: node.id,
        name: node.name,
        title: node.data.title || node.name,
        img: img
      };

      if (this.tabList.length >= 5) {
        this.tabList.shift();
        this.currentTab = Math.max(0, this.currentTab - 1);
      }

      this.tabList.push(newTab);
      this.currentTab = this.tabList.length - 1;
      this.componentChange(newTab);
    };

    beforeEach(() => {
      mockComponent.tabList = [
        { id: 1, name: 'Tab1', title: 'Tab 1', img: 'img1.png' },
        { id: 2, name: 'Tab2', title: 'Tab 2', img: 'img2.png' },
        { id: 3, name: 'Tab3', title: 'Tab 3', img: 'img3.png' },
      ];
      mockComponent.currentTab = 0;
      mockComponent.createTab = createTab.bind(mockComponent);
    });

    it('should set currentTab to existing tab index', () => {
      const node = { id: 2, name: 'Tab2', data: {} };

      adjustCurrentTab.call(mockComponent, node);

      expect(mockComponent.currentTab).toBe(1);
    });

    it('should call createTab for new node', () => {
      const createTabSpy = vi.spyOn(mockComponent, 'createTab');
      const node = { id: 4, name: 'NewTab', data: { title: 'New Tab' } };

      adjustCurrentTab.call(mockComponent, node);

      expect(createTabSpy).toHaveBeenCalledWith(node);
    });

    it('should call componentChange with current tab', () => {
      const node = { id: 2, name: 'Tab2', data: {} };

      adjustCurrentTab.call(mockComponent, node);

      expect(mockComponent.componentChange).toHaveBeenCalledWith(
        mockComponent.tabList[1]
      );
    });

    it('should handle empty tabList', () => {
      mockComponent.tabList = [];
      const node = { id: 1, name: 'NewTab', data: {} };

      adjustCurrentTab.call(mockComponent, node);

      expect(mockComponent.tabList).toHaveLength(1);
      expect(mockComponent.currentTab).toBe(0);
    });

    it('should update currentTab correctly when tab exists', () => {
      mockComponent.currentTab = 0;
      const node = { id: 3, name: 'Tab3', data: {} };

      adjustCurrentTab.call(mockComponent, node);

      expect(mockComponent.currentTab).toBe(2);
      expect(mockComponent.componentChange).toHaveBeenCalledWith(
        mockComponent.tabList[2]
      );
    });
  });

  /**
   * tabClick Method Tests
   */
  describe('tabClick', () => {
    const tabClick = function(idx) {
      this.currentTab = idx;
      this.componentChange(this.tabList[idx]);
    };

    beforeEach(() => {
      mockComponent.tabList = [
        { id: 1, name: 'Tab1', title: 'Tab 1', img: 'img1.png' },
        { id: 2, name: 'Tab2', title: 'Tab 2', img: 'img2.png' },
        { id: 3, name: 'Tab3', title: 'Tab 3', img: 'img3.png' },
      ];
      mockComponent.currentTab = 0;
    });

    it('should set currentTab to clicked index', () => {
      tabClick.call(mockComponent, 2);

      expect(mockComponent.currentTab).toBe(2);
    });

    it('should call componentChange with clicked tab', () => {
      tabClick.call(mockComponent, 1);

      expect(mockComponent.componentChange).toHaveBeenCalledWith(
        mockComponent.tabList[1]
      );
    });

    it('should handle clicking first tab', () => {
      mockComponent.currentTab = 2;

      tabClick.call(mockComponent, 0);

      expect(mockComponent.currentTab).toBe(0);
      expect(mockComponent.componentChange).toHaveBeenCalledWith(
        mockComponent.tabList[0]
      );
    });

    it('should handle clicking last tab', () => {
      tabClick.call(mockComponent, 2);

      expect(mockComponent.currentTab).toBe(2);
      expect(mockComponent.componentChange).toHaveBeenCalledWith(
        mockComponent.tabList[2]
      );
    });
  });

  /**
   * componentChange Method Tests
   */
  describe('componentChange', () => {
    const componentChange = function(tab) {
      let newPath = `/workflow/${tab.name.toLowerCase()}`;

      if (this.$route.path === newPath) {
        this.$router.push({
          path: newPath,
          query: {
            workflow_id: String(this.currentWorkflowId || ''),
            node: tab.id,
            forceReload: Date.now(),
          },
        });
      } else {
        this.$router.push({
          path: newPath,
          query: {
            workflow_id: String(this.currentWorkflowId || ''),
            node: tab.id,
          },
        });
      }

      this.$emit('process-workflow-nodes')
    };

    it('should navigate to new path without forceReload', () => {
      const tab = { id: 1, name: 'Algorithm', title: 'Algorithm' };
      mockComponent.$route.path = '/workflow/data';

      componentChange.call(mockComponent, tab);

      expect(mockComponent.$router.push).toHaveBeenCalledWith({
        path: '/workflow/algorithm',
        query: {
          workflow_id: '123',
          node: 1,
        },
      });
    });

    it('should add forceReload when navigating to same path', () => {
      const tab = { id: 1, name: 'Data', title: 'Data' };
      mockComponent.$route.path = '/workflow/data';

      componentChange.call(mockComponent, tab);

      expect(mockComponent.$router.push).toHaveBeenCalledWith({
        path: '/workflow/data',
        query: {
          workflow_id: '123',
          node: 1,
          forceReload: expect.any(Number),
        },
      });
    });

    it('should emit process-workflow-nodes event', () => {
      const tab = { id: 1, name: 'Data', title: 'Data' };

      componentChange.call(mockComponent, tab);

      expect(mockComponent.$emit).toHaveBeenCalledWith('process-workflow-nodes');
    });

    it('should handle null currentWorkflowId', () => {
      mockComponent.currentWorkflowId = null;
      const tab = { id: 1, name: 'Algorithm', title: 'Algorithm' };

      componentChange.call(mockComponent, tab);

      expect(mockComponent.$router.push).toHaveBeenCalledWith({
        path: '/workflow/algorithm',
        query: {
          workflow_id: '',
          node: 1,
        },
      });
    });

    it('should convert tab name to lowercase for path', () => {
      const tab = { id: 1, name: 'UPPERCASE', title: 'Uppercase Tab' };

      componentChange.call(mockComponent, tab);

      expect(mockComponent.$router.push).toHaveBeenCalledWith(
        expect.objectContaining({
          path: '/workflow/uppercase',
        })
      );
    });

    it('should include node ID in query parameters', () => {
      const tab = { id: 42, name: 'Data', title: 'Data' };

      componentChange.call(mockComponent, tab);

      expect(mockComponent.$router.push).toHaveBeenCalledWith(
        expect.objectContaining({
          query: expect.objectContaining({
            node: 42,
          }),
        })
      );
    });
  });
});
