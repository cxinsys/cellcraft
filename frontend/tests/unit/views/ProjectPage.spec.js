import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { shallowMount } from '@vue/test-utils';
import ProjectPage from '@/views/ProjectPage.vue';
import * as apiIndex from '@/api/index';

// Mock API functions
vi.mock('@/api/index', () => ({
  getUser: vi.fn(),
  getWorkflows: vi.fn(),
  deleteWorkflow: vi.fn(),
  getPlugins: vi.fn()
}));

describe('ProjectPage.vue', () => {
  let wrapper;
  let mockWorkflows;

  beforeEach(() => {
    // Setup mock workflows data
    mockWorkflows = [
      { id: 1, title: 'Workflow 1', updated_at: '2025-10-24T10:00:00', thumbnail: null },
      { id: 2, title: 'Workflow 2', updated_at: '2025-10-24T09:00:00', thumbnail: 'thumb2.png' },
      { id: 3, title: 'Workflow 3', updated_at: '2025-10-24T08:00:00', thumbnail: null }
    ];

    // Create wrapper with initial data
    wrapper = shallowMount(ProjectPage, {
      data() {
        return {
          workflows: [...mockWorkflows],
          workflow_id: 2,
          list_idx: 1,
          targetWorkflow: null,
          showDeleteModal: false,
          toggleMessage: false,
          messageStatus: '',
          messageContent: ''
        };
      },
      mocks: {
        $router: {
          push: vi.fn()
        }
      }
    });

    // Clear all mocks before each test
    vi.clearAllMocks();
  });

  afterEach(() => {
    vi.restoreAllMocks();
  });

  /**
   * Data Initialization
   */
  describe('Data Initialization', () => {
    it('should initialize with correct default data', () => {
      const freshWrapper = shallowMount(ProjectPage, {
        mocks: {
          $router: { push: vi.fn() }
        }
      });

      expect(freshWrapper.vm.toggleMenu).toBe(true);
      expect(freshWrapper.vm.workflows).toBeNull();
      expect(freshWrapper.vm.R_Mouse_isActive).toBe(false);
      expect(freshWrapper.vm.showDeleteModal).toBe(false);
      expect(freshWrapper.vm.isSelectModalVisible).toBe(false);
    });
  });

  /**
   * Modal Management
   */
  describe('Modal Management', () => {
    describe('createProject', () => {
      it('should set isSelectModalVisible to true', () => {
        wrapper.vm.isSelectModalVisible = false;

        wrapper.vm.createProject();

        expect(wrapper.vm.isSelectModalVisible).toBe(true);
      });
    });

    describe('closeSelectModal', () => {
      it('should set isSelectModalVisible to false', () => {
        wrapper.vm.isSelectModalVisible = true;

        wrapper.vm.closeSelectModal();

        expect(wrapper.vm.isSelectModalVisible).toBe(false);
      });
    });
  });

  /**
   * Context Menu
   */
  describe('Context Menu', () => {
    describe('RMouseClick', () => {
      it('should activate context menu with correct position', () => {
        const mockEvent = {
          clientX: 150,
          clientY: 200
        };

        wrapper.vm.RMouseClick(mockEvent, 5, 2);

        expect(wrapper.vm.R_Mouse_isActive).toBe(true);
        expect(wrapper.vm.xPosition).toBe('150px');
        expect(wrapper.vm.yPosition).toBe('200px');
        expect(wrapper.vm.workflow_id).toBe(5);
        expect(wrapper.vm.list_idx).toBe(2);
      });

      it('should reset menu state before showing', () => {
        wrapper.vm.R_Mouse_isActive = true;

        const mockEvent = { clientX: 100, clientY: 100 };
        wrapper.vm.RMouseClick(mockEvent, 1, 0);

        expect(wrapper.vm.R_Mouse_isActive).toBe(true);
      });
    });

    describe('ClickOut', () => {
      it('should deactivate context menu', () => {
        wrapper.vm.R_Mouse_isActive = true;

        wrapper.vm.ClickOut();

        expect(wrapper.vm.R_Mouse_isActive).toBe(false);
      });
    });
  });

  /**
   * Workflow Management
   */
  describe('Workflow Management', () => {
    describe('confirmDelete', () => {
      it('should set targetWorkflow and show delete modal', () => {
        wrapper.vm.workflows = mockWorkflows;
        wrapper.vm.list_idx = 1;
        wrapper.vm.showDeleteModal = false;

        wrapper.vm.confirmDelete();

        expect(wrapper.vm.targetWorkflow).toEqual(mockWorkflows[1]);
        expect(wrapper.vm.showDeleteModal).toBe(true);
      });
    });

    describe('removeWorkflow', () => {
      it('should delete workflow successfully', async () => {
        apiIndex.deleteWorkflow.mockResolvedValue({ success: true });

        wrapper.vm.workflows = [...mockWorkflows];
        wrapper.vm.workflow_id = 2;
        wrapper.vm.list_idx = 1;
        wrapper.vm.showDeleteModal = true;

        const initialLength = wrapper.vm.workflows.length;

        await wrapper.vm.removeWorkflow();

        expect(apiIndex.deleteWorkflow).toHaveBeenCalledWith({ id: 2 });
        expect(wrapper.vm.workflows).toHaveLength(initialLength - 1);
        expect(wrapper.vm.workflows.find(w => w.id === 2)).toBeUndefined();
        expect(wrapper.vm.showDeleteModal).toBe(false);
        expect(wrapper.vm.messageStatus).toBe('success');
        expect(wrapper.vm.messageContent).toContain('successfully');
      });

      it('should restore workflow on API failure', async () => {
        apiIndex.deleteWorkflow.mockRejectedValue(new Error('API error'));

        wrapper.vm.workflows = [...mockWorkflows];
        wrapper.vm.workflow_id = 2;
        wrapper.vm.list_idx = 1;

        const initialLength = wrapper.vm.workflows.length;
        const targetWorkflow = wrapper.vm.workflows[1];

        await wrapper.vm.removeWorkflow();

        // Workflow should be restored
        expect(wrapper.vm.workflows).toHaveLength(initialLength);
        expect(wrapper.vm.workflows[1]).toEqual(targetWorkflow);
        expect(wrapper.vm.messageStatus).toBe('error');
        expect(wrapper.vm.messageContent).toContain('Failed');
      });

      it('should restore workflow at correct index after error', async () => {
        apiIndex.deleteWorkflow.mockRejectedValue(new Error('Network error'));

        wrapper.vm.workflows = [...mockWorkflows];
        wrapper.vm.workflow_id = 3;
        wrapper.vm.list_idx = 2; // Last item

        await wrapper.vm.removeWorkflow();

        // Check that workflow 3 is back at index 2
        expect(wrapper.vm.workflows[2].id).toBe(3);
        expect(wrapper.vm.workflows[2].title).toBe('Workflow 3');
      });
    });

    describe('openWorkflow', () => {
      it('should navigate with workflow_id parameter when provided', () => {
        wrapper.vm.openWorkflow(5);

        expect(wrapper.vm.$router.push).toHaveBeenCalledWith({
          path: '/workflow',
          query: { workflow_id: '5' }
        });
      });

      it('should use stored workflow_id when not provided', () => {
        wrapper.vm.workflow_id = 10;

        wrapper.vm.openWorkflow();

        expect(wrapper.vm.$router.push).toHaveBeenCalledWith({
          path: '/workflow',
          query: { workflow_id: '10' }
        });
      });

      it('should convert workflow_id to string', () => {
        wrapper.vm.openWorkflow(42);

        expect(wrapper.vm.$router.push).toHaveBeenCalledWith({
          path: '/workflow',
          query: { workflow_id: '42' }
        });
      });
    });
  });

  /**
   * Plugin Selection
   */
  describe('Plugin Selection', () => {
    describe('selectPlugin', () => {
      it('should navigate with plugin_id and close modal', () => {
        const mockPlugin = { id: 'plugin-123', name: 'Test Plugin' };
        wrapper.vm.isSelectModalVisible = true;

        wrapper.vm.selectPlugin(mockPlugin);

        expect(wrapper.vm.$router.push).toHaveBeenCalledWith({
          path: '/workflow',
          query: { plugin_id: 'plugin-123' }
        });
        expect(wrapper.vm.isSelectModalVisible).toBe(false);
      });
    });

    describe('selectDefault', () => {
      it('should navigate to workflow without query and close modal', () => {
        wrapper.vm.isSelectModalVisible = true;

        wrapper.vm.selectDefault();

        expect(wrapper.vm.$router.push).toHaveBeenCalledWith({
          path: '/workflow'
        });
        expect(wrapper.vm.isSelectModalVisible).toBe(false);
      });
    });

    describe('selectPlugins', () => {
      it('should navigate with comma-separated plugin_ids', () => {
        const mockPlugins = [
          { id: 'plugin-1', name: 'Plugin 1' },
          { id: 'plugin-2', name: 'Plugin 2' },
          { id: 'plugin-3', name: 'Plugin 3' }
        ];
        wrapper.vm.isSelectModalVisible = true;

        wrapper.vm.selectPlugins(mockPlugins);

        expect(wrapper.vm.$router.push).toHaveBeenCalledWith({
          path: '/workflow',
          query: { plugin_ids: 'plugin-1,plugin-2,plugin-3' }
        });
        expect(wrapper.vm.isSelectModalVisible).toBe(false);
      });

      it('should handle single plugin in array', () => {
        const mockPlugins = [{ id: 'single-plugin', name: 'Single' }];

        wrapper.vm.selectPlugins(mockPlugins);

        expect(wrapper.vm.$router.push).toHaveBeenCalledWith({
          path: '/workflow',
          query: { plugin_ids: 'single-plugin' }
        });
      });
    });
  });

  /**
   * UI Feedback
   */
  describe('UI Feedback', () => {
    describe('setMessage', () => {
      beforeEach(() => {
        vi.useFakeTimers();
      });

      afterEach(() => {
        vi.restoreAllMocks();
      });

      it('should display success message', () => {
        wrapper.vm.setMessage('success', 'Operation completed');

        expect(wrapper.vm.toggleMessage).toBe(true);
        expect(wrapper.vm.messageStatus).toBe('success');
        expect(wrapper.vm.messageContent).toBe('Operation completed');
      });

      it('should display error message', () => {
        wrapper.vm.setMessage('error', 'Operation failed');

        expect(wrapper.vm.toggleMessage).toBe(true);
        expect(wrapper.vm.messageStatus).toBe('error');
        expect(wrapper.vm.messageContent).toBe('Operation failed');
      });

      it('should hide message after 5 seconds', () => {
        wrapper.vm.setMessage('success', 'Test message');

        expect(wrapper.vm.toggleMessage).toBe(true);

        vi.advanceTimersByTime(5000);

        expect(wrapper.vm.toggleMessage).toBe(false);
      });

      it('should not interfere with multiple messages', () => {
        wrapper.vm.setMessage('success', 'First message');
        vi.advanceTimersByTime(2000);

        wrapper.vm.setMessage('error', 'Second message');

        expect(wrapper.vm.messageStatus).toBe('error');
        expect(wrapper.vm.messageContent).toBe('Second message');

        vi.advanceTimersByTime(5000);
        expect(wrapper.vm.toggleMessage).toBe(false);
      });
    });
  });

  /**
   * Data Fetching
   */
  describe('Data Fetching', () => {
    describe('fetchUserProfile', () => {
      it('should fetch user profile and return username', async () => {
        const mockProfile = {
          data: {
            username: 'testuser',
            email: 'test@example.com'
          }
        };
        apiIndex.getUser.mockResolvedValue(mockProfile);

        const username = await wrapper.vm.fetchUserProfile();

        expect(apiIndex.getUser).toHaveBeenCalled();
        expect(wrapper.vm.profile).toEqual(mockProfile.data);
        expect(username).toBe('testuser');
      });

      it('should handle API response correctly', async () => {
        const mockProfile = {
          data: {
            username: 'admin',
            email: 'admin@cellcraft.com'
          }
        };
        apiIndex.getUser.mockResolvedValue(mockProfile);

        await wrapper.vm.fetchUserProfile();

        expect(wrapper.vm.profile.username).toBe('admin');
        expect(wrapper.vm.profile.email).toBe('admin@cellcraft.com');
      });
    });

    describe('fetchWorkflows', () => {
      it('should fetch and sort workflows by updated_at', async () => {
        const mockWorkflows = {
          data: [
            { id: 1, title: 'Old', updated_at: '2025-10-23T08:00:00' },
            { id: 2, title: 'New', updated_at: '2025-10-24T10:00:00' },
            { id: 3, title: 'Mid', updated_at: '2025-10-23T12:00:00' }
          ]
        };
        apiIndex.getWorkflows.mockResolvedValue(mockWorkflows);

        await wrapper.vm.fetchWorkflows();

        expect(apiIndex.getWorkflows).toHaveBeenCalled();
        expect(wrapper.vm.workflows).toHaveLength(3);
        expect(wrapper.vm.workflows[0].title).toBe('New');
        expect(wrapper.vm.workflows[1].title).toBe('Mid');
        expect(wrapper.vm.workflows[2].title).toBe('Old');
      });

      it('should handle empty workflows array', async () => {
        apiIndex.getWorkflows.mockResolvedValue({ data: [] });

        await wrapper.vm.fetchWorkflows();

        expect(wrapper.vm.workflows).toEqual([]);
      });

      it('should sort workflows in descending order', async () => {
        const mockWorkflows = {
          data: [
            { id: 1, updated_at: '2025-10-20T10:00:00' },
            { id: 2, updated_at: '2025-10-24T10:00:00' }
          ]
        };
        apiIndex.getWorkflows.mockResolvedValue(mockWorkflows);

        await wrapper.vm.fetchWorkflows();

        expect(wrapper.vm.workflows[0].id).toBe(2);
        expect(wrapper.vm.workflows[1].id).toBe(1);
      });
    });

    describe('fetchPlugins', () => {
      it('should fetch plugins and mark checked for current user', async () => {
        const mockPlugins = {
          data: {
            plugins: [
              { id: 1, name: 'Plugin 1', users: [{ username: 'testuser' }] },
              { id: 2, name: 'Plugin 2', users: [{ username: 'otheruser' }] }
            ]
          }
        };
        apiIndex.getPlugins.mockResolvedValue(mockPlugins);

        await wrapper.vm.fetchPlugins('testuser');

        expect(apiIndex.getPlugins).toHaveBeenCalled();
        expect(wrapper.vm.plugins).toHaveLength(2);
        expect(wrapper.vm.plugins[0].checked).toBe(true);
        expect(wrapper.vm.plugins[1].checked).toBe(false);
      });

      it('should handle plugins without matching users', async () => {
        const mockPlugins = {
          data: {
            plugins: [
              { id: 1, name: 'Plugin 1', users: [{ username: 'other' }] }
            ]
          }
        };
        apiIndex.getPlugins.mockResolvedValue(mockPlugins);

        await wrapper.vm.fetchPlugins('testuser');

        expect(wrapper.vm.plugins[0].checked).toBe(false);
      });

      it('should handle empty plugins array', async () => {
        apiIndex.getPlugins.mockResolvedValue({ data: { plugins: [] } });

        await wrapper.vm.fetchPlugins('testuser');

        expect(wrapper.vm.plugins).toEqual([]);
      });
    });

    describe('initializeData', () => {
      it('should initialize all data successfully', async () => {
        apiIndex.getUser.mockResolvedValue({
          data: { username: 'testuser', email: 'test@example.com' }
        });
        apiIndex.getWorkflows.mockResolvedValue({
          data: [{ id: 1, title: 'Test', updated_at: '2025-10-24T10:00:00' }]
        });
        apiIndex.getPlugins.mockResolvedValue({
          data: { plugins: [{ id: 1, name: 'Plugin 1', users: [{ username: 'testuser' }] }] }
        });

        await wrapper.vm.initializeData();

        expect(wrapper.vm.profile.username).toBe('testuser');
        expect(wrapper.vm.workflows).toHaveLength(1);
        expect(wrapper.vm.plugins).toHaveLength(1);
        expect(wrapper.vm.plugins[0].checked).toBe(true);
      });

      it('should handle user profile fetch failure', async () => {
        apiIndex.getUser.mockRejectedValue(new Error('User API error'));
        apiIndex.getWorkflows.mockResolvedValue({ data: [] });
        apiIndex.getPlugins.mockResolvedValue({ data: { plugins: [] } });

        const consoleSpy = vi.spyOn(console, 'error').mockImplementation();

        await wrapper.vm.initializeData();

        expect(consoleSpy).toHaveBeenCalledWith(
          'Failed to fetch user profile:',
          expect.any(Error)
        );

        consoleSpy.mockRestore();
      });

      it('should handle workflows fetch failure', async () => {
        apiIndex.getUser.mockResolvedValue({
          data: { username: 'testuser', email: 'test@example.com' }
        });
        apiIndex.getWorkflows.mockRejectedValue(new Error('Workflows API error'));
        apiIndex.getPlugins.mockResolvedValue({ data: { plugins: [] } });

        const consoleSpy = vi.spyOn(console, 'error').mockImplementation();

        await wrapper.vm.initializeData();

        expect(consoleSpy).toHaveBeenCalledWith(
          'Failed to fetch workflows:',
          expect.any(Error)
        );

        consoleSpy.mockRestore();
      });

      it('should handle plugins fetch failure', async () => {
        apiIndex.getUser.mockResolvedValue({
          data: { username: 'testuser', email: 'test@example.com' }
        });
        apiIndex.getWorkflows.mockResolvedValue({ data: [] });
        apiIndex.getPlugins.mockRejectedValue(new Error('Plugins API error'));

        const consoleSpy = vi.spyOn(console, 'error').mockImplementation();

        await wrapper.vm.initializeData();

        expect(consoleSpy).toHaveBeenCalledWith(
          'Failed to fetch plugins:',
          expect.any(Error)
        );

        consoleSpy.mockRestore();
      });
    });
  });

  /**
   * Computed Properties
   */
  describe('Computed Properties', () => {
    describe('filteredPlugins', () => {
      it('should return only checked analysis plugins', () => {
        wrapper.vm.plugins = [
          { id: 1, name: 'Plugin 1', checked: true, plugin_type: 'analysis' },
          { id: 2, name: 'Plugin 2', checked: false, plugin_type: 'analysis' },
          { id: 3, name: 'Plugin 3', checked: true, plugin_type: 'visualization' },
          { id: 4, name: 'Plugin 4', checked: true, plugin_type: 'analysis' }
        ];

        const filtered = wrapper.vm.filteredPlugins;

        expect(filtered).toHaveLength(2);
        expect(filtered[0].id).toBe(1);
        expect(filtered[1].id).toBe(4);
      });

      it('should return empty array when no plugins match', () => {
        wrapper.vm.plugins = [
          { id: 1, name: 'Plugin 1', checked: false, plugin_type: 'analysis' },
          { id: 2, name: 'Plugin 2', checked: true, plugin_type: 'visualization' }
        ];

        const filtered = wrapper.vm.filteredPlugins;

        expect(filtered).toHaveLength(0);
      });

      it('should return empty array when plugins is empty', () => {
        wrapper.vm.plugins = [];

        const filtered = wrapper.vm.filteredPlugins;

        expect(filtered).toHaveLength(0);
      });
    });

    describe('filteredMessageContent', () => {
      it('should split message by periods and filter empty strings', () => {
        wrapper.vm.messageContent = 'First sentence.Second sentence.Third sentence.';

        const filtered = wrapper.vm.filteredMessageContent;

        expect(filtered).toEqual(['First sentence', 'Second sentence', 'Third sentence']);
      });

      it('should handle message without periods', () => {
        wrapper.vm.messageContent = 'Single message';

        const filtered = wrapper.vm.filteredMessageContent;

        expect(filtered).toEqual(['Single message']);
      });

      it('should handle empty message', () => {
        wrapper.vm.messageContent = '';

        const filtered = wrapper.vm.filteredMessageContent;

        expect(filtered).toEqual([]);
      });

      it('should handle message with multiple consecutive periods', () => {
        wrapper.vm.messageContent = 'First..Second...Third.';

        const filtered = wrapper.vm.filteredMessageContent;

        expect(filtered).toEqual(['First', 'Second', 'Third']);
      });
    });
  });
});
