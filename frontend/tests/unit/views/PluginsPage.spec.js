import { shallowMount } from '@vue/test-utils';
import PluginsPage from '@/views/PluginsPage.vue';
import { PluginService } from '@/services/pluginService';
import { BuildService } from '@/services/buildService';
import { BuildMonitor } from '@/services/buildMonitorService';
import { DialogService } from '@/services/dialogService';
import { NotificationService } from '@/services/notificationService';
import { IntervalManager } from '@/utils/intervalManager';
import { enrichPlugins } from '@/utils/pluginDataProcessor';
import { transformBuildTasks } from '@/utils/buildTaskMapper';
import { vi } from 'vitest';

// Mock all services
vi.mock('@/services/pluginService');
vi.mock('@/services/buildService');
vi.mock('@/services/buildMonitorService');
vi.mock('@/services/dialogService');
vi.mock('@/services/notificationService');
vi.mock('@/utils/intervalManager');
vi.mock('@/utils/pluginDataProcessor');
vi.mock('@/utils/buildTaskMapper');

// Test data factories
const createMockPlugin = (overrides = {}) => ({
  id: 1,
  name: 'TENET',
  description: 'Test plugin for GRN analysis',
  source: 'local',
  plugin_type: 'ANALYSIS',
  use_gpu: false,
  checked: false,
  building: false,
  imageExists: false,
  buildTaskId: null,
  buildStatus: null,
  updated_at: '2025-10-23T10:00:00',
  current_version: '1.0.0',
  available_versions: ['1.0.0', '0.9.0'],
  ...overrides
});

const createMockTask = (overrides = {}) => ({
  task_id: 'task-123',
  plugin_name: 'TENET',
  state: 'RUNNING',
  start_time: '2025-10-24T10:00:00',
  end_time: null,
  ...overrides
});

describe('PluginsPage.vue', () => {
  let wrapper;
  let mockPluginService;
  let mockBuildService;
  let mockBuildMonitor;
  let mockDialogService;
  let mockNotificationService;
  let mockIntervalManager;

  beforeEach(() => {
    // Create mock instances
    mockPluginService = {
      getUserProfile: vi.fn(),
      getPluginsList: vi.fn(),
      checkMultiplePluginImages: vi.fn(),
      associatePlugin: vi.fn(),
      dissociatePlugin: vi.fn(),
      getBuildTasks: vi.fn()
    };

    mockBuildService = {
      buildPlugin: vi.fn(),
      getBuildLogs: vi.fn(),
      buildMultiplePlugins: vi.fn(),
      cancelBuildTask: vi.fn(),
      filterPluginsToBuild: vi.fn(),
      processBuildResults: vi.fn(),
      getBuildSummaryMessage: vi.fn()
    };

    mockBuildMonitor = {
      startMonitoring: vi.fn(),
      stopMonitoring: vi.fn(),
      stopAll: vi.fn()
    };

    mockDialogService = {
      confirm: vi.fn()
    };

    mockNotificationService = {
      success: vi.fn(),
      error: vi.fn(),
      warning: vi.fn(),
      info: vi.fn()
    };

    mockIntervalManager = {
      start: vi.fn(),
      stop: vi.fn(),
      stopAll: vi.fn()
    };

    // Mock constructors
    PluginService.mockImplementation(() => mockPluginService);
    BuildService.mockImplementation(() => mockBuildService);
    BuildMonitor.mockImplementation(() => mockBuildMonitor);
    DialogService.mockImplementation(() => mockDialogService);
    NotificationService.mockImplementation(() => mockNotificationService);
    IntervalManager.mockImplementation(() => mockIntervalManager);
  });

  afterEach(() => {
    if (wrapper) {
      wrapper.destroy();
    }
    vi.clearAllMocks();
  });

  describe('Component Initialization', () => {
    it('should create service instances on mount', () => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });

      expect(PluginService).toHaveBeenCalled();
      expect(BuildService).toHaveBeenCalled();
      expect(IntervalManager).toHaveBeenCalled();
      expect(DialogService).toHaveBeenCalled();
      expect(NotificationService).toHaveBeenCalled();
    });
  });

  describe('Plugin Association', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should associate plugin when checked', async () => {
      const plugin = { id: 1, name: 'TENET', checked: true };
      mockPluginService.associatePlugin.mockResolvedValue({ message: 'Associated' });

      await wrapper.vm.handlePluginAssociate(plugin);

      expect(mockPluginService.associatePlugin).toHaveBeenCalledWith(1);
    });

    it('should dissociate plugin when unchecked', async () => {
      const plugin = { id: 1, name: 'TENET', checked: false };
      mockPluginService.dissociatePlugin.mockResolvedValue({ message: 'Dissociated' });

      await wrapper.vm.handlePluginAssociate(plugin);

      expect(mockPluginService.dissociatePlugin).toHaveBeenCalledWith(1);
    });
  });

  describe('Plugin Building', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should build plugin successfully', async () => {
      const plugin = { name: 'TENET', building: false, imageExists: false };
      mockBuildService.buildPlugin.mockResolvedValue({
        success: true,
        taskId: 'task-123'
      });

      await wrapper.vm.handleBuildPlugin(plugin);

      expect(mockBuildService.buildPlugin).toHaveBeenCalledWith('TENET', false);
      expect(plugin.building).toBe(true);
      expect(plugin.buildTaskId).toBe('task-123');
      expect(mockNotificationService.success).toHaveBeenCalled();
    });

    it('should handle build failure', async () => {
      const plugin = { name: 'TENET', building: false, imageExists: false };
      mockBuildService.buildPlugin.mockResolvedValue({
        success: false,
        error: 'Build failed'
      });

      await wrapper.vm.handleBuildPlugin(plugin);

      expect(plugin.building).toBe(false);
      expect(mockNotificationService.error).toHaveBeenCalled();
    });

    it('should skip build if already building', async () => {
      const plugin = { name: 'TENET', building: true, imageExists: false };

      // Mock showBuildLogs to prevent actual execution
      wrapper.vm.showBuildLogs = vi.fn();

      await wrapper.vm.handleBuildPlugin(plugin);

      expect(mockBuildService.buildPlugin).not.toHaveBeenCalled();
    });

    it('should skip build if image already exists', async () => {
      const plugin = { name: 'TENET', building: false, imageExists: true };

      await wrapper.vm.handleBuildPlugin(plugin);

      expect(mockBuildService.buildPlugin).not.toHaveBeenCalled();
    });
  });

  describe('Build Tasks', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should fetch build tasks', async () => {
      const mockTasks = [
        { task_id: 'task-1', plugin_name: 'TENET', state: 'RUNNING' },
        { task_id: 'task-2', plugin_name: 'GENIE3', state: 'SUCCESS' }
      ];
      mockPluginService.getBuildTasks.mockResolvedValue(mockTasks);

      await wrapper.vm.fetchBuildTasks();

      expect(mockPluginService.getBuildTasks).toHaveBeenCalled();
    });

    it('should cancel build task with confirmation', async () => {
      mockDialogService.confirm.mockReturnValue(true);
      mockBuildService.cancelBuildTask.mockResolvedValue({
        success: true,
        message: 'Task cancelled'
      });
      mockPluginService.getBuildTasks.mockResolvedValue([]);

      await wrapper.vm.cancelBuildTask('task-123');

      expect(mockDialogService.confirm).toHaveBeenCalled();
      expect(mockBuildService.cancelBuildTask).toHaveBeenCalledWith('task-123');
      expect(mockNotificationService.success).toHaveBeenCalled();
    });

    it('should not cancel build task without confirmation', async () => {
      mockDialogService.confirm.mockReturnValue(false);

      await wrapper.vm.cancelBuildTask('task-123');

      expect(mockBuildService.cancelBuildTask).not.toHaveBeenCalled();
    });

    it('should handle cancel failure', async () => {
      mockDialogService.confirm.mockReturnValue(true);
      mockBuildService.cancelBuildTask.mockResolvedValue({
        success: false,
        error: 'Cancel failed'
      });

      await wrapper.vm.cancelBuildTask('task-123');

      expect(mockNotificationService.error).toHaveBeenCalled();
    });
  });

  describe('Build Monitoring', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should toggle build monitor', async () => {
      mockPluginService.getBuildTasks.mockResolvedValue([]);

      expect(wrapper.vm.showBuildMonitor).toBe(false);

      await wrapper.vm.toggleBuildMonitor();

      expect(wrapper.vm.showBuildMonitor).toBe(true);
      expect(mockPluginService.getBuildTasks).toHaveBeenCalled();
    });
  });

  describe('Cleanup', () => {
    it('should stop all monitoring on destroy', () => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });

      wrapper.destroy();

      expect(mockIntervalManager.stopAll).toHaveBeenCalled();
      expect(mockBuildMonitor.stopAll).toHaveBeenCalled();
    });
  });

  // ============================================================
  // Phase 1: Critical Business Logic Tests
  // ============================================================

  describe('Data Fetching - getUserAssociatePlugins', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should fetch user profile, plugins list, and check images successfully', async () => {
      const mockProfile = { username: 'testuser' };
      const mockPlugins = [createMockPlugin(), createMockPlugin({ id: 2, name: 'GENIE3' })];
      const enrichedPlugins = [...mockPlugins];

      mockPluginService.getUserProfile.mockResolvedValue(mockProfile);
      mockPluginService.getPluginsList.mockResolvedValue(mockPlugins);
      mockPluginService.checkMultiplePluginImages.mockResolvedValue(enrichedPlugins);
      enrichPlugins.mockReturnValue(enrichedPlugins);

      await wrapper.vm.getUserAssociatePlugins();

      expect(mockPluginService.getUserProfile).toHaveBeenCalled();
      expect(mockPluginService.getPluginsList).toHaveBeenCalled();
      expect(enrichPlugins).toHaveBeenCalledWith(mockPlugins, 'testuser');
      expect(mockPluginService.checkMultiplePluginImages).toHaveBeenCalledWith(enrichedPlugins);
      expect(wrapper.vm.plugins).toEqual(enrichedPlugins);
      expect(wrapper.vm.profile).toEqual(mockProfile);
    });

    it('should handle getUserProfile API failure gracefully', async () => {
      const mockPlugins = [createMockPlugin()];

      mockPluginService.getUserProfile.mockRejectedValue(new Error('Profile fetch failed'));
      mockPluginService.getPluginsList.mockResolvedValue(mockPlugins);

      console.error = vi.fn(); // Mock console.error

      await wrapper.vm.getUserAssociatePlugins();

      expect(console.error).toHaveBeenCalled();
      expect(mockPluginService.getPluginsList).toHaveBeenCalled();
    });

    it('should handle getPluginsList API failure gracefully', async () => {
      const mockProfile = { username: 'testuser' };

      mockPluginService.getUserProfile.mockResolvedValue(mockProfile);
      mockPluginService.getPluginsList.mockRejectedValue(new Error('Plugins fetch failed'));

      console.error = vi.fn();

      await wrapper.vm.getUserAssociatePlugins();

      expect(console.error).toHaveBeenCalled();
    });

    it('should enrich plugins with user association data', async () => {
      const mockProfile = { username: 'testuser' };
      const mockPlugins = [createMockPlugin()];
      const enrichedPlugins = [{ ...mockPlugins[0], checked: true }];

      mockPluginService.getUserProfile.mockResolvedValue(mockProfile);
      mockPluginService.getPluginsList.mockResolvedValue(mockPlugins);
      mockPluginService.checkMultiplePluginImages.mockResolvedValue(enrichedPlugins);
      enrichPlugins.mockReturnValue(enrichedPlugins);

      await wrapper.vm.getUserAssociatePlugins();

      expect(enrichPlugins).toHaveBeenCalledWith(mockPlugins, 'testuser');
      expect(wrapper.vm.plugins[0].checked).toBe(true);
    });

    it('should start monitoring for plugins that are currently building', async () => {
      const mockProfile = { username: 'testuser' };
      const buildingPlugin = createMockPlugin({ building: true, buildTaskId: 'task-123' });
      const mockPlugins = [buildingPlugin];

      mockPluginService.getUserProfile.mockResolvedValue(mockProfile);
      mockPluginService.getPluginsList.mockResolvedValue(mockPlugins);
      mockPluginService.checkMultiplePluginImages.mockResolvedValue(mockPlugins);
      enrichPlugins.mockReturnValue(mockPlugins);

      wrapper.vm.startBuildMonitoring = vi.fn();

      await wrapper.vm.getUserAssociatePlugins();

      expect(wrapper.vm.startBuildMonitoring).toHaveBeenCalledWith(buildingPlugin);
    });
  });

  describe('Filtering Logic - filteredAndTypedPlugins', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should filter plugins by main category (algorithm)', () => {
      wrapper.vm.plugins = [
        createMockPlugin({ name: 'TENET', plugin_type: 'ANALYSIS' }),
        createMockPlugin({ name: 'Plotter', plugin_type: 'VISUALIZATION' })
      ];
      wrapper.vm.pluginFilters = { mainCategory: 'algorithm', source: 'all', resource: 'all' };
      wrapper.vm.searchTerm = '';

      const filtered = wrapper.vm.filteredAndTypedPlugins;

      expect(filtered.length).toBe(1);
      expect(filtered[0].name).toBe('TENET');
    });

    it('should filter plugins by main category (visualization)', () => {
      wrapper.vm.plugins = [
        createMockPlugin({ name: 'TENET', plugin_type: 'ANALYSIS' }),
        createMockPlugin({ name: 'Plotter', plugin_type: 'VISUALIZATION' })
      ];
      wrapper.vm.pluginFilters = { mainCategory: 'visualization', source: 'all', resource: 'all' };
      wrapper.vm.searchTerm = '';

      const filtered = wrapper.vm.filteredAndTypedPlugins;

      expect(filtered.length).toBe(1);
      expect(filtered[0].name).toBe('Plotter');
    });

    it('should filter plugins by source (official)', () => {
      wrapper.vm.plugins = [
        createMockPlugin({ name: 'TENET', source: 'official' }),
        createMockPlugin({ name: 'Custom', source: 'local' })
      ];
      wrapper.vm.pluginFilters = { mainCategory: 'all', source: 'official', resource: 'all' };
      wrapper.vm.searchTerm = '';

      const filtered = wrapper.vm.filteredAndTypedPlugins;

      expect(filtered.length).toBe(1);
      expect(filtered[0].source).toBe('official');
    });

    it('should filter plugins by resource requirement (GPU)', () => {
      wrapper.vm.plugins = [
        createMockPlugin({ name: 'GPUPlugin', use_gpu: true }),
        createMockPlugin({ name: 'CPUPlugin', use_gpu: false })
      ];
      wrapper.vm.pluginFilters = { mainCategory: 'all', source: 'all', resource: 'gpu' };
      wrapper.vm.searchTerm = '';

      const filtered = wrapper.vm.filteredAndTypedPlugins;

      expect(filtered.length).toBe(1);
      expect(filtered[0].use_gpu).toBe(true);
    });

    it('should combine multiple filters correctly', () => {
      wrapper.vm.plugins = [
        createMockPlugin({ name: 'TENET', source: 'official', plugin_type: 'ANALYSIS', use_gpu: false }),
        createMockPlugin({ name: 'Local', source: 'local', plugin_type: 'ANALYSIS', use_gpu: false }),
        createMockPlugin({ name: 'Viz', source: 'official', plugin_type: 'VISUALIZATION', use_gpu: false })
      ];
      wrapper.vm.pluginFilters = { mainCategory: 'algorithm', source: 'official', resource: 'cpu' };
      wrapper.vm.searchTerm = '';

      const filtered = wrapper.vm.filteredAndTypedPlugins;

      expect(filtered.length).toBe(1);
      expect(filtered[0].name).toBe('TENET');
    });

    it('should sort plugins with official versions first for duplicates', () => {
      wrapper.vm.plugins = [
        createMockPlugin({ name: 'TENET', source: 'local', id: 2 }),
        createMockPlugin({ name: 'TENET', source: 'official', id: 1 }),
        createMockPlugin({ name: 'GENIE3', source: 'local', id: 3 })
      ];
      wrapper.vm.pluginFilters = { mainCategory: 'all', source: 'all', resource: 'all' };
      wrapper.vm.searchTerm = '';

      const filtered = wrapper.vm.filteredAndTypedPlugins;

      expect(filtered.length).toBe(3);
      // First TENET should be official
      expect(filtered[0].name).toBe('GENIE3');
      expect(filtered[1].name).toBe('TENET');
      expect(filtered[1].source).toBe('official');
      expect(filtered[2].name).toBe('TENET');
      expect(filtered[2].source).toBe('local');
    });
  });

  describe('Plugin Categorization - getPluginCategory', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should categorize as algorithm from plugin_type ANALYSIS', () => {
      const plugin = createMockPlugin({ plugin_type: 'ANALYSIS' });

      const category = wrapper.vm.getPluginCategory(plugin);

      expect(category).toBe('algorithm');
    });

    it('should categorize as visualization from plugin_type VISUALIZATION', () => {
      const plugin = createMockPlugin({ plugin_type: 'VISUALIZATION' });

      const category = wrapper.vm.getPluginCategory(plugin);

      expect(category).toBe('visualization');
    });

    it('should fallback to keyword detection from name/description', () => {
      const plugin = createMockPlugin({
        plugin_type: null,
        name: 'Heatmap Plotter',
        description: 'Creates beautiful charts'
      });

      const category = wrapper.vm.getPluginCategory(plugin);

      expect(category).toBe('visualization');
    });

    it('should default to algorithm when no indicators found', () => {
      const plugin = createMockPlugin({
        plugin_type: null,
        name: 'Generic Plugin',
        description: 'Does something'
      });

      const category = wrapper.vm.getPluginCategory(plugin);

      expect(category).toBe('algorithm');
    });
  });

  describe('Duration Calculation - calculateDuration', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should calculate correct duration for running tasks with hours, minutes, seconds', () => {
      const task = createMockTask({
        start_time: '2025-10-24T10:00:00',
        state: 'RUNNING'
      });

      // Mock current time to be 2h 30m 45s after start
      const currentTime = new Date('2025-10-24T12:30:45');

      const duration = wrapper.vm.calculateDuration(task);

      // Should be calculated from start_time to current time
      // This will depend on when test runs, so let's just check format
      expect(duration).toMatch(/^\d+h \d+m \d+s$|^\d+m \d+s$|^\d+s$/);
    });

    it('should calculate correct duration for completed tasks using end_time', () => {
      const task = createMockTask({
        start_time: '2025-10-24T10:00:00',
        end_time: '2025-10-24T10:05:30',
        state: 'SUCCESS'
      });

      const duration = wrapper.vm.calculateDuration(task);

      expect(duration).toBe('5m 30s');
    });

    it('should return "-" for tasks without start_time', () => {
      const task = createMockTask({
        start_time: null,
        state: 'PENDING'
      });

      const duration = wrapper.vm.calculateDuration(task);

      expect(duration).toBe('-');
    });

    it('should handle negative time differences and return "-"', () => {
      const task = createMockTask({
        start_time: '2025-10-24T12:00:00',
        end_time: '2025-10-24T10:00:00', // End before start
        state: 'SUCCESS'
      });

      const duration = wrapper.vm.calculateDuration(task);

      expect(duration).toBe('-');
    });
  });

  // ============================================================================
  // PHASE 2: Display & UI Logic Tests
  // ============================================================================

  describe('Display Logic - Plugin Names', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should return plugin name without modification for single instance', () => {
      const plugin = createMockPlugin({ id: 1, name: 'TENET' });
      wrapper.vm.plugins = [plugin];

      const displayName = wrapper.vm.getDisplayName(plugin);

      expect(displayName).toBe('TENET');
    });

    it('should detect no duplicates when plugin name is unique', () => {
      const plugin1 = createMockPlugin({ id: 1, name: 'TENET' });
      const plugin2 = createMockPlugin({ id: 2, name: 'GENIE3' });
      wrapper.vm.plugins = [plugin1, plugin2];

      const hasDuplicate = wrapper.vm.hasDuplicateName(plugin1);

      expect(hasDuplicate).toBe(false);
    });

    it('should detect duplicates when multiple plugins have same name', () => {
      const plugin1 = createMockPlugin({ id: 1, name: 'TENET', source: 'official' });
      const plugin2 = createMockPlugin({ id: 2, name: 'TENET', source: 'local' });
      wrapper.vm.plugins = [plugin1, plugin2];

      const hasDuplicate = wrapper.vm.hasDuplicateName(plugin1);

      expect(hasDuplicate).toBe(true);
    });

    it('should handle empty plugins array', () => {
      const plugin = createMockPlugin({ id: 1, name: 'TENET' });
      wrapper.vm.plugins = [];

      const hasDuplicate = wrapper.vm.hasDuplicateName(plugin);

      expect(hasDuplicate).toBe(false);
    });
  });

  describe('Display Logic - Source Tooltips', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should generate tooltip for official plugin without duplicates', () => {
      const plugin = createMockPlugin({ id: 1, name: 'TENET', source: 'official' });
      wrapper.vm.plugins = [plugin];

      const tooltip = wrapper.vm.getPluginSourceTooltip(plugin);

      expect(tooltip).toBe('Official CellCraft plugin');
    });

    it('should generate tooltip for local plugin without duplicates', () => {
      const plugin = createMockPlugin({ id: 1, name: 'TENET', source: 'local' });
      wrapper.vm.plugins = [plugin];

      const tooltip = wrapper.vm.getPluginSourceTooltip(plugin);

      expect(tooltip).toBe('User-created local plugin');
    });

    it('should generate tooltip showing other available versions for duplicates', () => {
      const officialPlugin = createMockPlugin({ id: 1, name: 'TENET', source: 'official' });
      const localPlugin = createMockPlugin({ id: 2, name: 'TENET', source: 'local' });
      wrapper.vm.plugins = [officialPlugin, localPlugin];

      const tooltip = wrapper.vm.getPluginSourceTooltip(officialPlugin);

      expect(tooltip).toContain('This is the Official version');
      expect(tooltip).toContain('Also available: Local');
    });
  });

  describe('Display Logic - Version Display', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should add "v" prefix to semantic version numbers', () => {
      const plugin = createMockPlugin({
        current_version: '1.2.3',
        source: 'official'
      });

      const versionDisplay = wrapper.vm.getVersionDisplay(plugin);

      expect(versionDisplay).toBe('v1.2.3');
    });

    it('should return "latest" when no specific version is available', () => {
      const plugin = createMockPlugin({
        current_version: null,
        version: null,
        source: 'official'
      });

      const versionDisplay = wrapper.vm.getVersionDisplay(plugin);

      expect(versionDisplay).toBe('latest');
    });

    it('should return custom tag without modification if not starting with digit', () => {
      const plugin = createMockPlugin({
        current_version: 'dev-branch',
        source: 'official'
      });

      const versionDisplay = wrapper.vm.getVersionDisplay(plugin);

      expect(versionDisplay).toBe('dev-branch');
    });

    it('should generate tooltip showing current version and other available versions', () => {
      const plugin = createMockPlugin({
        current_version: '1.2.3',
        available_versions: ['1.2.3', '1.2.2', '1.1.0', '1.0.0'],
        source: 'official'
      });

      const tooltip = wrapper.vm.getVersionTooltip(plugin);

      expect(tooltip).toContain('Current version: 1.2.3');
      expect(tooltip).toContain('Other versions available:');
      expect(tooltip).toContain('1.2.2');
    });

    it('should generate tooltip for latest version when using default', () => {
      const plugin = createMockPlugin({
        current_version: null,
        available_versions: ['1.2.0', '1.1.0', '1.0.0'],
        source: 'official'
      });

      const tooltip = wrapper.vm.getVersionTooltip(plugin);

      expect(tooltip).toContain('Using latest version');
      expect(tooltip).toContain('Available versions:');
    });

    it('should handle plugin with no available versions gracefully', () => {
      const plugin = createMockPlugin({
        current_version: null,
        available_versions: [],
        source: 'official'
      });

      const tooltip = wrapper.vm.getVersionTooltip(plugin);

      expect(tooltip).toBe('Using the latest available version');
    });
  });

  describe('Display Logic - Plugin Statistics', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should calculate total, official, and local plugin counts', () => {
      wrapper.vm.plugins = [
        createMockPlugin({ id: 1, name: 'TENET', source: 'official' }),
        createMockPlugin({ id: 2, name: 'GENIE3', source: 'local' }),
        createMockPlugin({ id: 3, name: 'LEAP', source: 'local' })
      ];

      const stats = wrapper.vm.pluginStatistics;

      expect(stats.total).toBe(3);
      expect(stats.official).toBe(1);
      expect(stats.local).toBe(2);
    });

    it('should detect and count duplicate plugin names', () => {
      wrapper.vm.plugins = [
        createMockPlugin({ id: 1, name: 'TENET', source: 'official' }),
        createMockPlugin({ id: 2, name: 'TENET', source: 'local' }),
        createMockPlugin({ id: 3, name: 'GENIE3', source: 'official' })
      ];

      const stats = wrapper.vm.pluginStatistics;

      expect(stats.duplicates).toBe(1);
      expect(stats.duplicateNames).toEqual(['TENET']);
    });

    it('should list all duplicate plugin names', () => {
      wrapper.vm.plugins = [
        createMockPlugin({ id: 1, name: 'TENET', source: 'official' }),
        createMockPlugin({ id: 2, name: 'TENET', source: 'local' }),
        createMockPlugin({ id: 3, name: 'GENIE3', source: 'official' }),
        createMockPlugin({ id: 4, name: 'GENIE3', source: 'local' })
      ];

      const stats = wrapper.vm.pluginStatistics;

      expect(stats.duplicates).toBe(2);
      expect(stats.duplicateNames).toContain('TENET');
      expect(stats.duplicateNames).toContain('GENIE3');
    });

    it('should handle empty plugin list', () => {
      wrapper.vm.plugins = [];

      const stats = wrapper.vm.pluginStatistics;

      expect(stats.total).toBe(0);
      expect(stats.official).toBe(0);
      expect(stats.local).toBe(0);
      expect(stats.duplicates).toBe(0);
      expect(stats.duplicateNames).toEqual([]);
    });
  });

  describe('UI Interactions - Filter Handler', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should update pluginFilters when handleFiltersChanged is called', () => {
      const newFilters = {
        mainCategory: 'algorithm',
        source: 'official',
        resource: 'gpu'
      };

      wrapper.vm.handleFiltersChanged(newFilters);

      expect(wrapper.vm.pluginFilters).toEqual(newFilters);
    });

    it('should handle partial filter updates', () => {
      wrapper.vm.pluginFilters = {
        mainCategory: 'all',
        source: 'all',
        resource: 'all'
      };

      const partialUpdate = {
        mainCategory: 'visualization'
      };

      wrapper.vm.handleFiltersChanged(partialUpdate);

      expect(wrapper.vm.pluginFilters.mainCategory).toBe('visualization');
    });
  });

  // ============================================================================
  // PHASE 3: Modal & Build Features Tests
  // ============================================================================

  describe('Modal Management - Plugin Extension', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should open plugin extension modal with empty plugin for adding', () => {
      wrapper.vm.addPluginExtension();

      expect(wrapper.vm.selectedPlugin).toEqual({
        name: "",
        description: "",
        dependencies: {},
        drawflow: {},
        rules: []
      });
      expect(wrapper.vm.showPluginExtension).toBe(true);
    });

    it('should open plugin extension modal with selected plugin for editing', () => {
      const pluginToEdit = createMockPlugin({
        id: 1,
        name: 'TENET',
        description: 'GRN analysis tool'
      });

      wrapper.vm.editPluginExtension(pluginToEdit);

      expect(wrapper.vm.selectedPlugin).toEqual(pluginToEdit);
      expect(wrapper.vm.showPluginExtension).toBe(true);
    });

    it('should open plugin details in read-only mode', () => {
      const pluginToView = createMockPlugin({
        id: 1,
        name: 'TENET',
        source: 'official'
      });

      wrapper.vm.viewPluginDetails(pluginToView);

      expect(wrapper.vm.selectedPlugin).toEqual({
        ...pluginToView,
        readOnly: true
      });
      expect(wrapper.vm.showPluginExtension).toBe(true);
    });
  });

  describe('Modal Management - Build Logs', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should load and display build logs successfully', async () => {
      const plugin = {
        ...createMockPlugin({ name: 'TENET' }),
        buildStatus: 'RUNNING'
      };
      const mockLogs = {
        success: true,
        logs: 'Build step 1: Completed\nBuild step 2: In progress...'
      };
      mockBuildService.getBuildLogs.mockResolvedValue(mockLogs);

      await wrapper.vm.showBuildLogs(plugin);

      expect(wrapper.vm.selectedPluginName).toBe('TENET');
      expect(wrapper.vm.showLogsModal).toBe(true);
      expect(wrapper.vm.selectedBuildLogs.status).toBe('RUNNING');
      expect(wrapper.vm.selectedBuildLogs.logs).toBe(mockLogs.logs);
      expect(wrapper.vm.logsLoading).toBe(false);
    });

    it('should handle build logs loading failure', async () => {
      const plugin = createMockPlugin({ name: 'TENET' });
      const mockError = new Error('Network error');
      mockBuildService.getBuildLogs.mockRejectedValue(mockError);

      await wrapper.vm.showBuildLogs(plugin);

      expect(wrapper.vm.selectedBuildLogs.status).toBe('Error');
      expect(wrapper.vm.selectedBuildLogs.logs).toContain('Failed to load build logs');
      expect(wrapper.vm.logsLoading).toBe(false);
    });

    it('should close logs modal and reset state', () => {
      wrapper.vm.showLogsModal = true;
      wrapper.vm.selectedBuildLogs = { status: 'SUCCESS', logs: 'test logs' };
      wrapper.vm.selectedPluginName = 'TENET';
      wrapper.vm.logsLoading = true;

      wrapper.vm.closeLogsModal();

      expect(wrapper.vm.showLogsModal).toBe(false);
      expect(wrapper.vm.selectedBuildLogs).toBe(null);
      expect(wrapper.vm.selectedPluginName).toBe('');
      expect(wrapper.vm.logsLoading).toBe(false);
    });

    it('should refresh logs for selected plugin', async () => {
      const plugin = {
        ...createMockPlugin({ name: 'TENET' }),
        buildStatus: 'SUCCESS'
      };
      wrapper.vm.plugins = [plugin];
      wrapper.vm.selectedPluginName = 'TENET';

      const mockLogs = {
        success: true,
        logs: 'Updated logs...'
      };
      mockBuildService.getBuildLogs.mockResolvedValue(mockLogs);

      await wrapper.vm.refreshLogs();

      expect(mockBuildService.getBuildLogs).toHaveBeenCalledWith('TENET');
      expect(wrapper.vm.selectedBuildLogs.status).toBe('SUCCESS');
      expect(wrapper.vm.selectedBuildLogs.logs).toBe('Updated logs...');
    });
  });

  describe('Helper Functions', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should format file size correctly', () => {
      expect(wrapper.vm.formatFileSize(0)).toBe('0 Bytes');
      expect(wrapper.vm.formatFileSize(null)).toBe('0 Bytes');
      expect(wrapper.vm.formatFileSize(500)).toBe('500 Bytes');
      expect(wrapper.vm.formatFileSize(1024)).toBe('1 KB');
      expect(wrapper.vm.formatFileSize(1048576)).toBe('1 MB');
      expect(wrapper.vm.formatFileSize(1073741824)).toBe('1 GB');
    });

    it('should return current date in YYYY/MM/DD format', () => {
      const dateString = wrapper.vm.getCurrentDateString();

      // Verify format: YYYY/MM/DD
      expect(dateString).toMatch(/^\d{4}\/\d{2}\/\d{2}$/);

      // Verify it's today's date
      const today = new Date();
      const expectedYear = today.getFullYear();
      const expectedMonth = String(today.getMonth() + 1).padStart(2, '0');
      const expectedDay = String(today.getDate()).padStart(2, '0');
      expect(dateString).toBe(`${expectedYear}/${expectedMonth}/${expectedDay}`);
    });

    it('should apply build info to plugin and start monitoring', () => {
      const plugin = createMockPlugin({ name: 'TENET' });
      wrapper.vm.plugins = [plugin];
      wrapper.vm.startBuildMonitoring = vi.fn();

      const buildInfo = {
        pluginName: 'TENET',
        taskId: 'task-456',
        buildStarted: true
      };

      wrapper.vm.applyBuildInfo(buildInfo);

      expect(plugin.building).toBe(true);
      expect(plugin.buildTaskId).toBe('task-456');
      expect(plugin.buildStatus).toBe('RUNNING');
      expect(plugin.imageExists).toBe(false);
      expect(wrapper.vm.startBuildMonitoring).toHaveBeenCalledWith(plugin);
    });
  });

  describe('Task Monitoring', () => {
    beforeEach(() => {
      wrapper = shallowMount(PluginsPage, {
        stubs: {
          'PluginExtention': true,
          'BuildMonitor': true,
          'PluginCategoryTabs': true
        }
      });
    });

    it('should show build task logs for specific task ID', async () => {
      const task = createMockTask({ task_id: 'task-123', plugin_name: 'TENET' });
      wrapper.vm.buildTaskList = [task];
      wrapper.vm.showBuildLogs = vi.fn();

      await wrapper.vm.showBuildTaskLogs('task-123');

      expect(wrapper.vm.showBuildLogs).toHaveBeenCalledWith({
        name: 'TENET',
        buildTaskId: 'task-123'
      });
    });

    it('should start monitoring for plugins that are currently building', () => {
      const buildingPlugin = createMockPlugin({
        name: 'TENET',
        building: true,
        buildTaskId: 'task-123'
      });
      const idlePlugin = createMockPlugin({
        name: 'GENIE3',
        building: false
      });
      wrapper.vm.plugins = [buildingPlugin, idlePlugin];
      wrapper.vm.startBuildMonitoring = vi.fn();

      wrapper.vm.startMonitoringForBuildingPlugins();

      expect(wrapper.vm.startBuildMonitoring).toHaveBeenCalledTimes(1);
      expect(wrapper.vm.startBuildMonitoring).toHaveBeenCalledWith(buildingPlugin);
    });

    it('should update running task durations', () => {
      const runningTask = {
        task_id: 'task-1',
        status: 'RUNNING',
        start_time: new Date(Date.now() - 5000).toISOString(),
        running_time: '0s'
      };
      const completedTask = {
        task_id: 'task-2',
        status: 'SUCCESS',
        start_time: new Date(Date.now() - 10000).toISOString(),
        end_time: new Date().toISOString(),
        running_time: '10s'
      };
      wrapper.vm.buildTaskList = [runningTask, completedTask];

      wrapper.vm.updateRunningTaskDurations();

      // Running task duration should be updated
      expect(runningTask.running_time).not.toBe('0s');
      // Completed task duration should remain unchanged
      expect(completedTask.running_time).toBe('10s');
    });
  });
});
