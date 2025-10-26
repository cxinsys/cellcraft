<template>
  <div class="layout">
    <div class="first-line">
      <div class="first-line__left">
        <div class="header__text">
          Plugins
          <!-- <div class="header__desc">
            Plugin for data analysis algorithm extensions
          </div> -->
        </div>
      </div>
      <div class="first-line__right">
        <button type="button" class="add__button" @click="addPluginExtension">
          <img class="add__button--icon" src="@/assets/add_circle.png" alt="Add Plugin" />
          <span>Add Plugin</span>
        </button>
        <button type="button" class="monitor__button" @click="toggleBuildMonitor">
          <img class="monitor__button--icon" src="@/assets/control_jobs.png" alt="Build Monitor" />
          <span>Build Monitor</span>
        </button>
        <!-- <div class="build-all__button" @click="buildAllPlugins">
          <img class="build-all__button--icon" src="@/assets/add_circle.png" />
          <h1>Build All</h1>
        </div> -->
        <div class="search">
          <label for="plugin-search" class="visually-hidden">Search plugins</label>
          <input type="text" id="plugin-search" v-model="searchTerm" placeholder="Search plugins..." />
        </div>
      </div>
    </div>
    <PluginExtention v-if="showPluginExtension" @close="closePluginExtension" :editName="selectedPlugin.name"
      :editDescription="selectedPlugin.description" :editDependencies="selectedPlugin.dependencies"
      :editDrawflow="selectedPlugin.drawflow" :editRules="selectedPlugin.rules" :readOnly="selectedPlugin.readOnly"
      :pluginSource="selectedPlugin.source" />
    <BuildMonitor v-if="showBuildMonitor" :show_monitor="showBuildMonitor" :buildTaskList="buildTaskList"
      @cancel-task="cancelBuildTask" @show-logs="showBuildTaskLogs" @close="toggleBuildMonitor" />
    <!-- Enhanced Plugin Category Tabs -->
    <PluginCategoryTabs :plugins="plugins" @filters-changed="handleFiltersChanged" />

    <table>
      <tbody>
        <tr v-for="plugin in filteredAndTypedPlugins" :key="plugin.id">
          <td>
            <div class="plugin-container">
              <div class="title-container">
                {{ getDisplayName(plugin) }}
                <!-- Official/Local 배지 추가 -->
                <span class="plugin-badge" :class="plugin.source === 'official' ? 'badge-official' : 'badge-local'"
                  :title="getPluginSourceTooltip(plugin)">
                  {{ plugin.source === 'official' ? 'OFFICIAL' : 'LOCAL' }}
                </span>
                <!-- Version badge (only for official plugins) -->
                <span v-if="plugin.source === 'official'" class="version-badge version-official"
                  :title="getVersionTooltip(plugin)">
                  {{ getVersionDisplay(plugin) }}
                </span>
                <!-- Duplicate plugin indicator -->
                <span v-if="hasDuplicateName(plugin)" class="duplicate-indicator"
                  :title="`Multiple versions of '${plugin.name}' are available`">
                  <img src="@/assets/view_gray.png" alt="Multiple versions" class="duplicate-icon" />
                </span>
                <!-- GPU indicator -->
                <span v-if="plugin.use_gpu" class="gpu-indicator"
                  :title="'This plugin requires GPU acceleration'">
                  <span class="gpu-text">GPU</span>
                </span>
              </div>
              <div class="description-container">
                {{ plugin.description }}
              </div>
              
              <!-- Removed version selector - moved to admin page -->
              
              <div class="lastUpdated-container">
                Last Updated: {{ plugin.updated_at.split('T')[0] }}
              </div>
            </div>
            <div class="option-container">
              <!-- 첫 번째 Row: 설정 버튼과 스위치 -->
              <div class="controls-row">
                <!-- Official 플러그인은 설정 버튼 숨김 -->
                <button v-if="plugin.source !== 'official'" type="button" class="setting"
                  @click="editPluginExtension(plugin)" :aria-label="`Edit ${plugin.name} plugin settings`">
                  <img class="setting__button--icon" src="@/assets/settings.png" alt="Settings" />
                </button>
                <!-- Official 플러그인은 읽기 전용 버튼 표시 -->
                <button v-else type="button" class="setting readonly" @click="viewPluginDetails(plugin)"
                  :aria-label="`View ${plugin.name} plugin details (read-only)`">
                  <img class="setting__button--icon" src="@/assets/view_gray.png" alt="View Details" />
                </button>
                <label class="switch">
                  <input v-model="plugin.checked" type="checkbox" @change="handlePluginAssociate(plugin)"
                    :disabled="isCheckboxDisabled" :aria-label="`Enable or disable ${getDisplayName(plugin)} plugin`" />
                  <span class="slider round"></span>
                </label>
              </div>

              <!-- 두 번째 Row: 빌드 버튼 -->
              <div class="build-row">
                <!-- Official 플러그인은 빌드 버튼 숨김 -->
                <button v-if="plugin.source !== 'official'" class="build-button" @click="handleBuildPlugin(plugin)"
                  :disabled="plugin.imageExists"
                  :class="{ 'building': plugin.building, 'image-exists': plugin.imageExists }">
                  <div v-if="plugin.building" class="building-content">
                    <div class="loading-spinner"></div>
                    <span>View Logs</span>
                  </div>
                  <span v-else-if="plugin.imageExists">Built</span>
                  <span v-else>Build</span>
                </button>
                <!-- Official 플러그인은 Pre-built 표시 -->
                <div v-else class="official-status">
                  <span class="pre-built">Pre-built</span>
                </div>
              </div>
            </div>
          </td>
        </tr>
      </tbody>
    </table>

    <!-- Removed version update modal - moved to admin page -->

    <!-- 빌드 로그 모달 -->
    <div v-if="showLogsModal" class="logs-modal-overlay" @click.self="closeLogsModal">
      <div class="logs-modal">
        <div class="logs-modal-header">
          <h3>Build Logs - {{ selectedPluginName }}</h3>
          <div class="logs-modal-controls">
            <button @click="refreshLogs" :disabled="logsLoading" class="refresh-btn">
              <img src="@/assets/refresh.png" alt="Refresh" class="refresh-icon" />
            </button>
            <button @click="closeLogsModal" class="close-btn">
              <img src="@/assets/close.png" alt="Close" class="close-icon" />
            </button>
          </div>
        </div>

        <div v-if="logsLoading" class="logs-loading">
          Loading build logs...
        </div>

        <div v-else-if="selectedBuildLogs" class="logs-content">
          <div class="logs-task-info">
            <p><strong>Plugin Name:</strong> {{ selectedPluginName }}</p>
            <p><strong>Build Status:</strong> {{ selectedBuildLogs.status || 'N/A' }}</p>
            <p><strong>Last Updated:</strong> {{ selectedBuildLogs.timestamp || 'N/A' }}</p>
          </div>

          <div v-if="!selectedBuildLogs.logs || selectedBuildLogs.logs.length === 0" class="no-logs">
            No build logs available yet.
          </div>

          <div v-else class="logs-files">
            <div class="log-file">
              <div class="log-file-header">
                <h4>Build Log</h4>
                <span class="log-file-size">{{ formatFileSize(selectedBuildLogs.logs.length) }}</span>
              </div>
              <pre class="log-file-content">{{ selectedBuildLogs.logs }}</pre>
            </div>
          </div>
        </div>

        <div v-else class="no-logs">
          No build logs available.
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import PluginExtention from "@/components/PluginExtention.vue";
import BuildMonitor from "@/components/pluginComponents/BuildMonitor.vue";
import PluginCategoryTabs from "@/components/pluginComponents/PluginCategoryTabs.vue";

// Services and utilities
import { PluginService } from "@/services/pluginService";
import { BuildMonitor as BuildMonitorService, BuildStatus } from "@/services/buildMonitorService";
import { DialogService } from "@/services/dialogService";
import { NotificationService } from "@/services/notificationService";
import { BuildService } from "@/services/buildService";
import { IntervalManager } from "@/utils/intervalManager";
import { enrichPlugins } from "@/utils/pluginDataProcessor";
import { transformBuildTasks } from "@/utils/buildTaskMapper";
import { analyzeError } from "@/utils/errorAnalyzer";

export default {
  components: {
    PluginExtention,
    BuildMonitor,
    PluginCategoryTabs,
  },
  data() {
    return {
      // Service instances
      pluginService: new PluginService(),
      intervalManager: new IntervalManager(),
      buildMonitor: null, // Will be initialized in created()
      dialogService: new DialogService(),
      notificationService: new NotificationService(),
      buildService: new BuildService(),

      showPluginExtension: false,
      showBuildMonitor: false,
      buildTaskList: [],
      isCheckboxDisabled: false,
      searchTerm: "",
      filterType: "all", // Legacy filter type for backward compatibility
      // Enhanced filtering
      pluginFilters: {
        mainCategory: 'all',
        source: 'all',
        resource: 'all'
      },
      plugins: [],
      profile: {},
      selectedPlugin: {
        name: "",
        description: "",
        dependencies: {},
        drawflow: {},
        rules: [],
      },
      // 빌드 로그 모달 관련
      showLogsModal: false,
      selectedBuildLogs: null,
      selectedPluginName: "",
      logsLoading: false,
      // Error and loading states
      isLoadingPlugins: false,
      pluginError: null,
      // Network status tracking
      isOffline: false,
      retryAttempts: 0,
      maxRetryAttempts: 3
    };
  },
  created() {
    // Initialize BuildMonitor with IntervalManager
    this.buildMonitor = new BuildMonitorService(this.intervalManager);
  },
  async mounted() {
    // Check network status
    this.isOffline = !navigator.onLine;
    window.addEventListener('online', this.handleOnline);
    window.addEventListener('offline', this.handleOffline);
    
    // Load plugins data
    try {
      await this.getUserAssociatePlugins();
    } catch (error) {
      console.error('Failed to load plugins on mount:', error);
      this.handleLoadingError(error);
    }
  },
  beforeDestroy() {
    // Clean up event listeners
    window.removeEventListener('online', this.handleOnline);
    window.removeEventListener('offline', this.handleOffline);

    // Stop all build monitoring
    if (this.buildMonitor) {
      this.buildMonitor.stopAll();
    }

    // Stop all intervals
    if (this.intervalManager) {
      this.intervalManager.stopAll();
    }
  },
  computed: {
    filteredPlugins() {
      return this.plugins.filter((plugin) =>
        plugin.name.toLowerCase().includes(this.searchTerm.toLowerCase())
      );
    },
    filteredAndTypedPlugins() {
      let filtered = this.filteredPlugins;

      // Apply enhanced filtering from PluginCategoryTabs
      if (this.pluginFilters.mainCategory !== 'all') {
        filtered = filtered.filter(plugin => {
          const category = this.getPluginCategory(plugin);
          return category === this.pluginFilters.mainCategory;
        });
      }

      // Apply source filtering
      if (this.pluginFilters.source !== 'all') {
        filtered = filtered.filter(plugin => {
          const source = plugin.source || 'local';
          return source === this.pluginFilters.source;
        });
      }

      // Apply resource filtering
      if (this.pluginFilters.resource !== 'all') {
        filtered = filtered.filter(plugin => {
          const useGpu = plugin.use_gpu || false;
          if (this.pluginFilters.resource === 'gpu') return useGpu;
          if (this.pluginFilters.resource === 'cpu') return !useGpu;
          return true;
        });
      }

      // Legacy filter type support for backward compatibility
      if (this.filterType !== 'all') {
        if (this.filterType === 'official') {
          filtered = filtered.filter(plugin => plugin.source === 'official');
        } else if (this.filterType === 'local') {
          filtered = filtered.filter(plugin => plugin.source === 'local' || !plugin.source);
        }
      }

      // Sort plugins to show official versions first when there are duplicates
      return filtered.sort((a, b) => {
        // If same name, prioritize official over local
        if (a.name === b.name) {
          if (a.source === 'official' && b.source !== 'official') return -1;
          if (b.source === 'official' && a.source !== 'official') return 1;
        }
        // Otherwise maintain alphabetical order by name
        return a.name.localeCompare(b.name);
      });
    },

    // Plugin statistics for better understanding
    pluginStatistics() {
      const stats = {
        total: this.plugins.length,
        official: this.plugins.filter(p => p.source === 'official').length,
        local: this.plugins.filter(p => p.source === 'local' || !p.source).length,
        duplicates: 0,
        duplicateNames: []
      };

      const nameGroups = this.plugins.reduce((acc, plugin) => {
        acc[plugin.name] = (acc[plugin.name] || 0) + 1;
        return acc;
      }, {});

      stats.duplicateNames = Object.keys(nameGroups).filter(name => nameGroups[name] > 1);
      stats.duplicates = stats.duplicateNames.length;

      return stats;
    },
  },
  watch: {
    plugins: {
      handler(newPlugins) {
        if (newPlugins.length > 0) {
          const stats = this.pluginStatistics;
          console.log('Plugin Statistics:', stats);

          if (stats.duplicates > 0) {
            console.log(`Found ${stats.duplicates} plugin names with multiple versions:`, stats.duplicateNames);

            stats.duplicateNames.forEach(name => {
              const duplicates = newPlugins.filter(p => p.name === name);
              console.log(`- "${name}":`, duplicates.map(d => `${d.source} (ID: ${d.id})`).join(', '));
            });
          }
          
        }
      },
      immediate: true
    }
  },
  methods: {
    // Enhanced filtering methods
    handleFiltersChanged(filters) {
      this.pluginFilters = { ...filters };
    },


    // Plugin display name and duplicate handling methods
    getDisplayName(plugin) {
      // const hasDuplicates = this.hasDuplicateName(plugin);
      // if (hasDuplicates) {
      //   return `${plugin.name} (${plugin.source === 'official' ? 'Official' : 'Local'})`;
      // }
      return plugin.name;
    },

    hasDuplicateName(plugin) {
      const pluginsWithSameName = this.plugins.filter(p => p.name === plugin.name);
      return pluginsWithSameName.length > 1;
    },

    getPluginSourceTooltip(plugin) {
      const hasDuplicates = this.hasDuplicateName(plugin);
      if (hasDuplicates) {
        const otherSources = this.plugins
          .filter(p => p.name === plugin.name && p.id !== plugin.id)
          .map(p => p.source === 'official' ? 'Official' : 'Local')
          .join(', ');
        return `This is the ${plugin.source === 'official' ? 'Official' : 'Local'} version. Also available: ${otherSources}`;
      }
      return plugin.source === 'official' ? 'Official CellCraft plugin' : 'User-created local plugin';
    },

    getVersionDisplay(plugin) {
      // This method is now only called for official plugins
      const currentVersion = plugin.current_version || plugin.version;
      
      if (currentVersion && currentVersion !== 'latest') {
        // If it looks like a semantic version (starts with digit), add 'v' prefix
        if (/^\d/.test(currentVersion)) {
          return `v${currentVersion}`;
        }
        // For other tag formats, show as-is
        return currentVersion;
      }
      
      // Fallback to 'latest' if no specific version is available
      return 'latest';
    },

    getVersionTooltip(plugin) {
      // This method is now only called for official plugins
      const currentVersion = plugin.current_version || plugin.version;
      const availableVersions = plugin.available_versions || [];
      
      if (currentVersion && currentVersion !== 'latest') {
        let tooltip = `Current version: ${currentVersion}`;
        if (availableVersions.length > 1) {
          const otherVersions = availableVersions.filter(v => v !== currentVersion).slice(0, 3);
          if (otherVersions.length > 0) {
            tooltip += `\nOther versions available: ${otherVersions.join(', ')}`;
            if (availableVersions.length > 4) {
              tooltip += ` and ${availableVersions.length - 4} more...`;
            }
          }
        }
        return tooltip;
      }
      
      if (availableVersions.length > 0) {
        return `Using latest version. Available versions: ${availableVersions.slice(0, 5).join(', ')}${availableVersions.length > 5 ? '...' : ''}`;
      }
      
      return 'Using the latest available version';
    },

    getPluginCategory(plugin) {
      // First check plugin_type field from backend (ANALYSIS/VISUALIZATION)
      if (plugin.plugin_type) {
        const type = plugin.plugin_type.toLowerCase();
        if (type === 'analysis' || type === 'algorithm') {
          return 'algorithm';
        } else if (type === 'visualization') {
          return 'visualization';
        }
      }

      // Legacy support: check category field
      if (plugin.category) {
        return plugin.category;
      }

      // Fallback: analyze plugin name/description for categorization
      const name = plugin.name.toLowerCase();
      const desc = (plugin.description || '').toLowerCase();

      const visualizationKeywords = ['plot', 'chart', 'graph', 'visual', 'heatmap', 'scatter'];
      const isVisualization = visualizationKeywords.some(keyword =>
        name.includes(keyword) || desc.includes(keyword)
      );

      return isVisualization ? 'visualization' : 'algorithm';
    },

    async closePluginExtension(buildInfo) {
      this.showPluginExtension = false;

      // Reload plugin list after extension is closed
      try {
        await this.getUserAssociatePlugins();

        // Apply build info if a new plugin build was started
        if (buildInfo && buildInfo.buildStarted) {
          this.$nextTick(() => {
            this.applyBuildInfo(buildInfo);
          });
        }
      } catch (error) {
        console.error(error);
      }
    },

    applyBuildInfo(buildInfo) {
      const plugin = this.plugins.find(p => p.name === buildInfo.pluginName);
      if (!plugin) return;

      // Force build state (override API data)
      plugin.building = true;
      plugin.buildTaskId = buildInfo.taskId;
      plugin.buildStatus = 'RUNNING';
      plugin.imageExists = false;

      console.log(`Plugin ${buildInfo.pluginName} build status forced to building`);

      // Start build monitoring
      this.startBuildMonitoring(plugin);
    },
    getCurrentDateString() {
      const today = new Date();
      const year = today.getFullYear();
      // 월은 0부터 시작하므로 1을 더해줍니다. 또한, 월과 일이 10보다 작을 때 앞에 '0'을 붙여줍니다.
      const month = String(today.getMonth() + 1).padStart(2, "0");
      const day = String(today.getDate()).padStart(2, "0");
      // YYYY/MM/DD 형식으로 문자열을 반환합니다.
      return `${year}/${month}/${day}`;
    },
    addPluginExtension() {
      this.selectedPlugin = {
        name: "",
        description: "",
        dependencies: {},
        drawflow: {},
        rules: [],
      };
      this.showPluginExtension = true;
    },
    editPluginExtension(plugin) {
      this.selectedPlugin = plugin;
      this.showPluginExtension = true;
    },
    viewPluginDetails(plugin) {
      // Official 플러그인은 읽기 전용 모드로 표시
      this.selectedPlugin = {
        ...plugin,
        readOnly: true // 읽기 전용 플래그 추가
      };
      this.showPluginExtension = true;
    },
    async getUserAssociatePlugins() {
      try {
        // Get user profile and plugins list using service
        this.profile = await this.pluginService.getUserProfile();
        const plugins = await this.pluginService.getPluginsList();

        // Enrich plugins with user association and build status
        this.plugins = enrichPlugins(plugins, this.profile.username);

        // Check plugin images
        this.plugins = await this.pluginService.checkMultiplePluginImages(this.plugins);

        // Start monitoring for plugins that are currently building
        this.startMonitoringForBuildingPlugins();
      } catch (error) {
        console.error(error);
      }
    },

    startMonitoringForBuildingPlugins() {
      this.plugins.forEach(plugin => {
        if (plugin.building && plugin.buildTaskId) {
          console.log(`Starting monitoring for plugin ${plugin.name} with task ${plugin.buildTaskId}`);
          this.startBuildMonitoring(plugin);
        }
      });
    },
    async handlePluginAssociate(plugin) {
      const pluginId = parseInt(plugin.id);

      try {
        let result;
        if (plugin.checked) {
          result = await this.pluginService.associatePlugin(pluginId);
        } else {
          result = await this.pluginService.dissociatePlugin(pluginId);
        }
        console.log(result);

        // 1초 동안 체크박스 비활성화
        this.isCheckboxDisabled = true;
        setTimeout(() => {
          this.isCheckboxDisabled = false;
        }, 1000);
      } catch (error) {
        console.error('Error associating/disassociating plugin:', error);
      }
    },
    async handleBuildPlugin(plugin) {
      if (plugin.building) {
        // 빌드 중인 경우 로그 보기
        await this.showBuildLogs(plugin);
        return;
      }

      if (plugin.imageExists) {
        return;
      }

      try {
        plugin.building = true;
        const result = await this.buildService.buildPlugin(plugin.name, false); // 기존 플러그인은 기본적으로 GPU 비활성화

        if (!result.success) {
          plugin.building = false;
          this.showErrorNotification(
            `Failed to start build for plugin ${plugin.name}`,
            'Build Failed'
          );
          return;
        }

        // 태스크 ID 저장
        plugin.buildTaskId = result.taskId;
        plugin.buildStatus = 'RUNNING';

        // 빌드 상태 주기적 모니터링 시작
        this.startBuildMonitoring(plugin);

        console.log('Build started:', result.taskId);
        this.showSuccessNotification(`Plugin ${plugin.name} build started!`, 'Build Started');
      } catch (error) {
        console.error(`Error starting build for plugin ${plugin.name}:`, error);
        plugin.building = false;
        this.showErrorNotification(
          `Failed to start build for plugin ${plugin.name}`,
          'Build Failed'
        );
      }
    },
    startBuildMonitoring(plugin) {
      this.buildMonitor.startMonitoring(
        plugin.name,
        plugin.buildTaskId,
        ({ pluginName, status, error }) => {
          // Find the plugin to update
          const pluginToUpdate = this.plugins.find(p => p.name === pluginName);
          if (!pluginToUpdate) return;

          // Update build status
          pluginToUpdate.buildStatus = status;

          console.log(`Build status for ${pluginName}: ${status}`);

          // Handle terminal states
          if (status === BuildStatus.SUCCESS) {
            pluginToUpdate.building = false;
            pluginToUpdate.imageExists = true;
            console.log(`Plugin ${pluginName} built successfully!`);
            this.showBuildNotification(pluginName, 'success');
          } else if (status === BuildStatus.FAILURE || status === BuildStatus.REVOKED) {
            pluginToUpdate.building = false;
            pluginToUpdate.imageExists = false;
            console.error(`Plugin ${pluginName} build failed with status: ${status}`);
            this.showBuildNotification(pluginName, 'failure');
          } else if (status === BuildStatus.ERROR) {
            pluginToUpdate.building = false;
            console.error(`Build monitoring error for ${pluginName}:`, error);
          } else if (status === BuildStatus.RUNNING || status === BuildStatus.PENDING) {
            pluginToUpdate.building = true;
          }
        },
        { interval: 2000 }
      );
    },
    showBuildNotification(pluginName, status) {
      if (status === 'success') {
        // 성공 알림 (브라우저 알림 또는 토스트 메시지)
        console.log(`Plugin ${pluginName} built successfully!`);
        // 필요시 토스트 라이브러리 사용 가능
      } else if (status === 'failure') {
        // 실패 알림
        console.error(`Plugin ${pluginName} build failed!`);
        // 필요시 토스트 라이브러리 사용 가능
      }
    },
    async showBuildLogs(plugin) {
      try {
        this.selectedPluginName = plugin.name;
        this.logsLoading = true;
        this.showLogsModal = true;
        this.selectedBuildLogs = null;

        const result = await this.buildService.getBuildLogs(plugin.name);

        this.selectedBuildLogs = {
          status: plugin.buildStatus || 'Unknown',
          timestamp: new Date().toLocaleString(),
          logs: result.success ? result.logs : result.error
        };

        console.log('Build logs loaded for plugin:', plugin.name);
      } catch (error) {
        console.error(`Error fetching build logs for plugin ${plugin.name}:`, error);
        this.selectedBuildLogs = {
          status: 'Error',
          timestamp: new Date().toLocaleString(),
          logs: 'Failed to load build logs: ' + (error.response?.data?.detail || error.message)
        };
      } finally {
        this.logsLoading = false;
      }
    },
    closeLogsModal() {
      this.showLogsModal = false;
      this.selectedBuildLogs = null;
      this.selectedPluginName = "";
      this.logsLoading = false;
    },
    async refreshLogs() {
      if (this.selectedPluginName) {
        const plugin = this.plugins.find(p => p.name === this.selectedPluginName);
        if (plugin) {
          await this.showBuildLogs(plugin);
        }
      }
    },
    formatFileSize(bytes) {
      if (!bytes || bytes === 0) return '0 Bytes';

      const k = 1024;
      const sizes = ['Bytes', 'KB', 'MB', 'GB'];
      const i = Math.floor(Math.log(bytes) / Math.log(k));

      return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
    },
    async toggleBuildMonitor() {
      this.showBuildMonitor = !this.showBuildMonitor;

      if (this.showBuildMonitor) {
        // 빌드 모니터가 열릴 때 태스크 목록 갱신
        await this.fetchBuildTasks();
        // 주기적으로 태스크 상태 업데이트
        this.startBuildTaskMonitoring();
      } else {
        // 모니터링 중지
        this.stopBuildTaskMonitoring();
      }
    },
    async fetchBuildTasks() {
      try {
        const tasks = await this.pluginService.getBuildTasks();

        // Transform tasks using utility function
        this.buildTaskList = transformBuildTasks(tasks, this.calculateDuration);

        console.log('Build tasks fetched:', this.buildTaskList);
      } catch (error) {
        console.error('Failed to fetch build tasks:', error);
      }
    },
    calculateDuration(task) {
      const startTime = task.start_time ? new Date(task.start_time) : null;
      const endTime = task.end_time ? new Date(task.end_time) : null;

      if (!startTime) return '-';

      let targetTime;
      if (task.state === 'RUNNING' || task.state === 'PENDING') {
        // 실행 중인 태스크는 현재 시간을 기준으로 계산
        targetTime = new Date();
      } else if (endTime) {
        // 완료된 태스크는 종료 시간을 기준으로 계산
        targetTime = endTime;
      } else {
        return '-';
      }

      const diff = targetTime - startTime;
      if (diff < 0) return '-';

      const totalSeconds = Math.floor(diff / 1000);
      const hours = Math.floor(totalSeconds / 3600);
      const minutes = Math.floor((totalSeconds % 3600) / 60);
      const seconds = totalSeconds % 60;

      if (hours > 0) {
        return `${hours}h ${minutes}m ${seconds}s`;
      } else if (minutes > 0) {
        return `${minutes}m ${seconds}s`;
      } else {
        return `${seconds}s`;
      }
    },
    startBuildTaskMonitoring() {
      // 5초마다 태스크 상태 업데이트
      this.intervalManager.start('buildTaskFetch', () => {
        if (this.showBuildMonitor) {
          this.fetchBuildTasks();
        }
      }, 5000);

      // 1초마다 실행 중인 태스크의 duration 업데이트
      this.intervalManager.start('durationUpdate', () => {
        if (this.showBuildMonitor) {
          this.updateRunningTaskDurations();
        }
      }, 1000);
    },
    updateRunningTaskDurations() {
      // RUNNING 또는 PENDING 상태인 태스크들의 duration만 업데이트
      this.buildTaskList.forEach(task => {
        if (task.status === 'RUNNING' || task.status === 'PENDING') {
          task.running_time = this.calculateDuration(task);
        }
      });
    },
    stopBuildTaskMonitoring() {
      this.intervalManager.stop('buildTaskFetch');
      this.intervalManager.stop('durationUpdate');
    },
    async cancelBuildTask(taskId) {
      const confirmed = this.dialogService.confirm('Are you sure you want to cancel this build task?');
      if (!confirmed) {
        return;
      }

      const result = await this.buildService.cancelBuildTask(taskId);

      if (result.success) {
        this.showSuccessNotification(result.message, 'Task Cancelled');
        await this.fetchBuildTasks();
      } else {
        this.showErrorNotification(result.error, 'Cancellation Failed');
      }
    },
    async showBuildTaskLogs(taskId) {
      // 태스크에서 plugin_name 찾기
      const task = this.buildTaskList.find(t => t.task_id === taskId);
      if (task && task.plugin_name) {
        await this.showBuildLogs({ name: task.plugin_name, buildTaskId: taskId });
      }
    },
    async buildAllPlugins() {
      try {
        // Filter plugins that need to be built
        const pluginsToBuild = this.buildService.filterPluginsToBuild(this.plugins);

        if (pluginsToBuild.length === 0) {
          this.notificationService.warning(
            'No plugins need to be built. All plugins are already built or currently building.',
            'Build Status'
          );
          return;
        }

        // Confirm with user
        const confirmed = this.dialogService.confirm(`Build ${pluginsToBuild.length} plugin(s)?`);
        if (!confirmed) {
          return;
        }

        // Mark plugins as building
        pluginsToBuild.forEach(plugin => {
          plugin.building = true;
        });

        // Build plugins in parallel
        const results = await this.buildService.buildMultiplePlugins(pluginsToBuild, false);

        // Update plugin states based on results
        results.forEach(result => {
          const plugin = this.plugins.find(p => p.name === result.plugin);
          if (plugin) {
            plugin.building = false;
            if (result.success) {
              plugin.imageExists = result.imageExists || false;
            }
          }
        });

        // Process results and show notification
        const summary = this.buildService.processBuildResults(results);
        const message = this.buildService.getBuildSummaryMessage(summary);

        switch (message.type) {
          case 'success':
            this.notificationService.success(message.text, 'Build Complete');
            break;
          case 'warning':
            this.notificationService.warning(message.text, 'Build Results');
            break;
          case 'error':
            this.notificationService.error(message.text, 'Build Failed');
            break;
          default:
            this.notificationService.info(message.text, 'Build Status');
        }
      } catch (error) {
        console.error('Error building all plugins:', error);
        this.notificationService.error('An error occurred while building plugins.', 'Build Error');
      }
    },


    // Error handling methods
    handleLoadingError(error) {
      this.pluginError = analyzeError(error);
    },


    // Network status handlers
    handleOnline() {
      this.isOffline = false;
      console.log('Connection restored');
      
      // Optionally reload plugins when connection is restored
      if (this.pluginError) {
        this.getUserAssociatePlugins().catch(error => {
          console.error('Failed to reload plugins after connection restored:', error);
        });
      }
    },

    handleOffline() {
      this.isOffline = true;
      console.log('Connection lost');

      this.notificationService.error(
        'Internet connection lost. Some features may not work properly.',
        'Network Status'
      );
    },

    // User notification methods
    showSuccessNotification(message, title = 'Success') {
      this.notificationService.success(message, title);
    },

    showErrorNotification(message, title = 'Error') {
      this.notificationService.error(message, title);
    },

    showWarningNotification(message, title = 'Warning') {
      this.notificationService.warning(message, title);
    },
  },
};
</script>

<style scoped>
.layout {
  padding: 10px 30px;
  overflow-y: auto;
}

table {
  width: 100%;
  height: 100%;
  border-collapse: separate;
  border-spacing: 5px;
  /* background-color: #c9c9c9; */
  transition: all 0.3s ease;
  border-radius: 15px;
  /* color: #ffffff; */
}

thead th,
td {
  padding: 10px;
  padding-left: 15px;
  text-align: left;
  border-radius: 10px;
  border: 1px solid #a8a8a8;
  /* box-shadow: 0px 4px 4px rgba(176, 169, 255, 0.25); */
}

th {
  text-transform: capitalize;
  background-color: #474747;
  color: #ffffff;
}

th:hover {
  background-color: #616161;
}

td {
  display: flex;
}

button {
  margin-right: 10px;
  color: black;
  padding: 5px;
  left: 10px;
  border-radius: 10px;
  background-color: #eaecff;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  text-transform: capitalize;
}

button:disabled {
  color: #ccc;
}

.sort-icon {
  color: rgb(199, 199, 199);
  font-weight: normal;
  font-size: small;
}

.first-line {
  width: calc(100% - 10px);
  margin: 1rem 5px;
  display: flex;
  align-items: center;
}

.first-line__left,
.first-line__right {
  width: calc(50% - 5px);
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.first-line__left {
  justify-content: left;
  align-items: end;
}

.header__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  line-height: 1rem;
  /* padding-left: 2rem; */
  color: rgba(0, 0, 0, 0.8);
}

.header__desc {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1rem;
  color: rgba(0, 0, 0, 0.5);
  display: flex;
  align-items: center;
  margin-top: 1rem;
}

.add__button {
  min-width: 8rem;
  height: 2rem;
  padding: 0.2rem 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  border: none;
  /* background: #ffffff; */
  border-radius: 1.2rem;
  margin-right: 1rem;
  box-shadow: rgba(0, 0, 0, 0.15) 0px 0px 4px;
}

.add__button:hover {
  cursor: pointer;
  box-shadow: rgba(0, 0, 0, 0.35) 0px 0px 4px;
}

.add__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
  opacity: 0.8;
  margin-right: 0.5rem;
}

.monitor__button {
  min-width: 8rem;
  height: 2rem;
  padding: 0.2rem 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  border: none;
  background: #2196f3;
  color: white;
  border-radius: 1.2rem;
  margin-right: 1rem;
  box-shadow: rgba(0, 0, 0, 0.15) 0px 0px 4px;
}

.monitor__button:hover {
  cursor: pointer;
  background: #1976d2;
  box-shadow: rgba(0, 0, 0, 0.35) 0px 0px 4px;
}

.monitor__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
  opacity: 0.9;
  margin-right: 0.5rem;
  filter: brightness(0) invert(1);
}

.monitor__button h1 {
  font-size: 0.9rem;
  font-weight: 500;
  margin: 0;
}

.setting {
  width: 2.5rem;
  height: 2.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: pointer;
  border-radius: 8px;
  transition: all 0.3s ease;
  background-color: #f5f5f5;
  border: 1px solid #e0e0e0;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

.setting:hover {
  background-color: #e8e8e8;
  transform: translateY(-1px);
  box-shadow: 0 3px 6px rgba(0, 0, 0, 0.15);
}

.setting__button--icon {
  width: 1.5rem;
  height: 1.5rem;
  object-fit: contain;
  opacity: 0.8;
  transition: all 0.3s ease;
}

.setting:hover .setting__button--icon {
  opacity: 1;
  transform: scale(1.1) rotate(90deg);
}

.search {
  display: flex;
  align-items: center;
}

.search input {
  width: 300px;
  height: 2.5rem;
  border: 1px solid #e1e1e1;
  border-radius: 1rem;
  padding: 0 2rem;
  outline-style: none;
  background: #f7f7f7;
}

.search input:focus {
  border: 1px solid #bcbcbc;
}

#pageSize {
  padding: 2px;
  border-radius: 5px;
  border: 1px solid #ccc;
  margin-bottom: 5px;
}

.pagination {
  display: flex;
  justify-content: center;
  margin: 20px 0px;
}

.pagination button {
  margin: -5px 10px 0px 10px;
}

.download-icon {
  margin: 0px 0px;
  width: 33px;
  height: 33px;
}

.plugin-container {
  width: calc(100% - 10rem);
}

.title-container {
  width: 100%;
  font-size: 1.4rem;
  align-items: center;
  display: flex;
  flex-wrap: wrap;
  font-weight: 600;
  margin-top: 5px;
  gap: 0.5rem;
}

.description-container {
  width: 100%;
  font-size: 1rem;
  font-weight: 400;
  color: #474747;
  margin: 5px 2px;
  margin-bottom: 1rem;
}

.lastUpdated-container {
  width: 100%;
  font-size: 1rem;
  font-weight: 400;
  color: #474747;
  display: inline-block;
}

.option-container {
  width: 10rem;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  gap: 0.75rem;
  padding: 1rem 0.5rem;
}

.controls-row {
  display: flex;
  justify-content: space-between;
  align-items: center;
  width: 100%;
  max-width: 8rem;
}

.build-row {
  display: flex;
  justify-content: center;
  align-items: center;
  width: 100%;
  max-width: 8rem;
}

.build-button {
  padding: 0.75rem 1rem;
  margin: 0;
  border: none;
  border-radius: 8px;
  font-size: 0.85rem;
  font-weight: 600;
  cursor: pointer;
  transition: all 0.3s ease;
  background-color: #2196f3;
  color: white;
  width: 100%;
  text-align: center;
  box-shadow: 0 2px 4px rgba(33, 150, 243, 0.2);
}

.build-button:hover:not(:disabled) {
  background-color: #1976d2;
  transform: translateY(-2px);
  box-shadow: 0 4px 8px rgba(33, 150, 243, 0.3);
}

.build-button:disabled {
  cursor: not-allowed;
  opacity: 0.6;
}

.build-button.building {
  background-color: #ff9800;
  cursor: pointer;
  box-shadow: 0 2px 4px rgba(255, 152, 0, 0.2);
}

.build-button.building:hover {
  background-color: #f57c00;
  transform: translateY(-2px);
  box-shadow: 0 4px 8px rgba(255, 152, 0, 0.3);
}

.build-button.image-exists {
  background-color: #4caf50;
  cursor: not-allowed;
  box-shadow: 0 2px 4px rgba(76, 175, 80, 0.2);
}

.building-content {
  display: flex;
  align-items: center;
  gap: 8px;
  justify-content: center;
}

.loading-spinner {
  width: 16px;
  height: 16px;
  border: 2px solid rgba(255, 255, 255, 0.3);
  border-top: 2px solid white;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  0% {
    transform: rotate(0deg);
  }

  100% {
    transform: rotate(360deg);
  }
}

.switch {
  position: relative;
  display: inline-block;
  width: 3.5rem;
  height: 2rem;
  margin-top: 0;
}

.switch input {
  opacity: 0;
  width: 0;
  height: 0;
}

.slider {
  position: absolute;
  cursor: pointer;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background-color: #e0e0e0;
  -webkit-transition: 0.4s;
  transition: 0.4s;
  border: 1px solid #d0d0d0;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

.slider:before {
  position: absolute;
  content: "";
  height: 1.4rem;
  width: 1.4rem;
  left: 2px;
  top: 50%;
  transform: translateY(-50%);
  background-color: white;
  -webkit-transition: 0.4s;
  transition: 0.4s;
  border-radius: 50%;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
}

input:checked+.slider.w-color {
  background-color: #ccc;
}

input:checked+.slider.icon {
  background-color: #a37eba;
}

.slider.icon:before {
  background-color: #ffe05d;
}

.slider.icon:after {
  background-color: #e2df23;
}

input:checked+.slider {
  background-color: #2196f3;
  border-color: #1976d2;
  box-shadow: 0 2px 4px rgba(33, 150, 243, 0.2);
}

input:checked+.slider:before {
  -webkit-transform: translateX(1.4rem) translateY(-50%);
  -ms-transform: translateX(1.4rem) translateY(-50%);
  transform: translateX(1.4rem) translateY(-50%);
}

/* Rounded sliders */
.slider.round {
  border-radius: 1rem;
}

.slider.round:before {
  border-radius: 50%;
}

.disabled-toggle {
  opacity: 0.5;
  cursor: default;
}

.build-all__button {
  min-width: 8rem;
  height: 2rem;
  padding: 0.2rem 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  border: none;
  background: #2196f3;
  color: white;
  border-radius: 1.2rem;
  margin-right: 1rem;
  box-shadow: rgba(0, 0, 0, 0.15) 0px 0px 4px;
}

.build-all__button:hover {
  cursor: pointer;
  background: #1976d2;
  box-shadow: rgba(0, 0, 0, 0.35) 0px 0px 4px;
}

.build-all__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
  opacity: 0.8;
  margin-right: 0.5rem;
  filter: brightness(0) invert(1);
}

/* 빌드 로그 모달 스타일 */
.logs-modal-overlay {
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: rgba(0, 0, 0, 0.8);
  display: flex;
  justify-content: center;
  align-items: center;
  z-index: 10000;
  backdrop-filter: blur(5px);
}

.logs-modal {
  background: #2c3e50;
  border-radius: 16px;
  max-width: 90vw;
  max-height: 90vh;
  width: 900px;
  height: 700px;
  display: flex;
  flex-direction: column;
  box-shadow: 0px 4px 20px rgba(0, 0, 0, 0.5);
  border: 1px solid rgba(255, 255, 255, 0.1);
}

.logs-modal-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 1.5rem;
  border-bottom: 1px solid rgba(255, 255, 255, 0.1);
  background: #1f2a38;
  border-radius: 16px 16px 0 0;
}

.logs-modal-header h3 {
  margin: 0;
  color: #ecf0f1;
  font-weight: 600;
  font-size: 1.2rem;
}

.logs-modal-controls {
  display: flex;
  gap: 0.75rem;
}

.close-btn,
.refresh-btn {
  background: #e74c3c;
  color: white;
  border: none;
  padding: 0.75rem 1rem;
  border-radius: 8px;
  cursor: pointer;
  font-size: 0.9rem;
  line-height: 1;
  font-weight: bold;
  transition: all 0.2s ease;
  display: flex;
  align-items: center;
  justify-content: center;
  width: 32px;
  height: 32px;
}

.refresh-btn {
  background: #007bff;
}

.close-btn:hover {
  background: #c0392b;
  transform: translateY(-1px);
}

.refresh-btn:hover:not(:disabled) {
  background: #0056b3;
  transform: translateY(-1px);
}

.refresh-btn:disabled {
  background: #576574;
  cursor: not-allowed;
  opacity: 0.6;
}

.refresh-icon,
.close-icon {
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  filter: invert(100%) sepia(0%) saturate(0%) hue-rotate(0deg) brightness(100%) contrast(100%);
}

.logs-loading {
  display: flex;
  justify-content: center;
  align-items: center;
  flex: 1;
  font-size: 1.2rem;
  color: #ecf0f1;
  background: #34495e;
}

.logs-content {
  flex: 1;
  overflow: auto;
  padding: 1.5rem;
  background: #34495e;
  border-radius: 0 0 16px 16px;
}

/* 스크롤바 스타일 */
.logs-content::-webkit-scrollbar {
  width: 8px;
}

.logs-content::-webkit-scrollbar-track {
  background: #2c3e50;
  border-radius: 8px;
}

.logs-content::-webkit-scrollbar-thumb {
  background: #576574;
  border-radius: 8px;
}

.logs-content::-webkit-scrollbar-thumb:hover {
  background: #5a6c7d;
}

.logs-task-info {
  background: rgba(31, 42, 56, 0.8);
  padding: 1.25rem;
  border-radius: 12px;
  margin-bottom: 1.5rem;
  border: 1px solid rgba(255, 255, 255, 0.1);
  backdrop-filter: blur(5px);
}

.logs-task-info p {
  margin: 0.75rem 0;
  font-size: 0.95rem;
  color: #ecf0f1;
  line-height: 1.5;
}

.logs-task-info strong {
  color: #3498db;
  font-weight: 600;
}

.no-logs {
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  flex: 1;
  color: #bdc3c7;
  font-size: 1.1rem;
  text-align: center;
  padding: 2rem;
}

.logs-files {
  display: flex;
  flex-direction: column;
  gap: 1.5rem;
}

.log-file {
  background: rgba(44, 62, 80, 0.6);
  border-radius: 12px;
  border: 1px solid rgba(255, 255, 255, 0.1);
  overflow: hidden;
}

.log-file-header {
  background: rgba(31, 42, 56, 0.9);
  padding: 1rem 1.5rem;
  display: flex;
  justify-content: space-between;
  align-items: center;
  border-bottom: 1px solid rgba(255, 255, 255, 0.1);
}

.log-file-header h4 {
  margin: 0;
  color: #ecf0f1;
  font-size: 1rem;
  font-weight: 600;
}

.log-file-size {
  background: #3498db;
  color: white;
  padding: 0.25rem 0.75rem;
  border-radius: 12px;
  font-size: 0.8rem;
  font-weight: 500;
}

.log-file-content {
  margin: 0;
  padding: 1.5rem;
  background: #2c3e50;
  color: #ecf0f1;
  font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
  font-size: 0.85rem;
  line-height: 1.6;
  white-space: pre-wrap;
  word-wrap: break-word;
  overflow-x: auto;
  max-height: 400px;
  overflow-y: auto;
}

/* 로그 파일 콘텐츠 스크롤바 스타일 */
.log-file-content::-webkit-scrollbar {
  width: 6px;
  height: 6px;
}

.log-file-content::-webkit-scrollbar-track {
  background: #34495e;
  border-radius: 6px;
}

.log-file-content::-webkit-scrollbar-thumb {
  background: #576574;
  border-radius: 6px;
}

.log-file-content::-webkit-scrollbar-thumb:hover {
  background: #5a6c7d;
}

/* Plugin category tabs are now handled by PluginCategoryTabs component */

/* 플러그인 배지 스타일 */
.plugin-badge {
  padding: 0.25rem 0.75rem;
  border-radius: 12px;
  font-size: 0.75rem;
  font-weight: 600;
  letter-spacing: 0.5px;
  transition: all 0.3s ease;
  cursor: help;
  white-space: nowrap;
}

.badge-official {
  background: #e3f2fd;
  color: #1565c0;
  border: 1px solid #bbdefb;
  box-shadow: 0 1px 3px rgba(21, 101, 192, 0.1);
}

.badge-official:hover {
  background: #d0e8fc;
  box-shadow: 0 2px 6px rgba(21, 101, 192, 0.2);
  transform: translateY(-1px);
}

.badge-local {
  background: #e8f5e9;
  color: #2e7d32;
  border: 1px solid #c8e6c9;
  box-shadow: 0 1px 3px rgba(46, 125, 50, 0.1);
}

.badge-local:hover {
  background: #dcedc8;
  box-shadow: 0 2px 6px rgba(46, 125, 50, 0.2);
  transform: translateY(-1px);
}

/* Version badge styles */
.version-badge {
  padding: 0.25rem 0.5rem;
  border-radius: 10px;
  font-size: 0.7rem;
  font-weight: 500;
  letter-spacing: 0.3px;
  transition: all 0.3s ease;
  cursor: help;
  white-space: nowrap;
  font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
}

.version-official {
  background: #f5f5f5;
  color: #666;
  border: 1px solid #e0e0e0;
  box-shadow: 0 1px 2px rgba(102, 102, 102, 0.1);
}

.version-official:hover {
  background: #eeeeee;
  color: #555;
  box-shadow: 0 2px 4px rgba(102, 102, 102, 0.15);
  transform: translateY(-1px);
}


/* Duplicate plugin indicator styles */
.duplicate-indicator {
  margin-left: 0.5rem;
  padding: 0.25rem;
  border-radius: 50%;
  background: #fff3cd;
  border: 1px solid #ffeaa7;
  cursor: help;
  transition: all 0.3s ease;
  display: inline-flex;
  align-items: center;
  justify-content: center;
}

.duplicate-indicator:hover {
  background: #ffeb3b;
  border-color: #fdd835;
  transform: scale(1.1);
}

.duplicate-icon {
  width: 12px;
  height: 12px;
  opacity: 0.7;
  filter: sepia(100%) saturate(200%) hue-rotate(30deg);
}

.duplicate-indicator:hover .duplicate-icon {
  opacity: 1;
}

/* GPU indicator styles */
.gpu-indicator {
  margin-left: 0.5rem;
  padding: 0.25rem 0.75rem;
  border-radius: 12px;
  background: linear-gradient(135deg, #ff9800 0%, #ff6f00 100%);
  color: white;
  font-size: 0.75rem;
  font-weight: 700;
  letter-spacing: 0.5px;
  cursor: help;
  transition: all 0.3s ease;
  box-shadow: 0 2px 4px rgba(255, 152, 0, 0.3);
  display: inline-flex;
  align-items: center;
  justify-content: center;
  position: relative;
  overflow: hidden;
}

.gpu-indicator::before {
  content: '';
  position: absolute;
  top: -50%;
  left: -50%;
  width: 200%;
  height: 200%;
  background: linear-gradient(45deg, 
    transparent 30%, 
    rgba(255, 255, 255, 0.3) 50%, 
    transparent 70%);
  transform: rotate(45deg);
  transition: all 0.6s;
  opacity: 0;
}

.gpu-indicator:hover {
  transform: translateY(-2px);
  box-shadow: 0 4px 8px rgba(255, 152, 0, 0.4);
}

.gpu-indicator:hover::before {
  animation: shimmer 0.6s ease;
}

@keyframes shimmer {
  0% {
    transform: translateX(-100%) translateY(-100%) rotate(45deg);
    opacity: 0;
  }
  50% {
    opacity: 1;
  }
  100% {
    transform: translateX(100%) translateY(100%) rotate(45deg);
    opacity: 0;
  }
}

.gpu-text {
  position: relative;
  z-index: 1;
}

/* 읽기 전용 설정 버튼 스타일 */
.setting.readonly {
  background-color: #e8f4fd;
  border-color: #bbdefb;
}

.setting.readonly:hover {
  background-color: #d6ebfa;
  transform: translateY(-1px);
  box-shadow: 0 3px 6px rgba(21, 101, 192, 0.15);
}

.setting.readonly:hover .setting__button--icon {
  transform: scale(1.1);
  transition: 0.3s;
}

/* Official 플러그인 상태 스타일 */
.official-status {
  display: flex;
  align-items: center;
  justify-content: center;
  margin: 0;
}

.pre-built {
  padding: 0.75rem 1rem;
  background: #f0f7ff;
  color: #1565c0;
  border-radius: 8px;
  font-size: 0.85rem;
  font-weight: 600;
  border: 1px solid #bbdefb;
  width: 100%;
  text-align: center;
  box-shadow: 0 2px 4px rgba(21, 101, 192, 0.1);
}


/* Accessibility helpers */
.visually-hidden {
  position: absolute !important;
  width: 1px !important;
  height: 1px !important;
  padding: 0 !important;
  margin: -1px !important;
  overflow: hidden !important;
  clip: rect(0, 0, 0, 0) !important;
  white-space: nowrap !important;
  border: 0 !important;
}
</style>
