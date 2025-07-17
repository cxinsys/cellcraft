<template>
  <div class="layout">
    <div class="first-line">
      <div class="first-line__left">
        <div class="header__text">
          Plugins
          <div class="header__desc">
            Plugin for data analysis algorithm extensions
          </div>
        </div>
      </div>
      <div class="first-line__right">
        <div class="add__button" @click="addPluginExtension">
          <img class="add__button--icon" src="@/assets/add_circle.png" />
          <h1>Add Plugin</h1>
        </div>
        <div class="monitor__button" @click="toggleBuildMonitor">
          <img class="monitor__button--icon" src="@/assets/control_jobs.png" />
          <h1>Build Monitor</h1>
        </div>
        <!-- <div class="build-all__button" @click="buildAllPlugins">
          <img class="build-all__button--icon" src="@/assets/add_circle.png" />
          <h1>Build All</h1>
        </div> -->
        <div class="search">
          <input type="text" v-model="searchTerm" placeholder="Search titles..." />
        </div>
      </div>
    </div>
    <PluginExtention v-if="showPluginExtension" @close="closePluginExtension" :editName="selectedPlugin.name"
      :editDescription="selectedPlugin.description" :editDependencies="selectedPlugin.dependencies"
      :editDrawflow="selectedPlugin.drawflow" :editRules="selectedPlugin.rules" />
    <BuildMonitor v-if="showBuildMonitor" :show_monitor="showBuildMonitor" :buildTaskList="buildTaskList"
      @cancel-task="cancelBuildTask" @show-logs="showBuildTaskLogs" @close="toggleBuildMonitor" />
    <table>
      <tbody>
        <tr v-for="plugin in filteredPlugins" :key="plugin.id">
          <td>
            <div class="plugin-container">
              <div class="title-container">
                {{ plugin.name }}
              </div>
              <div class="description-container">
                {{ plugin.description }}
              </div>
              <div class="lastUpdated-container">
                Last Updated: {{ plugin.updated_at.split('T')[0] }}
              </div>
            </div>
            <div class="option-container">
              <div class="setting" @click="editPluginExtension(plugin)">
                <img class="setting__button--icon" src="@/assets/settings.png" />
              </div>
              <button class="build-button" @click="handleBuildPlugin(plugin)" :disabled="plugin.imageExists"
                :class="{ 'building': plugin.building, 'image-exists': plugin.imageExists }">
                <div v-if="plugin.building" class="building-content">
                  <div class="loading-spinner"></div>
                  <span>View Logs</span>
                </div>
                <span v-else-if="plugin.imageExists">Built</span>
                <span v-else>Build</span>
              </button>
              <label class="switch">
                <input v-model="plugin.checked" type="checkbox" @change="handlePluginAssociate(plugin)"
                  :disabled="isCheckboxDisabled" />
                <span class="slider round"></span>
              </label>
            </div>
          </td>
        </tr>
      </tbody>
    </table>

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
import { getUser, getPlugins, associatePlugin, dissociatePlugin, buildPluginDocker, checkPluginImage, getBuildStatus, getBuildTasks, cancelBuildTask, getBuildLogs } from "@/api/index";
import PluginExtention from "@/components/PluginExtention.vue";
import BuildMonitor from "@/components/pluginComponents/BuildMonitor.vue";

export default {
  components: {
    PluginExtention,
    BuildMonitor,
  },
  data() {
    return {
      showPluginExtension: false,
      showBuildMonitor: false,
      buildTaskList: [],
      isCheckboxDisabled: false,
      searchTerm: "",
      plugins: [
        // {
        //   id: 1,
        //   name: "TENET",
        //   description:
        //     "A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data",
        //   lastUpdated: "2024/03/18",
        //   checked: true,
        // },
        // {
        //   id: 2,
        //   name: "TENET TF",
        //   description:
        //     "A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data",
        //   lastUpdated: "2024/03/18",
        //   checked: true,
        // },
      ],
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
    };
  },
  async mounted() {
    try {
      await this.getUserAssociatePlugins();
    } catch (error) {
      console.error(error);
    }
  },
  beforeDestroy() {
    // 컴포넌트가 파괴될 때 모니터링 중지
    this.stopBuildTaskMonitoring();
  },
  computed: {
    filteredPlugins() {
      return this.plugins.filter((plugin) =>
        plugin.name.toLowerCase().includes(this.searchTerm.toLowerCase())
      );
    },
  },
  methods: {
    async closePluginExtension(buildInfo) {
      this.showPluginExtension = false;

      // extension 완료했으니, 다시 plugin list를 불러옵니다.
      try {
        await this.getUserAssociatePlugins();

        // 빌드가 시작된 새 플러그인이 있으면 상태 업데이트 (API 데이터 로드 후에 덮어쓰기)
        if (buildInfo && buildInfo.buildStarted) {
          // 잠시 기다린 후 상태 업데이트 (API 응답 처리 완료 후)
          this.$nextTick(() => {
            const plugin = this.plugins.find(p => p.name === buildInfo.pluginName);
            if (plugin) {
              // 강제로 빌드 상태 설정 (API에서 받은 데이터 덮어쓰기)
              plugin.building = true;
              plugin.buildTaskId = buildInfo.taskId;
              plugin.buildStatus = 'RUNNING';
              plugin.imageExists = false;

              console.log(`Plugin ${buildInfo.pluginName} build status forced to building`);

              // 빌드 모니터링 시작
              this.startBuildMonitoring(plugin);
            }
          });
        }
      } catch (error) {
        console.error(error);
      }
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
    async getUserAssociatePlugins() {
      try {
        const profile = await getUser();
        this.profile = profile.data;

        const plugins = await getPlugins();
        console.log(plugins.data.plugins);
        const currentUser = this.profile.username;

        this.plugins = plugins.data.plugins.map(plugin => {
          const userIncluded = plugin.users.some(user => user.username === currentUser);
          const buildInfo = plugin.latest_build || {};
          const isBuilding = buildInfo.status === 'RUNNING' || buildInfo.status === 'PENDING';

          return {
            ...plugin,
            checked: userIncluded,
            building: isBuilding,
            imageExists: false,
            buildTaskId: buildInfo.task_id || null,
            buildStatus: buildInfo.status || null,
          };
        });

        // 각 플러그인의 이미지 존재 여부 확인
        await this.checkAllPluginImages();

        // 빌드 중인 플러그인들의 모니터링 시작
        this.plugins.forEach(plugin => {
          if (plugin.building && plugin.buildTaskId) {
            console.log(`Starting monitoring for plugin ${plugin.name} with task ${plugin.buildTaskId}`);
            this.startBuildMonitoring(plugin);
          }
        });
      } catch (error) {
        console.error(error);
      }
    },
    async handlePluginAssociate(plugin) {
      const pluginId = parseInt(plugin.id);

      try {
        let result;
        if (plugin.checked) {
          result = await associatePlugin(pluginId);
        } else {
          result = await dissociatePlugin(pluginId);
        }
        console.log(result.data);

        // 1초 동안 체크박스 비활성화
        this.isCheckboxDisabled = true;
        setTimeout(() => {
          this.isCheckboxDisabled = false;
        }, 1000);
      } catch (error) {
        console.error('Error associating/disassociating plugin:', error);
      }
    },
    async checkAllPluginImages() {
      for (let plugin of this.plugins) {
        try {
          const result = await checkPluginImage(plugin.name);
          plugin.imageExists = result.data.image_exists;
        } catch (error) {
          console.error(`Error checking image for plugin ${plugin.name}:`, error);
          plugin.imageExists = false;
        }
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
        const result = await buildPluginDocker(plugin.name, false); // 기존 플러그인은 기본적으로 GPU 비활성화

        // 태스크 ID 저장
        plugin.buildTaskId = result.data.task_id;
        plugin.buildStatus = 'RUNNING';

        // 빌드 상태 주기적 모니터링 시작
        this.startBuildMonitoring(plugin);

        console.log('Build started:', result.data);
        alert(`Plugin ${plugin.name} build started!`);
      } catch (error) {
        console.error(`Error starting build for plugin ${plugin.name}:`, error);
        plugin.building = false;

        // 에러 메시지 표시
        let errorMessage = `Failed to start build for plugin ${plugin.name}`;
        if (error.response && error.response.data && error.response.data.detail) {
          if (typeof error.response.data.detail === 'string') {
            errorMessage += `: ${error.response.data.detail}`;
          } else if (error.response.data.detail.message) {
            errorMessage += `: ${error.response.data.detail.message}`;
          }
        }
        alert(errorMessage);
      }
    },
    async startBuildMonitoring(plugin) {
      const checkInterval = setInterval(async () => {
        try {
          const result = await getBuildStatus(plugin.buildTaskId);
          const status = result.data.state; // 'state' not 'status'

          console.log(`Build status for ${plugin.name}: ${status}`);

          // 상태에 따라 플러그인 속성 업데이트
          plugin.buildStatus = status;

          if (status === 'SUCCESS') {
            plugin.building = false;
            plugin.imageExists = true;
            plugin.buildStatus = 'SUCCESS';
            clearInterval(checkInterval);
            console.log(`Plugin ${plugin.name} built successfully!`);

            // 성공 알림 (선택사항)
            this.showBuildNotification(plugin.name, 'success');
          } else if (status === 'FAILURE' || status === 'REVOKED') {
            plugin.building = false;
            plugin.buildStatus = status;
            plugin.imageExists = false; // 실패 시 이미지 존재하지 않음
            clearInterval(checkInterval);
            console.error(`Plugin ${plugin.name} build failed with status: ${status}`);

            // 실패 알림 (선택사항)
            this.showBuildNotification(plugin.name, 'failure');
          } else if (status === 'RUNNING' || status === 'PENDING') {
            // 빌드 중 상태 유지
            plugin.building = true;
            plugin.buildStatus = status;
          }
        } catch (error) {
          console.error(`Error checking build status for plugin ${plugin.name}:`, error);
          // 에러 발생 시 모니터링 중단
          clearInterval(checkInterval);
          plugin.building = false;
          plugin.buildStatus = 'ERROR';
        }
      }, 2000); // 2초마다 상태 확인
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

        const result = await getBuildLogs(plugin.name);
        this.selectedBuildLogs = {
          status: plugin.buildStatus || 'Unknown',
          timestamp: new Date().toLocaleString(),
          logs: result.data.log_content || 'No logs available'
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
        const response = await getBuildTasks();
        const tasks = response.data.tasks || [];

        // 태스크 데이터를 PopupJobTable 형식에 맞게 변환
        this.buildTaskList = tasks.map(task => {
          return {
            task_id: task.task_id,
            plugin_name: task.plugin_name,
            start_time: task.start_time,
            end_time: task.end_time,
            running_time: this.calculateDuration(task),
            status: task.state,
            error: task.error,
            info: task.info
          };
        });

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
      this.buildTaskInterval = setInterval(() => {
        if (this.showBuildMonitor) {
          this.fetchBuildTasks();
        }
      }, 5000);

      // 1초마다 실행 중인 태스크의 duration 업데이트
      this.durationUpdateInterval = setInterval(() => {
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
      if (this.buildTaskInterval) {
        clearInterval(this.buildTaskInterval);
        this.buildTaskInterval = null;
      }
      if (this.durationUpdateInterval) {
        clearInterval(this.durationUpdateInterval);
        this.durationUpdateInterval = null;
      }
    },
    async cancelBuildTask(taskId) {
      if (!confirm('정말로 이 빌드 작업을 취소하시겠습니까?')) {
        return;
      }

      try {
        await cancelBuildTask(taskId);
        alert('빌드 작업이 취소되었습니다.');
        await this.fetchBuildTasks();
      } catch (error) {
        console.error('Failed to cancel build task:', error);
        alert('빌드 작업 취소에 실패했습니다.');
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
        // 빌드가 필요한 플러그인들만 필터링 (building이 아니고 imageExists가 false인 것들)
        const pluginsToBuild = this.plugins.filter(plugin => !plugin.building && !plugin.imageExists);

        if (pluginsToBuild.length === 0) {
          alert('빌드할 플러그인이 없습니다. 모든 플러그인이 이미 빌드되었거나 빌드 중입니다.');
          return;
        }

        if (!confirm(`${pluginsToBuild.length}개의 플러그인을 빌드하시겠습니까?`)) {
          return;
        }

        const buildPromises = pluginsToBuild.map(plugin => {
          plugin.building = true;
          return buildPluginDocker(plugin.name, false) // 기존 플러그인은 기본적으로 GPU 비활성화
            .then(result => {
              console.log(`Build result for ${plugin.name}:`, result.data);
              setTimeout(async () => {
                try {
                  const checkResult = await checkPluginImage(plugin.name);
                  plugin.imageExists = checkResult.data.image_exists;
                } catch (error) {
                  console.error(`Error checking image after build for plugin ${plugin.name}:`, error);
                }
                plugin.building = false;
              }, 1000);
              return { success: true, plugin: plugin.name };
            })
            .catch(error => {
              console.error(`Error building plugin ${plugin.name}:`, error);
              plugin.building = false;
              return { success: false, plugin: plugin.name, error };
            });
        });

        const results = await Promise.all(buildPromises);
        const successful = results.filter(r => r.success).length;
        const failed = results.filter(r => !r.success).length;

        if (failed === 0) {
          alert(`모든 플러그인(${successful}개)이 성공적으로 빌드되었습니다!`);
        } else {
          alert(`${successful}개 플러그인은 성공, ${failed}개 플러그인은 실패했습니다.`);
        }
      } catch (error) {
        console.error('Error building all plugins:', error);
        alert('플러그인 빌드 중 오류가 발생했습니다.');
      }
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
  width: 2rem;
  height: 2rem;
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: pointer;
}

.setting__button--icon {
  width: 2rem;
  height: 2rem;
  object-fit: contain;
  opacity: 1;
}

.setting__button--icon:hover {
  /* 톱니바퀴 이미지 마우스 올렸을 때, 1.1배 커지고 rotate 애니메이션 */
  transition: 0.5s;
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
  width: calc(100% - 8rem);
}

.title-container {
  width: 100%;
  font-size: 1.4rem;
  align-items: center;
  display: flex;
  font-weight: 600;
  margin-top: 5px;
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
  width: 8rem;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
}

.build-button {
  padding: 8px 16px;
  margin: 8px 0;
  border: none;
  border-radius: 6px;
  font-size: 0.9rem;
  font-weight: 500;
  cursor: pointer;
  transition: all 0.3s ease;
  background-color: #2196f3;
  color: white;
  min-width: 80px;
}

.build-button:hover:not(:disabled) {
  background-color: #1976d2;
  transform: translateY(-1px);
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
}

.build-button:disabled {
  cursor: not-allowed;
  opacity: 0.6;
}

.build-button.building {
  background-color: #ff9800;
  cursor: pointer;
}

.build-button.building:hover {
  background-color: #f57c00;
}

.build-button.image-exists {
  background-color: #4caf50;
  cursor: not-allowed;
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
  width: 60px;
  height: 34px;
  margin-top: 1rem;
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
  background-color: #ccc;
  -webkit-transition: 0.4s;
  transition: 0.4s;
}

.slider:before {
  position: absolute;
  content: "";
  height: 26px;
  width: 26px;
  left: 4px;
  bottom: 4px;
  background-color: white;
  -webkit-transition: 0.4s;
  transition: 0.4s;
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
}

input:checked+.slider:before {
  -webkit-transform: translateX(26px);
  -ms-transform: translateX(26px);
  transform: translateX(26px);
}

/* Rounded sliders */
.slider.round {
  border-radius: 34px;
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
</style>
