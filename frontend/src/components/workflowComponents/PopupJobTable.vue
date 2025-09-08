<template>
  <div class="control-popup__jobs" v-if="show_jobs">
    <button class="close-button" @click="closePopup" aria-label="Close">
    </button>
    <div class="job-table-container">
      <table class="job-table">
        <thead>
          <tr>
            <th>No.</th>
            <th>Name</th>
            <th>Plugin</th>
            <th>Type</th>
            <th>Start</th>
            <th>End</th>
            <th>Time</th>
            <th>Status</th>
            <!-- <th>Action</th> -->
          </tr>
        </thead>
        <tbody>
          <tr v-for="(task, index) in taskList" :key="index" @click.right.prevent="RMouseClick($event, task)">
            <td>{{ index + 1 }}</td>
            <td>{{ task.workflow_title | titleNone }}</td>
            <td>{{ formatPluginInfo(task) }}</td>
            <td>{{ formatPluginType(task) }}</td>
            <td>{{ task.start_time | formatDateTime }}</td>
            <td>{{ task.end_time | formatDateTime }}</td>
            <td>{{ task.running_time }}</td>
            <td>
              <div class="task-status">
                <div class="status-indicator" :class="getStatusClass(task.status)"></div>
                {{ task.status }}
              </div>
            </td>
            <!-- <td>
              <img v-if="task.status === 'RUNNING' || task.status === 'PENDING' || task.status === 'INSTALLING'"
                @click="cancelTask(task.task_id)" class="cancel-icon" src="@/assets/multiply.png" alt="Cancel" />
            </td> -->
          </tr>
        </tbody>
      </table>
    </div>
    <ul class="toggle__menu" v-bind:class="{ open: R_Mouse_isActive }" :style="{ left: xPosition, top: yPosition }"
      @click.stop>
      <li @click="confirmDelete" v-if="isCompleted">Delete</li>
      <li @click="cancelTask" v-else>Cancle</li>
      <li @click="showLogs">View Logs</li>
      <li @click="viewProgress">View Progress</li>
      <li @click="downloadExecutionManifest" v-if="isCompleted && canDownloadManifest">Download Execution Manifest</li>
    </ul>
  </div>
</template>

<script>
import moment from "moment";

export default {
  props: {
    show_jobs: {
      type: Boolean,
      required: true
    },
    taskList: {
      type: Array,
      required: true
    }
  },
  data() {
    return {
      R_Mouse_isActive: false,
      xPosition: 0,
      yPosition: 0,
      isCompleted: false,
      currentTaskId: null,
      currentTask: null,
      canDownloadManifest: false
    };
  },
  created() {
    document.addEventListener('click', this.hideMenu);
  },
  beforeDestroy() {
    document.removeEventListener('click', this.hideMenu);
  },
  methods: {
    hideMenu() {
      this.R_Mouse_isActive = false;
    },
    cancelTask() {
      this.$emit('cancel-task', this.currentTaskId);
      this.R_Mouse_isActive = false;
    },
    confirmDelete() {
      this.$emit('confirm-delete', this.currentTaskId);
      this.R_Mouse_isActive = false;
    },
    showLogs() {
      this.$emit('show-logs', this.currentTaskId);
      this.R_Mouse_isActive = false;
    },
    viewProgress() {
      this.$emit('view-progress', this.currentTaskId);
      this.R_Mouse_isActive = false;
    },
    downloadExecutionManifest() {
      this.$emit('download-execution-manifest', this.currentTaskId);
      this.R_Mouse_isActive = false;
    },
    closePopup() {
      this.$emit('close-popup');
    },
    getStatusClass(status) {
      if (status === "SUCCESS") return "status-success";
      if (status === "FAILURE" || status === "REVOKED" || status === "RETRY") return "status-failure";
      if (status === "RUNNING" || status === "PENDING" || status === "INSTALLING") return "status-running";
    },
    RMouseClick(event, task) {
      // 기존 메뉴를 먼저 숨김
      this.R_Mouse_isActive = false;

      // Vue의 다음 틱에서 실행하여 DOM 업데이트 후 위치 계산
      this.$nextTick(() => {
        // 부모 컨테이너(.control-popup__jobs)의 위치 정보 가져오기
        const popupContainer = this.$el;
        const containerRect = popupContainer.getBoundingClientRect();

        // 마우스 클릭 위치를 부모 컨테이너 기준으로 변환
        const relativeX = event.clientX - containerRect.left;
        const relativeY = event.clientY - containerRect.top;

        // 메뉴 크기 계산을 위해 임시로 메뉴를 표시
        this.currentTaskId = task.task_id;
        this.currentTask = task;
        this.isCompleted = ["SUCCESS", "FAILURE", "REVOKED", "RETRY"].includes(task.status);
        this.canDownloadManifest = task.status === 'SUCCESS' && this.formatPluginType(task) === 'Analysis';
        this.R_Mouse_isActive = true;

        // 다음 틱에서 실제 메뉴 크기 측정
        this.$nextTick(() => {
          const menuElement = this.$el.querySelector('.toggle__menu');
          let menuWidth = 200; // 기본값
          let menuHeight = 80; // 기본값

          // 메뉴가 렌더링된 경우 실제 크기 측정
          if (menuElement) {
            const menuRect = menuElement.getBoundingClientRect();
            menuWidth = menuRect.width;
            menuHeight = menuRect.height;
          }

          let x = relativeX;
          let y = relativeY;

          // 컨테이너 경계를 고려한 위치 조정
          const containerWidth = containerRect.width;
          const containerHeight = containerRect.height;
          const margin = 10; // 안전 여백

          // 오른쪽 경계 체크 - 메뉴가 컨테이너 밖으로 나가면 왼쪽으로 이동
          if (x + menuWidth > containerWidth - margin) {
            x = Math.max(margin, containerWidth - menuWidth - margin);
          }

          // 하단 경계 체크 - 메뉴가 컨테이너 밖으로 나가면 위쪽으로 이동
          if (y + menuHeight > containerHeight - margin) {
            y = Math.max(margin, containerHeight - menuHeight - margin);
          }

          // 왼쪽 경계 체크 - 최소 여백 보장
          if (x < margin) {
            x = margin;
          }

          // 상단 경계 체크 - 최소 여백 보장
          if (y < margin) {
            y = margin;
          }

          // 최종 안전장치 - 메뉴가 너무 클 경우 컨테이너 내부에 강제 배치
          if (menuWidth > containerWidth - (margin * 2)) {
            x = margin;
          }

          if (menuHeight > containerHeight - (margin * 2)) {
            y = margin;
          }

          // 정수값으로 변환하여 픽셀 완전성 보장
          this.xPosition = Math.round(x) + 'px';
          this.yPosition = Math.round(y) + 'px';
        });
      });
    },
    formatPluginInfo(task) {
      if (!task.plugin_name) {
        return "N/A";
      }

      // Check if enhanced plugin data is available
      if (task.plugin && task.plugin.source && task.plugin.version) {
        const { source, version } = task.plugin;
        return `${source}/${task.plugin_name} : v${version}`;
      }

      // Legacy fallback for tasks without plugin relationship
      return task.plugin_name;
    },
    formatPluginType(task) {
      // Check if enhanced plugin data is available
      if (task.plugin && task.plugin.plugin_type) {
        const pluginType = task.plugin.plugin_type;
        if (pluginType === 'compile' || pluginType === 'analysis') {
          return 'Analysis';
        } else if (pluginType === 'visualization') {
          return 'Visualization';
        }
        return pluginType;
      }

      // Legacy fallback for tasks without plugin relationship
      let taskType = task.task_type;
      if (taskType === 'compile' || pluginType === 'analysis') {
        return 'Analysis';
      } else if (taskType === 'visualization') {
        return 'Visualization';
      }
      return taskType || 'unknown';
    },
  },
  filters: {
    formatDateTime(dateTime) {
      const date = moment(dateTime).format("YYYY.MM.DD-HH:mm");
      if (date === "Invalid date") return "Not Yet Completed"; // 날짜가 유효하지 않을 경우 처리
      return date;
    },
    titleNone(title) {
      if (title === null) return "Untitled";
      return title;
    }
  }
};
</script>

<style scoped>
/* ControlBar 위에 4px 공백을 두고 위치 */
.control-popup__jobs {
  height: 540px;
  max-height: 540px;
  /* width: 720px; */
  left: 50%;
  bottom: 92px;
  /* ControlBar(24px + 64px) + 4px 공백 */
  transform: translateX(-50%);
  overflow-y: auto;
  border-radius: 16px;
  position: absolute;
  z-index: 9998;
  opacity: 1;
  display: flex;
  /* align-items: center; */
  justify-content: center;
  background-color: #2c3e50;
  /* 다크 테마 */
  box-shadow: 0px 4px 10px rgba(0, 0, 0, 0.3);
  padding: 0.5rem;
}

/* 더 얇고 어두운 스크롤바 스타일 */
.job-table-container::-webkit-scrollbar {
  width: 6px;
  /* 기존 10px → 6px로 얇게 */
}

.job-table-container::-webkit-scrollbar-track {
  background: #2a3d55;
  /* 더 어두운 트랙 색상 */
  border-radius: 12px;
}

.job-table-container::-webkit-scrollbar-thumb {
  background: #444;
  /* 기본 스크롤바 색상 더 어둡게 */
  border-radius: 12px;
}

.job-table-container::-webkit-scrollbar-thumb:hover {
  background: #333;
  /* Hover 시 더 어두운 색 */
}

/* 테이블 스타일 */
.job-table-container {
  /* width: 95%; */
  height: auto;
  max-height: calc(540px - 1rem);
  overflow-y: auto;
  border-radius: 10px;
}

.job-table {
  width: 100%;
  border-collapse: collapse;
  background-color: #34495e;
  border-radius: 10px;
  overflow: hidden;
}

.job-table thead {
  background-color: #1f2a38;
  color: #ecf0f1;
}

.job-table th,
.job-table td {
  padding: 8px;
  text-align: center;
  border-bottom: 1px solid #576574;
  color: #ecf0f1;
}

.job-table td {
  font-size: 0.75rem;
  text-align: center !important;
  vertical-align: middle;
}

/* 내부 요소들도 중앙 정렬 */
.job-table td * {
  text-align: center;
}

.job-table tbody tr {
  transition: background-color 0.2s ease-in-out;
}

.job-table tbody tr:hover {
  background-color: #3d566e;
}


.toggle__menu {
  display: none;
  position: absolute;
  width: 200px;
  margin: 0;
  padding: 0;
  background: #ffffff;
  border-radius: 8px;
  list-style: none;
  box-shadow: 0 8px 25px rgba(0, 0, 0, 0.15), 0 3px 10px rgba(0, 0, 0, 0.1);
  overflow: hidden;
  z-index: 999999;
  text-transform: capitalize;
  pointer-events: none;
  /* 부모 컨테이너와 동일한 좌표계 사용 */
  transform: none;
}

.toggle__menu.open {
  display: block;
  opacity: 1;
  pointer-events: auto;
}

.toggle__menu>li {
  border-left: 3px solid transparent;
  transition: ease 0.2s;
  padding: 10px;
}

.toggle__menu>li:hover {
  background: #e5e5e5;
}

/* 상태 아이콘 */
.task-status {
  display: flex;
  align-items: center;
  justify-content: center;
}

.status-indicator {
  width: 8px;
  height: 8px;
  border-radius: 50%;
  margin-right: 8px;
}

.status-failure {
  background-color: #e74c3c;
}

.status-success {
  background-color: #2ecc71;
}

.status-running {
  background-color: #f39c12;
}

/* 취소 아이콘 */
.cancel-icon {
  width: 16px;
  height: 16px;
  cursor: pointer;
}

/* 닫기 버튼 스타일 */
.close-button {
  position: absolute;
  top: 10px;
  right: 10px;
  width: 24px;
  height: 24px;
  border-radius: 50%;
  background-color: rgba(231, 76, 60, 0.3);
  /* 30% 투명도 */
  border: none;
  cursor: pointer;
  display: flex;
  align-items: center;
  justify-content: center;
  transition: all 0.3s ease;
  z-index: 10;
  box-shadow: none;
}

.close-button:hover {
  background-color: rgba(231, 76, 60, 1);
  /* 호버 시 불투명 */
  transform: scale(1.1);
  box-shadow: 0 2px 6px rgba(0, 0, 0, 0.3);
}

.close-button:active {
  transform: scale(0.95);
}
</style>