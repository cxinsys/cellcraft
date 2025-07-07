<template>
  <div class="build-monitor-container" v-if="show_monitor">
    <div class="build-monitor">
      <div class="build-monitor-header">
        <h3>Plugin Build Monitor</h3>
        <button @click="$emit('close')" class="close-btn">
          <img src="@/assets/close.png" alt="Close" class="close-icon" />
        </button>
      </div>
      
      <div class="build-table-container">
        <table class="build-table">
          <thead>
            <tr>
              <th>Plugin Name</th>
              <th>Start Time</th>
              <th>End Time</th>
              <th>Duration</th>
              <th>Status</th>
            </tr>
          </thead>
          <tbody>
            <tr v-for="task in buildTaskList" :key="task.task_id" @click.right.prevent="RMouseClick($event, task)">
              <td>{{ task.plugin_name || 'Unknown' }}</td>
              <td>{{ task.start_time | formatDateTime }}</td>
              <td>{{ task.end_time | formatDateTime }}</td>
              <td>{{ task.running_time || '-' }}</td>
              <td>
                <div class="task-status">
                  <div class="status-indicator" :class="getStatusClass(task.status)"></div>
                  {{ task.status }}
                </div>
              </td>
            </tr>
          </tbody>
        </table>
        
        <div v-if="buildTaskList.length === 0" class="no-tasks">
          No build tasks found
        </div>
      </div>
    </div>
    
    <!-- Right-click context menu -->
    <ul class="toggle__menu" v-bind:class="{ open: R_Mouse_isActive }" :style="{ left: xPosition, top: yPosition }"
      @click.stop>
      <li @click="showLogs">View Logs</li>
      <li @click="cancelTask" v-if="!isCompleted">Cancel</li>
    </ul>
  </div>
</template>

<script>
import moment from "moment";

export default {
  props: {
    show_monitor: {
      type: Boolean,
      required: true
    },
    buildTaskList: {
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
      currentTaskId: null
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
    cancelTask(taskId) {
      if (taskId) {
        this.$emit('cancel-task', taskId);
      } else if (this.currentTaskId) {
        this.$emit('cancel-task', this.currentTaskId);
      }
      this.R_Mouse_isActive = false;
    },
    showLogs(taskId) {
      if (taskId) {
        this.$emit('show-logs', taskId);
      } else if (this.currentTaskId) {
        this.$emit('show-logs', this.currentTaskId);
      }
      this.R_Mouse_isActive = false;
    },
    getStatusClass(status) {
      if (status === "SUCCESS") return "status-success";
      if (status === "FAILURE" || status === "REVOKED" || status === "RETRY") return "status-failure";
      if (status === "RUNNING" || status === "PENDING" || status === "INSTALLING") return "status-running";
      return "";
    },
    RMouseClick(event, task) {
      this.R_Mouse_isActive = false;
      this.xPosition = Math.min(event.clientX, window.innerWidth - 210) + 'px';
      this.yPosition = Math.min(event.clientY, window.innerHeight - 60) + 'px';
      this.R_Mouse_isActive = true;
      this.currentTaskId = task.task_id;
      this.isCompleted = ["SUCCESS", "FAILURE", "REVOKED", "RETRY"].includes(task.status);
    },
  },
  filters: {
    formatDateTime(dateTime) {
      if (!dateTime) return "-";
      const date = moment(dateTime).format("YYYY.MM.DD-HH:mm");
      if (date === "Invalid date") return "-";
      return date;
    }
  }
};
</script>

<style scoped>
.build-monitor-container {
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: rgba(0, 0, 0, 0.5);
  display: flex;
  justify-content: center;
  align-items: center;
  z-index: 10000;
  backdrop-filter: blur(5px);
}

.build-monitor {
  background: #2c3e50;
  border-radius: 16px;
  width: 90%;
  max-width: 1200px;
  max-height: 80vh;
  display: flex;
  flex-direction: column;
  box-shadow: 0px 4px 20px rgba(0, 0, 0, 0.5);
  border: 1px solid rgba(255, 255, 255, 0.1);
}

.build-monitor-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 1.5rem;
  border-bottom: 1px solid rgba(255, 255, 255, 0.1);
  background: #1f2a38;
  border-radius: 16px 16px 0 0;
}

.build-monitor-header h3 {
  margin: 0;
  color: #ecf0f1;
  font-weight: 600;
  font-size: 1.2rem;
}

.close-btn {
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

.close-btn:hover {
  background: #c0392b;
  transform: translateY(-1px);
}

.close-icon {
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  filter: invert(100%) sepia(0%) saturate(0%) hue-rotate(0deg) brightness(100%) contrast(100%);
}

/* Table container */
.build-table-container {
  flex: 1;
  overflow-y: auto;
  padding: 1rem;
}

/* Scrollbar styling */
.build-table-container::-webkit-scrollbar {
  width: 6px;
}

.build-table-container::-webkit-scrollbar-track {
  background: #2a3d55;
  border-radius: 12px;
}

.build-table-container::-webkit-scrollbar-thumb {
  background: #444;
  border-radius: 12px;
}

.build-table-container::-webkit-scrollbar-thumb:hover {
  background: #333;
}

/* Table styling */
.build-table {
  width: 100%;
  border-collapse: collapse;
  background-color: #34495e;
  border-radius: 10px;
  overflow: hidden;
}

.build-table thead {
  background-color: #1f2a38;
  color: #ecf0f1;
}

.build-table th,
.build-table td {
  padding: 12px;
  text-align: center;
  border-bottom: 1px solid #576574;
  color: #ecf0f1;
}

.build-table td {
  font-size: 0.9rem;
}

.build-table tbody tr {
  transition: background-color 0.2s ease-in-out;
}

.build-table tbody tr:hover {
  background-color: #3d566e;
}

/* Status styling */
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
  animation: pulse 1.5s infinite;
}

@keyframes pulse {
  0% {
    opacity: 1;
  }
  50% {
    opacity: 0.5;
  }
  100% {
    opacity: 1;
  }
}


/* No tasks message */
.no-tasks {
  text-align: center;
  color: #bdc3c7;
  font-size: 1.1rem;
  padding: 3rem;
}

/* Context menu */
.toggle__menu {
  display: none;
  position: fixed;
  width: 200px;
  margin: 0;
  padding: 0;
  background: #ffffff;
  border-radius: 5px;
  list-style: none;
  box-shadow: 0 15px 35px rgba(50, 50, 90, 0.1), 0 5px 15px rgba(0, 0, 0, 0.07);
  overflow: hidden;
  z-index: 999999;
  text-transform: capitalize;
}

.toggle__menu.open {
  display: block;
  opacity: 1;
}

.toggle__menu > li {
  border-left: 3px solid transparent;
  transition: ease 0.2s;
  padding: 10px;
  cursor: pointer;
}

.toggle__menu > li:hover {
  background: #e5e5e5;
  border-left-color: #3498db;
}
</style>