<template>
  <div class="control-popup__jobs" v-if="show_jobs">
    <div class="job-table-container">
      <table class="job-table">
        <thead>
          <tr>
            <th>No.</th>
            <th>Name</th>
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
            <td>{{ task.title | titleNone }}</td>
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
    cancelTask() {
      this.$emit('cancel-task', this.currentTaskId);
      this.R_Mouse_isActive = false;
    },
    confirmDelete() {
      this.$emit('confirm-delete', this.currentTaskId);
      this.R_Mouse_isActive = false;
    },
    getStatusClass(status) {
      if (status === "SUCCESS") return "status-success";
      if (status === "FAILURE" || status === "REVOKED" || status === "RETRY") return "status-failure";
      if (status === "RUNNING" || status === "PENDING" || status === "INSTALLING") return "status-running";
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
/* 기존 포지션 유지 */
.control-popup__jobs {
  /* width: 720px;
  max-width: 720px; */
  height: 540px;
  max-height: 540px;
  left: calc(50%);
  overflow-y: auto;
  border-radius: 16px;
  position: absolute;
  bottom: 98px;
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
  font-size: 0.8rem;
}

.job-table tbody tr {
  transition: background-color 0.2s ease-in-out;
}

.job-table tbody tr:hover {
  background-color: #3d566e;
}


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
</style>