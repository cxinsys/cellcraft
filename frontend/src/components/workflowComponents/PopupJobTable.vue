<template>
  <div class="control-popup__jobs" v-if="show_jobs">
    <table class="control-popup__table">
      <thead>
        <tr>
          <th>No.</th>
          <th>Name</th>
          <th>Start</th>
          <th>End</th>
          <th>Time</th>
          <th>Status</th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="(task, index) in taskList" :key="index">
          <td>{{ index + 1 }}</td>
          <td>{{ task.title | titleNone }}</td>
          <td>{{ task.start_time | formatDateTime }}</td>
          <td>{{ task.end_time | formatDateTime }}</td>
          <td>{{ task.running_time }}</td>
          <td class="task-status">
            <div class="status-box__red" v-if="
              task.status === 'FAILURE' ||
              task.status === 'REVOKED' ||
              task.status === 'RETRY'
            "></div>
            <div class="status-box__yellow" v-if="
              task.status === 'RUNNING' ||
              task.status === 'PENDING'
            "></div>
            <div class="status-box__green" v-if="task.status === 'SUCCESS'"></div>
            {{ task.status }}
          </td>
          <td>
            <img v-if="task.status === 'RUNNING' || task.status === 'PENDING'" @click="cancelTask(task.task_id)"
              class="control-bar__icon" src="@/assets/multiply.png" />
          </td>
        </tr>
      </tbody>
    </table>
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
  methods: {
    cancelTask(taskId) {
      this.$emit('cancel-task', taskId);
    }
  },
  filters: {
    formatDateTime(dateTime) {
      const date = moment(dateTime).format("MMMM Do, HH:mm");
      if (date === "Invalid date") return "Not Yet Completed";
      return date;
    },
  },
};
</script>

<style scoped>
.control-popup__jobs {
  /* width: 8rem; */
  width: 40vw;
  max-width: 400px;
  /* height: 34rem; */
  height: 30vh;
  max-height: 300px;

  border-radius: 16px;
  background: rgba(244, 246, 251, 0.586);
  box-shadow: 0px 0px 5px 0px rgba(0, 0, 0, 1);
  position: absolute;
  bottom: 98px;
  z-index: 9998;
  opacity: 1;
  display: flex;
  align-items: center;
  justify-content: center;
}

.control-popup__jobs {
  max-width: 720px;
  max-height: 540px;
  width: 720px;
  height: 540px;
  left: calc(50% + 1vw);
  overflow-y: auto;
  border-radius: 16px;
  /* or whatever radius you prefer */
}

.control-popup__jobs::-webkit-scrollbar {
  width: 10px;
  /* width of the entire scrollbar */
}

.control-popup__jobs::-webkit-scrollbar-track {
  background: #f1f1f1;
  /* color of the tracking area */
  border-radius: 16px;
  /* keep the same radius as the container */
}

.control-popup__jobs::-webkit-scrollbar-thumb {
  background: #888;
  /* color of the scroll thumb */
  border-radius: 16px;
  /* keep the same radius as the container */
}

.control-popup__jobs::-webkit-scrollbar-thumb:hover {
  background: #555;
  /* color of the scroll thumb on hover */
}

.control-popup__table {
  width: 95%;
  height: auto;
  margin: auto;
  border-collapse: collapse;
  position: absolute;
  top: 20px;
}

.control-popup__table thead {
  height: 26px;
  font-weight: 500;
  color: rgb(49, 49, 49);
  border-bottom: 1px solid #6767678c;
}

.control-popup__table td.task-status {
  display: flex;
  /* align items horizontally */
  align-items: center;
  /* center items vertically */
  justify-content: center;
  /* center items horizontally */
}

.control-popup__table td {
  vertical-align: middle;
  font-weight: 400;
  text-align: center;
  color: rgb(68, 68, 68);
  padding: 0.7rem;
  margin: 1rem;
}

.control-popup__table__progress {
  width: 40%;
}

.status-box__red,
.status-box__green,
.status-box__yellow {
  width: 0.5rem;
  height: 0.5rem;
  border-radius: 50%;
  margin-right: 5px;
}

.status-box__red {
  background-color: #ff4444;
  /* 더 선명한 빨간색 */
}

.status-box__green {
  background-color: #00C851;
  /* 더 선명한 초록색 */
}

.status-box__yellow {
  background-color: #ffbb33;
  /* 더 선명한 노란색 */
}

.progress-bar {
  width: 100%;
  height: 5px;
  background-color: #eee;
  border-radius: 10px;
  overflow: hidden;
}

.progress {
  height: 100%;
  background-color: #3a98fc;
  transition: width 0.3s;
  border-radius: 10px;
}

.control-bar__icon {
  width: 1rem;
  height: 1rem;
  cursor: pointer;
}
</style>