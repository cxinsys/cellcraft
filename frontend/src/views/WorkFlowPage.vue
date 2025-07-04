<template>
  <div class="layout__workflow" ref="captureArea">
    <section class="node-bar">
      <ul class="node-bar__nodelist" draggable="false">
        <li class="node-bar__drag-drawflow" v-for="(node, idx) in listNodes.slice(0, 4)" :key="idx" draggable="true"
          :data-node="node.name" :nodeId="node.id" @dragstart="drag($event)">
          <img class="node-bar__img" :src="node.img" draggable="false" />
        </li>
      </ul>
      <ul class="node-bar__nodelist" draggable="false">
        <li class="node-bar__drag-drawflow" v-for="(node, idx) in listNodes.slice(4, 6)" :key="idx" draggable="true"
          :data-node="node.name" :nodeId="node.id" @dragstart="drag($event)">
          <img class="node-bar__img" :src="node.img" draggable="false" />
        </li>
      </ul>
    </section>
    <div id="drawflow" @drop="drop($event)" @dragover="allowDrop($event)"></div>
    <CompileCheck v-if="compile_check" @deactivate-compile-check="deactivateCompileCheck" @run-workflow="runWorkflow" />
    <!-- <FileTable :show_files="show_files" :files_list="files_list" /> -->
    <JobTable :show_jobs="show_jobs" :taskList="taskList" @cancel-task="cancelTask" @confirm-delete="confirmDelete"
      @show-logs="showLogs" />
    <ControlBar :on_progress="on_progress" :isTabView="isTabView" @toggle-file="toggleFile"
      @save-workflow-project="saveWorkflowProject" @activate-compile-check="activateCompileCheck"
      @toggle-task="toggleTask" @toggle-tab-view="toggleTabView" @zoom-in="zoomIn" @zoom-out="zoomOut" />
    <TabComponent ref="tabComponent" :initialTabList="initialTabList" :isTabView="isTabView"
      :currentWorkflowId="currentWorkflowId" @update:isTabView="updateIsTabView"
      @process-workflow-nodes="processWorkflowNodes" />
    <div class="message" v-bind:class="{ toggleMessage: !toggleMessage }">
      <!-- <p class="message__text">{{ messageContent }}</p> -->
      <img class="message__status" src="@/assets/succes.png" v-if="messageStatus === 'success'" />
      <img class="message__status" src="@/assets/error.png" v-else-if="messageStatus === 'error'" />
      <div class="message__box">
        <p class="message__text" v-for="(content, index) in filteredMessageContent" :key="index">
          {{ content }}
        </p>
      </div>
      <img class="message__close" @click="toggleMessage = !toggleMessage" src="@/assets/close.png" />
    </div>

    <!-- 로그 모달 -->
    <div v-if="showLogsModal" class="logs-modal-overlay" @click.self="closeLogsModal">
      <div class="logs-modal">
        <div class="logs-modal-header">
          <h3>Task Logs</h3>
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
          Loading logs...
        </div>

        <div v-else-if="selectedTaskLogs" class="logs-content">
          <div class="logs-task-info">
            <p><strong>Task ID:</strong> {{ selectedTaskLogs.task_info.task_id }}</p>
            <p><strong>Status:</strong> {{ selectedTaskLogs.task_info.status }}</p>
            <p><strong>Workflow ID:</strong> {{ selectedTaskLogs.task_info.workflow_id }}</p>
            <p><strong>Algorithm ID:</strong> {{ selectedTaskLogs.task_info.algorithm_id }}</p>
          </div>

          <div v-if="selectedTaskLogs.logs.length === 0" class="no-logs">
            No log files have been generated yet.
          </div>

          <div v-else class="logs-files">
            <div v-for="(logFile, index) in selectedTaskLogs.logs" :key="index" class="log-file">
              <div class="log-file-header">
                <h4>{{ logFile.filename }}</h4>
                <span class="log-file-size">{{ formatFileSize(logFile.size) }}</span>
              </div>
              <pre class="log-file-content">{{ logFile.content }}</pre>
            </div>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import Vue from "vue";
import moment from "moment";
/* eslint-disable */

// import FileTable from "@/components/workflowComponents/PopupFileTable.vue"
import JobTable from "@/components/workflowComponents/PopupJobTable.vue"
import ControlBar from "@/components/workflowComponents/ControlBar.vue"
import TabComponent from "@/components/workflowComponents/TabComponent.vue"
import CompileCheck from "@/components/workflowComponents/CompileCheck.vue"

//노드 import (3번)
import InputFile from "@/components/nodes/InputFileNode.vue";
import DataTable from "@/components/nodes/DataTableNode.vue";
import ScatterPlot from "@/components/nodes/ScatterPlotNode.vue";
import Algorithm from "@/components/nodes/AlgorithmNode.vue";
import ResultFile from "@/components/nodes/ResultFileNode.vue";
import Visualization from "@/components/nodes/VisualizationNode.vue";

import {
  exportData,
  findWorkflow,
  saveWorkflow,
  userTaskMonitoring,
  findFolder,
  revokeTask,
  deleteTask,
  getTaskLogs,
  getPluginTemplate,
  createTaskEventSource,
} from "@/api/index";

export default {
  components: {
    // FileTable,
    JobTable,
    ControlBar,
    TabComponent,
    CompileCheck,
  },
  data() {
    return {
      exportValue: null,
      isTabView: true,
      // 왼쪽에 보여지는 노드 목록 (1번)
      listNodes: [
        // {
        //   name: "File",
        //   img: require("@/assets/file-upload2.png"),
        //   input: 0,
        //   output: 1,
        // },
        {
          name: "InputFile",
          img: require("@/assets/InputFile_logo.png"),
          input: 0,
          output: 1,
        },
        {
          name: "DataTable",
          img: require("@/assets/DataTable_logo.png"),
          input: 1,
          output: 1,
        },
        {
          name: "ScatterPlot",
          img: require("@/assets/ScatterPlot_logo.png"),
          input: 1,
          output: 1,
        },
        {
          name: "Algorithm",
          img: require("@/assets/Algorithm_logo.png"),
          input: 1,
          output: 1,
        },
        {
          name: "ResultFile",
          img: require("@/assets/ResultFile_logo.png"),
          input: 1,
          output: 1,
        },
        {
          name: "Visualization",
          img: require("@/assets/Visualization_logo.png"),
          input: 1,
          output: 0,
        },
      ],
      initialTabList: [],
      currentTab: 0,
      show_files: false,
      show_jobs: false,
      on_progress: false,
      compile_check: false,
      eventSources: {}, // Use an object to manage multiple event sources
      taskList: [],
      currentTime: new Date(),
      timeInterval: null,
      files_list: [],
      // workflow로 넘어왔을 때, 쿼리 데이터 관리
      currentWorkflowId: this.$route.query.workflow_id,
      basedPluginId: this.$route.query.plugin_id,
      toggleMessage: false,
      messageContent: "",
      messageStatus: "",
      notRemoveConnectionOutputId: "",
      showLogsModal: false,
      selectedTaskLogs: null,
      logsLoading: false,
    };
  },
  async mounted() {
    const id = document.getElementById("drawflow");
    Vue.prototype.$df = new Drawflow(id, Vue, this);
    this.$df.start();

    //노드 등록 (2번)
    this.$df.registerNode("InputFile", InputFile, {}, {});
    this.$df.registerNode("DataTable", DataTable, {}, {});
    this.$df.registerNode("ScatterPlot", ScatterPlot, {}, {});
    this.$df.registerNode("Algorithm", Algorithm, {}, {});
    this.$df.registerNode("ResultFile", ResultFile, {}, {});
    this.$df.registerNode("Visualization", Visualization, {}, {});

    // 노드 수직 연결선
    this.$df.curvature = 0.5;
    this.$df.reroute_curvature_start_end = 0;
    this.$df.reroute_curvature = 0;

    this.$df.on("nodeCreated", async (ev) => {
      // 노드 생성시 탭 생성
      const node = this.$df.getNodeFromId(ev);
      console.log(node);
      this.createNewTab(node);
      this.setCurrentWorkflowInfo();
    });

    this.$df.on("nodeRemoved", async (ev) => {
      this.removeTab(parseInt(ev));
      this.$store.commit("removeWorkflowFile", parseInt(ev));
      this.setCurrentWorkflowInfo();
    });
    this.$df.on("connectionCreated", async (ev) => {
      const input_node = this.$df.getNodeFromId(ev.input_id);
      const output_node = this.$df.getNodeFromId(ev.output_id);
      console.log(input_node, output_node);

      // InputFile 노드 : DataTable, ScatterPlot, Algorithm 제외하고 연결 불가능
      if (output_node.name === "InputFile") {
        if (input_node.name !== "DataTable" && input_node.name !== "ScatterPlot" && input_node.name !== "Algorithm") {
          this.$df.removeSingleConnection(ev.output_id, ev.input_id, ev.output_class, ev.input_class);
          this.setMessage("error", "InputFile node must be connected to DataTable, ScatterPlot, Algorithm node");
          return;
        }
      }
      // DataTable 노드 - ScatterPlot, Algorithm 제외하고 연결 불가능
      if (output_node.name === "DataTable") {
        if (input_node.name !== "ScatterPlot" && input_node.name !== "Algorithm") {
          this.$df.removeSingleConnection(ev.output_id, ev.input_id, ev.output_class, ev.input_class);
          this.setMessage("error", "DataTable node must be connected to ScatterPlot, Algorithm node");
          return;
        }
      }
      // ScatterPlot 노드 - DataTable, Algorithm 제외하고 연결 불가능
      if (output_node.name === "ScatterPlot") {
        if (input_node.name !== "DataTable" && input_node.name !== "Algorithm") {
          this.$df.removeSingleConnection(ev.output_id, ev.input_id, ev.output_class, ev.input_class);
          this.setMessage("error", "ScatterPlot node must be connected to DataTable, Algorithm node");
          return;
        }
      }
      // Algorithm 노드 - ResultFile, Visualization 제외하고 연결 불가능
      if (output_node.name === "Algorithm") {
        if (input_node.name !== "ResultFile" && input_node.name !== "Visualization") {
          this.$df.removeSingleConnection(ev.output_id, ev.input_id, ev.output_class, ev.input_class);
          this.setMessage("error", "Algorithm node must be connected to ResultFile, Visualization node");
          return;
        }
      }
      // ResultFile 노드 - Visualization 제외하고 연결 불가능
      if (output_node.name === "ResultFile") {
        if (input_node.name !== "Visualization") {
          this.$df.removeSingleConnection(ev.output_id, ev.input_id, ev.output_class, ev.input_class);
          this.setMessage("error", "ResultFile node must be connected to Visualization node");
          return;
        }
      }
      // // input_node의 name이 Algorithm, Visualization 아닌 경우, 다중 연결 검토
      // if (input_node.name !== "Algorithm" && input_node.name !== "Visualization") {
      //   // 다중 연결 시, 이전 연결 끊어주기지만,, 일단 현재 연결 끊어주는 것으로 대체하기
      //   // this.checkAndRemoveConnection(input_node, output_node);
      //   this.$df.removeSingleConnection(ev.output_id, ev.input_id, ev.output_class, ev.input_class);
      //   this.setMessage("error", `${input_node.name} node is multiple connections are not allowed`);
      // }
      console.log(ev);
      this.setCurrentWorkflowInfo();
    });
    this.$df.on("connectionRemoved", async (ev) => {
      // const input_id = this.$df.getNodeFromId(ev.input_id);
      // const output_id = this.$df.getNodeFromId(ev.output_id);

      // this.notRemoveConnectionOutputId가 ev.output_id와 같은 경우, 연결 끊어주지 않기
      if (this.notRemoveConnectionOutputId === ev.output_id) {
        this.$df.addConnection(ev.output_id, ev.input_id, ev.output_class, ev.input_classs)
        this.notRemoveConnectionOutputId = "";
        return;
      }

      this.$store.commit("removeWorkflowFile", parseInt(ev.output_id));

      this.setCurrentWorkflowInfo();
    });
    this.$df.on("nodeDataChanged", (ev) => {
      const node = this.$df.getNodeFromId(ev);
      console.log(node);

      this.setCurrentWorkflowInfo();
    });
    this.$df.on("clickEnd", (ev) => {
      // ev 값에 따라 기능 구분
      console.dir(ev.target.className);
      if (ev.detail === 2 && this.$df.node_selected) {
        this.isTabView = true;

        // 드래그 상태 해제
        this.$df.drag = false;
        this.$df.drag_point = false;
        this.$df.editor_selected = false;

        const node_id = this.$df.node_selected.id.replace(/node-/g, "");
        const node = this.$df.getNodeFromId(node_id);
        console.log(node);
        this.adjustCurrentTab(node);
      }
    });

    try {
      if (this.currentWorkflowId) {
        const workflowInfo = {
          id: this.currentWorkflowId,
        };
        const workflow_data = await findWorkflow(workflowInfo);
        console.log(workflow_data.data.workflow_info);
        this.$df.import(workflow_data.data.workflow_info);
        this.$store.commit("setTitle", workflow_data.data.title);
        this.$store.commit("setThumbnail", workflow_data.data.thumbnail);
      }
    } catch (error) {
      console.error(error);
    }

    const storeItem = localStorage.getItem("vuex");
    if (storeItem) {
      const workflow_data = JSON.parse(storeItem);
      // console.log(workflow_data);
      if (workflow_data.workflow.workflow_info) {
        this.$df.import(workflow_data.workflow.workflow_info);
        this.$store.commit("setTitle", workflow_data.workflow.title);
        this.$store.commit("setThumbnail", workflow_data.workflow.thumbnail);
      }
    }

    // Task 모니터링
    try {
      const userTasks = await userTaskMonitoring();
      console.log(userTasks.data);
      if (userTasks.data.some(task => task.status === "RUNNING" || task.status === "PENDING" || task.status === "INSTALLING")) {
        this.on_progress = true;
      }
    } catch (error) {
      console.error(error);
    }

    // plugin_id 기반으로 drawflow 템플릿 불러오기
    if (this.basedPluginId) {
      const pluginTemplate = await getPluginTemplate(this.basedPluginId);
      console.log(pluginTemplate.data);
      const drawflow_template = pluginTemplate.data.drawflow
      console.log("Drawflow Template :" + drawflow_template);
      this.$df.import(drawflow_template);
    }

    // workflow 들어오자마자 저장
    const currentWorkflow = await this.setCurrentWorkflow();
    console.log(currentWorkflow);
  },
  methods: {
    createEventSource(task_id) {
      this.on_progress = true;

      this.eventSources[task_id] = createTaskEventSource(task_id, {
        onMessage: (event) => {
          console.log("Received update: ", event.data);
          // PENDING 상태도 진행 중으로 처리
          if (event.data === "RUNNING" || event.data === "PENDING" || event.data === "INSTALLING") {
            this.on_progress = true;
          }
          if (event.data === "INSTALLING") {
            // taskList에서 해당 task에 대해서 status를 INSTALLING으로 변경
            const taskIndex = this.taskList.findIndex(task => task.task_id === task_id);
            this.$set(this.taskList, taskIndex, {
              ...this.taskList[taskIndex],
              status: "INSTALLING",
            });
            console.log(this.taskList);
          }
        },
        onComplete: (status) => {
          console.log("close");
          this.on_progress = false;
          this.closeEventSource(task_id);
          clearInterval(this.timeInterval);
        },
        onError: (error) => {
          console.error("SSE Error:", error);
          this.on_progress = false;
          this.closeEventSource(task_id);
          clearInterval(this.timeInterval);
        }
      });
    },
    closeEventSource(task_id) {
      // If the event source exists, close it
      if (this.eventSources[task_id]) {
        this.eventSources[task_id].close();
        delete this.eventSources[task_id];
      }
    },
    activateCompileCheck() {
      this.compile_check = true;
      this.isTabView = false;
    },
    deactivateCompileCheck() {
      this.compile_check = false;
    },
    async runWorkflow() {
      try {
        this.updateWorkflowInfo();
        this.exportValue = this.$df.export();
        const title = this.$store.getters.getTitle;
        const thumbnail = this.$store.getters.getThumbnail;
        // console.log(JSON.stringify(this.exportValue));
        console.log(this.$df.drawflow.drawflow[this.$df.module]);
        const workflow = {
          id: this.currentWorkflowId,
          title: title,
          thumbnail: thumbnail,
          workflow_info: this.exportValue,
        };
        const workflow_data = await exportData(workflow);

        // task_ids가 배열이므로 각 태스크 ID에 대해 createEventSource 실행
        if (workflow_data.data.task_ids && Array.isArray(workflow_data.data.task_ids)) {
          workflow_data.data.task_ids.forEach(task_id => {
            this.createEventSource(task_id);
          });
        }

        if (this.show_jobs) {
          this.show_jobs = false;
          this.toggleTask();
        }
      } catch (error) {
        console.error(error);
      }
    },
    updateWorkflowInfo() {
      const workflow_info = this.$store.getters.getWorkflowInfo
      this.$df.import(workflow_info);
    },
    processWorkflowNodes() {
      const workflow_info = this.$store.getters.getWorkflowInfo;
      const nodes = workflow_info.drawflow.Home.data;

      // console.log("workflow_info from store :", JSON.stringify(workflow_info, null, 2));
      // console.log("nodes from drawflow :", JSON.stringify(nodes, null, 2));

      for (const nodeId in nodes) {
        if (nodes.hasOwnProperty(nodeId)) {
          const node = nodes[nodeId];
          // console.log("node :", JSON.stringify(node, null, 2));
          this.$df.updateNodeDataFromId(node.id, node.data);
        }
      }
    },
    drag(event) {
      event.dataTransfer.setData(
        "node",
        event.target.getAttribute("data-node")
      );
    },
    drop(event) {
      event.preventDefault();
      const data = event.dataTransfer.getData("node");
      this.addNodeToDrawFlow(data, event.clientX, event.clientY);
    },
    allowDrop(event) {
      event.preventDefault();
    },
    zoomIn() {
      this.$df.zoom_in();
    },
    zoomOut() {
      this.$df.zoom_out();
    },
    addNodeToDrawFlow(name, pos_x, pos_y) {
      pos_x =
        pos_x *
        (this.$df.precanvas.clientWidth /
          (this.$df.precanvas.clientWidth * this.$df.zoom)) -
        this.$df.precanvas.getBoundingClientRect().x *
        (this.$df.precanvas.clientWidth /
          (this.$df.precanvas.clientWidth * this.$df.zoom));
      pos_y =
        pos_y *
        (this.$df.precanvas.clientHeight /
          (this.$df.precanvas.clientHeight * this.$df.zoom)) -
        this.$df.precanvas.getBoundingClientRect().y *
        (this.$df.precanvas.clientHeight /
          (this.$df.precanvas.clientHeight * this.$df.zoom));

      const nodeSelected = this.listNodes.find((ele) => ele.name === name);
      this.$df.addNode(
        name,
        nodeSelected.input,
        nodeSelected.output,
        pos_x,
        pos_y,
        name,
        { "title": name },
        name,
        "vue"
      );
    },

    async saveWorkflowProject() {
      try {
        const currentWorkflow = await this.setCurrentWorkflow();
        console.log(currentWorkflow);
        this.setMessage("success", "Save workflow successfully!");
      } catch (error) {
        console.error(error);
      }
    },
    async captureWorkflow() {
      // console.log("이미지 캡쳐");
      try {
        await html2canvas(this.$refs.captureArea).then(canvas => {
          const originalCanvasWidth = canvas.width;
          const originalCanvasHeight = canvas.height;

          const newCanvasWidth = originalCanvasWidth / 8;
          const newCanvasHeight = originalCanvasHeight / 8;

          const resizedCanvas = document.createElement('canvas');
          resizedCanvas.width = newCanvasWidth;
          resizedCanvas.height = newCanvasHeight;

          const ctx = resizedCanvas.getContext('2d');
          ctx.drawImage(canvas, 0, 0, newCanvasWidth, newCanvasHeight);

          const thumbnailURL = resizedCanvas.toDataURL("image/png");
          this.$store.commit("setThumbnail", thumbnailURL);
          // console.log(thumbnailURL);
        });
      } catch (error) {
        console.error(error);
      }
    },
    async toggleTask() {
      try {
        if (!this.show_jobs) {
          const user_tasks = await userTaskMonitoring();
          console.log(user_tasks);
          this.taskList = user_tasks.data;

          this.taskList.forEach((task, idx) => {
            if (task.status === "SUCCESS" ||
              task.status === "FAILURE" ||
              task.status === "REVOKED" ||
              task.status === "RETRY") {
              this.taskList[idx].running_time = this.getTimeDifference(
                task.start_time,
                task.end_time
              );
            } else if (task.status === "RUNNING" || task.status === "PENDING" || task.status === "INSTALLING") {
              // RUNNING 또는 PENDING 상태일 때 타이머 시작
              this.timeInterval = this.startTimer(idx);
              // Task가 실행 중이므로 on_progress를 true로 설정
              this.on_progress = true;
            }
          });
        } else {
          clearInterval(this.timeInterval);
        }

        setTimeout(() => {
          this.show_jobs = !this.show_jobs;
        }, 300);
      } catch (error) {
        console.error(error);
        this.setMessage(
          "error",
          "No tasks have been executed yet. Please run workflow"
        );
      }
    },
    async cancelTask(task_id) {
      try {
        const revoke_task = await revokeTask(task_id);
        console.log(revoke_task);
        this.toggleTask();
        this.setMessage("success", "Cancel task successfully!");
      } catch (error) {
        console.error(error);
      }
    },
    async confirmDelete(task_id) {
      try {
        const delete_task = await deleteTask(task_id);
        console.log(delete_task);
        this.toggleTask();
        this.setMessage("success", "Delete task successfully!");
      } catch (error) {
        console.error(error);
      }
    },
    async showLogs(task_id) {
      try {
        this.logsLoading = true;
        this.showLogsModal = true;

        const response = await getTaskLogs(task_id);
        console.log(response.data);

        this.selectedTaskLogs = response.data;
        this.setMessage("success", "Logs loaded successfully!");
      } catch (error) {
        console.error(error);
        this.setMessage("error", "Failed to load logs: " + error.message);
        this.showLogsModal = false;
      } finally {
        this.logsLoading = false;
      }
    },
    async refreshLogs() {
      if (this.selectedTaskLogs && this.selectedTaskLogs.task_info) {
        await this.showLogs(this.selectedTaskLogs.task_info.task_id);
      }
    },
    closeLogsModal() {
      this.showLogsModal = false;
      this.selectedTaskLogs = null;
      this.logsLoading = false;
    },
    formatFileSize(bytes) {
      if (bytes === 0) return '0 Bytes';
      const k = 1024;
      const sizes = ['Bytes', 'KB', 'MB', 'GB'];
      const i = Math.floor(Math.log(bytes) / Math.log(k));
      return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
    },
    async toggleFile() {
      if (!this.show_files) {
        try {
          const filesList = await findFolder({
            folder_name: "data",
          });
          console.log(filesList.data);
          this.files_list = filesList.data;
        } catch (error) {
          this.show_files = !this.show_files;
          console.error(error);
          this.setMessage(
            "error",
            "No files have been uploaded yet. Please upload files"
          );
        }
      }
      this.show_files = !this.show_files;
    },
    startTimer(idx) {
      const interval = setInterval(() => {
        if (!this.on_progress) {
          this.show_jobs = false;
          clearInterval(interval);
        }

        // RUNNING 또는 PENDING 상태일 때 running_time 계산
        if (this.taskList[idx].status === "RUNNING" ||
          this.taskList[idx].status === "PENDING" || this.taskList[idx].status === "INSTALLING") {
          let currentTime = new Date();
          let running_time = this.getRunningTime(
            this.taskList[idx].start_time,
            currentTime
          );

          this.$set(this.taskList, idx, {
            ...this.taskList[idx],
            running_time: running_time,
          });
        }
      }, 1000);
      return interval;
    },
    getTimeDifference(start_time, end_time) {
      const start = new Date(start_time);
      const end = new Date(end_time);
      let diff = Math.abs(end - start); // Difference in milliseconds
      let hours = Math.floor(diff / 3600000);
      diff -= hours * 3600000;

      const minutes = Math.floor(diff / 60000);
      diff -= minutes * 60000;

      const seconds = Math.floor(diff / 1000);

      return `${hours.toString().padStart(2, "0")}:${minutes
        .toString()
        .padStart(2, "0")}:${seconds.toString().padStart(2, "0")}`;
    },
    getRunningTime(startTime, currentTime) {
      const start = new Date(startTime);
      let diff = Math.abs(currentTime - start); // Difference in milliseconds
      let hours = Math.floor(diff / 3600000);
      diff -= hours * 3600000;

      const minutes = Math.floor(diff / 60000);
      diff -= minutes * 60000;

      const seconds = Math.floor(diff / 1000);

      hours = hours;
      console.log(
        `${hours.toString().padStart(2, "0")}:${minutes
          .toString()
          .padStart(2, "0")}:${seconds.toString().padStart(2, "0")}`
      );

      return `${hours.toString().padStart(2, "0")}:${minutes
        .toString()
        .padStart(2, "0")}:${seconds.toString().padStart(2, "0")}`;
    },
    async getWorkflowTitle(id) {
      const workflow_id = {
        id: id,
      };
      const user_workflow = await findWorkflow(workflow_id);
      console.log(user_workflow.data.title);
      return user_workflow.data.title;
    },
    setMessage(status, content) {
      this.toggleMessage = true;
      this.messageStatus = status;
      this.messageContent = content;
      setTimeout(() => {
        this.toggleMessage = false;
      }, 5000);
    },
    async setCurrentWorkflow() {
      try {
        this.setCurrentWorkflowInfo();
        await this.captureWorkflow();
        const title = this.$store.getters.getTitle;
        const thumbnail = this.$store.getters.getThumbnail;
        const workflow = {
          id: this.currentWorkflowId,
          title: title,
          thumbnail: thumbnail,
          workflow_info: this.exportValue,
        };
        console.log("currentWorkflowId : " + workflow.id + "type : " + typeof workflow.id);
        const workflow_data = await saveWorkflow(workflow);
        this.currentWorkflowId = workflow_data.data.id;
        return workflow_data.data;
      } catch (error) {
        console.error(error);
      }
    },
    setCurrentWorkflowInfo() {
      this.exportValue = this.$df.export();
      this.$store.commit("setWorkflow", this.exportValue);
    },
    getNodeTitleById(id) {
      const node = this.$df.getNodeFromId(id);
      return node.data.title;
    },
    toggleTabView() {
      this.isTabView = !this.isTabView;
    },
    downloadDrawflow() {
      const drwaflow = this.$df.export();
      const blob = new Blob([JSON.stringify(drwaflow)], {
        type: "application/json",
      });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "workflow.json";
      a.click();
    },
    updateIsTabView(newIsTabView) {
      this.isTabView = newIsTabView;
    },
    createNewTab(node) {
      this.$refs.tabComponent.createTab(node);
    },
    adjustCurrentTab(node) {
      this.$refs.tabComponent.adjustCurrentTab(node);
    },
    removeTab(id) {
      this.$refs.tabComponent.removeTab(id);
    },
    checkAndRemoveConnection(input_node, output_node) {
      // Check if input_node.inputs.input_1.connections exists and has a length of 2 or more
      const input_node_id = String(input_node.id);
      const output_node_id = String(output_node.id);
      if (input_node.inputs && input_node.inputs.input_1 && input_node.inputs.input_1.connections.length >= 2) {
        // Find the connection with output_node.id
        const connectionIndex = input_node.inputs.input_1.connections.findIndex(connection => connection.node === output_node_id);
        console.log(connectionIndex);
        // If the connection with output_node.id exists
        if (connectionIndex !== -1) {
          const otherConnection = input_node.inputs.input_1.connections.find((_, index) => index !== connectionIndex);
          if (otherConnection) {
            this.$df.removeConnectionNodeId("node-" + otherConnection.node);
            this.notRemoveConnectionOutputId = otherConnection.node;
            // this.$df.removeSingleConnection(otherConnection.node, input_node_id, "output_1", "input_1");
            console.log(otherConnection.node, input_node_id, "output_1", "input_1");
          }
        }
      }
    }

  },
  beforeDestroy() {
    // Close all the event source connections before the component is destroyed
    for (let task_id in this.eventSources) {
      this.closeEventSource(task_id);
    }
  },
  filters: {
    titleNone(value) {
      if (value === "" || !value) return "Untitled";
      return value;
    },
    formatBytes(a, b) {
      if (a === 0) return "0 Bytes";
      const c = 1024;
      const d = b || 2;
      const e = ["Bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"];
      const f = Math.floor(Math.log(a) / Math.log(c));

      return parseFloat((a / Math.pow(c, f)).toFixed(d)) + " " + e[f];
    },
    cutFromT(value) {
      return value.split("T")[0];
    },
    cutFromDotName(value) {
      return value.split(".")[0];
    },
    cutFromDotType(value) {
      return value.split(".")[1];
    },
    formatDateTime(dateTime) {
      const date = moment(dateTime).format("MMMM Do, HH:mm");
      if (date === "Invalid date") return "Not Yet Completed";
      return date;
    },
  },
  computed: {
    filteredMessageContent() {
      // The split() function is used to divide the messageContent string into an array of substrings,
      // using the dot character as the delimiter.
      return this.messageContent.split(".").filter(Boolean);
    },
  },
};
</script>

<style>
.layout__workflow {
  width: 100%;
  height: 100%;
  position: relative;
  overflow: hidden;
  background: #191f26;
  /* background-image: url("@/assets/fantastic_background3.png"); */
  background-size: cover;
  background-position: center center;
  background-repeat: no-repeat;
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 1rem;
}

.node-bar {
  width: 5rem;
  border-radius: 16px;
  background-color: rgba(255, 255, 255, 0.1);
  box-shadow: 0px 0px 0px 1px rgba(255, 255, 255, 0.1);
  z-index: 9998;
  opacity: 1;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}

/* .node-bar_output {
  width: 80px;
  height: 400px;
  border-radius: 16px;
  background-color: rgba(69, 69, 69, 0.5);
  box-shadow: 0px 0px 0px 1px rgba(255, 255, 255, 0.1);
  position: absolute;
  top: calc(50% - 208px);
  right: 12px;
  z-index: 9998;
  opacity: 1;
  display: flex;
  align-items: center;
  justify-content: center;
  backdrop-filter: blur(10px);
} */

.node-bar__nodelist {
  width: 80%;
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 1.5rem;
  padding: 1rem 0;
  border-bottom: 0.1rem solid rgba(255, 255, 255, 0.1);
}

.node-bar__nodelist:last-child {
  border-bottom: none;
}

.node-bar__drag-drawflow {
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: move;
  /* background: rgba(69, 69, 69); */
  background: rgb(228, 228, 228);
  /* color: rgb(245, 245, 245); */
  border-radius: 1rem;
  /* width: 5rem;
  height: 5rem; */
  width: 4rem;
  height: 4rem;
  box-shadow: 0px 0px 0px 1px rgba(255, 255, 255, 0.3);
}

.node-bar__img {
  /* width: 3rem;
  height: 3rem; */
  width: 2.2rem;
  height: 2.2rem;
  object-fit: contain;
  -webkit-user-drag: none;
  -khtml-user-drag: none;
  -moz-user-drag: none;
  -o-user-drag: none;
}

.white {
  filter: invert(100%) sepia(75%) saturate(0%) hue-rotate(51deg) brightness(115%) contrast(101%);
}

.control-bar__icon:hover,
.loader.drawflow-node:hover,
.loader_done:hover {
  opacity: 0.7;
  /* You can adjust this value to your liking */
  transform: scale(1.1);
}

#drawflow {
  width: calc(100% - 8rem);
  height: calc(100% - 2rem);
  border-radius: 0.8%;
  box-shadow: 0px 0px 1px 1px rgba(255, 255, 255, 0.2);
  backdrop-filter: blur(10px);

  background: var(--dfBackgroundColor);
  background-size: var(--dfBackgroundSize) var(--dfBackgroundSize);
  background-image: var(--dfBackgroundImage);
}

#drawflow .drawflow .drawflow-node {
  display: var(--dfNodeType);
  background: var(--dfNodeBackgroundColor);
  backdrop-filter: blur(5px);
  color: var(--dfNodeTextColor);
  border: var(--dfNodeBorderSize) solid var(--dfNodeBorderColor);
  border-radius: var(--dfNodeBorderRadius);
  min-height: var(--dfNodeMinHeight);
  width: 4rem;
  height: 4rem;
  min-width: var(--dfNodeMinWidth);
  padding-top: var(--dfNodePaddingTop);
  padding-bottom: var(--dfNodePaddingBottom);
  padding-left: var(--dfNodePaddingLeft);
  padding-right: var(--dfNodePaddingRight);
  -webkit-box-shadow: var(--dfNodeBoxShadowHL) var(--dfNodeBoxShadowVL) var(--dfNodeBoxShadowBR) var(--dfNodeBoxShadowS) var(--dfNodeBoxShadowColor);
  box-shadow: var(--dfNodeBoxShadowHL) var(--dfNodeBoxShadowVL) var(--dfNodeBoxShadowBR) var(--dfNodeBoxShadowS) var(--dfNodeBoxShadowColor);
}

#drawflow .drawflow .drawflow-node:hover {
  background: var(--dfNodeHoverBackgroundColor);
  color: var(--dfNodeHoverTextColor);
  border: var(--dfNodeHoverBorderSize) solid var(--dfNodeHoverBorderColor);
  border-radius: var(--dfNodeHoverBorderRadius);
  -webkit-box-shadow: var(--dfNodeHoverBoxShadowHL) var(--dfNodeHoverBoxShadowVL) var(--dfNodeHoverBoxShadowBR) var(--dfNodeHoverBoxShadowS) var(--dfNodeHoverBoxShadowColor);
  box-shadow: var(--dfNodeHoverBoxShadowHL) var(--dfNodeHoverBoxShadowVL) var(--dfNodeHoverBoxShadowBR) var(--dfNodeHoverBoxShadowS) var(--dfNodeHoverBoxShadowColor);
}

#drawflow .drawflow .drawflow-node.selected {
  background: var(--dfNodeSelectedBackgroundColor);
  color: var(--dfNodeSelectedTextColor);
  border: var(--dfNodeSelectedBorderSize) solid var(--dfNodeSelectedBorderColor);
  border-radius: var(--dfNodeSelectedBorderRadius);
  -webkit-box-shadow: var(--dfNodeSelectedBoxShadowHL) var(--dfNodeSelectedBoxShadowVL) var(--dfNodeSelectedBoxShadowBR) var(--dfNodeSelectedBoxShadowS) var(--dfNodeSelectedBoxShadowColor);
  box-shadow: var(--dfNodeSelectedBoxShadowHL) var(--dfNodeSelectedBoxShadowVL) var(--dfNodeSelectedBoxShadowBR) var(--dfNodeSelectedBoxShadowS) var(--dfNodeSelectedBoxShadowColor);
}

#drawflow .drawflow .drawflow-node .input {
  left: var(--dfInputLeft);
  background: var(--dfInputBackgroundColor);
  border: var(--dfInputBorderSize) solid var(--dfInputBorderColor);
  border-radius: var(--dfInputBorderRadius);
  height: var(--dfInputHeight);
  width: var(--dfInputWidth);
  box-shadow: 0px 0px 1px 1px rgba(255, 255, 255, 0.15);
}

#drawflow .drawflow .drawflow-node .input:hover {
  background: var(--dfInputHoverBackgroundColor);
  border: var(--dfInputHoverBorderSize) solid var(--dfInputHoverBorderColor);
  border-radius: var(--dfInputHoverBorderRadius);
}

#drawflow .drawflow .drawflow-node .outputs {
  float: var(--dfNodeTypeFloat);
}

#drawflow .drawflow .drawflow-node .output {
  right: var(--dfOutputRight);
  background: var(--dfOutputBackgroundColor);
  border: var(--dfOutputBorderSize) solid var(--dfOutputBorderColor);
  border-radius: var(--dfOutputBorderRadius);
  height: var(--dfOutputHeight);
  width: var(--dfOutputWidth);
  box-shadow: 0px 0px 1px 1px rgba(255, 255, 255, 0.15);
}

#drawflow .drawflow .drawflow-node .output:hover {
  background: var(--dfOutputHoverBackgroundColor);
  border: var(--dfOutputHoverBorderSize) solid var(--dfOutputHoverBorderColor);
  border-radius: var(--dfOutputHoverBorderRadius);
}

#drawflow .drawflow .connection .main-path {
  stroke-width: var(--dfLineWidth);
  stroke: var(--dfLineColor);
}

#drawflow .drawflow .connection .main-path:hover {
  stroke: var(--dfLineHoverColor);
}

#drawflow .drawflow .connection .main-path.selected {
  stroke: var(--dfLineSelectedColor);
}

#drawflow .drawflow .connection .point {
  stroke: var(--dfRerouteBorderColor);
  stroke-width: var(--dfRerouteBorderWidth);
  fill: var(--dfRerouteBackgroundColor);
}

#drawflow .drawflow .connection .point:hover {
  stroke: var(--dfRerouteHoverBorderColor);
  stroke-width: var(--dfRerouteHoverBorderWidth);
  fill: var(--dfRerouteHoverBackgroundColor);
}

#drawflow .drawflow-delete {
  content: "";
  color: rgba(0, 0, 0, 0);
  display: var(--dfDeleteDisplay);
  background: var(--dfDeleteBackgroundColor);
  border: var(--dfDeleteBorderSize) solid var(--dfDeleteBorderColor);
  border-radius: var(--dfDeleteBorderRadius);
  width: 15px;
  height: 15px;
}

#drawflow .drawflow-delete::before,
#drawflow .drawflow-delete::after {
  font-size: x-large;
  color: var(--dfDeleteColor);
  content: '-';
  position: absolute;
  /* left: 0px; */
  top: -6px;
  right: 0;
  bottom: 0;
  width: 15px;
  height: 15px;
}

#drawflow .drawflow-delete::before {
  transform: rotate(45deg);
  left: 6px;
}

#drawflow .drawflow-delete::after {
  transform: rotate(-45deg);
  left: -6px;
}

#drawflow .parent-node .drawflow-delete {
  top: var(--dfDeleteTop);
  right: var(--dfDeleteRight);
  border-radius: var(--dfDeleteHoverBorderRadius);
}

#drawflow .drawflow-delete:hover {
  color: var(--dfDeleteHoverColor);
  background: var(--dfDeleteHoverBackgroundColor);
  border: var(--dfDeleteHoverBorderSize) solid var(--dfDeleteHoverBorderColor);
  border-radius: var(--dfDeleteHoverBorderRadius);
}

.vdr-stick {
  opacity: 0;
}

.vdr.active:before {
  border: none;
  outline: none !important;
}

.message {
  width: 30rem;
  height: 6rem;
  display: flex;
  align-items: center;
  justify-content: space-around;
  position: absolute;
  bottom: 7rem;
  left: calc(50% - 16rem);
  background: rgba(0, 0, 0, 0.8);
  border-radius: 1rem;
  padding: 0 1rem;
}

.message__status {
  width: 2rem;
  height: 2rem;
  object-fit: contain;
  margin: 0 1rem;
}

.message__box {
  display: flex;
  flex-direction: column;
  justify-content: center;
  width: 100%;
  height: 100%;
}

.message__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.1rem;
  line-height: 1.4rem;
  color: #ffffff;
}

.message__close {
  cursor: pointer;
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  margin: 0 0.5rem;
  opacity: 0.5;
  filter: invert(100%) sepia(0%) saturate(0%) hue-rotate(0deg) brightness(100%) contrast(100%);
}

.toggleMessage {
  display: none;
}

.node-zoom-buttons {
  position: absolute;
  bottom: 2rem;
  right: 7rem;
  display: flex;
  gap: 0.5rem;
  z-index: 9999;
}

.node-zoom-button {
  width: 3rem;
  height: 3rem;
  background-color: #007BFF;
  color: white;
  border: none;
  border-radius: 6px;
  cursor: pointer;
}

.node-zoom-button img {
  width: 100%;
  height: 100%;
  object-fit: contain;
}

.margin__top-4 {
  margin-top: 4px;
}

/* 로그 모달 스타일 */
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

.refresh-icon {
  width: 1.2rem;
  height: 1.2rem;
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
  text-align: center;
  color: #bdc3c7;
  font-style: italic;
  padding: 3rem;
  font-size: 1.1rem;
  background: rgba(31, 42, 56, 0.5);
  border-radius: 12px;
  border: 2px dashed rgba(255, 255, 255, 0.2);
}

.logs-files {
  max-height: 450px;
  overflow-y: auto;
}

/* 로그 파일 스크롤바 */
.logs-files::-webkit-scrollbar {
  width: 6px;
}

.logs-files::-webkit-scrollbar-track {
  background: #2c3e50;
  border-radius: 6px;
}

.logs-files::-webkit-scrollbar-thumb {
  background: #576574;
  border-radius: 6px;
}

.logs-files::-webkit-scrollbar-thumb:hover {
  background: #5a6c7d;
}

.log-file {
  margin-bottom: 1.5rem;
  border: 1px solid rgba(255, 255, 255, 0.15);
  border-radius: 12px;
  overflow: hidden;
  background: #2c3e50;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.2);
}

.log-file-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 1rem 1.25rem;
  background: #1f2a38;
  border-bottom: 1px solid rgba(255, 255, 255, 0.1);
}

.log-file-header h4 {
  margin: 0;
  font-size: 1.05rem;
  color: #ecf0f1;
  font-weight: 600;
  font-family: 'Courier New', monospace;
}

.log-file-size {
  font-size: 0.85rem;
  color: #95a5a6;
  background: rgba(52, 73, 94, 0.7);
  padding: 0.25rem 0.5rem;
  border-radius: 6px;
  font-family: 'Courier New', monospace;
}

.log-file-content {
  margin: 0;
  padding: 1.25rem;
  background: #1e2833;
  white-space: pre-wrap;
  word-wrap: break-word;
  max-height: 350px;
  overflow-y: auto;
  font-family: 'Courier New', monospace;
  font-size: 0.9rem;
  line-height: 1.6;
  color: #e8f4f8;
  border: none;
  font-weight: 400;
}

/* 로그 내용 스크롤바 */
.log-file-content::-webkit-scrollbar {
  width: 6px;
}

.log-file-content::-webkit-scrollbar-track {
  background: #1a252f;
  border-radius: 6px;
}

.log-file-content::-webkit-scrollbar-thumb {
  background: #34495e;
  border-radius: 6px;
}

.log-file-content::-webkit-scrollbar-thumb:hover {
  background: #4a6378;
}

/* @media (prefers-color-scheme: dark) {
  .node-bar__img {
    filter: invert(100%) sepia(3%) saturate(2008%) hue-rotate(348deg)
      brightness(125%) contrast(111%);
  }
  .tab__item {
    background: rgb(32, 34, 39);
    color: rgb(255, 255, 255);
    box-shadow: 0px -6px 5px 0px rgba(0, 0, 0, 0.5);
  }
  .currentTab {
    background: rgb(53, 55, 60);
    box-shadow: 0px -6px 5px 0px rgba(0, 0, 0, 0.5);
  }
  .content-view {
    background: rgb(53, 55, 60);
    box-shadow: 0px -5px 5px 0px rgba(0, 0, 0, 0.5);
  }
  .node-bar {
    background: rgb(53, 55, 60);
    box-shadow: 0px 0px 5px 0px rgba(0, 0, 0, 1);
  }
  .node-bar__drag-drawflow {
    background: rgb(32, 34, 39);
    color: rgb(245, 245, 245);
  }
} */
</style>
