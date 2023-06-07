<template>
  <div class="layout__workflow">
    <div id="drawflow" @drop="drop($event)" @dragover="allowDrop($event)"></div>
    <section class="node-bar">
      <ul class="node-bar__nodelist" draggable="false">
        <li
          class="node-bar__drag-drawflow"
          v-for="(node, idx) in listNodes"
          :key="idx"
          draggable="true"
          :data-node="node.name"
          @dragstart="drag($event)"
        >
          <img class="node-bar__img" :src="node.img" draggable="false" />
        </li>
      </ul>
    </section>
    <div class="control-popup__files" v-if="show_files">
      <table class="control-popup__table">
        <thead>
          <tr>
            <th>Name</th>
            <th>Date</th>
            <th>Type</th>
            <th>Size</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>{{ "name" }}</td>
            <td>{{ "date" }}</td>
            <td>{{ "type" }}</td>
            <td>{{ "size" }}</td>
          </tr>
        </tbody>
      </table>
    </div>
    <div class="control-popup__jobs" v-if="show_jobs">
      <table class="control-popup__table">
        <thead>
          <tr>
            <th>No.</th>
            <th>Name</th>
            <th>Running time</th>
            <th>Status</th>
          </tr>
        </thead>
        <tbody>
          <tr v-for="(task, index) in taskList" :key="index">
            <td>{{ index + 1 }}</td>
            <td>{{ task.title }}</td>
            <td>{{ task.running_time }}</td>
            <td>{{ task.status }}</td>
          </tr>
        </tbody>
      </table>
    </div>
    <section class="control-bar">
      <ul class="control-bar__btnList">
        <li class="control-bar__button" @click="show_files = !show_files">
          <img class="control-bar__icon" src="@/assets/control_files.png" />
        </li>
        <li class="control-bar__button" @click="saveWorkflowProject">
          <img class="control-bar__icon" src="@/assets/control_save.png" />
        </li>
        <li class="control-bar__button">
          <button class="run_button" @click="exportdf">
            <img class="control-bar__icon" src="@/assets/control_run.png" />
          </button>
        </li>
        <li>
          <div
            class="loader"
            @click="toggleTask"
            v-if="on_progress == true"
          ></div>
          <div class="loader_done" @click="toggleTask" v-else></div>
        </li>
        <li class="control-bar__button">
          <img class="control-bar__icon" src="@/assets/control_export.png" />
        </li>
      </ul>
    </section>
    <VueDragResize
      contentClass="content-component"
      v-if="tabList.length != 0 && isTabView"
      :isActive="true"
      :x="600"
      :y="64"
      :w="880"
      :h="672"
      :minw="820"
      :minh="540"
      :stickSize="14"
      :sticks="['tl', 'ml', 'tr', 'bl', 'br']"
    >
      <ul class="content-tab" v-if="tabList.length != 0 && isTabView">
        <li
          class="tab__item"
          v-for="(tab, idx) in tabList"
          :key="idx"
          v-bind:class="{ currentTab: currentTab === idx }"
          @click="tabClick(idx)"
        >
          <div class="tab__name">
            <img class="tab__icon" :src="tab.img" />
            <p class="tab__text">{{ tab.name }}</p>
          </div>
          <img class="tab__close" @click="closeTab" src="@/assets/close.png" />
        </li>
      </ul>
      <div class="content-view" v-if="tabList.length != 0 && isTabView">
        <fileuploadModal
          v-show="tabList[currentTab].name === 'File'"
        ></fileuploadModal>
        <dataTableModal
          v-show="tabList[currentTab].name === 'DataTable'"
        ></dataTableModal>
        <scatterPlotModal
          v-show="tabList[currentTab].name === 'scatterPlot'"
        ></scatterPlotModal>
        <heatMapModal
          v-show="tabList[currentTab].name === 'heatMap'"
        ></heatMapModal>
        <!-- <div class="content__handle" @click="hideContent">
          <img
            class="handle--img"
            src="@/assets/arrow-right.png"
            draggable="false"
          />
        </div> -->
      </div>
    </VueDragResize>
  </div>
</template>

<script>
import Vue from "vue";
import VueDragResize from "vue-drag-resize";
/* eslint-disable */
// import Drawflow from 'drawflow'
// import styleDrawflow from 'drawflow/dist/drawflow.min.css' // eslint-disable-line no-use-before-define

//노드 import (3번)
import scatterPlot from "@/components/nodes/scatterPlotNode.vue";
import fileUpload from "@/components/nodes/fileUploadNode.vue";
import dataTable from "@/components/nodes/dataTableNode.vue";
import heatMap from "@/components/nodes/heatMapNode.vue";

//노드 모달 컴포넌트 import (4번)
import dataTableModal from "@/components/modals/datatable.vue";
import fileuploadModal from "@/components/modals/fileupload.vue";
import scatterPlotModal from "@/components/modals/scatterPlot.vue";
import heatMapModal from "@/components/modals/heatMap.vue";
import Fileupload from "../components/modals/fileupload.vue";
import {
  exportData,
  findWorkflow,
  saveWorkflow,
  userTaskMonitoring,
} from "@/api/index";

export default {
  components: {
    dataTableModal,
    fileuploadModal,
    scatterPlotModal,
    Fileupload,
    heatMapModal,
    VueDragResize,
  },
  data() {
    return {
      editor: null,
      exportValue: null,
      isTabView: true,
      // 왼쪽에 보여지는 노드 목록 (1번)
      listNodes: [
        {
          name: "File",
          // name2: "File",
          img: require("@/assets/file-upload2.png"),
          input: 0,
          output: 1,
        },
        {
          name: "DataTable",
          // name2: "DataTable",
          img: require("@/assets/table2.png"),
          input: 1,
          output: 1,
        },
        {
          name: "scatterPlot",
          // name2: "Plot",
          img: require("@/assets/scatter-plot2.png"),
          input: 1,
          output: 0,
        },
        {
          name: "heatMap",
          img: require("@/assets/heatMap2.png"),
          input: 1,
          output: 0,
        },
        {
          name: "Algorithm",
          // name2: "Algorithm",
          img: require("@/assets/algorithm2.png"),
          input: 1,
          output: 0,
        },
      ],
      tabList: [],
      is_show_modal: false,
      is_show_info: false,
      show_modal: null,
      compile_check: null,
      node_info: {
        name: null,
        desc: null,
        input: null,
        output: null,
        content: null,
      },
      node_connection: [],
      selected_file: null,
      file_name: null,
      currentTab: 0,
      isHide: false,
      show_files: false,
      show_jobs: false,
      on_progress: false, // 나중에 false로 바꾸기
      eventSources: {}, // Use an object to manage multiple event sources
      taskList: [],
      taskTitleList: [],
      currentTime: new Date(),
      timeInterval: null,
    };
  },
  async mounted() {
    const id = document.getElementById("drawflow");
    Vue.prototype.$df = new Drawflow(id, Vue, this);
    //this.$df == editor
    this.$df.start();

    //노드 등록 (2번)
    this.$df.registerNode("File", fileUpload, {}, {});
    this.$df.registerNode("DataTable", dataTable, {}, {});
    this.$df.registerNode("scatterPlot", scatterPlot, {}, {});
    this.$df.registerNode("heatMap", heatMap, {}, {});

    // 노드 수직 연결선
    this.$df.curvature = 0.5;
    this.$df.reroute_curvature_start_end = 0;
    this.$df.reroute_curvature = 0;
    this.$df.on("nodeCreated", (ev) => {
      const node = this.$df.getNodeFromId(ev);
      this.tabList.push({
        id: node.id,
        name: node.name,
        img: require(`@/assets/${node.name}.png`),
      });
      console.log(node);
      this.$store.commit("createNode", {
        id: node.id,
        name: node.name,
        file: "",
      });
      this.$store.commit("changeNode", node.id);
      console.log(this.currentTab);
      if (this.tabList.length != 1) {
        this.currentTab = this.tabList.length - 1;
      }
      this.exportValue = this.$df.export();
      this.$store.commit("setWorkflow", this.exportValue);
    });
    this.$df.on("nodeRemoved", (ev) => {
      this.tabList.forEach((ele, idx) => {
        if (ele.id === parseInt(ev)) {
          this.tabList.splice(idx, 1);
          //현재 탭이 지워지는 탭일 때
          if (this.currentTab === idx) {
            this.currentTab -= 1;
          } else if (this.currentTab > idx) {
            this.currentTab -= 1;
          } else if (this.currentTab < idx) {
            this.currentTab = this.currentTab;
          }
        }
      });
      this.$store.commit("deleteNode", {
        id: parseInt(ev),
      });
      this.exportValue = this.$df.export();
      this.$store.commit("setWorkflow", this.exportValue);
    });
    this.$df.on("connectionCreated", (ev) => {
      // ev 값에 따라 기능 구분
      console.log(ev);
      // const input_id = this.$df.getNodeFromId(ev.input_id);
      // const output_id = this.$df.getNodeFromId(ev.output_id);
      const lastNodeInfo = this.$store.getters.getNodeInfo(
        parseInt(ev.input_id)
      );
      // console.log(lastNodeInfo.name);
      this.$store.commit("createConnection", {
        connection: [parseInt(ev.output_id), parseInt(ev.input_id)],
        file: "",
        lastNode: lastNodeInfo.name,
      });
      this.$store.commit("shareConnectionFile");
      this.exportValue = this.$df.export();
      this.$store.commit("setWorkflow", this.exportValue);
    });
    this.$df.on("connectionRemoved", (ev) => {
      // ev 값에 따라 기능 구분
      console.log(ev);
      // const input_id = this.$df.getNodeFromId(ev.input_id);
      // const output_id = this.$df.getNodeFromId(ev.output_id);
      this.$store.commit("deleteConnection", [
        parseInt(ev.output_id),
        parseInt(ev.input_id),
      ]);
      this.exportValue = this.$df.export();
      this.$store.commit("setWorkflow", this.exportValue);
    });
    this.$df.on("nodeDataChanged", (ev) => {
      // nodeData 바뀌게 되면 Connection Update
      // console.log(ev)
      const node = this.$df.getNodeFromId(ev);
      console.log(node);
      this.$df.updateConnectionNodes(ev);
    });
    this.$df.on("clickEnd", (ev) => {
      // ev 값에 따라 기능 구분
      // console.log(ev);
      if (ev.detail === 2 && this.$df.node_selected) {
        // 해당 노드와 연결되어 있는 File 정보 추출
        this.isTabView = true;

        const node_id = this.$df.node_selected.id.replace(/node-/g, "");
        console.log(typeof node_id);
        //Vuex 에서의 Current Node 변경
        this.$store.commit("changeNode", parseInt(node_id));
        //this.tabList에 추가
        const node = this.$store.getters.getNodeInfo(parseInt(node_id));
        console.log(node);
        this.tabList.push({
          id: node.id,
          name: node.name,
          img: require(`@/assets/${node.name}.png`),
        });
        //this.currentTab 바꾸기
        this.tabList.forEach((ele, idx) => {
          if (ele.id === parseInt(node_id)) {
            this.currentTab = idx;
          }
        });
        console.log(this.tabList);
      }
    });
    try {
      if (this.$route.query.id) {
        const workflowInfo = {
          id: this.$route.query.id,
        };
        const workflow_data = await findWorkflow(workflowInfo);
        console.log(workflow_data.data);
        this.$df.import(workflow_data.data.workflow_info);
        this.$store.commit("setNodes", workflow_data.data.nodes);
        this.$store.commit("setLinkedNodes", workflow_data.data.linked_nodes);
        this.$store.commit("setTitle", workflow_data.data.title);
      }
    } catch (error) {
      console.error(error);
    }

    const storeItem = localStorage.getItem("vuex");
    if (storeItem) {
      const workflow_data = JSON.parse(storeItem);
      console.log(workflow_data);
      if (workflow_data.workflow.workflow_info) {
        this.$df.import(workflow_data.workflow.workflow_info);
        this.$store.commit("setNodes", workflow_data.workflow.nodes);
        this.$store.commit(
          "setLinkedNodes",
          workflow_data.workflow.linked_nodes
        );
        this.$store.commit("setTitle", workflow_data.workflow.title);
      }
    }
  },
  methods: {
    // Define a function that returns a task_id and creates a new EventSource instance
    createEventSource(task_id) {
      this.on_progress = true;
      // Initialize a new event source with the provided URL
      this.eventSources[task_id] = new EventSource(
        `http://127.0.0.1:8000/routes/workflow/task/${task_id}`
      );

      // Event handler for the 'onmessage' event
      this.eventSources[task_id].onmessage = (event) => {
        // Log the received event data
        console.log("Received update: ", event.data);
        if (event.data === "SUCCESS" || event.data === "FAILURE") {
          console.log("close");
          this.on_progress = false;
          this.closeEventSource(task_id);
        }
      };
    },
    // Define a function to close a specific event source
    closeEventSource(task_id) {
      // If the event source exists, close it
      if (this.eventSources[task_id]) {
        this.eventSources[task_id].close();
        delete this.eventSources[task_id];
      }
    },
    async exportdf() {
      try {
        this.exportValue = this.$df.export();
        const nodes = this.$store.getters.getNodes;
        const linked_nodes = this.$store.getters.getLinkedNodes;
        const title = this.$store.getters.getTitle;
        // console.log(JSON.stringify(this.exportValue));
        console.log(this.$df.drawflow.drawflow[this.$df.module]);
        const workflow = {
          id: this.$route.query.id,
          title: title,
          workflow_info: this.exportValue,
          nodes: nodes,
          linked_nodes: linked_nodes,
        };
        console.log(workflow);
        // this.compile_check = "loading";
        const workflow_data = await exportData(workflow);
        this.createEventSource(workflow_data.data.task_id);
        if (this.show_jobs) {
          this.show_jobs = false;
          this.toggleTask();
        }
      } catch (error) {
        console.error(error);
      }
    },
    handleCompletedWorkflow(workflowData) {
      console.log("Workflow completed!", workflowData);
    },
    importdf() {
      this.$df.import(this.exportValue);
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
      console.log(nodeSelected);
      this.$df.addNode(
        name,
        nodeSelected.input,
        nodeSelected.output,
        pos_x,
        pos_y,
        name,
        {},
        name,
        "vue"
      );
    },
    tabClick(idx) {
      this.currentTab = idx;
      this.$store.commit("changeNode", this.tabList[idx].id);
    },
    hideContent(event) {
      this.isHide = !this.isHide;
      console.log(event);
    },
    moveResizing(event) {
      console.log(event);
    },
    closeTab(event) {
      setTimeout(() => {
        console.log(this.tabList);
        const currentNodeId = this.$store.getters.getCurrentNode;
        this.tabList.splice(currentNodeId - 1, 1);
        this.currentTab -= 1;
        console.log(this.tabList);
      }, "100");
    },
    async saveWorkflowProject() {
      try {
        this.exportValue = this.$df.export();
        const nodes = this.$store.getters.getNodes;
        const linked_nodes = this.$store.getters.getLinkedNodes;
        const title = this.$store.getters.getTitle;
        console.log(title);
        const workflow = {
          id: this.$route.query.id,
          title: title,
          workflow_info: this.exportValue,
          nodes: nodes,
          linked_nodes: linked_nodes,
        };
        console.log(workflow);
        const workflow_data = await saveWorkflow(workflow);
        console.log(workflow_data);
      } catch (error) {
        console.error(error);
      }
    },
    async toggleTask() {
      if (!this.show_jobs) {
        const user_tasks = await userTaskMonitoring();
        console.log(user_tasks);
        this.taskList = user_tasks.data;
        this.taskList.forEach(async (task, idx) => {
          if (task.status === "SUCCESS" || task.status === "FAILURE") {
            this.taskList[idx].running_time = this.getTimeDifference(
              task.start_time,
              task.end_time
            );
          } else {
            this.timeInterval = this.startTimer(idx);
          }
          const workflow = await findWorkflow({
            id: task.workflow_id,
          });
          this.taskList[idx].title = workflow.data.title;
        });
      }
      //stop interval
      else {
        clearInterval(this.timeInterval);
      }
      console.log(this.taskList);
      setTimeout(() => {
        this.show_jobs = !this.show_jobs;
      }, 100);
    },
    startTimer(idx) {
      const interval = setInterval(() => {
        if (!this.on_progress) {
          this.show_jobs = false;
          clearInterval(interval);
          this.toggleTask();
        }
        let currentTime = new Date();
        let running_time = this.getRunningTime(
          this.taskList[idx].start_time,
          currentTime
        );
        // this.taskList[idx].running_time = running_time;
        this.$set(this.taskList, idx, {
          ...this.taskList[idx],
          running_time: running_time,
        });
      }, 1000);
      return interval;
    },
    getTimeDifference(start_time, end_time) {
      const start = new Date(start_time);
      const end = new Date(end_time);
      let diff = Math.abs(end - start); // Difference in milliseconds

      const hours = Math.floor(diff / 3600000);
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

      const hours = Math.floor(diff / 3600000);
      diff -= hours * 3600000;

      const minutes = Math.floor(diff / 60000);
      diff -= minutes * 60000;

      const seconds = Math.floor(diff / 1000);

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
  },
  beforeDestroy() {
    // Close all the event source connections before the component is destroyed
    for (let task_id in this.eventSources) {
      this.closeEventSource(task_id);
    }
  },
};
</script>

<style>
@keyframes rotate-loading {
  0% {
    transform: rotate(0deg);
    -ms-transform: rotate(0deg);
    -webkit-transform: rotate(0deg);
    -o-transform: rotate(0deg);
    -moz-transform: rotate(0deg);
  }
  100% {
    transform: rotate(360deg);
    -ms-transform: rotate(360deg);
    -webkit-transform: rotate(360deg);
    -o-transform: rotate(360deg);
    -moz-transform: rotate(360deg);
  }
}

@keyframes loading-text-opacity {
  0% {
    opacity: 0;
  }
  20% {
    opacity: 0;
  }
  50% {
    opacity: 1;
  }
  100% {
    opacity: 0;
  }
}

@keyframes complete-text-opacity {
  0% {
    opacity: 0;
  }
  5% {
    opacity: 1;
  }
  95% {
    opacity: 1;
  }
  100% {
    opacity: 0;
  }
}

.layout__workflow {
  width: 100%;
  height: 100%;
  position: relative;
  overflow: hidden;
  background-image: url("@/assets/fantastic_background2.png");
  background-size: cover;
  background-position: center center;
  background-repeat: no-repeat;
}
.main__bg-video {
  position: fixed;
  width: 100%;
  height: 100%;
  overflow: hidden;
  z-index: -1;
}
.main-video {
  position: absolute;
  top: 0;
  right: 17px;
  width: 100%;
  /* max-width: 98vw; */
  height: 120vh;
  object-fit: cover;
  /* width: 100rem; */
  /* max-width: calc(99%); */
  /* height: 20rem; */
  /* object-fit: cover;
  object-position: left; */
  /* border-radius: 1.5rem; */
}
.content-component {
  width: 55rem;
  height: 42rem;
  position: absolute;
  right: -55rem;
  top: calc(50% - 21rem);
}
.tab_actvie {
  right: 0;
}
.tab_hide {
  left: 100%;
}
.content-tab {
  width: 100%;
  height: 2.5rem;
  display: flex;
  z-index: 9998;
  background: rgba(223, 225, 229, 0.3);
  /* border-radius: 0.5rem 0.5rem 0 0; */
}
.tab__item {
  cursor: pointer;
  width: 10rem;
  height: 2.2rem;
  top: 0.3rem;
  border-radius: 0.5rem 0.5rem 0 0;
  border-right: 1px solid #7f7f7f;
  display: flex;
  align-items: center;
  background: rgba(149, 151, 154, 0.6);
  color: rgb(255, 255, 255);
  position: relative;
  opacity: 1;
  box-shadow: inset 0 -5px 10px -5px rgba(0, 0, 0, 0.3);
}
.tab__item:last-child {
  border-right: none;
}
.currentTab {
  background: rgba(244, 246, 251, 0.5);
  color: rgb(51, 51, 51);
  box-shadow: inset 0 -5px 10px -5px rgba(0, 0, 0, 0);
}
.tab__name {
  width: 8rem;
  height: 100%;
  display: flex;
  align-items: center;
  position: absolute;
  left: 0.5rem;
}
.tab__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 0.8rem;
  line-height: 1rem;
  display: flex;
  align-items: center;
  justify-content: flex-start;
  overflow: hidden;
}
.tab__icon {
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  margin-right: 0.5rem;
}
.tab__close {
  width: 0.7rem;
  height: 0.7rem;
  object-fit: contain;
  position: absolute;
  right: 1rem;
}
.content-view {
  width: 100%;
  height: calc(100% - 2rem);
  background: rgba(244, 248, 251, 0.6);
  border-radius: 0 0 0.5rem 0.5rem;
}
.content__handle {
  position: absolute;
  left: -0.75rem;
  top: calc(50% - 2rem);

  cursor: pointer;
  width: 1.5rem;
  height: 2rem;
  border-radius: 3px;
  box-shadow: rgba(6, 24, 44, 0.4) 0px 0px 0px 2px,
    rgba(6, 24, 44, 0.65) 0px 4px 6px -1px,
    rgba(255, 255, 255, 0.08) 0px 1px 0px inset;
  background: rgb(255, 255, 255);
  z-index: 9998;

  display: flex;
  align-items: center;
}
.handle--img {
  width: 1.5rem;
  height: 1.5rem;
  object-fit: contain;
}
.isResizing {
  right: -55rem;
}
.node-bar {
  /* width: 8rem; */
  width: 80px;
  /* height: 34rem; */
  height: 400px;
  border-radius: 16px;
  background: rgba(244, 246, 251, 0.586);
  box-shadow: 0px 0px 5px 0px rgba(0, 0, 0, 1);
  position: absolute;
  /* top: calc(50% - 17rem); */
  top: calc(50% - 208px);
  left: 12px;
  z-index: 9998;
  opacity: 1;
  display: flex;
  align-items: center;
  justify-content: center;
}
.node-bar__nodelist {
  width: 80%;
  height: 100%;
  display: flex;
  margin-left: 3px;
  flex-direction: column;
  justify-content: space-evenly;
  align-items: center;
}
.node-bar__drag-drawflow {
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: move;
  background: rgba(15, 19, 70, 0.868);
  color: rgb(245, 245, 245);
  border-radius: 1rem;
  /* width: 5rem;
  height: 5rem; */
  width: 4rem;
  height: 4rem;
  box-shadow: 0px 0px 6px 0px rgb(99, 99, 99);
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
.control-popup__files,
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
.control-popup__files {
  right: calc(50% + 1vw);
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
.control-popup__table td {
  vertical-align: middle;
  font-weight: 400;
  text-align: center;
  color: rgb(68, 68, 68);
  padding: 0.7rem;
  margin: 1rem;
}
.control-popup__jobs {
  left: calc(50% + 1vw);
  overflow-y: scroll;
  border-radius: 16px; /* or whatever radius you prefer */
}

.control-popup__jobs::-webkit-scrollbar {
  width: 10px; /* width of the entire scrollbar */
}

.control-popup__jobs::-webkit-scrollbar-track {
  background: #f1f1f1; /* color of the tracking area */
  border-radius: 16px; /* keep the same radius as the container */
}

.control-popup__jobs::-webkit-scrollbar-thumb {
  background: #888; /* color of the scroll thumb */
  border-radius: 16px; /* keep the same radius as the container */
}

.control-popup__jobs::-webkit-scrollbar-thumb:hover {
  background: #555; /* color of the scroll thumb on hover */
}

.control-popup__table__progress {
  width: 40%;
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

.control-bar {
  height: 50px;
  width: 300px;
  border-radius: 10px;
  background: rgba(0, 0, 0, 0.776);
  position: absolute;
  bottom: 24px;
  left: calc(50% - 150px);
  display: flex;
  align-items: center;
  justify-content: center;
}
.control-bar__btnList {
  width: 100%;
  height: 100%;
  display: flex;
  justify-content: center;
  align-items: center;
}
.control-bar__button {
  width: 24px;
  height: 24px;
  padding: 8px;
}

.control-bar__icon {
  max-width: 24px;
  max-height: 24px;
  object-fit: cover;
  opacity: 0.72;
}
.loader,
.loader_done {
  border: 4px solid #f3f3f3bf;
  border-radius: 50%;
  margin-left: 8px;
  margin-right: 6px;
  width: 20px;
  height: 20px;
}
.loader {
  border-top: 4px solid #41b3ff;
  animation: spin 3s linear infinite;
}

@keyframes spin {
  0% {
    transform: rotate(0deg);
  }
  100% {
    transform: rotate(360deg);
  }
}

.run_button {
  width: 100%;
  height: 100%;
  background: none;
  border: none;
  cursor: pointer;
}
/* .run_button__icon {
  width: 1.5rem;
  height: 1.5rem;
  object-fit: contain;
  filter: invert(100%) sepia(3%) saturate(2008%) hue-rotate(348deg)
    brightness(125%) contrast(111%);
} */
#drawflow {
  width: calc(100% - 180px);
  height: calc(100% - 50px);
  top: 30px;
  left: 101px;
  position: relative;
  border-radius: 0.8%;

  background: var(--dfBackgroundColor);
  background-size: var(--dfBackgroundSize) var(--dfBackgroundSize);
  background-image: var(--dfBackgroundImage);
}

.drawflow .drawflow-node {
  display: var(--dfNodeType);
  background: var(--dfNodeBackgroundColor);
  color: var(--dfNodeTextColor);
  border: var(--dfNodeBorderSize) solid var(--dfNodeBorderColor);
  border-radius: var(--dfNodeBorderRadius);
  min-height: var(--dfNodeMinHeight);
  width: auto;
  min-width: var(--dfNodeMinWidth);
  padding-top: var(--dfNodePaddingTop);
  padding-bottom: var(--dfNodePaddingBottom);
  -webkit-box-shadow: var(--dfNodeBoxShadowHL) var(--dfNodeBoxShadowVL)
    var(--dfNodeBoxShadowBR) var(--dfNodeBoxShadowS) var(--dfNodeBoxShadowColor);
  box-shadow: var(--dfNodeBoxShadowHL) var(--dfNodeBoxShadowVL)
    var(--dfNodeBoxShadowBR) var(--dfNodeBoxShadowS) var(--dfNodeBoxShadowColor);
}

.drawflow .drawflow-node:hover {
  background: var(--dfNodeHoverBackgroundColor);
  color: var(--dfNodeHoverTextColor);
  border: var(--dfNodeHoverBorderSize) solid var(--dfNodeHoverBorderColor);
  border-radius: var(--dfNodeHoverBorderRadius);
  -webkit-box-shadow: var(--dfNodeHoverBoxShadowHL)
    var(--dfNodeHoverBoxShadowVL) var(--dfNodeHoverBoxShadowBR)
    var(--dfNodeHoverBoxShadowS) var(--dfNodeHoverBoxShadowColor);
  box-shadow: var(--dfNodeHoverBoxShadowHL) var(--dfNodeHoverBoxShadowVL)
    var(--dfNodeHoverBoxShadowBR) var(--dfNodeHoverBoxShadowS)
    var(--dfNodeHoverBoxShadowColor);
}

.drawflow .drawflow-node.selected {
  background: var(--dfNodeSelectedBackgroundColor);
  color: var(--dfNodeSelectedTextColor);
  border: var(--dfNodeSelectedBorderSize) solid var(--dfNodeSelectedBorderColor);
  border-radius: var(--dfNodeSelectedBorderRadius);
  -webkit-box-shadow: var(--dfNodeSelectedBoxShadowHL)
    var(--dfNodeSelectedBoxShadowVL) var(--dfNodeSelectedBoxShadowBR)
    var(--dfNodeSelectedBoxShadowS) var(--dfNodeSelectedBoxShadowColor);
  box-shadow: var(--dfNodeSelectedBoxShadowHL) var(--dfNodeSelectedBoxShadowVL)
    var(--dfNodeSelectedBoxShadowBR) var(--dfNodeSelectedBoxShadowS)
    var(--dfNodeSelectedBoxShadowColor);
}

.drawflow .drawflow-node .input {
  left: var(--dfInputLeft);
  background: var(--dfInputBackgroundColor);
  border: var(--dfInputBorderSize) solid var(--dfInputBorderColor);
  border-radius: var(--dfInputBorderRadius);
  height: var(--dfInputHeight);
  width: var(--dfInputWidth);
}

.drawflow .drawflow-node .input:hover {
  background: var(--dfInputHoverBackgroundColor);
  border: var(--dfInputHoverBorderSize) solid var(--dfInputHoverBorderColor);
  border-radius: var(--dfInputHoverBorderRadius);
}

.drawflow .drawflow-node .outputs {
  float: var(--dfNodeTypeFloat);
}

.drawflow .drawflow-node .output {
  right: var(--dfOutputRight);
  background: var(--dfOutputBackgroundColor);
  border: var(--dfOutputBorderSize) solid var(--dfOutputBorderColor);
  border-radius: var(--dfOutputBorderRadius);
  height: var(--dfOutputHeight);
  width: var(--dfOutputWidth);
}

.drawflow .drawflow-node .output:hover {
  background: var(--dfOutputHoverBackgroundColor);
  border: var(--dfOutputHoverBorderSize) solid var(--dfOutputHoverBorderColor);
  border-radius: var(--dfOutputHoverBorderRadius);
}

.drawflow .connection .main-path {
  stroke-width: var(--dfLineWidth);
  stroke: var(--dfLineColor);
}

.drawflow .connection .main-path:hover {
  stroke: var(--dfLineHoverColor);
}

.drawflow .connection .main-path.selected {
  stroke: var(--dfLineSelectedColor);
}

.drawflow .connection .point {
  stroke: var(--dfRerouteBorderColor);
  stroke-width: var(--dfRerouteBorderWidth);
  fill: var(--dfRerouteBackgroundColor);
}

.drawflow .connection .point:hover {
  stroke: var(--dfRerouteHoverBorderColor);
  stroke-width: var(--dfRerouteHoverBorderWidth);
  fill: var(--dfRerouteHoverBackgroundColor);
}

.drawflow-delete {
  display: var(--dfDeleteDisplay);
  color: var(--dfDeleteColor);
  background: var(--dfDeleteBackgroundColor);
  border: var(--dfDeleteBorderSize) solid var(--dfDeleteBorderColor);
  border-radius: var(--dfDeleteBorderRadius);
}

.parent-node .drawflow-delete {
  top: var(--dfDeleteTop);
}

.drawflow-delete:hover {
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
  outline: none;
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
