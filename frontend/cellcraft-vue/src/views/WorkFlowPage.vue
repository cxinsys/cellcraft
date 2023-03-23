<template>
  <div class="layout__workflow">
    <!-- <div class="main__bg-video">
      <video class="main-video" autoplay loop muted>
        <source src="@/assets/main_bg_fantastic.mp4" type="video/mp4" />
      </video>
    </div> -->
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
    <section class="control-bar">
      <ul class="control-bar__btnList">
        <li class="control-bar__button">
          <img class="control-bar__icon" src="@/assets/control_files.png" />
        </li>
        <li class="control-bar__button">
          <img class="control-bar__icon" src="@/assets/control_zoom.png" />
        </li>
        <li class="control-bar__button">
          <img class="control-bar__icon" src="@/assets/control_save.png" />
        </li>
        <li class="control-bar__button">
          <button class="run_button" @click="exportdf">
            <img class="control-bar__icon" src="@/assets/control_run.png" />
          </button>
        </li>
        <li class="control-bar__button">
          <img class="control-bar__icon" src="@/assets/control_jobs.png" />
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
      :sticks="['tl']"
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
        <!-- <div class="content__handle" @mouseup="resizingContent">
          <img class="handle--img" src="@/assets/lines.png" draggable="false" />
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
import { exportData, findWorkflow } from "@/api/index";

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

    const workflowInfo = {
      id: this.$route.query.id,
    };
    const workflow_data = await findWorkflow(workflowInfo);
    console.log(workflow_data.data);
    this.$df.import(workflow_data.data.workflow_info);
    this.$store.commit("setNodes", workflow_data.data.nodes);
    this.$store.commit("setLinkedNodes", workflow_data.data.linked_nodes);
  },
  methods: {
    connectionParsing(IO) {
      if (Object.keys(IO)[0] == null) {
        return null;
      } else {
        if (Object.values(Object.values(Object.values(IO)[0])[0])[0]) {
          const node_id = Object.entries(
            Object.values(Object.values(Object.values(IO)[0])[0])[0]
          )[0][1];
          return node_id;
        }
        return null;
      }
    },
    async exportdf() {
      try {
        //원래 코드
        // const result = await this.$store.dispatch("compileNodes");
        // console.log(result);
        this.exportValue = this.$df.export();
        const nodes = this.$store.getters.getNodes;
        const linked_nodes = this.$store.getters.getLinkedNodes;
        // console.log(JSON.stringify(this.exportValue));
        console.log(this.$df.drawflow.drawflow[this.$df.module]);
        const workflow = {
          title: "Untitled",
          workflow_info: this.exportValue,
          nodes: nodes,
          linked_nodes: linked_nodes,
        };
        console.log(workflow);
        // this.compile_check = "loading";
        const workflow_data = await exportData(workflow);
        console.log(workflow_data);
        // console.log(JSON.stringify(this.exportValue.drawflow.Home.data))
        // this.compile_check = "complete";
        // console.log(
        //   typeof JsonData.data.recived_data,
        //   JsonData.data.recived_data
        // );
        // this.node_connection = JsonData.data.recived_data;
      } catch (error) {
        console.error(error);
      }
    },
    importdf() {
      this.$df.import(this.exportValue);
    },
    drag(event) {
      event.dataTransfer.setData(
        "node",
        event.target.getAttribute("data-node")
      );
      // console.log(event.dataTransfer, event.target);
      // 모바일
      // if (event.type === 'touchstart') {
      //   mobile_item_selec = event.target.closest('node-bar__drag-drawflow').getAttribute('data-node')
      // }
    },
    drop(event) {
      // 모바일
      // if (event.type === 'touchend') {
      //   var parentdrawflow = document.elementFromPoint(mobile_last_move.touches[0].clientX, mobile_last_move.touches[0].clientY).closest('#drawflow')
      //   if (parentdrawflow != null) {
      //     addNodeToDrawFlow(mobile_item_selec, mobile_last_move.touches[0].clientX, mobile_last_move.touches[0].clientY)
      //   }
      //   mobile_item_selec = ''
      // }
      // console.log(event);
      event.preventDefault();
      const data = event.dataTransfer.getData("node");
      this.addNodeToDrawFlow(data, event.clientX, event.clientY);
      // console.log(data);
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
    resizingContent(event) {
      this.isResizing = !this.isResizing;
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

  cursor: col-resize;
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
