<template>
  <div class="layout__workflow">
    <div id="drawflow" @drop="drop($event)" @dragover="allowDrop($event)">
      <button class="run_button" @click="exportdf">Run</button>
    </div>
    <section class="node-bar">
      <ul class="node-bar__nodelist">
        <li
          class="node-bar__drag-drawflow"
          v-for="(node, idx) in listNodes"
          :key="idx"
          draggable="true"
          :data-node="node.name"
          @dragstart="drag($event)"
        >
          <img class="node-bar__img" :src="node.img" />
        </li>
      </ul>
      <!-- <div class="node-bar__row">
          <div class="loading_bg">
            <h2 v-if="compile_check == 'loading'" class="loading_container">
              <div class="loading"></div>
              <div class="loading_text">Loading...</div>
            </h2>
            <div class="complete_container">
              <h2 v-if="compile_check == 'complete'" class="complete">√</h2>
            </div>
          </div>
        </div> -->
    </section>
    <main
      class="content-component"
      v-bind:class="{ tab_actvie: tabList.length != 0 }"
    >
      <ul class="content-tab">
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
        </li>
      </ul>
      <div class="content-view" v-if="tabList.length != 0">
        <fileuploadModal
          v-show="tabList[currentTab].name === 'File'"
        ></fileuploadModal>
        <dataTableModal
          v-show="tabList[currentTab].name === 'DataTable'"
        ></dataTableModal>
        <PlotModal v-show="tabList[currentTab].name === 'Plot'"></PlotModal>
      </div>
    </main>
  </div>
</template>

<script>
import Vue from "vue";
/* eslint-disable */
// import Drawflow from 'drawflow'
// import styleDrawflow from 'drawflow/dist/drawflow.min.css' // eslint-disable-line no-use-before-define
import Plot from "@/components/nodes/PlotNode.vue";
import fileUpload from "@/components/nodes/fileUploadNode.vue";
import dataTable from "@/components/nodes/dataTableNode.vue";
import dataTableModal from "@/components/modals/datatable.vue";
import fileuploadModal from "@/components/modals/fileupload.vue";
import PlotModal from "@/components/modals/Plot.vue";
import scatterPlotModal from "@/components/modals/scatterPlot.vue";

import { exportData, getCheckCompile } from "@/api/index";
import Fileupload from "../components/modals/fileupload.vue";

export default {
  components: {
    dataTableModal,
    fileuploadModal,
    PlotModal,
    scatterPlotModal,
    Fileupload,
  },
  data() {
    return {
      rightSidebar_isActive: true,
      editor: null,
      exportValue: null,
      listNodes: [
        {
          name: "File",
          name2: "File",
          img: require("@/assets/file-upload.png"),
          input: 0,
          output: 1,
        },
        {
          name: "DataTable",
          name2: "DataTable",
          img: require("@/assets/table.png"),
          input: 1,
          output: 1,
        },
        {
          name: "Plot",
          name2: "Plot",
          img: require("@/assets/scatter-plot.png"),
          input: 1,
          output: 0,
        },
        {
          name: "Algorithm",
          name2: "Algorithm",
          img: require("@/assets/algorithm.png"),
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
  mounted() {
    const id = document.getElementById("drawflow");
    Vue.prototype.$df = new Drawflow(id, Vue, this);
    this.$df.start();

    this.$df.registerNode("Plot", Plot, {}, {});
    this.$df.registerNode("File", fileUpload, {}, {});
    this.$df.registerNode("DataTable", dataTable, {}, {});
    // 노드 수직 연결선
    this.$df.curvature = 0.5;
    this.$df.reroute_curvature_start_end = 0;
    this.$df.reroute_curvature = 0;
    //노드 연결선 화살표 추가
    // this.$df.createCurvature = function(start_pos_x, start_pos_y, end_pos_x, end_pos_y, curvature_value, type) {
    //   var line_x = start_pos_x;
    //   var line_y = start_pos_y;
    //   var x = end_pos_x;
    //   var y = end_pos_y;
    //   var curvature = curvature_value;

    //   //type openclose open close other
    //   switch (type) {
    //     case 'open':
    //       if(start_pos_x >= end_pos_x) {
    //         var hx1 = line_x + Math.abs(x - line_x) * curvature;
    //         var hx2 = x - Math.abs(x - line_x) * (curvature*-1);
    //       } else {
    //         var hx1 = line_x + Math.abs(x - line_x) * curvature;
    //         var hx2 = x - Math.abs(x - line_x) * curvature;
    //       }
    //       return ' M '+ line_x +' '+ line_y +' C '+ hx1 +' '+ line_y +' '+ hx2 +' ' + y +' ' + x +'  ' + y;

    //       break
    //     case 'close':
    //       if(start_pos_x >= end_pos_x) {
    //         var hx1 = line_x + Math.abs(x - line_x) * (curvature*-1);
    //         var hx2 = x - Math.abs(x - line_x) * curvature;
    //       } else {
    //         var hx1 = line_x + Math.abs(x - line_x) * curvature;
    //         var hx2 = x - Math.abs(x - line_x) * curvature;
    //       }                                                                                                                  //M0 75H10L5 80L0 75Z
    //       return ' M '+ line_x +' '+ line_y +' C '+ hx1 +' '+ line_y +' '+ hx2 +' ' + y +' ' + x +'  ' + y +' M '+ (x-11)  + ' ' + y + ' L'+(x-20)+' '+ (y-5)+'  L'+(x-20)+' '+ (y+5)+' Z' +' M '+ (x-11)  + ' ' + y + ' L'+(x-20)+' '+ (y-3)+'  L'+(x-20)+' '+ (y+3)+' Z' +' M '+ (x-11)  + ' ' + y + ' L'+(x-20)+' '+ (y-1)+'  L'+(x-20)+' '+ (y+1)+' Z';
    //       // return ' M '+ line_x +' '+ line_y +' C '+ hx1 +' '+ line_y +' '+ hx2 +' ' + y +' ' + x +'  ' + y +' M '+ (x-11)  + ' ' + y + ' L'+(x-20)+' '+ (y-5)+'  L'+(x-20)+' '+ (y+5)+'Z';
    //       break;
    //     case 'other':
    //       if(start_pos_x >= end_pos_x) {
    //         var hx1 = line_x + Math.abs(x - line_x) * (curvature*-1);
    //         var hx2 = x - Math.abs(x - line_x) * (curvature*-1);
    //       } else {
    //         var hx1 = line_x + Math.abs(x - line_x) * curvature;
    //         var hx2 = x - Math.abs(x - line_x) * curvature;
    //       }
    //       return ' M '+ line_x +' '+ line_y +' C '+ hx1 +' '+ line_y +' '+ hx2 +' ' + y +' ' + x +'  ' + y;
    //       break;
    //     default:

    //       var hx1 = line_x + Math.abs(x - line_x) * curvature;
    //       var hx2 = x - Math.abs(x - line_x) * curvature;

    //       //return ' M '+ line_x +' '+ line_y +' C '+ hx1 +' '+ line_y +' '+ hx2 +' ' + y +' ' + x +'  ' + y;
    //       return ' M '+ line_x +' '+ line_y +' C '+ hx1 +' '+ line_y +' '+ hx2 +' ' + y +' ' + x +'  ' + y +' M '+ (x-11)  + ' ' + y + ' L'+(x-20)+' '+ (y-5)+'  L'+(x-20)+' '+ (y+5)+'Z';
    //   }
    // }
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
      const lastNodeInfo = this.$store.getters.getNodeInfo(parseInt(ev.input_id));
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
    this.$df.on("nodeSelected", (ev) => {
      // ev 값에 따라 기능 구분
      const node = this.$df.getNodeFromId(ev);
      console.log(node.inputs, node.outputs);
      // this.node_info.name = node.name
      //Plot 임시 코드
      if (node.name != "Plot") {
        this.node_info.name = node.name;
      } else {
        this.node_info.name = "Plot";
      }

      if (node.name == "File") {
        this.node_info.desc = "Read data from an input file";
      } else if (node.name == "DataTable") {
        this.node_info.desc = "View the dataset in a spreadsheet";
      } else if (node.name == "Plot") {
        // this.node_info.desc = 'Interactive Plot visualization'
        this.node_info.desc = "Interactive plot visualization";
      }
      if (this.connectionParsing(node.inputs)) {
        const input_node = this.$df.getNodeFromId(
          this.connectionParsing(node.inputs)
        );
        console.log(input_node.name);
        this.node_info.input = `${input_node.name} -> `;
      } else {
        this.node_info.input = null;
      }
      if (this.connectionParsing(node.outputs)) {
        const output_node = this.$df.getNodeFromId(
          this.connectionParsing(node.outputs)
        );
        console.log(output_node.name);
        // this.node_info.output = ` -> ${output_node.name}`
        //Plot 임시 코드
        if (output_node.name != "Plot") {
          this.node_info.output = ` -> ${output_node.name}`;
        } else {
          this.node_info.output = " -> Plot";
        }
      } else {
        this.node_info.output = null;
      }
    });
    this.$df.on("nodeUnselected", (ev) => {
      this.is_show_info = false;
    });
    this.$df.on("clickEnd", (ev) => {
      // ev 값에 따라 기능 구분
      // console.log(ev);
      if (ev.detail === 2 && this.$df.node_selected) {
        // 해당 노드와 연결되어 있는 File 정보 추출
        const node_id = this.$df.node_selected.id.replace(/node-/g, "");
        console.log(node_id);
        this.$store.commit("changeNode", parseInt(node_id));
        this.tabList.forEach((ele, idx) => {
          if (ele.id === parseInt(node_id)) {
            this.currentTab = idx;
          }
        });
        this.node_connection.forEach((connection) => {
          // console.log(connection)
          connection.forEach((node) => {
            // console.log(node)
            if (node == node_id) {
              const node_info = this.$df.getNodeFromId(connection[0]);
              const file_name = node_info.data.file
                .replace(/C:\\fakepath\\/, "")
                .replace(/.csv/, "");
              this.file_name = file_name;
              console.log(connection[0], file_name);
            }
          });
        });
        // 해당 노드 File 노드 아니면 modal 보여줌
        console.dir(this.$df.node_selected.innerText.replace(/(\s*)/g, ""));
        if (
          (this.$df.node_selected.innerText.replace(/(\s*)/g, "") ==
            "DataTable") |
          (this.$df.node_selected.innerText.replace(/(\s*)/g, "") == "Plot")
        ) {
          this.is_show_modal = true;
          this.show_modal = this.$df.node_selected.innerText.replace(
            /(\s*)/g,
            ""
          );
        }
      }
    });
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
        const result = await this.$store.dispatch("compileNodes");
        console.log(result);
        // this.exportValue = this.$df.export();
        // this.compile_check = "loading";
        // console.log(this.exportValue);
        // const JsonData = await exportData(
        //   JSON.stringify(this.exportValue.drawflow.Home.data)
        // );
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
}
.content-component {
  width: 55rem;
  height: 42rem;
  position: absolute;
  right: -55rem;
  top: calc(50% - 21rem);
  transition: all 0.8s;
}
.tab_actvie {
  right: 0;
}
.content-tab {
  width: 100%;
  height: 2rem;
  display: flex;
}
.tab__item {
  cursor: pointer;
  width: 10rem;
  height: 100%;
  border-radius: 0.5rem 0.5rem 0 0;
  display: flex;
  align-items: center;
  background: rgb(32, 33, 36);
  color: rgb(255, 255, 255);
  position: relative;
  opacity: 1;
  box-shadow: 0px -6px 5px 0px rgba(0, 0, 0, 0.5);
}
.currentTab {
  background: rgb(53, 54, 58);
  box-shadow: 0px -6px 5px 0px rgba(0, 0, 0, 0.5);
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
  width: 0.8rem;
  height: 0.8rem;
  object-fit: contain;
  position: absolute;
  right: 1rem;
}
.content-view {
  width: 100%;
  height: 39rem;
  background: rgb(53, 54, 58);
  border-radius: 0 0 0 0.5rem;
  box-shadow: 0px -5px 5px 0px rgba(0, 0, 0, 0.5);
}
.node-bar {
  /* width: 8rem; */
  width: 4.5rem;
  /* height: 34rem; */
  height: 20rem;
  border-radius: 1rem;
  background: rgb(53, 54, 58);
  box-shadow: 0px 0px 5px 0px rgba(0, 0, 0, 1);
  position: absolute;
  /* top: calc(50% - 17rem); */
  top: calc(50% - 13rem);
  left: 1rem;
  z-index: 9998;
  opacity: 1;
  display: flex;
  align-items: center;
  justify-content: center;
}
.node-bar__nodelist {
  width: 80%;
  height: 105%;
  display: flex;
  flex-direction: column;
  justify-content: space-evenly;
  align-items: center;
}
.node-bar__drag-drawflow {
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: move;
  background: rgb(31, 31, 31);
  color: rgb(245, 245, 245);
  border-radius: 1rem;
  /* width: 5rem;
  height: 5rem; */
  width: 4rem;
  height: 4rem;
}
.node-bar__img {
  /* width: 3rem;
  height: 3rem; */
  width: 2.2rem;
  height: 2.2rem;
  object-fit: contain;
  filter: invert(100%) sepia(3%) saturate(2008%) hue-rotate(348deg)
    brightness(125%) contrast(111%);
  -webkit-user-drag: none;
  -khtml-user-drag: none;
  -moz-user-drag: none;
  -o-user-drag: none;
}
.run_button {
  width: 8rem;
  height: 5rem;
  background: rgb(170, 193, 240);
  border-radius: 1rem;
  border: none;

  position: absolute;
  bottom: 2rem;
  left: 1rem;

  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1.3rem;
  line-height: 1.3rem;

  z-index: 9997;
  cursor: pointer;
}
#drawflow {
  width: 100%;
  height: 100%;
  position: relative;

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
</style>
