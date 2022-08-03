<template>
<div>
  <div id="drawflow" @drop="drop($event)" @dragover="allowDrop($event)"></div>
  <section class="right-sidebar">
    <section class="right-sidebar__main" v-bind:class="{open: rightSidebar_isActive}">
      <div class="right-sidebar__row">
        <ul>
          <li class="right-sidebar__drag-drawflow"  v-for="(node, idx) in listNodes" :key="idx" draggable="true" :data-node="node.name" @dragstart="drag($event)">
            <div class="right-sidebar__node">{{ node.name }}
              <img :src="node.img" alt="" class="right-sidebar__img">
            </div>
          </li>
        </ul>
      </div>
      <div class="right-sidebar__row">
        <button class="right-sidebar__button" @click="exportdf">complie</button>
      </div>

    </section>
    <div class="popBtn" @click="openRightsidebar">
      <div class="popBtn__txt">></div>
    </div>
  </section>

</div>
</template>

<script>
import Vue from 'vue'
/* eslint-disable */
// import Drawflow from 'drawflow'
// import styleDrawflow from 'drawflow/dist/drawflow.min.css' // eslint-disable-line no-use-before-define
import scatterPlot from '@/components/nodes/scatterPlotNode.vue'
import fileUpload from '@/components/nodes/fileUploadNode.vue'
import dataTable from '@/components/nodes/dataTableNode.vue'

import { exportData } from '@/api/index'
import { requireFileAsExpression } from 'webpack/lib/ParserHelpers'


export default {
  data () {
    return {
      rightSidebar_isActive: true,
      editor: null,
      exportValue: null,
      listNodes: [
        {
          name: 'fileUpload',
          img: require('@/assets/file-upload.png'),
          input: 0,
          output: 1
        },
        {
          name: 'dataTable',
          img: require('@/assets/table.png'),
          input: 1,
          output: 1
        },
        {
          name: 'scatterPlot',
          img: require('@/assets/scatter-plot.png'),
          input: 1,
          output: 0
        }
      ],
    }
  },
  mounted () {
    const id = document.getElementById('drawflow')
    Vue.prototype.$df = new Drawflow(id, Vue, this)
    this.$df.start()

    this.$df.registerNode('scatterPlot', scatterPlot, {}, {})
    this.$df.registerNode('fileUpload', fileUpload, {}, {})
    this.$df.registerNode('dataTable', dataTable, {}, {})
    this.$df.on('nodeDataChanged', (ev) => {
      //nodeData 바뀌게 되면 Connection Update
      // console.log(ev)
      const node = this.$df.getNodeFromId(ev)
      console.log(node);
      this.$df.updateConnectionNodes(ev)
    })
    this.$df.on('connectionCreated', (ev) => {
      //ev 값에 따라 기능 구분
      console.log(ev);
    })
  },
  methods: {
    async exportdf () {
      try {
       this.exportValue = this.$df.export()
        const JsonData = await exportData(JSON.stringify(this.exportValue.drawflow.Home.data))
        // console.log(JSON.stringify(this.exportValue.drawflow.Home.data))
        console.log(typeof (JsonData), JsonData) 
      } catch (error) {
        console.error(error)
      }
    },
    importdf () {
      this.$df.import(this.exportValue)
    },
    drag (event) {
      event.dataTransfer.setData('node', event.target.getAttribute('data-node'))
    // console.log(event.dataTransfer, event.target);
    // 모바일
    // if (event.type === 'touchstart') {
    //   mobile_item_selec = event.target.closest('right-sidebar__drag-drawflow').getAttribute('data-node')
    // }
    },
    drop (event) {
      // 모바일
      // if (event.type === 'touchend') {
      //   var parentdrawflow = document.elementFromPoint(mobile_last_move.touches[0].clientX, mobile_last_move.touches[0].clientY).closest('#drawflow')
      //   if (parentdrawflow != null) {
      //     addNodeToDrawFlow(mobile_item_selec, mobile_last_move.touches[0].clientX, mobile_last_move.touches[0].clientY)
      //   }
      //   mobile_item_selec = ''
      // }
      // console.log(event);
      event.preventDefault()
      const data = event.dataTransfer.getData('node')
      this.addNodeToDrawFlow(data, event.clientX, event.clientY)
      // console.log(data);
    },
    allowDrop (event) {
      event.preventDefault()
    },
    addNodeToDrawFlow (name, pos_x, pos_y) {
      pos_x = pos_x * (this.$df.precanvas.clientWidth / (this.$df.precanvas.clientWidth * this.$df.zoom)) - (this.$df.precanvas.getBoundingClientRect().x * (this.$df.precanvas.clientWidth / (this.$df.precanvas.clientWidth * this.$df.zoom)))
      pos_y = pos_y * (this.$df.precanvas.clientHeight / (this.$df.precanvas.clientHeight * this.$df.zoom)) - (this.$df.precanvas.getBoundingClientRect().y * (this.$df.precanvas.clientHeight / (this.$df.precanvas.clientHeight * this.$df.zoom)))

      const nodeSelected = this.listNodes.find(ele => ele.name === name)
      console.log(nodeSelected)
      this.$df.addNode(name, nodeSelected.input, nodeSelected.output, pos_x, pos_y, name, {}, name, 'vue')
    },
    openRightsidebar () {
      this.rightSidebar_isActive = !this.rightSidebar_isActive
    },

  }
}
</script>

<style>
#drawflow {
  width: 100vw;
  height: 94vh;

  background: rgba(0, 0, 0, 1);
  background-size: 30px 30px;
  background-image: radial-gradient(rgba(111, 109, 109, 0.6) 1px, transparent 1px);
}
:root {
  --dfBackgroundColor: rgba(255, 255, 255, 1);
  --dfBackgroundSize: 16px;
  --dfBackgroundImage: radial-gradient(rgba(0, 29, 110, 1) 1px, transparent 1px);

  --dfNodeType: flex;
  --dfNodeTypeFloat: none;
  --dfNodeBackgroundColor: rgba(219, 223, 253, 1);
  --dfNodeTextColor: #000000;
  --dfNodeBorderSize: 2px;
  --dfNodeBorderColor: rgba(0, 0, 0, 1);
  --dfNodeBorderRadius: 50px;
  --dfNodeMinHeight: 50px;
  --dfNodeMinWidth: 50px;
  --dfNodePaddingTop: 15px;
  --dfNodePaddingBottom: 15px;
  --dfNodeBoxShadowHL: 0px;
  --dfNodeBoxShadowVL: 0px;
  --dfNodeBoxShadowBR: 0px;
  --dfNodeBoxShadowS: 0px;
  --dfNodeBoxShadowColor: rgba(0, 0, 0, 1);

  --dfNodeHoverBackgroundColor: rgba(36, 47, 155, 1);
  --dfNodeHoverTextColor: rgba(255, 255, 255, 1);
  --dfNodeHoverBorderSize: 2px;
  --dfNodeHoverBorderColor: #000000;
  --dfNodeHoverBorderRadius: 50px;

  --dfNodeHoverBoxShadowHL: 2px;
  --dfNodeHoverBoxShadowVL: 3px;
  --dfNodeHoverBoxShadowBR: 15px;
  --dfNodeHoverBoxShadowS: 2px;
  --dfNodeHoverBoxShadowColor: rgba(100, 111, 212, 1);

  --dfNodeSelectedBackgroundColor: rgba(36, 47, 155, 1);
  --dfNodeSelectedTextColor: #ffffff;
  --dfNodeSelectedBorderSize: 2px;
  --dfNodeSelectedBorderColor: #000000;
  --dfNodeSelectedBorderRadius: 50px;

  --dfNodeSelectedBoxShadowHL: 2px;
  --dfNodeSelectedBoxShadowVL: 2px;
  --dfNodeSelectedBoxShadowBR: 15px;
  --dfNodeSelectedBoxShadowS: 2px;
  --dfNodeSelectedBoxShadowColor: rgba(100, 111, 212, 1);

  --dfInputBackgroundColor: #ffffff;
  --dfInputBorderSize: 2px;
  --dfInputBorderColor: #000000;
  --dfInputBorderRadius: 50px;
  --dfInputLeft: -23px;
  --dfInputHeight: 10px;
  --dfInputWidth: 10px;

  --dfInputHoverBackgroundColor: rgba(0, 0, 0, 1);
  --dfInputHoverBorderSize: 2px;
  --dfInputHoverBorderColor: rgba(0, 0, 0, 1);
  --dfInputHoverBorderRadius: 50px;

  --dfOutputBackgroundColor: #ffffff;
  --dfOutputBorderSize: 2px;
  --dfOutputBorderColor: #000000;
  --dfOutputBorderRadius: 50px;
  --dfOutputRight: -8px;
  --dfOutputHeight: 10px;
  --dfOutputWidth: 10px;

  --dfOutputHoverBackgroundColor: rgba(0, 0, 0, 1);
  --dfOutputHoverBorderSize: 2px;
  --dfOutputHoverBorderColor: rgba(0, 0, 0, 1);
  --dfOutputHoverBorderRadius: 50px;

  --dfLineWidth: 4px;
  --dfLineColor: rgba(100, 111, 212, 1);
  --dfLineHoverColor: rgba(88, 0, 255, 1);
  --dfLineSelectedColor: rgba(88, 0, 255, 1);

  --dfRerouteBorderWidth: 2px;
  --dfRerouteBorderColor: #000000;
  --dfRerouteBackgroundColor: #ffffff;

  --dfRerouteHoverBorderWidth: 2px;
  --dfRerouteHoverBorderColor: #000000;
  --dfRerouteHoverBackgroundColor: #ffffff;

  --dfDeleteDisplay: block;
  --dfDeleteColor: #ffffff;
  --dfDeleteBackgroundColor: #000000;
  --dfDeleteBorderSize: 2px;
  --dfDeleteBorderColor: #ffffff;
  --dfDeleteBorderRadius: 50px;
  --dfDeleteTop: -15px;

  --dfDeleteHoverColor: #000000;
  --dfDeleteHoverBackgroundColor: #ffffff;
  --dfDeleteHoverBorderSize: 2px;
  --dfDeleteHoverBorderColor: #000000;
  --dfDeleteHoverBorderRadius: 50px;

}

#drawflow {
  background: var(--dfBackgroundColor);
  background-size: var(--dfBackgroundSize) var(--dfBackgroundSize);
  background-image: var(--dfBackgroundImage);
}

.drawflow .drawflow-node {
  display: var(--dfNodeType);
  background: var(--dfNodeBackgroundColor);
  color: var(--dfNodeTextColor);
  border: var(--dfNodeBorderSize)  solid var(--dfNodeBorderColor);
  border-radius: var(--dfNodeBorderRadius);
  min-height: var(--dfNodeMinHeight);
  width: auto;
  min-width: var(--dfNodeMinWidth);
  padding-top: var(--dfNodePaddingTop);
  padding-bottom: var(--dfNodePaddingBottom);
  -webkit-box-shadow: var(--dfNodeBoxShadowHL) var(--dfNodeBoxShadowVL) var(--dfNodeBoxShadowBR) var(--dfNodeBoxShadowS) var(--dfNodeBoxShadowColor);
  box-shadow:  var(--dfNodeBoxShadowHL) var(--dfNodeBoxShadowVL) var(--dfNodeBoxShadowBR) var(--dfNodeBoxShadowS) var(--dfNodeBoxShadowColor);
}

.drawflow .drawflow-node:hover {
  background: var(--dfNodeHoverBackgroundColor);
  color: var(--dfNodeHoverTextColor);
  border: var(--dfNodeHoverBorderSize)  solid var(--dfNodeHoverBorderColor);
  border-radius: var(--dfNodeHoverBorderRadius);
  -webkit-box-shadow: var(--dfNodeHoverBoxShadowHL) var(--dfNodeHoverBoxShadowVL) var(--dfNodeHoverBoxShadowBR) var(--dfNodeHoverBoxShadowS) var(--dfNodeHoverBoxShadowColor);
  box-shadow:  var(--dfNodeHoverBoxShadowHL) var(--dfNodeHoverBoxShadowVL) var(--dfNodeHoverBoxShadowBR) var(--dfNodeHoverBoxShadowS) var(--dfNodeHoverBoxShadowColor);
}

.drawflow .drawflow-node.selected {
  background: var(--dfNodeSelectedBackgroundColor);
  color: var(--dfNodeSelectedTextColor);
  border: var(--dfNodeSelectedBorderSize)  solid var(--dfNodeSelectedBorderColor);
  border-radius: var(--dfNodeSelectedBorderRadius);
  -webkit-box-shadow: var(--dfNodeSelectedBoxShadowHL) var(--dfNodeSelectedBoxShadowVL) var(--dfNodeSelectedBoxShadowBR) var(--dfNodeSelectedBoxShadowS) var(--dfNodeSelectedBoxShadowColor);
  box-shadow:  var(--dfNodeSelectedBoxShadowHL) var(--dfNodeSelectedBoxShadowVL) var(--dfNodeSelectedBoxShadowBR) var(--dfNodeSelectedBoxShadowS) var(--dfNodeSelectedBoxShadowColor);
}

.drawflow .drawflow-node .input {
  left: var(--dfInputLeft);
  background: var(--dfInputBackgroundColor);
  border: var(--dfInputBorderSize)  solid var(--dfInputBorderColor);
  border-radius: var(--dfInputBorderRadius);
  height: var(--dfInputHeight);
  width: var(--dfInputWidth);
}

.drawflow .drawflow-node .input:hover {
  background: var(--dfInputHoverBackgroundColor);
  border: var(--dfInputHoverBorderSize)  solid var(--dfInputHoverBorderColor);
  border-radius: var(--dfInputHoverBorderRadius);
}

.drawflow .drawflow-node .outputs {
  float: var(--dfNodeTypeFloat);
}

.drawflow .drawflow-node .output {
  right: var(--dfOutputRight);
  background: var(--dfOutputBackgroundColor);
  border: var(--dfOutputBorderSize)  solid var(--dfOutputBorderColor);
  border-radius: var(--dfOutputBorderRadius);
  height: var(--dfOutputHeight);
  width: var(--dfOutputWidth);
}

.drawflow .drawflow-node .output:hover {
  background: var(--dfOutputHoverBackgroundColor);
  border: var(--dfOutputHoverBorderSize)  solid var(--dfOutputHoverBorderColor);
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

.right-sidebar{
  margin-top: 50px;
  position: fixed;
  right: 0;
  top: 0;
  display: flex;
  flex-direction: row-reverse;
  align-items: center;
}

.right-sidebar__main{
  width: 0;
  height: 94vh;
  background: #DBDFFD;
  /* position: fixed; */
  right: 0;
  top: 0;
  /* margin-top: 6vh; */
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content:start;
}



.right-sidebar__main.open {
  width: 9vw;
}

.right-sidebar__main.open {
  width: 9vw;
}

.right-sidebar__row{
  width: 100%;
  height: 50%;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
  margin-top: 13vh;
}
.right-sidebar__button{
  display: flex;
  align-items: center;
  justify-content: center;
  width: 8vw;
  height: 4vh;
  border-radius: 30px;
  color: white;
  background-color: #242F9B;
  padding: 10px;
  margin: 10px 0px;
  cursor: pointer;
  font-size: 1.2vw;
  border: none;
  margin-top: 7vh;
}
.right-sidebar__node{
  display: flex;
  flex-direction: column-reverse;
  align-items: center;
  justify-content: center;
  text-align: center;
  border-radius: 8px;
  /* border: 2px solid #494949; */
  height: 60px;
  line-height:40px;
  padding: 20px;
  margin: 50px 0px;
  cursor: move;
  font-weight: bold;
  font-size: 1vw;
  white-space: nowrap;
}
.right-sidebar__img{
  width: 2vw;
}

.right-sidebar__main > * {
  visibility: hidden;
}


.right-sidebar__main.open > * {
  visibility: visible;
  opacity: 1;
}


.popBtn{
  width: 40px;
  height: 8vh;
  background-color: #242F9B;
  border-radius: 30px;
  margin-right: -20px;
  /* position: fixed; */
  right: 0;
  top: 0;
  margin-top: 4vh;
  text-align: center;
  color:white;
  z-index: -1;
}
.popBtn__txt{
  margin-top: 25px;
  margin-right: 20px;
  cursor: default;
}

.right-sidebar__main.open ~ .popBtn > .popBtn__txt{
  margin-top: 27px;
  transform: rotate( 180deg );
}
</style>
