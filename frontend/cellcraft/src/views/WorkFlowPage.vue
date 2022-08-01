<template>
<div>
  <div id="drawflow" @drop="drop($event)" @dragover="allowDrop($event)"></div>
  <section class="right-sidebar">
    <div class="right-sidebar__row">
      <button class="right-sidebar__button" @click="exportdf">complie</button>
    </div>
    <div class="right-sidebar__row">
      <ul>
        <li class="right-sidebar__drag-drawflow" v-for="(node, idx) in listNodes" :key="idx" draggable="true" :data-node="node.name" @dragstart="drag($event)">
          <div class="right-sidebar__node" :style="`background: ${node.color}`" >{{ node.name }}</div>
        </li>
      </ul>
    </div>
  </section>
</div>
</template>

<script>
import Vue from 'vue'
/* eslint-disable */
// import Drawflow from 'drawflow'
// import styleDrawflow from 'drawflow/dist/drawflow.min.css' // eslint-disable-line no-use-before-define
import selectReqest from '@/components/nodes/selectNode.vue'
import inputNum from '@/components/nodes/inputNumNode.vue'
import operator from '@/components/nodes/operatorNode.vue'
import result from '@/components/nodes/resultNode.vue'
import fileUpload from '@/components/nodes/fileUploadNode.vue'
import dataTable from '@/components/nodes/dataTableNode.vue'


export default {
  data () {
    return {
      editor: null,
      exportValue: null,
      listNodes: [
        {
          name: 'selectReqest',
          color: 'white',
          input: 0,
          output: 1
        },
        {
          name: 'operator',
          color: 'white',
          input: 2,
          output: 1
        },
        {
          name: 'fileUpload',
          color: 'white',
          input: 0,
          output: 1
        },
        {
          name: 'dataTable',
          color: 'white',
          input: 1,
          output: 1
        }
      ]
    }
  },
  mounted () {
    const id = document.getElementById('drawflow')
    Vue.prototype.$df = new Drawflow(id, Vue, this)
    this.$df.start()

    this.$df.registerNode('selectReqest', selectReqest, {}, {})
    this.$df.registerNode('inputNum', inputNum, {}, {})
    this.$df.registerNode('operator', operator, {}, {})
    this.$df.registerNode('result', result, {}, {})
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
    exportdf () {
      this.exportValue = this.$df.export()
      console.log(typeof (this.exportValue), this.exportValue)
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
    }
  }
}
</script>

<style>
#drawflow {
  width: 90vw;
  height: 95vh;

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
  width: 15vw;
  height: 100vh;
  background: rgb(216, 223, 222);
  position: fixed;
  right: 0;
  top: 0;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content:start;
}
.right-sidebar__row{
  width: 100%;
  height: 50%;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
}
.right-sidebar__button{
  display: flex;
  align-items: center;
  justify-content: center;
  width: 100px;
  height: 50px;
  border-radius: 8px;
  border: 2px solid #494949;
  line-height:40px;
  padding: 10px;
  margin: 10px 0px;
  cursor: pointer;
}
.right-sidebar__node{
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
  border-radius: 8px;
  border: 2px solid #494949;
  height: 60px;
  line-height:40px;
  padding: 10px;
  margin: 10px 0px;
  cursor: move;
}
</style>
