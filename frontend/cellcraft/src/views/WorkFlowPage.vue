<template>
<div>
  <div id="drawflow" @drop="drop($event)" @dragover="allowDrop($event)"></div>
  <section class="right-sidebar">
    <div class="right-sidebar__row">
      <button class="right-sidebar__button" @click="exportdf">Export</button>
      <button class="right-sidebar__button" @click="importdf">Import</button>
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
