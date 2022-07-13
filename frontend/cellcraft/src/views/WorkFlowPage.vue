<template>
<div>
  <div id="drawflow"></div>
  <section class="right-sidebar">
    <div class="right-sidebar__row">
      <button class="right-sidebar__button" @click="exportdf">Export</button>
      <button class="right-sidebar__button" @click="importdf">Import</button>
    </div>
    <div class="right-sidebar__row">
      <button class="right-sidebar__button" @click="addSelectNodedf">Select Type</button>
      <button class="right-sidebar__button" @click="addInputNodedf">Input Num</button>
    </div>
  </section>
</div>
</template>

<script>
import Vue from 'vue'
/* eslint-disable */
// import Drawflow from 'drawflow'
// import styleDrawflow from 'drawflow/dist/drawflow.min.css' // eslint-disable-line no-use-before-define
import selectType from '@/components/nodes/selectNode.vue'
import inputNum from '@/components/nodes/inputNumNode.vue'

export default {
  components: {
    selectType
  },
  data () {
    return {
      editor: null,
      exportValue: null,
    }
  },
  mounted () {
    const id = document.getElementById('drawflow')
    Vue.prototype.$df = new Drawflow(id, Vue, this);
    this.$df.start()

    this.$df.registerNode('selectType', selectType, {}, {})
    this.$df.registerNode('inputNum', inputNum, {}, {})
    this.$df.addNode('selectType', 0, 1, 150, 300, 'selectType', {}, 'selectType', 'vue')
    this.$df.addNode('inputNum', 0, 1, 150, 300, 'inputNum', {}, 'inputNum', 'vue')
  },
  methods: {
    exportdf () {
      this.exportValue = this.$df.export()
      console.log(this.exportValue)
    },
    importdf () {
      this.$df.import(this.exportValue)
    },
    addSelectNodedf () {
      this.$df.addNode('selectType', 0, 1, 150, 300, 'selectType', {}, 'selectType', 'vue')
    },
    addInputNodedf () {
      this.$df.addNode('inputNum', 0, 1, 150, 300, 'inputNum', {}, 'inputNum', 'vue')
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
  width: 10vw;
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
  cursor: move;
}
</style>
