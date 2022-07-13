<template>
<div>
    <button @click="exportdf">Export</button>
    <button @click="importdf">Import</button>
    <button @click="addnodedf">Add Node</button>
    <div id="drawflow"></div>
</div>
</template>

<script>
import Vue from 'vue'
/* eslint-disable */
// import Drawflow from 'drawflow'
// import styleDrawflow from 'drawflow/dist/drawflow.min.css' // eslint-disable-line no-use-before-define
import node from '@/components/nodes/node.vue'

export default {
  components: {
    node
  },
  data () {
    return {
      editor: null
    }
  },
  mounted () {
    const id = document.getElementById('drawflow')
    this.editor = new Drawflow(id, Vue)
    this.editor.start()

    this.editor.registerNode('node', node, {}, {})
    this.editor.addNode('node', 0, 1, 150, 300, 'node', {}, 'node', 'vue')
  },
  methods: {
    exportdf () {
      this.exportValue = this.editor.export()
      console.log(this.exportValue)
    },
    importdf () {
      this.editor.import(this.exportValue)
    },
    addnodedf () {
      this.editor.addNode('node', 0, 1, 150, 300, 'node', {select_type: '1'}, 'node', 'vue')
    }
  }
}
</script>

<style>
#drawflow {
  width: 100%;
  height: 800px;
  border: 1px solid red;
}
</style>
