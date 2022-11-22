<template>
  <div class="layout">
    <h1>Data Table</h1>
    <div class="dataTable" v-html="dataTable"></div>
  </div>
</template>

<script>
import { getResult } from '@/api/index'
export default {
  props: {
    file_name: null
  },
  data () {
    return {
      node_name: 'DataTable',
      dataTable: null
    }
  },
  async mounted () {
    const filename = {'filename': `${this.node_name}_${this.file_name}`}
    console.log(filename)
    const dataTableResult = await getResult(filename)
    console.log(dataTableResult.data)
    this.dataTable = dataTableResult.data
  }
}
</script>

<style>
.dataTable{
    display: flex;
    align-items: center;
    justify-content: center;
    width: 100%;
    height: 400px;
    overflow-x: auto;
    overflow-y: auto;
    margin: 10px;
}

table {
  border-collapse: separate;
  text-align: center;
  line-height: 1.5;
  border-top: 1px solid #ccc;
  margin: 20px 10px;
}
th {
  padding: 10px;
  font-weight: bold;
  vertical-align: top;
  border-bottom: 1px solid #ccc;
  background: #efefef;
}
td {
  width: 350px;
  padding: 10px;
  vertical-align: top;
  border-bottom: 1px solid #ccc;
}
</style>
