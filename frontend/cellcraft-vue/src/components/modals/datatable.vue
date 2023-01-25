<template>
  <div class="layout">
    <div class="dataTable" v-html="dataTable"></div>
  </div>
</template>

<script>
import { getResult } from "@/api/index";
export default {
  props: {
    file_name: null,
  },
  data() {
    return {
      node_name: "DataTable",
      dataTable: null,
    };
  },
  async mounted() {
    const filename = { filename: `${this.node_name}_${this.file_name}` };
    console.log(filename);
    const dataTableResult = await getResult(filename);
    console.log(dataTableResult.data);
    this.dataTable = dataTableResult.data;
  },
};
</script>

<style scoped>
.dataTable {
  display: flex;
  align-items: center;
  justify-content: flex-start;
  align-content: flex-start;
  flex-wrap: wrap;
  /* width: 100%; */
  height: 400px;
  overflow-x: scroll;
  overflow-y: scroll;
  margin-top: 15px;
  margin-right: 10px;
  margin-left: 0px;
}
table {
  border-collapse: separate;
  text-align: center;
  line-height: 10px;
  border-top: 1px solid #ccc;
  margin: 20px 10px;
}
th {
  width: auto;
  max-height: 10px;
  padding: 10px;
  font-weight: bold;
  vertical-align: top;
  border-bottom: 1px solid #ccc;
  background: #efefef;
  word-wrap: normal;
  text-align: center;
  position: sticky;
  top: 0;
}
td {
  width: auto;
  text-align: start;
  max-height: 10px;
  /* padding: 10px; */
  vertical-align: middle;
  border: 1px solid #ccc;
  word-wrap: break-word;
  word-break: keep-all;
  padding: 0px 3px;
  overflow: hidden;
  white-space: nowrap;
  text-overflow: ellipsis;
}
</style>
