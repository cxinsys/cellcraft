<template>
  <div class="layout">
    <div class="table-layout">
      <vue-good-table
        class="dataTable"
        :columns="columns"
        :rows="rows"
        theme="polar-bear"
      />
    </div>
  </div>
</template>

<script>
import { getResult } from "@/api/index";
import "vue-good-table/dist/vue-good-table.css";
import { VueGoodTable } from "vue-good-table";

export default {
  props: {
    file_name: null,
  },
  components: {
    VueGoodTable,
  },
  data() {
    return {
      node_name: "DataTable",
      dataTable: null,
      current_file: null,
      columns: [
        // {
        //   label: "Name",
        //   field: "name",
        // }
      ],
      rows: [
        // { id: 1, name: "John", age: 20, createdAt: "", score: 0.03343 }
      ],
    };
  },
  mounted() {},
  computed: {
    checkCurrentNode() {
      return this.$store.getters.getCurrentNode;
    },
  },
  watch: {
    async checkCurrentNode(val) {
      const current_node = this.$store.getters.getNodeInfo(val);
      this.current_file = this.$store.getters.getCurrentFile.file;
      // console.log(current_node);
      // console.log(this.current_file.file);
      if (current_node.name === "DataTable") {
        const filename = {
          filename: `${this.node_name}_${this.current_file.replace(
            ".csv",
            ""
          )}`,
        };
        console.log(filename);
        const dataTableResult = await getResult(filename);
        console.log(dataTableResult.data);
        this.dataTable = dataTableResult.data;

        //dataTable 데이터 추가
        this.columns = [];
        this.rows = [];
        Object.keys(this.dataTable).forEach((column) => {
          this.columns.push({
            label: column,
            field: column,
          });
        });
        // console.log(Object.entries(Object.values(this.dataTable)[0]).length);
        const rowLength = Object.entries(
          Object.values(this.dataTable)[0]
        ).length;
        for (let i = 0; i < rowLength; i++) {
          let row = [];
          Object.keys(this.dataTable).forEach((column) => {
            row.push([column, this.dataTable[column][i]]);
          });
          this.rows.push(Object.fromEntries(row));
        }
      }
    },
  },
};
</script>

<style scoped>
.layout {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}
.table-layout {
  width: 90%;
  height: 90%;
  background: #fff;
  box-shadow: rgba(100, 100, 111, 0.2) 0px 7px 29px 0px;
  border-radius: 0.5rem;
  padding: 0.5rem;
}
.dataTable {
  width: 100%;
  height: 100%;
  overflow: scroll;
}
</style>
