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
        // },
        // {
        //   label: "Age",
        //   field: "age",
        //   type: "number",
        // },
        // {
        //   label: "Created On",
        //   field: "createdAt",
        //   type: "date",
        //   dateInputFormat: "yyyy-MM-dd",
        //   dateOutputFormat: "MMM do yy",
        // },
        // {
        //   label: "Percent",
        //   field: "score",
        //   type: "percentage",
        // },
      ],
      rows: [
        // { id: 1, name: "John", age: 20, createdAt: "", score: 0.03343 },
        // {
        //   id: 2,
        //   name: "Jane",
        //   age: 24,
        //   createdAt: "2011-10-31",
        //   score: 0.03343,
        // },
        // {
        //   id: 3,
        //   name: "Susan",
        //   age: 16,
        //   createdAt: "2011-10-30",
        //   score: 0.03343,
        // },
        // {
        //   id: 4,
        //   name: "Chris",
        //   age: 55,
        //   createdAt: "2011-10-11",
        //   score: 0.03343,
        // },
        // {
        //   id: 5,
        //   name: "Dan",
        //   age: 40,
        //   createdAt: "2011-10-21",
        //   score: 0.03343,
        // },
        // {
        //   id: 6,
        //   name: "John",
        //   age: 20,
        //   createdAt: "2011-10-31",
        //   score: 0.03343,
        // },
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

        Object.keys(this.dataTable).forEach((column) => {
          this.columns.push({
            label: column,
            field: column,
          });
        });
        // console.log(Object.values(this.dataTable));

        // for (const [key, value] of Object.entries(
        //   Object.values(this.dataTable)
        // )) {
        //   console.log(`${key}: ${Object.values(value).length}`);
        // }
        console.log(Object.entries(Object.values(this.dataTable)[0]).length);
        const rowLength = Object.entries(
          Object.values(this.dataTable)[0]
        ).length;
        for (let i = 0; i < rowLength; i++) {
          let row = [];
          Object.keys(this.dataTable).forEach((column) => {
            row.push([column, this.dataTable[column][i]]);
          });
          console.log(row);
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
