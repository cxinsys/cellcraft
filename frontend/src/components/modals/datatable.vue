<template>
  <div class="layout">
    <div v-if="isLoading" class="loading-layout">
      <span v-if="current_file !== null"> </span>
    </div>
    <div v-else class="table-layout">
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
      lines: null,
      firstLine: null,
      isLoading: true,
      columns: [
        // columns 데이터 형식
        // {
        //   label: "Name",
        //   field: "name",
        // }
      ],
      rows: [
        // rows 데이터 형식
        // { id: 1, name: "John", age: 20, createdAt: "", score: 0.03343 }
      ],
    };
  },
  async mounted() {
    this.current_file = this.$store.getters.getCurrentFile.file;
    console.log(this.current_file);
    if (this.current_file !== "") {
      this.isLoading = true;
      try {
        console.log(this.current_file);
        const dataTableResult = await getResult({
          filename: this.current_file,
        });
        console.log(dataTableResult.data);

        //백엔드에서 넘겨준 dataTable 데이터
        this.lines = dataTableResult.data.split("\n").map((x) => x.split(","));
        this.firstLine = this.lines.splice(0, 1)[0];
        this.columns = this.firstLine.slice(1).map((x) => {
          return { label: x, field: x };
        });
        this.rows = this.lines.map((x) => {
          return Object.assign(
            ...this.firstLine.map((k, i) => ({ [k]: x[i] }))
          );
        });
        this.isLoading = false;
      } catch (error) {
        console.error(error);
        this.isLoading = false;
      }
    }
  },
  computed: {
    checkCurrentNode() {
      return this.$store.getters.getCurrentNode;
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

.loading-layout {
  display: flex;
  align-items: center;
  justify-content: center;
}

.loading-layout span {
  border: 4px solid #f3f3f3;
  border-top: 4px solid #3498db;
  border-radius: 50%;
  width: 40px;
  height: 40px;
  animation: spin 2s linear infinite;
}

@keyframes spin {
  0% {
    transform: rotate(0deg);
  }
  100% {
    transform: rotate(360deg);
  }
}
</style>
