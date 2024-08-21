<template>
  <div class="layout">
    <div v-if="current_file === null">
      <span> NO DATA FOR TABLE</span>
    </div>
    <!-- <div v-else-if="isLoading" class="loading-layout">
      <span v-if="current_files !== null"></span>
    </div> -->
    <div v-else class="table-layout">
      <vue-good-table class="dataTable" :columns="columns" :rows="rows" mode="remote" @on-page-change="onPageChange"
        @on-sort-change="onSortChange" @on-column-filter="onColumnFilter" @on-per-page-change="onPerPageChange"
        :line-numbers="true" theme="polar-bear" max-height="31.5rem" :isLoading.sync="isLoading"
        :totalRows="totalRecords" :pagination-options="{
          enabled: true,
          perPageDropdownEnabled: false,
          jumpFirstOrLast : true,
          firstLabel: 'First Page',
          lastLabel: 'Last Page'
        }"></vue-good-table>
    </div>
  </div>
</template>

<script>
import { getDataTableFile } from "@/api/index";
import "vue-good-table/dist/vue-good-table.css";
import { VueGoodTable } from "vue-good-table";

export default {
  components: {
    VueGoodTable,
  },
  data() {
    return {
      node_name: "DataTable",
      dataTable: null,
      current_file: null,
      current_files: null,
      lines: null,
      firstLine: null,
      isLoading: true,
      workflowId: this.$route.query.workflow_id,
      nodeId: this.$route.query.node,
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
      serverParams: {
        file_name: "",
        columnFilters: {
        },
        sort: [
          {
            field: '',
            type: ''
          }
        ],
        page: 1,
        perPage: 10
      },
      totalRecords: 0,
    };
  },
  methods: {
    updateParams(newProps) {
      this.serverParams = Object.assign({}, this.serverParams, newProps);
    },
    async onPageChange(params) {
      try {
        console.log(params);
        this.updateParams({ page: params.currentPage });
        console.log(this.serverParams);
        const dataTableResult = await getDataTableFile(this.serverParams);
        const { columns, rows, totalRecords } = dataTableResult.data;
        console.log(columns);
        this.rows = rows;
        this.totalRecords = totalRecords;
      } catch (error) {
        console.error(error);
      }
    },
    async onSortChange(params) {
      try {
        console.log(params);
        this.updateParams({
          sort: params,
        });
        const dataTableResult = await getDataTableFile(this.serverParams);
        const { columns, rows, totalRecords } = dataTableResult.data;
        console.log(columns);
        this.rows = rows;
        this.totalRecords = totalRecords;
      } catch (error) {
        console.error(error);
      }
    },
    async onColumnFilter(params) {
      try {
        console.log(params);
        this.updateParams(params);
        const dataTableResult = await getDataTableFile(this.serverParams);
        const { columns, rows, totalRecords } = dataTableResult.data;
        console.log(columns);
        this.rows = rows;
        this.totalRecords = totalRecords;
      } catch (error) {
        console.error(error);
      }
    },
    async onPerPageChange(params) {
      try {
        this.updateParams({ per_page: params.currentPerPage });
        const dataTableResult = await getDataTableFile(this.serverParams);
        const { columns, rows, totalRecords } = dataTableResult.data;
        console.log(columns);
        this.rows = rows;
        this.totalRecords = totalRecords;
      } catch (error) {
        console.error(error);
      }
    },
  },
  async mounted() {
    this.current_file = this.$store.getters.getWorkflowNodeFileInfo(this.nodeId);
    console.log(this.current_file);
    // const initial_file_ids = Object.keys(this.current_files);
    if (this.current_file) {
      this.isLoading = true;
      try {
        this.serverParams.file_name = this.current_file
        const dataTableResult = await getDataTableFile(this.serverParams);

        console.log(dataTableResult.data);

        const { columns, rows, totalRecords } = dataTableResult.data;
        console.log(columns, rows);
        this.columns = columns;
        this.rows = rows;
        this.totalRecords = totalRecords;

        this.isLoading = false;
      } catch (error) {
        console.error(error);
        this.isLoading = false;
      }
    }
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
  overflow: auto;
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
</style>
