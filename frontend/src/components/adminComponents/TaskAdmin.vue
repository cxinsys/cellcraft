<template>
  <div class="layout_admin">
    <div class="first-line">
      <div class="header__text">Tasks</div>
      <div class="search">
        <input type="text" v-model="searchTerm" placeholder="Search by username..." />
      </div>
      <div class="page-size">
        <label for="pageSize">Page Size : </label>
        <select id="pageSize" v-model="pageSize" @change="updatePage">
          <option value="5">5</option>
          <option value="10">10</option>
          <option value="15">15</option>
          <!-- <option value="20">20</option>
          <option value="50">50</option> -->
        </select>
      </div>
    </div>
    <table>
      <thead>
        <tr>
          <th @click="sortTable('no')" style="width: 70px">
            No. <span class="sort-icon">{{ sortIcon("no") }}</span>
          </th>
          <th @click="sortTable('username')">
            Username <span class="sort-icon">{{ sortIcon("username") }}</span>
          </th>
          <th @click="sortTable('workflowTitle')">
            Workflow Title <span class="sort-icon">{{ sortIcon("workflowTitle") }}</span>
          </th>
          <th @click="sortTable('status')">
            status <span class="sort-icon">{{ sortIcon("status") }}</span>
          </th>
          <th @click="sortTable('time')">
            submitted time <span class="sort-icon">{{ sortIcon("time") }}</span>
          </th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="job in displayedJobs" :key="job.no">
          <td>{{ job.no }}</td>
          <td>{{ job.username }}</td>
          <td>{{ job.workflowTitle }}</td>
          <td :class="job.status.replace(' ', '-').toLowerCase()">
            {{ job.status }}
          </td>
          <td>{{ job.time }}</td>
          <!-- <td> -->
          <button class="table-button" @click="dismissJob(job)">
            Dismiss Job
          </button>
          <!-- </td> -->
        </tr>
      </tbody>
    </table>
    <div class="pagination">
      <button :disabled="currentPage === 1" @click="prevPage">Prev</button>
      <span>{{ currentPage }} / {{ totalPages }}</span>
      <button :disabled="currentPage === totalPages" @click="nextPage">Next</button>
    </div>
  </div>
</template>

<script>
import { getFilteredTasks } from '@/api';

export default {
  data() {
    return {
      jobs: [],
      sortKey: "id",
      sortDirection: "dsc",
      pageSize: 15,
      currentPage: 1,
      searchTerm: "",
      searchProject: "",
      totalCount: 0
    };
  },
  async created() {
    await this.fetchTasks();
  },
  computed: {
    sortedJobs() {
      return this.jobs;
    },
    totalPages() {
      return Math.ceil(this.totalCount / this.pageSize);
    },
    displayedJobs() {
      return this.jobs;
    }
  },
  methods: {
    async fetchTasks() {
      try {
        const conditions = {
          amount: this.pageSize,
          page_num: this.currentPage,
          sort: this.sortKey,
          order: this.sortDirection === 'asc' ? 'asc' : 'desc',
          searchTerm: this.searchTerm
        };

        const response = await getFilteredTasks(conditions);

        this.jobs = response.data.map((task, index) => ({
          no: index + 1,
          id: task.id,
          userId: task.user_id,
          username: task.username,
          workflowId: task.workflow_id,
          workflowTitle: task.workflow_title,
          status: task.status,
          time: new Date(task.start_time).toLocaleDateString()
        }));

        this.totalCount = response.total_count;
      } catch (error) {
        console.error('Error fetching tasks:', error);
      }
    },
    async sortTable(key) {
      if (this.sortKey === key) {
        this.sortDirection = this.sortDirection === "asc" ? "desc" : "asc";
      } else {
        this.sortKey = key;
        this.sortDirection = "asc";
      }
      this.currentPage = 1; // 정렬 시 첫 페이지로 이동
      await this.fetchTasks();
    },
    sortIcon(key) {
      if (this.sortKey === key) {
        return this.sortDirection === "asc" ? "▽▲" : "▼△";
      }
      return "▽△";
    },
    async dismissJob(job) {
      try {
        // TODO: Implement dismiss API call
        await this.fetchTasks();
      } catch (error) {
        console.error('Error dismissing job:', error);
      }
    },
    async prevPage() {
      if (this.currentPage > 1) {
        this.currentPage--;
        await this.fetchTasks();
      }
    },
    async nextPage() {
      if (this.currentPage < this.totalPages) {
        this.currentPage++;
        await this.fetchTasks();
      }
    },
    async updatePage() {
      this.currentPage = 1;
      await this.fetchTasks();
    }
  },
  watch: {
    searchTerm: {
      handler: 'updatePage',
      immediate: false
    },
    pageSize: {
      handler: 'updatePage',
      immediate: false
    }
  }
};
</script>

<style scoped>
table {
  width: 100%;
  height: 100%;
  border-collapse: separate;
  border-spacing: 5px;
  /* background-color: #c9c9c9; */
  transition: all 0.3s ease;
  border-radius: 15px;
  /* color: #ffffff; */
}

thead th,
td {
  padding: 10px;
  padding-left: 25px;
  text-align: left;
  border-radius: 10px;
  border: 1px solid #a8a8a8;
  /* box-shadow: 0px 8px 20px rgba(176, 169, 255, 0.25); */
}

th {
  text-transform: capitalize;
  background-color: #474747;
  color: #ffffff;
}

td {
  /* background-color: #535353; */
  transition: all 0.3s ease;
}

th:hover {
  background-color: #616161;
}

button {
  margin-right: 10px;
  color: black;
  padding: 5px;
  left: 10px;
  border-radius: 10px;
  background-color: #eaecff;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  text-transform: capitalize;
}

button:disabled {
  color: #ccc;
}

.table-button {
  color: rgb(255, 255, 255);
  width: 100%;
  height: 100%;
  background-color: #474747;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  text-transform: capitalize;
}

.table-button:hover {
  background-color: #616161;
}

.sort-icon {
  color: rgb(199, 199, 199);
  font-weight: normal;
  font-size: small;
}

.first-line {
  height: 40px;
  width: calc(100% - 10px);
  padding: 5px 5px 0px 5px;
  display: flex;
  justify-content: space-between;
  flex-direction: row;
  align-items: center;
}

.reset-button {
  /* margin-top: -7px; */
  width: 1.5rem;
  height: 1.5rem;
  opacity: 0.7;
}

.reset-button:hover {
  opacity: 0.5;
  cursor: pointer;
}

#pageSize {
  padding: 2px;
  border-radius: 5px;
  border: 1px solid #ccc;
  margin-bottom: 5px;
}

.first-line {
  width: calc(100% - 10px);
  padding: 5px 5px 10px 5px;
  display: flex;
  justify-content: space-between;
  flex-direction: row;
}

.second-line {
  width: calc(100% - 10px);
  padding: 5px 5px 5px 5px;
  display: flex;
  justify-content: space-between;
  flex-direction: row;
}

.search {
  display: flex;
  align-items: center;
}

.search input {
  width: 300px;
  height: 2.5rem;
  border: 1px solid #e1e1e1;
  border-radius: 1rem;
  padding: 0 2rem;
  outline-style: none;
  background: #f7f7f7;
}

.search input:focus {
  border: 1px solid #bcbcbc;
}

.pagination {
  display: flex;
  justify-content: center;
  margin: 20px 0px;
}

.pagination button {
  margin: -5px 10px 0px 10px;
}

.success {
  color: rgb(0, 201, 0);
}

.running {
  color: rgb(254, 151, 49);
}

.failure {
  color: rgb(255, 57, 57);
}

.wait {
  color: gray;
}

.layout_admin {
  padding: 0 2rem 0 1rem;
}

.header__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  line-height: 1rem;
  /* padding-left: 2rem; */
  color: rgba(0, 0, 0, 0.8);
}
</style>
