<template>
  <div>
    <div class="first-line">
      <div class="search">
        <input type="text" v-model="searchTerm" placeholder="Search by id..." />
        <img
          class="reset-button"
          src="@/assets/reset.png"
          alt="reset"
          @click="resetSearch"
        />
      </div>
      <div class="page-size">
        <label for="pageSize">Page Size : </label>
        <select id="pageSize" v-model="pageSize" @change="updatePage">
          <option value="5">5</option>
          <option value="10">10</option>
          <option value="15">15</option>
          <option value="20">20</option>
          <option value="50">50</option>
        </select>
      </div>
    </div>
    <div class="second-line">
      <div class="search">
        <input
          type="text"
          v-model="searchProject"
          placeholder="Search by projects..."
        />
        <img
          class="reset-button"
          src="@/assets/reset.png"
          alt="reset"
          @click="resetProjectSearch"
        />
      </div>
    </div>
    <table>
      <thead>
        <tr>
          <th @click="sortTable('no')">
            No. <span class="sort-icon">{{ sortIcon("no") }}</span>
          </th>
          <th @click="sortTable('userId')">
            id <span class="sort-icon">{{ sortIcon("userId") }}</span>
          </th>
          <th>project</th>
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
          <td>{{ job.userId }}</td>
          <td>{{ job.project }}</td>
          <td :class="job.status.replace(' ', '-').toLowerCase()">
            {{ job.status }}
          </td>
          <td>{{ job.time }}</td>
          <td>
            <button @click="dismissJob(job)">Dismiss Job</button>
          </td>
        </tr>
      </tbody>
    </table>
    <div class="pagination">
      <button :disabled="currentPage === 1" @click="currentPage--">Prev</button>
      <span>{{ currentPage }}</span>
      <button :disabled="currentPage === totalPages" @click="currentPage++">
        Next
      </button>
    </div>
  </div>
</template>

<script>
export default {
  data() {
    return {
      jobs: [
        {
          no: 1,
          userId: "johndoe",
          project: "project1",
          status: "success",
          time: "2023-06-29",
        },
        {
          no: 2,
          userId: "janesmith",
          project: "project2",
          status: "running",
          time: "2023-06-30",
        },
        {
          no: 3,
          userId: "bobjohnson",
          project: "project12",
          status: "failure",
          time: "2023-07-01",
        },
        {
          no: 4,
          userId: "alicebrown",
          project: "project14",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 5,
          userId: "samwilson",
          project: "project19",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 6,
          userId: "emilydavis",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 7,
          userId: "michaelwilson",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 8,
          userId: "oliviajohnson",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 9,
          userId: "sophiamiller",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 10,
          userId: "williamanderson",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 11,
          userId: "benjamingarcia",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 12,
          userId: "avamartinez",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 13,
          userId: "miathompson",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 14,
          userId: "ethanlopez",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 15,
          userId: "jameswilson",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 16,
          userId: "liamwhite",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 17,
          userId: "sophiabrown",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 18,
          userId: "charlottedavis",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 19,
          userId: "alexanderjohnson",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
        {
          no: 20,
          userId: "emmamiller",
          project: "project1",
          status: "wait",
          time: "2023-07-01",
        },
      ],
      sortKey: "no", // Set initial sort key to 'id'
      sortDirection: "dsc", // Set initial sort direction to 'asc'
      pageSize: 20,
      currentPage: 1,
      searchTerm: "",
      searchProject: "",
    };
  },
  computed: {
    sortedJobs() {
      const jobsCopy = [...this.jobs];
      if (this.sortKey) {
        jobsCopy.sort((a, b) => {
          const aValue = a[this.sortKey];
          const bValue = b[this.sortKey];
          if (aValue < bValue) return this.sortDirection === "asc" ? -1 : 1;
          if (aValue > bValue) return this.sortDirection === "asc" ? 1 : -1;
          return 0;
        });
      }
      return jobsCopy;
    },

    totalPages() {
      return Math.ceil(this.filteredJobs.length / this.pageSize);
    },
    displayedJobs() {
      const startIndex = (this.currentPage - 1) * this.pageSize;
      const endIndex = startIndex + this.pageSize;
      return this.filteredJobs.slice(startIndex, endIndex);
    },
    filteredJobs() {
      if (this.searchTerm || this.searchProject) {
        const searchTermLower = this.searchTerm.toLowerCase();
        const searchProjectLower = this.searchProject.toLowerCase();
        return this.sortedJobs.filter((job) => {
          const userIdMatch = job.userId
            .toLowerCase()
            .includes(searchTermLower);
          const projectMatch = job.project
            .toLowerCase()
            .includes(searchProjectLower);
          return userIdMatch && projectMatch;
        });
      } else {
        return this.sortedJobs;
      }
    },
  },
  methods: {
    sortTable(key) {
      if (this.sortKey === key) {
        this.sortDirection = this.sortDirection === "asc" ? "desc" : "asc";
      } else {
        this.sortKey = key;
        this.sortDirection = "asc";
      }
    },
    sortIcon(key) {
      if (this.sortKey === key) {
        return this.sortDirection === "asc" ? "▽▲" : "▼△";
      }
      return "▽△";
    },
    resetSearch() {
      this.searchTerm = "";
    },
    resetProjectSearch() {
      this.searchProject = "";
    },
    dismissJob(job) {
      const index = this.jobs.findIndex((j) => j.id === job.id);
      if (index !== -1) {
        this.jobs.splice(index, 1);
      }
    },
    updatePage() {
      this.currentPage = 1; // Reset to first page when page size changes
    },
  },
};
</script>

<style scoped>
table {
  width: 100%;
  height: 100%;
  border-collapse: collapse;
}

thead th {
  background-color: #f5f5f5;
  font-weight: bold;
  text-align: left;
  padding: 10px;
  border-bottom: 1px solid #ccc;
  cursor: pointer;
  text-transform: capitalize;
}

tbody td {
  max-width: 30px;
  padding: 10px;
  white-space: nowrap; /* 텍스트 줄 바꿈 비활성화 */
  overflow: hidden; /* 텍스트가 넘칠 경우 숨김 처리 */
  text-overflow: ellipsis; /* 텍스트가 넘칠 경우 ...으로 표시 */
  border-bottom: 1px solid #ccc;
}
button {
  margin-right: 10px;
  color: black;
  padding: 2px;
  left: 10px;
  border-radius: 5px;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  margin-bottom: 5px;
  text-transform: capitalize;
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
.sort-icon {
  color: rgb(34, 34, 34);
  font-weight: normal;
  font-size: small;
}
.reset-button {
  margin-top: -7px;

  width: 1.5rem;
  height: 1.5rem;
  opacity: 1;
}
.reset-button:hover {
  opacity: 0.8;
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
  padding: 5px 5px 0px 5px;
  display: flex;
  justify-content: space-between;
  flex-direction: row;
}
.second-line {
  width: calc(100% - 10px);
  padding: 0px 5px 10px 5px;
  display: flex;
  justify-content: space-between;
  flex-direction: row;
}
.search {
  display: flex;
  align-items: center;
}

.search input {
  margin-right: 10px;
  color: black;
  padding: 2px;
  left: 10px;
  border-radius: 5px;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  margin-bottom: 5px;
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
  color: green;
}
.running {
  color: rgb(255, 142, 29);
}
.failure {
  color: red;
}
.wait {
  color: gray;
}
</style>
