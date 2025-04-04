<template>
  <div class="layout_admin">
    <div class="first-line">
      <div class="header__text">Users</div>
      <div class="search">
        <input type="text" v-model="searchTerm" placeholder="Search by keyword..." />
      </div>
      <div class="page-size">
        <label for="pageSize">Page Size : </label>
        <select id="pageSize" v-model="pageSize" @change="resetPageNum">
          <option value="5">5</option>
          <option value="10">10</option>
          <option value="15">15</option>
          <option value="20">20</option>
          <option value="50">50</option>
        </select>
      </div>
    </div>
    <table>
      <thead>
        <tr>
          <th @click="sortTable('id')" style="width: 70px">
            id <span class="sort-icon">{{ sortIcon("id") }}</span>
          </th>
          <th @click="sortTable('username')">
            name <span class="sort-icon">{{ sortIcon("username") }}</span>
          </th>
          <th @click="sortTable('email')">
            e-mail <span class="sort-icon">{{ sortIcon("email") }}</span>
          </th>
          <!-- <th>password</th> -->
        </tr>
      </thead>
      <tbody>
        <tr v-for="user in users" :key="user.id">
          <td>{{ user.id }}</td>
          <td>{{ user.username }}</td>
          <td>{{ user.email }}</td>
          <!-- <td>
            <span class="blind-password">****</span>
          </td> -->
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
import { getFilteredUsers } from "@/api/index";

export default {
  data() {
    return {
      users: [],
      sortKey: "id", // Set initial sort key to 'id'
      sortDirection: "asc", // Set initial sort direction to 'asc'
      pageSize: 20,
      currentPage: 1,
      searchTerm: "",
      usersCount: 0,
    };
  },
  async mounted() {
    await this.updateUsers();
  },
  computed: {
    totalPages() {
      return Math.ceil(this.usersCount / this.pageSize);
    },
  },
  methods: {
    async updateUsers() {
      const conditions = {
        amount: this.pageSize,
        page_num: this.currentPage,
        sort: this.sortKey,
        order: this.sortDirection,
        searchTerm: this.searchTerm,
      };
      const response = await getFilteredUsers(conditions);
      this.users = response.data;
      this.usersCount = response.total_count;
    },
    async sortTable(key) {
      if (this.sortKey === key) {
        this.sortDirection = this.sortDirection === "asc" ? "desc" : "asc";
      } else {
        this.sortKey = key;
        this.sortDirection = "asc";
      }
      await this.updateUsers();
    },
    sortIcon(key) {
      if (this.sortKey === key) {
        return this.sortDirection === "asc" ? "▽▲" : "▼△";
      }
      return "▽△";
    },
    resetPageNum() {
      this.currentPage = 1;
      this.updateUsers();
    },
  },
  watch: {
    searchTerm: {
      handler: function () {
        this.currentPage = 1;
        this.updateUsers();
      },
      immediate: true,
    },
    pageSize: {
      handler: function () {
        this.currentPage = 1;
        this.updateUsers();
      },
      immediate: true,
    },
    currentPage: {
      handler: function () {
        this.updateUsers();
      },
      immediate: true,
    },
  },
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
  /* box-shadow: 0px 4px 4px rgba(176, 169, 255, 0.25); */
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

.sort-icon {
  color: rgb(199, 199, 199);
  font-weight: normal;
  font-size: small;
}

.first-line {
  height: 40px;
  margin-bottom: 10px;
  width: calc(100% - 10px);
  padding: 5px 5px 0px 5px;
  display: flex;
  justify-content: space-between;
  flex-direction: row;
  align-items: center;
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

#pageSize {
  padding: 2px;
  border-radius: 5px;
  border: 1px solid #ccc;
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

.header__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  line-height: 1rem;
  /* padding-left: 2rem; */
  color: rgba(0, 0, 0, 0.8);
}

.layout_admin {
  padding: 0 2rem 0 1rem;
}
</style>
