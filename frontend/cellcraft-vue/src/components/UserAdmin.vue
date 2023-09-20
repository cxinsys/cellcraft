<template>
  <div>
    <div class="first-line">
      <div class="search">
        <input
          type="text"
          v-model="searchTerm"
          placeholder="Enter keyword for search..."
        />
        <button @click="updateUsers">search</button>
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
          <th @click="sortTable('id')">
            id <span class="sort-icon">{{ sortIcon("id") }}</span>
          </th>
          <th @click="sortTable('username')">
            name <span class="sort-icon">{{ sortIcon("username") }}</span>
          </th>
          <th @click="sortTable('email')">
            e-mail <span class="sort-icon">{{ sortIcon("email") }}</span>
          </th>
          <th>password</th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="user in users" :key="user.id">
          <td>{{ user.id }}</td>
          <td>{{ user.username }}</td>
          <td>{{ user.email }}</td>
          <td>
            <span class="blind-password">****</span>
            <!-- {{ user.hashed_password }} -->
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
import { getUsersCount, getFilteredUsers } from "@/api/index";

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
      const response = await getUsersCount();
      this.usersCount = response.data;
      console.log(this.usersCount);

      const conditions = {
        amount: this.pageSize,
        page_num: this.currentPage,
        sort: this.sortKey,
        order: this.sortDirection,
        searchTerm: this.searchTerm,
      };
      const filteredUsers = await getFilteredUsers(conditions);
      console.log(filteredUsers.data);
      this.users = filteredUsers.data;
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
      this.currentPage = 1; // Reset to first page when page size changes
    },
  },
};
</script>

<style scoped>
/* table {
  width: 100%;
  height: 100%;
  border-collapse: collapse;
  border-spacing: 0;
  background-color: #f9f9f9;
  transition: all 0.3s ease;
}

thead th,
td {
  padding: 20px;
  text-align: center;
  border-radius: 12px;
}

th {
  background-color: #eee;
  color: #333;
}

td {
  background-color: #fff;
  box-shadow: 0px 4px 12px rgba(0, 0, 0, 0.1);
  transition: all 0.3s ease;
}

tr:hover td {
  background-color: #f2f2f2;
  box-shadow: 0px 6px 16px rgba(0, 0, 0, 0.2);
}

td:hover {
  background-color: #e6e6e6;
  box-shadow: 0px 8px 20px rgba(0, 0, 0, 0.25);
} */

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
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
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
button:disabled {
  color: #ccc;
  border-color: #ccc;
}
.sort-icon {
  color: rgb(34, 34, 34);
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
.search {
  display: flex;
  align-items: center;
}

.search input {
  margin-right: 10px;
  width: 200px;
  color: black;
  padding: 5px;
  left: 10px;
  border-radius: 10px;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  margin-top: 10px;
  margin-bottom: 15px;
}
.blind-password {
  position: absolute;
  padding-left: 10px;
  padding-right: 100px;
  background: white;
}
.blind-password:hover {
  opacity: 0;
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
</style>
