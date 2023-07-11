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
    <table>
      <thead>
        <tr>
          <th @click="sortTable('id')">
            id <span class="sort-icon">{{ sortIcon("id") }}</span>
          </th>
          <th @click="sortTable('name')">
            name <span class="sort-icon">{{ sortIcon("name") }}</span>
          </th>
          <th @click="sortTable('email')">
            e-mail <span class="sort-icon">{{ sortIcon("email") }}</span>
          </th>
          <th>password</th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="user in displayedUsers" :key="user.id">
          <td>{{ user.id }}</td>
          <td>{{ user.name }}</td>
          <td>{{ user.email }}</td>
          <td><span class="blind-password">****</span>{{ user.password }}</td>
          <td>
            <button @click="resetPassword(user)">Reset Password</button>
          </td>
          <td>
            <button @click="deleteUser(user)">Delete User</button>
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
      users: [
        {
          id: "johndoe",
          name: "John Doe",
          email: "johndoe@example.com",
          password: "1234",
        },
        {
          id: "janesmith",
          name: "Jane Smith",
          email: "janesmith@example.com",
          password: "1234",
        },
        {
          id: "bobjohnson",
          name: "Bob Johnson",
          email: "bobjohnson@example.com",
          password: "1234",
        },
        {
          id: "alicebrown",
          name: "Alice Brown",
          email: "alicebrown@example.com",
          password: "1234",
        },
        {
          id: "samwilson",
          name: "Sam Wilson",
          email: "samwilson@example.com",
          password: "1234",
        },
        {
          id: "emilydavis",
          name: "Emily Davis",
          email: "emilydavis@example.com",
          password: "1234",
        },
        {
          id: "michaelwilson",
          name: "Michael Wilson",
          email: "michaelwilson@example.com",
          password: "1234",
        },
        {
          id: "oliviajohnson",
          name: "Olivia Johnson",
          email: "oliviajohnson@example.com",
          password: "1234",
        },
        {
          id: "sophiamiller",
          name: "Sophia Miller",
          email: "sophiamiller@example.com",
          password: "1234",
        },
        {
          id: "williamanderson",
          name: "William Anderson",
          email: "williamanderson@example.com",
          password: "1234",
        },
        {
          id: "benjamingarcia",
          name: "Benjamin Garcia",
          email: "benjamingarcia@example.com",
          password: "1234",
        },
        {
          id: "avamartinez",
          name: "Ava Martinez",
          email: "avamartinez@example.com",
          password: "1234",
        },
        {
          id: "miathompson",
          name: "Mia Thompson",
          email: "miathompson@example.com",
          password: "1234",
        },
        {
          id: "ethanlopez",
          name: "Ethan Lopez",
          email: "ethanlopez@example.com",
          password: "1234",
        },
        {
          id: "jameswilson",
          name: "James Wilson",
          email: "jameswilson@example.com",
          password: "1234",
        },
        {
          id: "liamwhite",
          name: "Liam White",
          email: "liamwhite@example.com",
          password: "1234",
        },
        {
          id: "sophiabrown",
          name: "Sophia Brown",
          email: "sophiabrown@example.com",
          password: "1234",
        },
        {
          id: "charlottedavis",
          name: "Charlotte Davis",
          email: "charlottedavis@example.com",
          password: "1234",
        },
        {
          id: "alexanderjohnson",
          name: "Alexander Johnson",
          email: "alexanderjohnson@example.com",
          password: "1234",
        },
        {
          id: "emmamiller",
          name: "Emma Miller",
          email: "emmamiller@example.com",
          password: "1234",
        },
      ],
      sortKey: "id", // Set initial sort key to 'id'
      sortDirection: "asc", // Set initial sort direction to 'asc'
      pageSize: 20,
      currentPage: 1,
      searchTerm: "",
    };
  },
  computed: {
    sortedUsers() {
      const usersCopy = [...this.users];
      if (this.sortKey) {
        usersCopy.sort((a, b) => {
          const aValue = a[this.sortKey];
          const bValue = b[this.sortKey];
          if (aValue < bValue) return this.sortDirection === "asc" ? -1 : 1;
          if (aValue > bValue) return this.sortDirection === "asc" ? 1 : -1;
          return 0;
        });
      }
      return usersCopy;
    },
    totalPages() {
      return Math.ceil(this.filteredUsers.length / this.pageSize);
    },
    displayedUsers() {
      const startIndex = (this.currentPage - 1) * this.pageSize;
      const endIndex = startIndex + this.pageSize;
      return this.filteredUsers.slice(startIndex, endIndex);
    },
    filteredUsers() {
      if (this.searchTerm) {
        const searchTermLower = this.searchTerm.toLowerCase();
        return this.sortedUsers.filter((user) =>
          user.id.toLowerCase().includes(searchTermLower)
        );
      } else {
        return this.sortedUsers;
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
    resetPassword(user) {
      user.password = "0000";
    },
    deleteUser(user) {
      const index = this.users.findIndex((u) => u.id === user.id);
      if (index !== -1) {
        this.users.splice(index, 1);
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
  color: black;
  padding: 2px;
  left: 10px;
  border-radius: 5px;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  margin-bottom: 5px;
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
