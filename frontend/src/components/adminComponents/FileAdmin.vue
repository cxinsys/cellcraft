<template>
  <div class="layout_admin">
    <div class="first-line">
      <div class="header__text">Files</div>
      <div class="search">
        <input type="text" v-model="searchTerm" placeholder="Search by title..." />
        <!-- <img
          class="reset-button"
          src="@/assets/reset.png"
          alt="reset"
          @click="resetSearch"
        /> -->
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
    <!-- <div class="second-line">
      <a class="upload-button" href="https://github.com/chxhyxn/TmpCellcraftBoard">
        ⇪ upload new dataset
      </a>
    </div> -->
    <table>
      <thead>
        <tr>
          <th>types</th>
          <th @click="sortTable('name')">
            title <span class="sort-icon">{{ sortIcon("name") }}</span>
          </th>
          <th @click="sortTable('size')">
            size <span class="sort-icon">{{ sortIcon("size") }}</span>
          </th>
          <th @click="sortTable('username')">
            username <span class="sort-icon">{{ sortIcon("username") }}</span>
          </th>
          <th>File URL</th>
          <th>Actions</th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="dataset in displayeddatasets" :key="dataset.no">
          <td>{{ dataset.type }}</td>
          <td>{{ dataset.name }}</td>
          <td>{{ dataset.size }}</td>
          <td>{{ dataset.username }}</td>
          <td><a :href="dataset.url">View File</a></td>
          <td>
            <button @click="deleteFile(dataset)" class="table-button delete">
              Delete
            </button>
          </td>
        </tr>
      </tbody>
    </table>
    <div class="pagination">
      <button :disabled="currentPage === 1" @click="prevPage">Prev</button>
      <span>{{ currentPage }} / {{ totalPages }}</span>
      <button :disabled="currentPage === totalPages" @click="nextPage">
        Next
      </button>
    </div>
  </div>
</template>

<script>
import { getFilteredFiles, deleteFile, getFilesCount } from '@/api';

function formatBytes(bytes, decimals = 2) {
  if (bytes === 0) return "0 Bytes";

  const k = 1024;
  const dm = decimals < 0 ? 0 : decimals;
  const sizes = ["Bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"];

  const i = Math.floor(Math.log(bytes) / Math.log(k));

  return parseFloat((bytes / Math.pow(k, i)).toFixed(dm)) + " " + sizes[i];
}

export default {
  data() {
    return {
      datasets: [],
      sortKey: "id",
      sortDirection: "dsc",
      pageSize: 15,
      currentPage: 1,
      searchTerm: "",
      totalCount: 0
    };
  },
  async created() {
    await this.fetchFiles();
  },
  computed: {
    sorteddatasets() {
      const datasetsCopy = [...this.datasets];
      if (this.sortKey) {
        datasetsCopy.sort((a, b) => {
          let aValue = a[this.sortKey];
          let bValue = b[this.sortKey];

          // 파일 크기 정렬을 위한 처리
          if (this.sortKey === 'size') {
            aValue = this.parseFileSize(aValue);
            bValue = this.parseFileSize(bValue);
          }

          if (aValue < bValue) return this.sortDirection === "asc" ? -1 : 1;
          if (aValue > bValue) return this.sortDirection === "asc" ? 1 : -1;
          return 0;
        });
      }
      return datasetsCopy;
    },
    totalPages() {
      return Math.ceil(this.totalCount / this.pageSize);
    },
    displayeddatasets() {
      return this.sorteddatasets;
    },
    filtereddatasets() {
      if (this.searchTerm) {
        const searchTermLower = this.searchTerm.toLowerCase();
        return this.sorteddatasets.filter((dataset) =>
          dataset.name.toLowerCase().includes(searchTermLower)
        );
      } else {
        return this.sorteddatasets;
      }
    }
  },
  methods: {
    async fetchFiles() {
      try {
        const [filesResponse, countResponse] = await Promise.all([
          getFilteredFiles({
            amount: this.pageSize,
            page_num: this.currentPage,
            sort: this.sortKey,
            order: this.sortDirection,
            searchTerm: this.searchTerm,
          }),
          getFilesCount(),
        ]);

        this.datasets = filesResponse.data.data.map((file) => ({
          id: file.id,
          type: file.file_name.split('.').pop(),
          name: file.file_name,
          size: formatBytes(file.file_size),
          url: file.file_path,
          username: file.username
        }));
        this.totalCount = countResponse.data;
      } catch (error) {
        console.error("Error fetching files:", error);
      }
    },
    async sortTable(key) {
      if (this.sortKey === key) {
        this.sortDirection = this.sortDirection === "asc" ? "desc" : "asc";
      } else {
        this.sortKey = key;
        this.sortDirection = "asc";
      }
      await this.fetchFiles();
    },
    sortIcon(key) {
      if (this.sortKey === key) {
        return this.sortDirection === "asc" ? "▽▲" : "▼△";
      }
      return "▽△";
    },
    async updatePage() {
      this.currentPage = 1;
      await this.fetchFiles();
    },
    async deleteFile(dataset) {
      if (confirm(`Are you sure you want to delete file ${dataset.name}?`)) {
        try {
          await deleteFile(dataset.id);
          await this.fetchFiles();
        } catch (error) {
          console.error('Error deleting file:', error);
        }
      }
    },
    async prevPage() {
      if (this.currentPage > 1) {
        this.currentPage--;
        await this.fetchFiles();
      }
    },
    async nextPage() {
      if (this.currentPage < this.totalPages) {
        this.currentPage++;
        await this.fetchFiles();
      }
    },
    // 파일 크기 문자열을 바이트 단위로 변환하는 함수
    parseFileSize(sizeString) {
      const units = {
        'B': 1,
        'KB': 1024,
        'MB': 1024 * 1024,
        'GB': 1024 * 1024 * 1024,
        'TB': 1024 * 1024 * 1024 * 1024
      };

      const match = sizeString.match(/^([\d.]+)\s*([A-Za-z]+)$/);
      if (match) {
        const value = parseFloat(match[1]);
        const unit = match[2].toUpperCase();
        return value * (units[unit] || 1);
      }
      return 0;
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
  border-collapse: separate;
  border-spacing: 5px;
  transition: all 0.3s ease;
  border-radius: 15px;
  table-layout: fixed;
}

thead th,
td {
  padding: 10px;
  padding-left: 15px;
  text-align: left;
  border-radius: 10px;
  border: 1px solid #a8a8a8;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}

th {
  text-transform: capitalize;
  background-color: #474747;
  color: #ffffff;
  position: sticky;
  top: 0;
}

td {
  transition: all 0.3s ease;
  background-color: #ffffff;
}

th:hover {
  background-color: #616161;
}

/* 컬럼 너비 설정 */
th:nth-child(1) {
  width: 10%;
}

/* types */
th:nth-child(2) {
  width: 25%;
}

/* title */
th:nth-child(3) {
  width: 10%;
}

/* size */
th:nth-child(4) {
  width: 20%;
}

/* username */
th:nth-child(5) {
  width: 15%;
}

/* File URL */
th:nth-child(6) {
  width: 20%;
}

/* Actions */

/* 반응형 스타일 */
@media screen and (max-width: 1200px) {
  .layout_admin {
    padding: 0 1rem;
  }

  .search input {
    width: 200px;
  }
}

@media screen and (max-width: 768px) {
  .first-line {
    flex-direction: column;
    height: auto;
    gap: 10px;
  }

  .search input {
    width: 100%;
  }

  .page-size {
    width: 100%;
  }

  th,
  td {
    padding: 8px;
    font-size: 0.9rem;
  }
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
  background-color: #323232;
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

a {
  /* color: #f0f1fb; */
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

.second-line {
  width: calc(100% - 10px);
  padding: 0px 5px 5px 5px;
  margin-bottom: 10px;
  display: flex;
  justify-content: flex-end;
  flex-direction: row;
}

#pageSize {
  padding: 2px;
  border-radius: 5px;
  border: 1px solid #ccc;
  margin-bottom: 5px;
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

.upload-button {
  border-radius: 5px;
  padding: 7px 10px;
  font-weight: 500;
  /* background-color: #eaecff; */
  background-color: #5a5a5a;
  /* border-color: #e7eaff; */
  border: 2px solid #e7eaff;
  border-radius: 10px;
  text-decoration: none;
  color: rgb(255, 255, 255);
  text-transform: capitalize;
}

.upload-button:hover {
  cursor: pointer;
  background-color: #7d7d7d;
}

.reset-button {
  margin-top: -7px;
  width: 1.5rem;
  height: 1.5rem;
  opacity: 0.7;
}

.reset-button:hover {
  opacity: 0.5;
  cursor: pointer;
}

.pagination {
  display: flex;
  justify-content: center;
  margin: 20px 0px;
}

.pagination button {
  margin: -5px 10px 0px 10px;
}

button:disabled {
  color: #ccc;
  border-color: #ccc;
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

.table-button.delete {
  background-color: #ff4444;
}

.table-button.delete:hover {
  background-color: #cc0000;
}
</style>
