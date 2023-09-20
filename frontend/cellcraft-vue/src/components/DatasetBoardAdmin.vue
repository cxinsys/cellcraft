<template>
  <div>
    <div class="first-line">
      <div class="search">
        <input
          type="text"
          v-model="searchTerm"
          placeholder="Search by title..."
        />
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
      <a
        class="upload-button"
        href="https://github.com/chxhyxn/TmpCellcraftBoard"
      >
        ⇪ UPLOAD NEW DATASET
      </a>
    </div>
    <table>
      <thead>
        <tr>
          <th>types</th>
          <th @click="sortTable('title')">
            title <span class="sort-icon">{{ sortIcon("title") }}</span>
          </th>
          <th>size</th>
          <th>File URL</th>
          <th>Post URL</th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="dataset in displayeddatasets" :key="dataset.no">
          <td>{{ dataset.type }}</td>
          <td>{{ dataset.path }}</td>
          <td>{{ dataset.size }}</td>
          <td><a :href="dataset.url">View File</a></td>
          <td><a :href="dataset.posturl">View Post</a></td>
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
import axios from "axios";

function formatBytes(bytes, decimals = 2) {
  if (bytes === 0) return "0 Bytes";

  const k = 1024;
  const dm = decimals < 0 ? 0 : decimals;
  const sizes = ["Bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"];

  const i = Math.floor(Math.log(bytes) / Math.log(k));

  return parseFloat((bytes / Math.pow(k, i)).toFixed(dm)) + " " + sizes[i];
}

async function checkFileExists(repo, path) {
  try {
    const { data } = await axios.get(
      `https://api.github.com/repos/${repo}/contents/${path}`
    );
    if (data && data.sha) {
      return true;
    }
  } catch (e) {
    return false;
  }
  return false;
}

export default {
  data() {
    return {
      datasets: [],
      sortKey: "no",
      sortDirection: "dsc", // Set initial sort direction to 'asc'
      pageSize: 20,
      currentPage: 1,
      searchTerm: "",
    };
  },
  async mounted() {
    try {
      const { data: repoData } = await axios.get(
        "https://api.github.com/repos/chxhyxn/TmpCellcraftBoard/git/trees/main?recursive=1"
      );
      const dirSha = repoData.tree.find(
        (item) => item.path === "datasets" && item.type === "tree"
      ).sha;

      const { data: dirData } = await axios.get(
        `https://api.github.com/repos/chxhyxn/TmpCellcraftBoard/git/trees/${dirSha}`
      );
      this.datasets = dirData.tree.filter((item) => item.type === "blob");

      const repo = "chxhyxn/TmpCellcraftBoard";

      for (const dataset of this.datasets) {
        const lastIndex = dataset.path.lastIndexOf(".");
        dataset.type =
          lastIndex !== -1 ? dataset.path.substring(lastIndex + 1) : "";
        dataset.size = formatBytes(dataset.size);
        dataset.url = `https://github.com/${repo}/blob/main/datasets/${dataset.path}`;

        const baseName = dataset.path.substring(0, lastIndex);
        const possibleExtensions = ["ipynb", "md", "txt"];

        for (const ext of possibleExtensions) {
          const exists = await checkFileExists(
            repo,
            `posts/${baseName}.${ext}`
          );
          if (exists) {
            dataset.posturl = `https://github.com/${repo}/blob/main/posts/${baseName}.${ext}`;
            console.log(dataset.posturl);
            break;
          }
        }
      }
      this.sortTable();
    } catch (error) {
      console.error("Error fetching data:", error);
    }
  },
  computed: {
    sorteddatasets() {
      const datasetsCopy = [...this.datasets];
      if (this.sortKey) {
        datasetsCopy.sort((a, b) => {
          const aValue = a[this.sortKey];
          const bValue = b[this.sortKey];
          if (aValue < bValue) return this.sortDirection === "asc" ? -1 : 1;
          if (aValue > bValue) return this.sortDirection === "asc" ? 1 : -1;
          return 0;
        });
      }
      return datasetsCopy;
    },
    totalPages() {
      return Math.ceil(this.filtereddatasets.length / this.pageSize);
    },
    displayeddatasets() {
      const startIndex = (this.currentPage - 1) * this.pageSize;
      const endIndex = startIndex + this.pageSize;
      return this.filtereddatasets.slice(startIndex, endIndex);
    },
    filtereddatasets() {
      // 이전 filtereddatasets 메서드 내용에 태그 검색 기능 추가
      if (this.searchTerm) {
        const searchTermLower = this.searchTerm.toLowerCase();
        return this.sorteddatasets.filter((dataset) => {
          return dataset.path.toLowerCase().includes(searchTermLower);
        });
      } else {
        return this.sorteddatasets;
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
    resetTypeSearch() {
      this.searchType = ""; // 태그 검색어 초기화
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
}
.second-line {
  width: calc(100% - 10px);
  padding: 0px 5px 5px 5px;
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
.upload-button {
  border-radius: 5px;
  padding: 5px 10px;
  font-weight: 500;
  border: 1px solid #ccc;
  background-color: #f0f1fd;
  text-decoration: none;
  color: black;
}
.upload-button:hover {
  cursor: pointer;
  padding: 4px;
  background-color: #fbfcff;
  border: 2px solid #ccc;
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
</style>
