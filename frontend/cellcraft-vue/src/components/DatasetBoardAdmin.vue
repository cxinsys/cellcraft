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
      <div class="search">
        <input
          type="text"
          v-model="searchTag"
          placeholder="Search by tags..."
        />
        <img
          class="reset-button"
          src="@/assets/reset.png"
          alt="reset"
          @click="resetTagSearch"
        />
      </div>
      <label class="upload-button">
        ⇪ UPLOAD NEW DATASET
        <input
          type="file"
          ref="fileInput"
          style="display: none"
          @change="uploadFile"
        />
      </label>
    </div>
    <table>
      <thead>
        <tr>
          <th @click="sortTable('no')">
            No. <span class="sort-icon">{{ sortIcon("no") }}</span>
          </th>
          <th>tags</th>
          <th @click="sortTable('title')">
            title <span class="sort-icon">{{ sortIcon("title") }}</span>
          </th>
          <th>description</th>
          <th @click="sortTable('userId')">
            Uploader id <span class="sort-icon">{{ sortIcon("userId") }}</span>
          </th>
          <th @click="sortTable('time')">
            uploaded date <span class="sort-icon">{{ sortIcon("time") }}</span>
          </th>
          <th>view post</th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="dataset in displayeddatasets" :key="dataset.no">
          <td>{{ dataset.no }}</td>
          <td class="td_tag">
            <span
              v-for="tag in dataset.tags"
              :key="tag"
              :class="[getTagClass(tag), 'tags']"
              >{{ tag }}</span
            >
          </td>
          <td>{{ dataset.title }}</td>
          <td class="description-cell">{{ dataset.description }}</td>
          <td>{{ dataset.userId }}</td>
          <td>{{ dataset.time }}</td>
          <td><a :href="getLinkUrl(dataset)">Link</a></td>
          <td>
            <button @click="deleteDataset(dataset)">delete dataset</button>
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
      datasets: [
        {
          no: 1,
          userId: "johndoe",
          title: "dataset1",
          time: "2023-06-29",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 2,
          userId: "janesmith",
          title: "dataset1",
          time: "2023-06-30",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 3,
          userId: "bobjohnson",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 4,
          userId: "alicebrown",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 5,
          userId: "samwilson",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 6,
          userId: "emilydavis",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 7,
          userId: "michaelwilson",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 8,
          userId: "oliviajohnson",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 9,
          userId: "sophiamiller",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 10,
          userId: "williamanderson",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["xml"],
          description: "temporary description",
        },
        {
          no: 11,
          userId: "benjamingarcia",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv"],
          description: "temporary description",
        },
        {
          no: 12,
          userId: "avamartinez",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv"],
          description: "temporary description",
        },
        {
          no: 13,
          userId: "miathompson",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv"],
          description: "temporary description",
        },
        {
          no: 14,
          userId: "ethanlopez",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv"],
          description: "temporary description",
        },
        {
          no: 15,
          userId: "jameswilson",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv"],
          description: "temporary description",
        },
        {
          no: 16,
          userId: "liamwhite",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv"],
          description: "temporary description",
        },
        {
          no: 17,
          userId: "sophiabrown",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv"],
          description: "temporary description",
        },
        {
          no: 18,
          userId: "charlottedavis",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv"],
          description: "temporary description",
        },
        {
          no: 19,
          userId: "alexanderjohnson",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv"],
          description: "temporary description",
        },
        {
          no: 20,
          userId: "emmamiller",
          title: "dataset1",
          time: "2023-07-01",
          tags: ["csv", "xml"],
          description: "temporary description",
        },
      ],
      sortKey: "no",
      sortDirection: "dsc", // Set initial sort direction to 'asc'
      pageSize: 20,
      currentPage: 1,
      searchTerm: "",
      searchTag: "",
    };
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
      if (this.searchTerm || this.searchTag) {
        const searchTermLower = this.searchTerm.toLowerCase();
        const searchTagLower = this.searchTag.toLowerCase();
        return this.sorteddatasets.filter((dataset) => {
          const titleMatch = dataset.title
            .toLowerCase()
            .includes(searchTermLower);
          const tagMatch = dataset.tags.some((tag) =>
            tag.toLowerCase().includes(searchTagLower)
          );
          return titleMatch && tagMatch;
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
    getTagClass(tag) {
      if (tag === "xml") {
        return "tag-red";
      } else if (tag === "csv") {
        return "tag-green";
      } else {
        return "tag-blue";
      }
    },
    resetSearch() {
      this.searchTerm = "";
    },
    resetTagSearch() {
      this.searchTag = ""; // 태그 검색어 초기화
    },
    deleteDataset(dataset) {
      const index = this.datasets.findIndex((j) => j.id === dataset.id);
      if (index !== -1) {
        this.datasets.splice(index, 1);
      }
    },
    updatePage() {
      this.currentPage = 1; // Reset to first page when page size changes
    },
    getLinkUrl(dataset) {
      console.log(dataset);
      // 데이터셋에 따른 링크 URL을 반환하는 메서드
      // 예시로 데이터셋의 no 값을 링크 URL에 사용
      return `https://singlecell.broadinstitute.org/single_cell/study/SCP2221/primary-nasal-viral-infection-rewires-the-tissue-scale-memory-response-rechallenge-rm`;
    },
    uploadFile(event) {
      const file = event.target.files[0]; // Get the selected file
      // Perform the necessary operations with the file, such as uploading to a server or processing it
      // You can access the file using the 'file' variable
      console.log(file);
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
  justify-content: space-between;
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
  color: black;
  padding: 2px;
  left: 10px;
  border-radius: 5px;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  margin-bottom: 5px;
}
.upload-button {
  border-radius: 5px;
  padding: 5px 10px;
  font-weight: 500;
  border: 1px solid #ccc;
  background-color: #f0f1fd;
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

.description-cell {
  max-width: 200px; /* 텍스트의 최대 너비 설정 */
}

.pagination {
  display: flex;
  justify-content: center;
  margin: 20px 0px;
}

.pagination button {
  margin: -5px 10px 0px 10px;
}
.td_tag {
  overflow: scroll;
  text-overflow: clip;
  max-width: 60px;
}
.tags {
  color: white;
  font-weight: bold;
  border-radius: 5px;
  padding: 2px 8px;
  margin: 0px 1px;
}
.tag-red {
  background-color: rgb(255, 57, 57);
}

.tag-green {
  background-color: rgb(6, 152, 6);
}

.tag-blue {
  background-color: blue;
}
</style>
