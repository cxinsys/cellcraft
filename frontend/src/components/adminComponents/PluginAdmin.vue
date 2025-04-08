<template>
  <div class="layout_admin">
    <div class="first-line">
      <div class="header__text">Algorithms</div>
      <div class="search">
        <input type="text" v-model="searchTerm" placeholder="Search by name..." />
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
      <div></div>
      <label class="upload-button">
        ⇪ upload new algorithm
        <input type="file" ref="fileInput" style="display: none" @change="uploadFile" />
      </label>
    </div> -->

    <table>
      <thead>
        <tr>
          <th @click="sortTable('no')" style="width: 70px">
            No. <span class="sort-icon">{{ sortIcon("no") }}</span>
          </th>
          <th @click="sortTable('name')">
            name <span class="sort-icon">{{ sortIcon("name") }}</span>
          </th>
          <th>description</th>
          <th @click="sortTable('userId')">
            uploader id <span class="sort-icon">{{ sortIcon("userId") }}</span>
          </th>
          <th @click="sortTable('time')">
            uploaded date <span class="sort-icon">{{ sortIcon("time") }}</span>
          </th>
          <th>Actions</th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="algorithm in displayedalgorithms" :key="algorithm.no">
          <td>{{ algorithm.no }}</td>
          <td>{{ algorithm.name }}</td>
          <td class="description-cell">{{ algorithm.description }}</td>
          <td>{{ algorithm.userId }}</td>
          <td>{{ algorithm.time }}</td>
          <td>
            <button @click="deleteAlgorithm(algorithm)" class="table-button delete">
              Delete
            </button>
          </td>
        </tr>
      </tbody>
    </table>
    <div class="pagination">
      <button :disabled="currentPage === 1" @click="currentPage--">Prev</button>
      <span>{{ currentPage }} / {{ totalPages }}</span>
      <button :disabled="currentPage === totalPages" @click="currentPage++">
        Next
      </button>
    </div>
    <div class="layout-edit-setting" v-if="this.editingAlgorithm">
      <!-- <div class="layout-edit-setting" v-if="true"> -->
      <div id="layout">
        <div class="input-layout">
          <div class="input-title">input node</div>
          <div v-for="inputDescription in editingAlgorithm['inputDescriptions']" :key="inputDescription.id"
            class="input-description" :class="{ linked: inputDescription.linked }">
            <input type="text" class="description-id" :value="inputDescription.id" :v-model="inputDescription.value" />
            <input type="text" class="description-tooltip" :value="inputDescription.description"
              :v-model="inputDescription.value" />
          </div>
        </div>
        <div class="algorithm-layout">
          <img class="algorithm-logo" :src="this.editingAlgorithm.logo" alt="Fast Tenet" />
          <div class="algorithm-parts">
            <div v-for="part in editingAlgorithm['algorithmParts']" :key="part.id">
              <div class="part-title">{{ part.id }}</div>
              <div v-for="parameter in part.parameters" :key="parameter.id">
                <div class="parameters">
                  <input type="text" class="description-id" :value="parameter.id" :v-model="parameter.value" />
                  <input type="text" :placeholder="parameter.default" class="parameter__textInput"
                    :v-model="parameter.value" :class="{ 'red-text': !parameter.value }"
                    :disabled="parameter.disabled" />
                  <label class="switch" :v-model="parameter.value"
                    v-on:click="parameter.disabled = !parameter.disabled">
                    <input type="checkbox" />
                    <span class="slider_button round"></span>
                  </label>
                  <input type="text" class="description-tooltip" :value="parameter.description"
                    :v-model="parameter.value" />
                </div>
              </div>
            </div>
          </div>
        </div>
        <div class="output-layout">
          <div class="output-title">output node</div>
          <div v-for="outputDescription in editingAlgorithm['outputDescriptions']" :key="outputDescription.id"
            class="output-description" :class="{ linked: outputDescription.linked }">
            <input type="text" class="description-id" :value="outputDescription.id"
              :v-model="outputDescription.value" />
            <input type="text" class="description-tooltip" :value="outputDescription.description"
              :v-model="outputDescription.value" />
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import { getFilteredPlugins, getPluginsCount } from '@/api';

export default {
  data() {
    return {
      editingAlgorithm: null,
      algorithms: [],
      sortKey: "id",
      sortDirection: "dsc",
      pageSize: 15,
      currentPage: 1,
      searchTerm: "",
      totalCount: 0
    };
  },
  async created() {
    await this.fetchPlugins();
  },
  computed: {
    sortedAlgorithms() {
      const algorithmsCopy = [...this.algorithms];
      if (this.sortKey) {
        algorithmsCopy.sort((a, b) => {
          const aValue = a[this.sortKey];
          const bValue = b[this.sortKey];
          if (aValue < bValue) return this.sortDirection === "asc" ? -1 : 1;
          if (aValue > bValue) return this.sortDirection === "asc" ? 1 : -1;
          return 0;
        });
      }
      return algorithmsCopy;
    },
    totalPages() {
      return Math.ceil(this.totalCount / this.pageSize);
    },
    displayedalgorithms() {
      return this.sortedAlgorithms;
    },
    filteredalgorithms() {
      if (this.searchTerm) {
        const searchTermLower = this.searchTerm.toLowerCase();
        return this.sortedAlgorithms.filter((algorithm) =>
          algorithm.name.toLowerCase().includes(searchTermLower)
        );
      } else {
        return this.sortedAlgorithms;
      }
    },
    showEditAlgorithmView() {
      if (this.editingAlgorithm == []) {
        return true;
      } else {
        return false;
      }
    },
  },
  methods: {
    async fetchPlugins() {
      try {
        const conditions = {
          amount: this.pageSize,
          page_num: this.currentPage,
          sort: this.sortKey,
          order: this.sortDirection === 'asc' ? 'asc' : 'desc',
          searchTerm: this.searchTerm
        };

        const [pluginsResponse, countResponse] = await Promise.all([
          getFilteredPlugins(conditions),
          getPluginsCount()
        ]);

        this.algorithms = pluginsResponse.data.map((plugin, index) => ({
          no: index + 1,
          id: plugin.id,
          userId: plugin.author,
          name: plugin.name,
          time: new Date(plugin.created_at).toLocaleDateString(),
          description: plugin.description || 'No description available',
          logo: plugin.logo_url || require("@/assets/fasttenet.png"),
          inputDescriptions: plugin.input_descriptions || [],
          algorithmParts: plugin.algorithm_parts || [],
          outputDescriptions: plugin.output_descriptions || []
        }));

        this.totalCount = countResponse.data;
      } catch (error) {
        console.error('Error fetching plugins:', error);
      }
    },
    async sortTable(key) {
      if (this.sortKey === key) {
        this.sortDirection = this.sortDirection === "asc" ? "desc" : "asc";
      } else {
        this.sortKey = key;
        this.sortDirection = "asc";
      }
      await this.fetchPlugins();
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
    resetTagSearch() {
      this.searchTag = ""; // 태그 검색어 초기화
    },
    async deleteAlgorithm(algorithm) {
      if (confirm(`Are you sure you want to delete algorithm ${algorithm.name}?`)) {
        try {
          // TODO: Implement delete API call
          await this.fetchPlugins();
        } catch (error) {
          console.error('Error deleting plugin:', error);
        }
      }
    },
    async updatePage() {
      this.currentPage = 1;
      await this.fetchPlugins();
    },
    getLinkUrl() {
      // 데이터셋에 따른 링크 URL을 반환하는 메서드
      // 예시로 데이터셋의 no 값을 링크 URL에 사용
      return `https://singlecell.broadinstitute.org/single_cell/study/SCP2221/primary-nasal-viral-infection-rewires-the-tissue-scale-memory-response-rechallenge-rm`;
    },
    uploadFile() {
      // Perform the necessary operations with the file, such as uploading to a server or processing it
      // You can access the file using the 'file' variable
    },
    editAlgorithmSetting(algorithm) {
      this.editingAlgorithm = algorithm;
    },
    async prevPage() {
      if (this.currentPage > 1) {
        this.currentPage--;
        await this.fetchPlugins();
      }
    },
    async nextPage() {
      if (this.currentPage < this.totalPages) {
        this.currentPage++;
        await this.fetchPlugins();
      }
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
  width: 8%;
}

/* No. */
th:nth-child(2) {
  width: 20%;
}

/* name */
th:nth-child(3) {
  width: 30%;
}

/* description */
th:nth-child(4) {
  width: 15%;
}

/* uploader id */
th:nth-child(5) {
  width: 15%;
}

/* uploaded date */
th:nth-child(6) {
  width: 12%;
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

  .description-cell {
    max-height: 100px;
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

a {
  color: #f0f1fb;
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
  justify-content: space-between;
  flex-direction: row;
}

.layout-edit-setting {
  position: fixed;
  top: calc(50% - 400px);
  left: calc(50% - 440px);
  width: 880px;
  height: 800px;
  background-color: blue;
  border-radius: 15px;
}

#layout {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: row;
}

.input-layout,
.output-layout {
  width: 25%;
  height: 85%;
  margin-left: 1rem;
  display: flex;
  align-items: center;
  /* justify-content: center; */
  flex-direction: column;
  padding: 1rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
}

.input-title,
.output-title {
  text-transform: capitalize;
  font-weight: bold;
  font-size: larger;
  color: #494949;
  margin: 50px 0px 20px 0px;
}

.input-description,
.output-description {
  width: 100%;
  height: 5rem;
  text-align: center;
  justify-content: center;
  display: flex;
  flex-direction: column;
  color: #353535;
  background-color: rgb(224, 224, 224);
  border-radius: 1rem;
  margin: 0.5rem 0rem;
  border: 2px solid #e7e7e7;
  opacity: 0.9;
}

.input-description.linked,
.output-description.linked {
  opacity: 1;
  background-color: rgb(202, 214, 255);
  border: 2px solid #ecebff;
}

.algorithm-layout {
  width: 50%;
  height: 95%;
  /* background-color: blue; */
  margin: 1rem;
  display: flex;
  align-items: center;
  /* justify-content: center; */
  flex-direction: column;
  padding: 1rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
  overflow-y: scroll;
}

.algorithm-logo {
  width: 100%;
  height: 10%;
  top: 1rem;
  object-fit: contain;
}

.algorithm-parts {
  align-content: start;
  width: 100%;
  height: 90%;
}

.part-title {
  text-transform: capitalize;
  font-weight: bold;
  font-size: large;
  color: #353535;

  margin: 20px 0px 10px 0px;
}

.parameters {
  display: flex;
  direction: row;
  justify-content: space-between;
  width: 80%;
  color: #353535;
  margin: 2px;
  flex-wrap: wrap;
}

.parameter__textInput {
  width: 80%;
  color: black;
  padding: 5px;
  /* right: 10px; */
  border-radius: 3px;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  margin-bottom: 0px;
}

.parameter__textInput:disabled {
  background-color: lightgray;
}

.parameter__textInput:disabled::placeholder {
  color: black;
}

.output-layout {
  margin-left: 0rem;
  margin-right: 1rem;
}

.description-id {
  width: 100%;
}

.description-tooltip {
  width: 100%;
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
  background-color: #5a5a5a;
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

.description-cell {
  max-height: 200px;
  overflow-y: auto;
}

.pagination {
  display: flex;
  justify-content: center;
  margin: 20px 0px;
}

.pagination button {
  margin: -5px 10px 0px 10px;
}

.switch {
  position: absolute;
  display: inline-block;
  width: 50px;
  height: 24px;
  margin-top: -5px;
  right: 10px;
}

/* Hide default HTML checkbox */
.switch input {
  opacity: 0;
  width: 0;
  height: 0;
}

/* The slider_button */
.slider_button {
  position: absolute;
  cursor: pointer;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background-color: #ccc;
  -webkit-transition: 0.4s;
  transition: 0.4s;
}

.slider_button:before {
  position: absolute;
  content: "";
  height: 20px;
  width: 20px;
  left: 2px;
  bottom: 2px;
  background-color: white;
  -webkit-transition: 0.4s;
  transition: 0.4s;
}

input:checked+.slider_button {
  background-color: #53b2ff;
}

input:focus+.slider_button {
  box-shadow: 0 0 1px #53b2ff;
}

input:checked+.slider_button:before {
  -webkit-transform: translateX(26px);
  -ms-transform: translateX(26px);
  transform: translateX(26px);
}

/* Rounded slider_buttons */
.slider_button.round {
  border-radius: 34px;
}

.slider_button.round:before {
  border-radius: 50%;
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
