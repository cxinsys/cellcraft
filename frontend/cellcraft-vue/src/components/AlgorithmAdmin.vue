<template>
  <div>
    <div class="first-line">
      <div class="search">
        <input
          type="text"
          v-model="searchTerm"
          placeholder="Search by name..."
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
      <div></div>
      <label class="upload-button">
        ⇪ UPLOAD NEW ALGORITHM
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
          <th>setting</th>
          <th></th>
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
            <button @click="editAlgorithmSetting(algorithm)">edit</button>
          </td>
          <td>
            <button @click="deleteAlgorithm(algorithm)">
              delete algorithm
            </button>
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
    <div class="layout-edit-setting" v-if="this.editingAlgorithm">
      <!-- <div class="layout-edit-setting" v-if="true"> -->
      <div id="layout">
        <div class="input-layout">
          <div class="input-title">input node</div>
          <div
            v-for="inputDescription in editingAlgorithm['inputDescriptions']"
            :key="inputDescription.id"
            class="input-description"
            :class="{ linked: inputDescription.linked }"
          >
            <input
              type="text"
              class="description-id"
              :value="inputDescription.id"
              :v-model="inputDescription.value"
            />
            <input
              type="text"
              class="description-tooltip"
              :value="inputDescription.description"
              :v-model="inputDescription.value"
            />
          </div>
        </div>
        <div class="algorithm-layout">
          <img
            class="algorithm-logo"
            :src="this.editingAlgorithm.logo"
            alt="Fast Tenet"
          />
          <div class="algorithm-parts">
            <div
              v-for="part in editingAlgorithm['algorithmParts']"
              :key="part.id"
            >
              <div class="part-title">{{ part.id }}</div>
              <div v-for="parameter in part.parameters" :key="parameter.id">
                <div class="parameters">
                  <input
                    type="text"
                    class="description-id"
                    :value="parameter.id"
                    :v-model="parameter.value"
                  />
                  <input
                    type="text"
                    :placeholder="parameter.default"
                    class="parameter__textInput"
                    :v-model="parameter.value"
                    :class="{ 'red-text': !parameter.value }"
                    :disabled="parameter.disabled"
                  />
                  <label
                    class="switch"
                    :v-model="parameter.value"
                    v-on:click="parameter.disabled = !parameter.disabled"
                  >
                    <input type="checkbox" />
                    <span class="slider_button round"></span>
                  </label>
                  <input
                    type="text"
                    class="description-tooltip"
                    :value="parameter.description"
                    :v-model="parameter.value"
                  />
                </div>
              </div>
            </div>
          </div>
        </div>
        <div class="output-layout">
          <div class="output-title">output node</div>
          <div
            v-for="outputDescription in editingAlgorithm['outputDescriptions']"
            :key="outputDescription.id"
            class="output-description"
            :class="{ linked: outputDescription.linked }"
          >
            <input
              type="text"
              class="description-id"
              :value="outputDescription.id"
              :v-model="outputDescription.value"
            />
            <input
              type="text"
              class="description-tooltip"
              :value="outputDescription.description"
              :v-model="outputDescription.value"
            />
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  data() {
    return {
      editingAlgorithm: null,
      algorithms: [
        {
          no: 0,
          userId: "cislab",
          name: "FastTenet",
          time: "2023-06-30",
          description: "temporary description",
          logo: require("@/assets/fasttenet.png"),
          inputDescriptions: [
            {
              id: "Expression Data (dpath_exp_data)",
              description:
                "This is the path to the expression data file. It contains gene expression values for different samples or cells.",
            },
            {
              id: "Trajectory Data (dpath_trj_data)",
              description:
                "This is the path to the trajectory data file. It represents the progression or trajectory of cells or samples in a biological process.",
            },
            {
              id: "Branch Data (dpath_branch_data)",
              description:
                "This is the path to the branch or cell select data file. It specifies the branches or groups of cells in the trajectory.",
            },
            {
              id: "TF Data (dpath_tf_data)",
              description:
                "This is the path to the transcription factor (TF) data file. It contains information about the transcription factors and their regulatory relationships.",
            },
          ],
          algorithmParts: [
            {
              id: "fasttenet parameters",
              parameters: [
                {
                  id: "dpath_exp_data",
                  description: "expression data path",
                  default: "",
                  value: "",
                  required: true,
                  disabled: true,
                },
                {
                  id: "dpath_trj_data",
                  description: "trajectory data path",
                  default: "",
                  value: "",
                  required: true,
                  disabled: true,
                },
                {
                  id: "dpath_branch_data",
                  description: "branch(cell select) data path",
                  default: "",
                  value: "",
                  required: true,
                  disabled: true,
                },
                {
                  id: "dpath_tf_data",
                  description: "tf data path",
                  default: "",
                  value: "",
                  required: true,
                  disabled: true,
                },
                {
                  id: "spath_result_matrix",
                  description: "spath_result_matrix",
                  default: "None",
                  value: "None",
                  required: false,
                  disabled: true,
                },
                {
                  id: "make_binary",
                  description:
                    "if True, make binary expression and node name file",
                  default: "False",
                  value: "False",
                  required: false,
                  disabled: false,
                },
              ],
            },
            {
              id: "worker run parameters",
              parameters: [
                {
                  id: "device",
                  description: "cpu or gpu",
                  default: "gpu",
                  value: "gpu",
                  required: false,
                  disabled: true,
                },
                {
                  id: "device_ids",
                  description: "[0](cpu) or [list of whole gpu devices](gpu)",
                  default: "[0, 1, 2, 3, 4, 5, 6, 7]",
                  value: "[0, 1, 2, 3, 4, 5, 6, 7]",
                  required: false,
                  disabled: true,
                },
                {
                  id: "batch_size",
                  description: "batch size",
                  default: "2 ** 16",
                  value: "2 ** 16",
                  required: true,
                  disabled: true,
                },
                {
                  id: "kp",
                  description: "kernel percentail",
                  default: "0.5",
                  value: "0.5",
                  required: false,
                  disabled: false,
                },
                {
                  id: "percentile",
                  description: "data crop percentile",
                  default: "0",
                  value: "0",
                  required: false,
                  disabled: false,
                },
                {
                  id: "win_length",
                  description: "smoothe func window length parameter",
                  default: "10",
                  value: "10",
                  required: false,
                  disabled: false,
                },
                {
                  id: "polyorder",
                  description: "smoothe func polyorder parameter",
                  default: "3",
                  value: "3",
                  required: false,
                  disabled: false,
                },
              ],
            },
          ],
          outputDescriptions: [
            {
              id: "Result Matrix (spath_result_matrix)",
              description:
                "This is the path to the result matrix data file. It stores the results of the FastTENET algorithm, which includes the inferred regulatory relationships between genes.",
            },
            {
              id: "Binary Expression and Node Name Files (make_binary)",
              description:
                "If the make_binary parameter is set to True, FastTENET generates binary expression and node name files. The binary expression file contains a binary representation of the gene expression data, while the node name file contains the names or identifiers of the genes.",
            },
            {
              id: "GRN (Gene Regulatory Network) Files:",
              description:
                "The GRN files are generated by running the make_grn.py script. They include the inferred gene regulatory network based on the result matrix, node name file, and TF data. The output includes a file with a '.sif' extension, which represents the network structure, and a file with a '.outdegrees.txt' extension, which contains information about the outdegrees (number of outgoing connections) of each gene in the network.",
            },
          ],
        },
      ],
      sortKey: "no",
      sortDirection: "dsc", // Set initial sort direction to 'asc'
      pageSize: 20,
      currentPage: 1,
      searchTerm: "",
    };
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
      return Math.ceil(this.filteredalgorithms.length / this.pageSize);
    },
    displayedalgorithms() {
      const startIndex = (this.currentPage - 1) * this.pageSize;
      const endIndex = startIndex + this.pageSize;
      return this.filteredalgorithms.slice(startIndex, endIndex);
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
    resetTagSearch() {
      this.searchTag = ""; // 태그 검색어 초기화
    },
    deleteAlgorithm(algorithm) {
      const index = this.algorithms.findIndex((j) => j.id === algorithm.id);
      if (index !== -1) {
        this.algorithms.splice(index, 1);
      }
    },
    updatePage() {
      this.currentPage = 1; // Reset to first page when page size changes
    },
    getLinkUrl(algorithm) {
      console.log(algorithm);
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
    editAlgorithmSetting(algorithm) {
      this.editingAlgorithm = algorithm;
      console.log(this.editingAlgorithm);
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

input:checked + .slider_button {
  background-color: #53b2ff;
}

input:focus + .slider_button {
  box-shadow: 0 0 1px #53b2ff;
}

input:checked + .slider_button:before {
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
</style>
