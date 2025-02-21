<template>
  <div class="layout">
    <div class="first-line">
      <div class="header__text">Datasets</div>
      <div class="search">
        <input type="text" v-model="searchTerm" placeholder="Search titles..." />
      </div>
      <div class="pagination">
        <button @click="prevPage" :disabled="currentPage <= 1">Prev</button>
        <span>{{ currentPage }}</span>
        <button @click="nextPage" :disabled="currentPage >= totalPages">Next</button>
      </div>
    </div>
    <table>
      <tbody>
        <tr v-for="item in paginatedStudies" :key="item.id">
          <td>
            <div>
              <div class="title-container">
                {{ item.title }}
                <button @click="downloadFile(item.title)" class="download-button">
                  <img src="@/assets/floppy-disk.png" alt="Download" class="download-icon" />
                </button>
              </div>
              <div class="description-container">
                {{ item.description }}
              </div>
              <div v-if="item.cells !== ''" class="cells-container">{{ item.cells }}</div>
            </div>
          </td>
        </tr>
      </tbody>
    </table>
  </div>
</template>

<script>
import { getTutorialFileDownload } from "@/api/index"; // API 함수 가져오기

export default {
  data() {
    return {
      currentPage: 1,
      pageSize: 5,
      searchTerm: "",
      studies: [
        {
          id: 1,
          title: "pbmc_light_1000.h5ad",
          cells: "32738 Cells",
          description:
            "A tutorial dataset in CellCraft containing 1,000 PBMC cells from a healthy donor. The data is structured as a 1000 × 32,738 matrix of cells by genes and is designed for testing built-in plugins.",
        },
        {
          id: 2,
          title: "unique_genes.txt",
          cells: "",
          description:
            "A tutorial dataset in CellCraft containing 971 unique differentially expressed genes (DEGs), aggregated from the top 150 DEGs of each group.",
        },
        {
          id: 3,
          title: "human_TFs.txt",
          cells: "",
          description:
            "A tutorial dataset in CellCraft containing a list of transcription factors (TFs) related to the regulation of transcription and sequence-specific DNA binding.",
        },
      ],
    };
  },
  computed: {
    totalPages() {
      return Math.ceil(this.filteredStudies.length / this.pageSize);
    },
    filteredStudies() {
      return this.studies.filter((study) =>
        study.title.toLowerCase().includes(this.searchTerm.toLowerCase())
      );
    },
    paginatedStudies() {
      const start = (this.currentPage - 1) * this.pageSize;
      const end = start + this.pageSize;
      return this.filteredStudies.slice(start, end);
    },
  },
  methods: {
    nextPage() {
      if (this.currentPage < this.totalPages) this.currentPage++;
    },
    prevPage() {
      if (this.currentPage > 1) this.currentPage--;
    },
    async downloadFile(filename) {
      try {
        const response = await getTutorialFileDownload(filename); // API 호출
        const blob = new Blob([response.data]); // Blob 생성
        const link = document.createElement("a");
        link.href = URL.createObjectURL(blob);
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
      } catch (error) {
        console.error("Download failed:", error);
      }
    },
  },
};
</script>

<style scoped>
.layout {
  padding: 10px 30px;
}

table {
  width: 100%;
  height: 100%;
  border-collapse: separate;
  border-spacing: 5px;
  transition: all 0.3s ease;
  border-radius: 15px;
}

thead th,
td {
  padding: 10px;
  padding-left: 15px;
  text-align: left;
  border-radius: 10px;
  border: 1px solid #a8a8a8;
}

th {
  text-transform: capitalize;
  background-color: #474747;
  color: #ffffff;
}

th:hover {
  background-color: #616161;
}

button {
  margin-right: 10px;
  color: black;
  padding: 5px;
  border-radius: 10px;
  background-color: #eaecff;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
}

button:disabled {
  color: #ccc;
}

.first-line {
  height: 40px;
  margin-bottom: 10px;
  width: calc(100% - 10px);
  padding: 5px;
  display: flex;
  justify-content: space-between;
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

.pagination {
  display: flex;
  justify-content: center;
  margin: 20px 0px;
}

.pagination button {
  margin: -5px 10px 0px 10px;
}

.download-icon {
  margin: 0px 0px;
  width: 33px;
  height: 33px;
  cursor: pointer;
}

.download-button {
  background: none;
  border: none;
  cursor: pointer;
  padding: 0;
  display: inline-block;
}

.title-container {
  font-size: 20px;
  align-items: center;
  display: flex;
  font-weight: 600;
  margin-top: 5px;
}

.description-container {
  font-size: 16px;
  color: #474747;
  margin: 5px 2px;
}

.cells-container {
  font-size: 14px;
  font-weight: 600;
  padding: 6px 10px;
  margin: 2px 0px;
  border-radius: 10px;
  color: white;
  background: rgb(40, 84, 197);
  display: inline-block;
}

.header__text {
  font-family: "Montserrat", sans-serif;
  font-weight: 600;
  font-size: 2rem;
  color: rgba(0, 0, 0, 0.8);
}
</style>
