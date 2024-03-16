<template>
  <div class="layout">
    <div class="first-line">
      <div class="header__text">Datasets</div>
      <div class="search">
        <input
          type="text"
          v-model="searchTerm"
          placeholder="Search titles..."
        />
      </div>
      <div class="pagination">
        <button @click="prevPage" :disabled="currentPage <= 1">Prev</button>
        <span>{{ currentPage }}</span>
        <button @click="nextPage" :disabled="currentPage >= totalPages">
          Next
        </button>
      </div>
    </div>
    <table>
      <tbody>
        <tr v-for="item in filteredStudies" :key="item.id">
          <td>
            <div>
              <div class="title-container">
                {{ item.title }}
                <a :href="item.link" target="_blank">
                  <img
                    src="@/assets/floppy-disk.png"
                    alt="Download"
                    class="download-icon"
                  />
                </a>
              </div>
              <div class="description-contatiner">
                {{ item.description }}
              </div>
              <div class="cells-contatiner">{{ item.cells }}</div>
            </div>
          </td>
        </tr>
      </tbody>
    </table>
  </div>
</template>

<script>
export default {
  data() {
    return {
      currentPage: 1,
      pageSize: 5,
      searchTerm: "",
      studies: [
        {
          id: 1,
          title: "Tuck_PAGA3281genes.h5ad",
          cells: "3281 Cells",
          description:
            "This dataset is scRNA-seq data obtained from mouse embryonic stem cells (mESC), and pseudo-time analysis has been conducted using PAGA. The data is structured as a 459 x 3281 Matrix of cells by genes.",
          link: "https://github.com/neocaleb/TENET/raw/master/Data.Tuck/Tuck_PAGA3281genes.h5ad",
        },
        {
          id: 2,
          title: "Tuck_PAGA510genes.h5ad",
          cells: "510 Cells",
          description:
            "This dataset is scRNA-seq data obtained from mouse embryonic stem cells (mESC), and pseudo-time analysis has been conducted using PAGA. The data is structured as a 459 x 510 Matrix of cells by genes.",
          link: "https://github.com/neocaleb/TENET/raw/master/Data.Tuck/Tuck_PAGA510genes.h5ad",
        },
        {
          id: 3,
          title: "pbmc_raw.h5ad",
          cells: "32738 Cells",
          description:
            "The data consists of PBMCs obtained from a healthy donor and is available through 10x Genomics. The data is structured as a 2700 x 32738 Matrix of cells by genes.",
          link: "https://github.com/mindongdong/h5adDatabase/raw/main/pbmc_raw.h5ad",
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
  /* background-color: #c9c9c9; */
  transition: all 0.3s ease;
  border-radius: 15px;
  /* color: #ffffff; */
}

thead th,
td {
  padding: 10px;
  padding-left: 15px;
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

.download-icon {
  margin: 0px 0px;
  width: 33px;
  height: 33px;
}

.title-container {
  font-size: 20px;
  align-items: center;
  display: flex;
  font-weight: 600;
  margin-top: 5px;
}

.description-contatiner {
  font-size: 16px;
  color: #474747;
  margin: 5px 2px;
}

.cells-contatiner {
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
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  line-height: 1rem;
  /* padding-left: 2rem; */
  color: rgba(0, 0, 0, 0.8);
}
</style>
