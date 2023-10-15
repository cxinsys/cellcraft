<template>
  <div id="layout">
    <div class="input-layout">
      <div class="input-title">input node</div>
      <div
        v-for="inputDescription in inputDescriptions"
        :key="inputDescription.id"
        class="input-description"
        :class="{ linked: inputDescription.linked }"
      >
        <span class="description-id">
          {{ inputDescription.id }}
          <span class="description-tooltip">{{
            inputDescription.description
          }}</span>
        </span>
      </div>
    </div>
    <div class="algorithm-layout">
      <img class="algorithm-logo" :src="this.logo" alt="Fast Tenet" />
      <div class="algorithm-parts">
        <div>
          <div class="part-title">Select Annotation Column</div>
          <div class="parameters">
            <span class="parameter-id"
              >File
              <span class="parameter-tooltip">Selected File Name</span>
            </span>
            <input
              type="text"
              :placeholder="this.current_file"
              class="parameter__textInput"
              :v-model="this.current_file"
              :class="{ 'red-text': !this.current_file }"
              :disabled="true"
            />
          </div>
          <div class="parameters" style="height: 1.8rem">
            <span class="parameter-id"
              >Annotation Column
              <span class="parameter-tooltip">Selected File Name</span>
            </span>
            <select
              class="parameter__dropdown"
              v-model="annotationColumn"
              @change="selectColumns($event)"
            >
              <option class="parameter__menu" disabled value="">
                Select annotation
              </option>
              <option
                v-for="(column, index) in annotations"
                :key="index"
                :value="column"
              >
                {{ column }}
              </option>
            </select>
          </div>
          <div class="parameters" style="height: 1.8rem">
            <span class="parameter-id">
              Pseudotime Column
              <span class="parameter-tooltip">Selected File Name</span>
            </span>
            <select
              class="parameter__dropdown"
              v-model="pseudotimeColumn"
              @change="selectColumns($event)"
            >
              <option class="parameter__menu" disabled value="">
                Select Pseudotime
              </option>
              <option
                class="parameter__menu"
                v-for="(column, index) in pseudotime"
                :key="index"
                :value="column"
              >
                {{ column }}
              </option>
            </select>
          </div>
          <div class="parameters" style="height: 1.8rem">
            <span class="parameter-id"
              >Clusters
              <span class="parameter-tooltip">Selected File Name</span>
            </span>
            <!-- <select
              class="parameter__dropdown"
              v-model="clusterColumn"
              @change="selectColumns($event)"
            >
              <option class="parameter__menu" disabled value="">
                Select Cluster
              </option>
              <option
                class="parameter__menu"
                v-for="(column, index) in clusters"
                :key="index"
                :value="column"
              >
                {{ column }}
              </option>
            </select> -->
            <div class="parameter__checkbox">
              <div v-for="(column, index) in clusters" :key="index">
                <input
                  type="checkbox"
                  class="parameter__check"
                  :value="column"
                  v-model="clusterColumn"
                />
                {{ column }}
              </div>
            </div>
          </div>
        </div>
        <div v-for="part in algorithmParts" :key="part.id">
          <div class="part-title">{{ part.id }}</div>
          <div v-for="parameter in part.parameters" :key="parameter.id">
            <div class="parameters">
              <span class="parameter-id"
                >{{ parameter.id }}
                <span class="parameter-tooltip">{{
                  parameter.description
                }}</span>
              </span>
              <input
                type="text"
                :placeholder="parameter.default"
                class="parameter__textInput"
                :v-model="parameter.value"
                :class="{ 'red-text': !parameter.value }"
                :disabled="parameter.disabled"
              />
            </div>
          </div>
        </div>
      </div>
      <!-- run button -->
      <div>
        <div class="options__item">
          <label class="form__button--setup" v-bind:class="{ apply: activate }">
            Set up
            <input
              @click="setOption"
              class="form__input"
              type="submit"
              value="업로드"
            />
          </label>
        </div>
      </div>
    </div>
    <div class="output-layout">
      <div class="output-title">output node</div>
      <div
        v-for="outputDescription in outputDescriptions"
        :key="outputDescription.id"
        class="output-description"
        :class="{ linked: outputDescription.linked }"
      >
        <span class="description-id">
          {{ outputDescription.id }}
          <span class="description-tooltip">{{
            outputDescription.description
          }}</span>
        </span>
      </div>
    </div>
  </div>
</template>

<script>
import { getColumns, getClusters, setupAlgorithm } from "@/api/index";

export default {
  data() {
    return {
      logo: require("@/assets/fasttenet.png"),
      current_file: null,
      inputDescriptions: [
        {
          id: "Expression Data",
          description:
            "This is the path to the expression data file. It contains gene expression values for different samples or cells.",
          linked: false,
        },
        {
          id: "Trajectory Data",
          description:
            "This is the path to the trajectory data file. It represents the progression or trajectory of cells or samples in a biological process.",
          linked: false,
        },
        {
          id: "Branch Data",
          description:
            "This is the path to the branch or cell select data file. It specifies the branches or groups of cells in the trajectory.",
          linked: false,
        },
        // {
        //   id: "TF Data",
        //   description:
        //     "This is the path to the transcription factor (TF) data file. It contains information about the transcription factors and their regulatory relationships.",
        //   linked: false,
        // },
      ],
      algorithmParts: [
        {
          id: "fasttenet parameters",
          parameters: [
            // {
            //   id: "dpath_exp_data",
            //   description: "expression data path",
            //   default: "",
            //   value: "",
            //   required: true,
            //   disabled: true,
            // },
            // {
            //   id: "dpath_trj_data",
            //   description: "trajectory data path",
            //   default: "",
            //   value: "",
            //   required: true,
            //   disabled: true,
            // },
            // {
            //   id: "dpath_branch_data",
            //   description: "branch(cell select) data path",
            //   default: "",
            //   value: "",
            //   required: true,
            //   disabled: true,
            // },
            // {
            //   id: "dpath_tf_data",
            //   description: "tf data path",
            //   default: "",
            //   value: "",
            //   required: true,
            //   disabled: true,
            // },
            // {
            //   id: "spath_result_matrix",
            //   description: "spath_result_matrix",
            //   default: "None",
            //   value: "None",
            //   required: false,
            //   disabled: true,
            // },
            {
              id: "make_binary",
              description: "if True, make binary expression and node name file",
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
              id: "Device",
              description: "cpu or gpu",
              default: "gpu",
              value: "gpu",
              required: false,
              disabled: false,
            },
            {
              id: "Device Ids",
              description: "[0](cpu) or [list of whole gpu devices](gpu)",
              default: "[0, 1, 2, 3, 4, 5, 6, 7]",
              value: [0, 1, 2, 3, 4, 5, 6, 7],
              required: false,
              disabled: false,
            },
            {
              id: "Batch Size",
              description: "batch size",
              default: "2 ** 16",
              value: 2 ** 16,
              required: true,
              disabled: false,
            },
            {
              id: "Kp",
              description: "kernel percentail",
              default: "0.5",
              value: 0.5,
              required: false,
              disabled: false,
            },
            {
              id: "Percentile",
              description: "data crop percentile",
              default: "0",
              value: 0,
              required: false,
              disabled: false,
            },
            {
              id: "Win length",
              description: "smoothe func window length parameter",
              default: "10",
              value: 10,
              required: false,
              disabled: false,
            },
            {
              id: "Polyorder",
              description: "smoothe func polyorder parameter",
              default: "3",
              value: 3,
              required: false,
              disabled: false,
            },
          ],
        },
      ],
      outputDescriptions: [
        {
          id: "Result Matrix",
          description:
            "This is the path to the result matrix data file. It stores the results of the FastTENET algorithm, which includes the inferred regulatory relationships between genes.",
          linked: false,
        },
      ],
      annotations: [],
      pseudotime: [],
      clusters: [],
      annotationColumn: "",
      pseudotimeColumn: "",
      clusterColumn: [],
      activate: false,
    };
  },
  async mounted() {
    this.current_file = this.filterAndAddSuffix(
      this.$store.getters.getCurrentFile.file
    );
    console.log(this.current_file);
    if (this.current_file !== "") {
      try {
        console.log(this.$store.getters.getCurrentFile.file);
        const result = await getColumns({
          file_name: this.current_file,
        });
        console.log(result.data);
        this.annotations = result.data.anno_columns;
        this.pseudotime = result.data.pseudo_columns;
      } catch (error) {
        console.error(error);
      }
    }
  },
  computed: {
    checkCurrentNode() {
      return this.$store.getters.getCurrentNode;
    },
  },
  methods: {
    filterAndAddSuffix(inputString) {
      // Check if the inputString contains an underscore
      if (inputString.includes("_")) {
        // Find the position of the first underscore
        let underscorePosition = inputString.indexOf("_");
        // Add ".h5ad" before the first underscore and exclude everything after underscore
        let modifiedString =
          inputString.substring(0, underscorePosition) + ".h5ad";
        // Return the modified string
        return modifiedString;
      }
      // If no underscore found, return the original string
      return inputString;
    },
    async selectColumns(event) {
      // this.current_file = event.target.value;
      // this.$store.dispatch("setCurrentFile", {
      //   file: this.current_file,
      // });
      // this.$store.dispatch("setCurrentNode", {
      //   node: "algorithm",
      // });
      console.log(event.target.value);
      if (this.annotationColumn !== "") {
        try {
          const result = await getClusters({
            file_name: this.current_file,
            anno_column: this.annotationColumn,
          });
          console.log(result.data);
          this.clusters = result.data.clusters;
        } catch (error) {
          console.error(error);
        }
      }
    },
    async setOption() {
      const make_binary = this.algorithmParts[0].parameters[0].value
        ? "True"
        : "False";
      console.log(this.clusterColumn);
      const optionData = {
        file_name: this.current_file,
        anno_of_interest: this.annotationColumn,
        pseudo_of_interest: this.pseudotimeColumn,
        clusters_of_interest: this.clusterColumn,
        make_binary: make_binary,
        device: this.algorithmParts[1].parameters[0].value,
        device_ids: this.algorithmParts[1].parameters[1].value,
        batch_size: this.algorithmParts[1].parameters[2].value,
        kp: this.algorithmParts[1].parameters[3].value,
        percentile: this.algorithmParts[1].parameters[4].value,
        win_length: this.algorithmParts[1].parameters[5].value,
        polyorder: this.algorithmParts[1].parameters[6].value,
      };
      console.log(optionData);
      try {
        const result = await setupAlgorithm(optionData);
        console.log(result.data);
        this.activate = true;
      } catch (error) {
        console.error(error);
      }
    },
  },
};
</script>

<style scoped>
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
.description-id {
  position: relative;
  display: inline-block;
  margin: 10px;
}

.description-id .description-tooltip {
  visibility: hidden;
  width: 200px;
  background-color: #555;
  color: #fff;
  text-align: center;
  padding: 15px 10px;
  border-radius: 10px;

  position: absolute;
  z-index: 1;
  bottom: 125%;
  left: 50%;
  margin-left: -60px;

  opacity: 0;
  transition: opacity 0.3s;
}

.description-id .description-tooltip::after {
  content: "";
  position: absolute;
  top: 100%;
  left: 30%;
  margin-left: -5px;
  border-width: 5px;
  border-style: solid;
  border-color: #555 transparent transparent transparent;
}

/* Show the tooltip text when you mouse over the tooltip container */
.description-id:hover .description-tooltip {
  visibility: visible;
  opacity: 1;
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
  scrollbar-width: thin;
  scrollbar-color: #888 #f5f5f5;
}

/* 웹킷 브라우저(크롬, 사파리 등)에 대한 스크롤바 스타일 */
.algorithm-layout::-webkit-scrollbar {
  width: 5px; /* 스크롤바 너비 */
}

.algorithm-layout::-webkit-scrollbar-thumb {
  background-color: #cbcbcb; /* 스크롤바 색상 */
  border-radius: 1rem; /* 스크롤바 border-radius */
}

.algorithm-layout::-webkit-scrollbar-thumb:hover {
  background-color: #8d8b8b; /* 스크롤바 호버 시 색상 */
}

.algorithm-logo {
  width: 100%;
  height: 10%;
  /* background-color: blue; */
  top: 1rem;
  object-fit: contain;
}
.algorithm-parts {
  align-content: start;
  width: 100%;
  height: 90%;
  margin-bottom: 2rem;
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
  width: 100%;
  color: #353535;

  margin: 2px;
}
.parameter-id {
  position: relative;
  display: inline-block;
}

.parameter-id .parameter-tooltip {
  visibility: hidden;
  width: 200px;
  background-color: #555;
  color: #fff;
  text-align: center;
  padding: 15px 10px;
  border-radius: 10px;

  position: absolute;
  z-index: 1;
  bottom: 125%;
  left: 50%;
  margin-left: -60px;

  opacity: 0;
  transition: opacity 0.3s;
}

.parameter-id .parameter-tooltip::after {
  content: "";
  position: absolute;
  top: 100%;
  left: 30%;
  margin-left: -5px;
  border-width: 5px;
  border-style: solid;
  border-color: #555 transparent transparent transparent;
}

/* Show the tooltip text when you mouse over the tooltip container */
.parameter-id:hover .parameter-tooltip {
  visibility: visible;
  opacity: 1;
}

.parameter__dropdown {
  width: 9.3rem;
  height: 100%;
  color: black;
  padding: 5px;
  /* right: 10px; */
  border: 0.5px solid black;
  border-radius: 3px;
  border-color: black;
  font-size: small;
  text-align: center;
  margin-bottom: 0px;
}
.parameter__checkbox {
  width: 8.6rem;
  height: 115%;
  color: black;
  padding: 5px;
  /* right: 10px; */
  border: 0.5px solid black;
  border-radius: 3px;
  border-color: black;
  font-size: 1rem;
  font-weight: bolder;
  text-align: center;
  margin-bottom: 0px;
  overflow-y: auto; /* 세로 방향으로 오버플로우가 발생하면 스크롤 표시 */
  display: flex;
  flex-direction: column;
  align-items: center;
}
.parameter__check {
  transform: scale(1.5);
  margin-bottom: 1rem;
  margin-right: 0.5rem;
}
.parameter__check:nth-child(2) {
  margin-left: 0rem;
}
.parameter__menu {
  width: 100%;
  height: 100%;
  color: black;
  padding: 5px;
  /* right: 10px; */
  border-radius: 3px;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  margin-bottom: 0px;
}
.parameter__textInput {
  color: black;
  padding: 5px;
  /* right: 10px; */
  border-radius: 3px;
  border-color: #f1f2fc;
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
.options-layout {
  width: 25%;
  height: 95%;
  /* background-color: red; */
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
}
.options__item {
  /* margin: auto; */
  color: black;
  /* padding: 2% 0; */
  margin: 2%;
}
.form__input {
  position: absolute;
  width: 0;
  height: 0;
  padding: 0;
  overflow: hidden;
  border: 0;
}

.form__button--setup {
  cursor: pointer;
  width: 15rem;
  height: 5rem;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  background: rgb(202, 214, 255);
  /* background: rgb(40, 84, 197); */
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.2rem;
  line-height: 1rem;
  color: black;
}

.apply {
  background: rgb(40, 84, 197);
  color: white;
}

@media (prefers-color-scheme: dark) {
  /* .plotly-layout {
    background-color: rgb(41, 43, 48);
  } */
  .options__item {
    /* color: rgb(255, 255, 255); */
  }
}
</style>
