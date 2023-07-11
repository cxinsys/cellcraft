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
import { getResult } from "@/api/index";

export default {
  data() {
    return {
      logo: require("@/assets/fasttenet.png"),
      current_file: null,
      inputDescriptions: [
        {
          id: "Expression Data (dpath_exp_data)",
          description:
            "This is the path to the expression data file. It contains gene expression values for different samples or cells.",
          linked: true,
        },
        {
          id: "Trajectory Data (dpath_trj_data)",
          description:
            "This is the path to the trajectory data file. It represents the progression or trajectory of cells or samples in a biological process.",
          linked: false,
        },
        {
          id: "Branch Data (dpath_branch_data)",
          description:
            "This is the path to the branch or cell select data file. It specifies the branches or groups of cells in the trajectory.",
          linked: false,
        },
        {
          id: "TF Data (dpath_tf_data)",
          description:
            "This is the path to the transcription factor (TF) data file. It contains information about the transcription factors and their regulatory relationships.",
          linked: false,
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
          linked: false,
        },
        {
          id: "Binary Expression and Node Name Files (make_binary)",
          description:
            "If the make_binary parameter is set to True, FastTENET generates binary expression and node name files. The binary expression file contains a binary representation of the gene expression data, while the node name file contains the names or identifiers of the genes.",
          linked: false,
        },
        {
          id: "GRN (Gene Regulatory Network) Files:",
          description:
            "The GRN files are generated by running the make_grn.py script. They include the inferred gene regulatory network based on the result matrix, node name file, and TF data. The output includes a file with a '.sif' extension, which represents the network structure, and a file with a '.outdegrees.txt' extension, which contains information about the outdegrees (number of outgoing connections) of each gene in the network.",
          linked: true,
        },
      ],
    };
  },
  computed: {
    checkCurrentNode() {
      return this.$store.getters.getCurrentNode;
    },
  },
  watch: {
    async checkCurrentNode(val) {
      const current_node = this.$store.getters.getNodeInfo(val);
      this.current_file = this.$store.getters.getCurrentFile.file;
      if (current_node.name === "algorithm") {
        const filename = {
          filename: `${current_node.name}_${this.current_file.replace(
            ".csv",
            ""
          )}`,
        };
        console.log(filename);
        const dataTableResult = await getResult(filename);
        console.log(dataTableResult);
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
.parameter__textInput {
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

@media (prefers-color-scheme: dark) {
  /* .plotly-layout {
    background-color: rgb(41, 43, 48);
  } */
  .options__item {
    /* color: rgb(255, 255, 255); */
  }
}
</style>
