<template>
  <div id="layout" @click="uncheckCheckbox">
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
      <div class="algorithm-select">
        <!-- arrow left / right -->
        <div class="arrow__left" @click="algorithmSelect('left')"></div>
        <img
          v-if="current_algorithm == 'fasttenet'"
          class="algorithm-logo"
          :src="fasttentLogo"
          alt="fasttenet"
        />
        <img
          v-else-if="current_algorithm == 'tenet'"
          class="algorithm-logo"
          :src="tenetLogo"
          alt="tenet"
        />
        <div class="arrow__right" @click="algorithmSelect('right')"></div>
      </div>
      <div class="algorithm-select__tenet" v-if="current_algorithm === 'tenet'">
        <input
          type="radio"
          name="select"
          id="option-1"
          value="TENET"
          v-model="selectedTenetOption"
          checked
        />
        <input
          type="radio"
          name="select"
          id="option-2"
          value="TENET_TF"
          v-model="selectedTenetOption"
        />
        <label for="option-1" class="option option-1">
          <div class="dot"></div>
          <span>TENET</span>
        </label>
        <label for="option-2" class="option option-2">
          <div class="dot"></div>
          <span>TENET_TF</span>
        </label>
      </div>
      <div class="algorithm-parts">
        <div>
          <div class="part-title">Select Input</div>
          <div class="parameters">
            <span class="parameter-id"
              >File
              <span class="parameter-tooltip">Selected File Name</span>
            </span>
            <select
              class="parameter__dropdown"
              :class="{
                wiggle: checkSelectedIndices && isFileSelected === false,
              }"
              v-if="checkSelectedIndices"
              v-model="isSelectedIndices"
              @change="isFileSelected = true"
            >
              <option class="parameter__menu" disabled value="">
                Select File
              </option>
              <option
                class="parameter__menu"
                :value="this.$store.getters.getCurrentFile.file"
              >
                {{ this.current_file }}
              </option>
              <option class="parameter__menu" value="selected">
                {{ this.current_file }} (selected)
              </option>
            </select>
            <input
              v-else
              type="text"
              :placeholder="this.current_file"
              class="parameter__textInput"
              :v-model="this.current_file"
              :class="{ 'red-text': !this.current_file }"
              :disabled="true"
            />
          </div>
          <div class="parameters">
            <span class="parameter-id"
              >Cell Group
              <span class="parameter-tooltip">Selected Cell Group</span>
            </span>
            <select
              v-if="checkFileSelected || !isFileSelected"
              class="parameter__dropdown"
              :class="{ wiggle: this.commonOptions.annotationColumn === '' }"
              v-model="commonOptions.annotationColumn"
              @change="selectColumns($event)"
            >
              <option class="parameter__menu" disabled value="">
                Select Cell Group
              </option>
              <option
                v-for="(column, index) in annotations"
                :key="index"
                :value="column"
                class="parameter__menu"
              >
                {{ column }}
              </option>
            </select>
            <input
              v-else
              type="text"
              :placeholder="this.$store.getters.getCurrentLinkedNodes[0].group"
              class="parameter__textInput"
              :disabled="true"
            />
          </div>
          <div class="parameters">
            <span class="parameter-id"
              >Clusters
              <span class="parameter-tooltip">Selected File Name</span>
            </span>
            <label
              v-if="checkFileSelected || !isFileSelected"
              class="parameter__dropdown"
              for="touch"
              :class="{
                wiggle:
                  this.commonOptions.annotationColumn !== '' &&
                  this.commonOptions.clusterColumn.length === 0,
              }"
            >
              <span>Select Clusters</span>
              <img
                v-if="commonOptions.clusterColumn.length !== 0"
                class="parameter__icon"
                src="@/assets/check.png"
              />
              <img v-else class="parameter__icon" src="@/assets/checkbox.png" />
            </label>
            <input
              v-else
              type="text"
              placeholder="already selected"
              class="parameter__textInput"
              :disabled="true"
            />
            <input
              ref="checkbox"
              class="parameter__drop"
              type="checkbox"
              id="touch"
              :disabled="this.commonOptions.annotationColumn === ''"
            />
            <div ref="checkbox__area" class="parameter__checkbox">
              <div
                class="parameter__menu"
                v-for="(column, index) in clusters"
                :key="index"
                @click="clusterToggle(column)"
              >
                {{ column }}
                <img
                  v-if="commonOptions.clusterColumn.includes(column)"
                  class="parameter__icon"
                  src="@/assets/check.png"
                  alt="check"
                />
                <img
                  v-else
                  class="parameter__icon"
                  src="@/assets/checkbox.png"
                  alt="check"
                />
              </div>
            </div>
          </div>
          <div class="parameters">
            <span class="parameter-id">
              Pseudotime Column
              <span class="parameter-tooltip">Selected File Name</span>
            </span>
            <select
              class="parameter__dropdown"
              :class="{ wiggle: this.commonOptions.pseudotimeColumn === '' }"
              v-model="commonOptions.pseudotimeColumn"
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
        </div>
        <div v-for="part in currentAlgorithmData" :key="part.id">
          <div class="part-title">{{ part.id }}</div>
          <div
            v-for="parameter in part.parameters"
            :key="parameter.id"
            v-show="
              parameter.id !== 'species' || selectedTenetOption === 'TENET_TF'
            "
          >
            <div class="parameters">
              <span class="parameter-id"
                >{{ parameter.id }}
                <span class="parameter-tooltip">{{
                  parameter.description
                }}</span>
              </span>
              <input
                type="number"
                v-if="parameter.type === 'number'"
                class="parameter__textInput"
                v-model="parameter.value"
                :class="{ 'red-text': !parameter.value }"
                :disabled="parameter.disabled"
              />
              <div
                v-else-if="parameter.type === 'radio'"
                class="parameter__radioBox"
              >
                <div
                  v-for="(option, idx) in Object.keys(parameter.value)"
                  :class="{ selected: parameter.value[option] }"
                  class="parameter__radio"
                  :key="idx"
                  @click="selectOption(option, parameter.value)"
                >
                  {{ option }}
                </div>
              </div>
              <input
                type="text"
                v-else
                class="parameter__textInput"
                v-model="parameter.value"
                :class="{ 'red-text': !parameter.value }"
                :disabled="parameter.disabled"
              />
            </div>
          </div>
        </div>
      </div>
      <!-- run button -->
      <!-- <div>
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
      </div> -->
    </div>
    <div class="right-layout">
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
      <div>
        <div>
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
  </div>
</template>

<script>
import { getColumns, getClusters, setupAlgorithm } from "@/api/index";

export default {
  data() {
    return {
      current_algorithm: "fasttenet",
      fasttentLogo: require("@/assets/fasttenet.png"),
      tenetLogo: require("@/assets/tenet.png"),
      algorithms: ["fasttenet", "tenet"],
      isAnimated: false,
      isFileSelected: false,
      isSelectedIndices: "",
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
      fasttenetAlgorithm: [
        {
          id: "fasttenet parameters",
          parameters: [
            {
              id: "make_binary",
              type: "boolean",
              description: "if True, make binary expression and node name file",
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
              type: "string",
              description: "cpu or gpu",
              value: "gpu",
              required: false,
              disabled: false,
            },
            {
              id: "Device Ids",
              type: "list",
              description: "[0](cpu) or [list of whole gpu devices](gpu)",
              value: [0, 1, 2, 3, 4, 5, 6, 7],
              required: false,
              disabled: false,
            },
            {
              id: "Batch Size",
              type: "number",
              description: "batch size",
              value: 2 ** 16,
              required: true,
              disabled: false,
            },
            {
              id: "Kp",
              type: "number",
              description: "kernel percentail",
              value: 0.5,
              required: false,
              disabled: false,
            },
            {
              id: "Percentile",
              type: "number",
              description: "data crop percentile",
              value: 0,
              required: false,
              disabled: false,
            },
            {
              id: "Win length",
              type: "number",
              description: "smoothe func window length parameter",
              value: 10,
              required: false,
              disabled: false,
            },
            {
              id: "Polyorder",
              type: "number",
              description: "smoothe func polyorder parameter",
              value: 3,
              required: false,
              disabled: false,
            },
          ],
        },
      ],
      tenetAlgorithm: [
        {
          id: "tenet parameters",
          parameters: [
            {
              id: "number of threads",
              type: "number",
              description: "number of threads",
              value: 10,
              required: true,
              disabled: false,
            },
            {
              id: "history length",
              type: "number",
              description: "history length",
              value: 1,
              required: true,
              disabled: false,
            },
            {
              id: "species",
              type: "radio",
              description: "human or mouse",
              value: {
                human: true,
                mouse: false,
              },
              required: true,
              disabled: true,
            },
          ],
        },
        {
          id: "filter parameters",
          parameters: [
            {
              id: "cutoff for FDR",
              type: "number",
              description: "A cutoff value for FDR by z-test",
              value: 0.01,
              required: true,
              disabled: false,
            },
            {
              id: "number of links",
              type: "number",
              description: "The number of links of the GRN",
              value: 1000,
              required: true,
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
      current_file: null,
      commonOptions: {
        annotationColumn: "",
        pseudotimeColumn: "",
        clusterColumn: [],
      },
      activate: false,
      optionData: {},
      selectedTenetOption: "TENET",
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
    // Vuex에서 algorithmOptions 가져오기
    const algorithmOptions = this.$store.getters.getAlgorithmOptions;
    console.log(algorithmOptions);
    if (
      algorithmOptions.algorithm !== null ||
      algorithmOptions.algorithm !== undefined
    )
      this.current_algorithm = algorithmOptions.algorithm;
    if (
      algorithmOptions.commonOptions !== null ||
      algorithmOptions.commonOptions !== undefined
    )
      this.commonOptions = algorithmOptions.commonOptions;
    if (algorithmOptions.fasttenetOptions !== null)
      this.fasttenetAlgorithm = algorithmOptions.fasttenetOptions;
    if (algorithmOptions.tenetOptions !== null)
      this.tenetAlgorithm = algorithmOptions.tenetOptions;

    if (this.commonOptions.annotationColumn !== "") {
      try {
        const result = await getClusters({
          file_name: this.current_file,
          anno_column: this.commonOptions.annotationColumn,
        });
        console.log(result.data);
        this.clusters = result.data.clusters;
      } catch (error) {
        console.error(error);
      }
    }
  },
  computed: {
    checkCurrentNode() {
      return this.$store.getters.getCurrentNode;
    },
    checkSelectedIndices() {
      return Object.keys(this.$store.getters.getCurrentLinkedNodes[0]).includes(
        "selectedIndices"
      );
    },
    checkFileSelected() {
      return this.isSelectedIndices === this.$store.getters.getCurrentFile.file;
    },
    currentAlgorithmData() {
      switch (this.current_algorithm) {
        case "fasttenet":
          return this.fasttenetAlgorithm;
        case "tenet":
          return this.tenetAlgorithm;
        default:
          return [];
      }
    },
  },
  watch: {
    current_algorithm: {
      handler(newVal) {
        this.$store.commit("setAlgorithm", newVal);
      },
      deep: true,
    },
    commonOptions: {
      handler(newVal) {
        this.$store.commit("setCommonOptions", newVal);
      },
      deep: true,
    },
    fasttenetAlgorithm: {
      handler(newVal) {
        this.$store.commit("setFasttenetOptions", newVal);
      },
      deep: true,
    },
    tenetAlgorithm: {
      handler(newVal) {
        this.$store.commit("setTenetOptions", newVal);
      },
      deep: true,
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
    algorithmSelect(direction) {
      let index = this.algorithms.indexOf(this.current_algorithm);
      if (direction === "left") {
        if (index === 0) {
          index = this.algorithms.length - 1;
        } else {
          index -= 1;
        }
      } else if (direction === "right") {
        if (index === this.algorithms.length - 1) {
          index = 0;
        } else {
          index += 1;
        }
      }
      this.current_algorithm = this.algorithms[index];
      console.log(this.current_algorithm);

      this.isAnimated = true;
      setTimeout(() => {
        this.isAnimated = false;
      }, 1000);
    },
    uncheckCheckbox(event) {
      // 체크박스 자체가 클릭된 경우에는 아무 것도 하지 않습니다.
      if (
        this.$refs.checkbox.contains(event.target) ||
        this.$refs.checkbox__area.contains(event.target)
      ) {
        return;
      }
      // 체크박스 외의 영역이 클릭되면 체크박스의 상태를 해제합니다.
      this.$refs.checkbox__area.checked = false;
    },
    clusterToggle(column) {
      if (this.commonOptions.clusterColumn.includes(column)) {
        this.commonOptions.clusterColumn =
          this.commonOptions.clusterColumn.filter((item) => {
            return item !== column;
          });
      } else {
        this.commonOptions.clusterColumn.push(column);
      }
      console.log(this.commonOptions.clusterColumn);
    },
    selectOption(selectedKey, values) {
      console.log(selectedKey, values);
      Object.keys(values).forEach((key) => {
        values[key] = key === selectedKey;
      });
    },
    async selectColumns(event) {
      console.log(event.target.value);
      if (this.commonOptions.annotationColumn !== "") {
        try {
          const result = await getClusters({
            file_name: this.current_file,
            anno_column: this.commonOptions.annotationColumn,
          });
          console.log(result.data);
          this.clusters = result.data.clusters;
        } catch (error) {
          console.error(error);
        }
      }
    },
    async setOption() {
      // 알고리즘 종류에 따라 this.optionData 설정
      if (this.current_algorithm == "fasttenet") {
        const make_binary = this.fasttenetAlgorithm[0].parameters[0].value
          ? "True"
          : "False";
        this.optionData = {
          algorithm: this.current_algorithm,
          file_name: this.current_file,
          anno_of_interest: this.commonOptions.annotationColumn,
          pseudo_of_interest: this.commonOptions.pseudotimeColumn,
          clusters_of_interest: this.commonOptions.clusterColumn,
          make_binary: make_binary,
          device: this.fasttenetAlgorithm[1].parameters[0].value,
          device_ids: this.fasttenetAlgorithm[1].parameters[1].value,
          batch_size: this.fasttenetAlgorithm[1].parameters[2].value,
          kp: this.fasttenetAlgorithm[1].parameters[3].value,
          percentile: this.fasttenetAlgorithm[1].parameters[4].value,
          win_length: this.fasttenetAlgorithm[1].parameters[5].value,
          polyorder: this.fasttenetAlgorithm[1].parameters[6].value,
        };
        console.log(this.optionData);
      } else if (this.current_algorithm == "tenet") {
        this.optionData = {
          algorithm: this.current_algorithm,
          selceted_tenet: this.selectedTenetOption,
          file_name: this.current_file,
          anno_of_interest: this.commonOptions.annotationColumn,
          pseudo_of_interest: this.commonOptions.pseudotimeColumn,
          clusters_of_interest: this.commonOptions.clusterColumn,
          num_of_threads: this.tenetAlgorithm[0].parameters[0].value,
          history_length: this.tenetAlgorithm[0].parameters[1].value,
          species: this.tenetAlgorithm[0].parameters[2].value["human"]
            ? "human"
            : "mouse",
          cutoff_for_fdr: this.tenetAlgorithm[1].parameters[0].value,
          num_of_links: this.tenetAlgorithm[1].parameters[1].value,
        };
        console.log(this.optionData);
      }
      // 현재 링크된 노드 안에 selectedIndices가 있는지 확인 후, 있으면 optionData에 추가
      const linkedNodes = this.$store.getters.getCurrentLinkedNodes;
      if (
        Object.keys(linkedNodes[0]).includes("selectedIndices") &&
        this.isSelectedIndices !== this.$store.getters.getCurrentFile.file
      ) {
        console.log(typeof linkedNodes[0].selectedIndices[0]);
        this.optionData["selected_indices"] = linkedNodes[0].selectedIndices;
        this.optionData["anno_of_interest"] = linkedNodes[0].group;
        this.optionData["clusters_of_interest"] = [];
      }
      try {
        const result = await setupAlgorithm(this.optionData);
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

.algorithm-select {
  width: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}

.algorithm-logo {
  /* width: 50%; */
  height: 5rem;
  /* background-color: blue; */
  top: 1rem;
  object-fit: contain;
}

.arrow__left::after {
  content: "◀";
  font-size: 2rem;
  color: rgb(48, 48, 48);
  cursor: pointer;
  margin-right: 1rem;
}

.arrow__right::after {
  content: "▶";
  font-size: 2rem;
  color: rgb(48, 48, 48);
  cursor: pointer;
  margin-left: 1rem;
}

.algorithm-parts {
  align-content: start;
  width: 100%;
  margin-bottom: 2rem;
}
.algorithm-select__tenet {
  width: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 1rem;
  box-sizing: border-box;
}
.option {
  background: #fff;
  height: 100%;
  width: 100%;
  display: flex;
  align-items: center;
  justify-content: space-evenly;
  margin: 0 10px;
  border-radius: 5px;
  cursor: pointer;
  padding: 0 10px;
  border: 2px solid lightgrey;
  transition: all 0.3s ease;
}
.dot {
  height: 20px;
  width: 20px;
  background: #d9d9d9;
  border-radius: 50%;
  position: relative;
}
.dot::before {
  position: absolute;
  content: "";
  top: 4px;
  left: 4px;
  width: 12px;
  height: 12px;
  background: #0069d9;
  border-radius: 50%;
  opacity: 0;
  transform: scale(1.5);
  transition: all 0.3s ease;
}
#option-1:checked:checked ~ .option-1,
#option-2:checked:checked ~ .option-2 {
  border-color: #0069d9;
  background: #0069d9;
}
#option-1:checked:checked ~ .option-1 .dot,
#option-2:checked:checked ~ .option-2 .dot {
  background: #fff;
}
#option-1:checked:checked ~ .option-1 .dot::before,
#option-2:checked:checked ~ .option-2 .dot::before {
  opacity: 1;
  transform: scale(1);
}
.wrapper .option span {
  font-size: 20px;
  color: #808080;
}
#option-1:checked:checked ~ .option-1 span,
#option-2:checked:checked ~ .option-2 span {
  color: #fff;
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
  height: 2rem;
  color: #353535;
  position: relative;
  margin: 2px;
  margin-bottom: 0.5rem;
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
  color: black;
  border: 0.5px solid black;
  border-radius: 3px;
  border-color: black;
  font-size: small;
  text-align: left;
  margin-bottom: 0px;
  display: flex;
  align-items: center;
  justify-content: left;
  position: relative;
}

.parameter__dropdown:focus {
  outline: none;
  border-color: #aaa;
}

.parameter__dropdown span {
  color: black;
  font-size: 13px;
  font-weight: 400;
  cursor: pointer;
  display: block;
  font-family: Arial, Helvetica, sans-serif;
  margin-left: 0.8rem;
}

/* .parameter__dropdown span::after {
  position: absolute;
  right: 0.1rem;
  top: calc(50% - 0.8rem);
  content: "☑";
  font-size: 1.5rem;
} */

.parameter__drop {
  position: absolute;
  opacity: 0;
  height: 0px;
}

.parameter__menu:checked::before {
  content: "✔";
}
.parameter__radioBox {
  width: 9.3rem;
  display: flex;
  align-items: center;
  justify-content: center;
}
.parameter__radio {
  width: 50%;
  border: 1px solid #8b8887;
  padding: 0.5rem 0;
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: pointer;
}
.parameter__checkbox {
  display: none;
  position: absolute;
  clear: both;
  width: 9.3rem;
  max-height: 300px;
  overflow-y: auto;
  overflow-x: hidden;
  text-align: center;
  background: #6c6867;
  transition: height 0.4s ease;
  border-radius: 0.5rem;
}
/* 웹킷 브라우저(크롬, 사파리 등)에 대한 스크롤바 스타일 */
.parameter__checkbox::-webkit-scrollbar {
  width: 5px; /* 스크롤바 너비 */
}

.parameter__checkbox::-webkit-scrollbar-thumb {
  background-color: #cbcbcb; /* 스크롤바 색상 */
  border-radius: 1rem; /* 스크롤바 border-radius */
}

.parameter__checkbox::-webkit-scrollbar-thumb:hover {
  background-color: #8d8b8b; /* 스크롤바 호버 시 색상 */
}
.parameter__checkbox div {
  height: 26px;
  line-height: 1.2rem;
  color: white;
  display: flex;
  align-items: center;
  justify-content: center;
}
#touch:checked + .parameter__checkbox {
  border: 1px solid #8b8887;
  display: block;
  top: calc(100%);
  right: 0;
  z-index: 9999;
}
.parameter__check {
  transform: scale(1.5);
  margin-bottom: 1rem;
  margin-right: 0.5rem;
}
.parameter__check:nth-child(2) {
  margin-left: 0rem;
}
.parameter__icon {
  position: absolute;
  right: 0.6rem;
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  padding-left: 2rem;
}
.parameter__menu {
  width: 85%;
  height: 100%;
  color: black;
  padding: 3px 7.5%;
  /* right: 10px; */
  border-radius: 3px;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  margin-bottom: 0px;
  cursor: pointer;
  position: relative;
  word-wrap: break-word; /* 긴 단어가 div의 너비를 넘어갈 때 줄바꿈 */
  overflow-wrap: break-word; /* CSS3에서 word-wrap 대신 사용 */
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
.right-layout {
  width: 25%;
  height: 95%;
  margin-right: 1rem;
  /* background-color: green; */
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
}
.output-layout {
  width: 100%;
  margin: 1rem;
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
  width: 12rem;
  height: 3rem;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.5);
  border-radius: 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  /* background: rgb(163, 163, 163); */
  background: rgb(59, 108, 221);
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.2rem;
  line-height: 1rem;
  color: rgb(202, 202, 202);
}

.apply {
  cursor: pointer;
  background: rgb(40, 84, 197);
  color: white;
}
.wiggle {
  animation: wiggle 2s linear infinite;
}
.selected {
  /* 선택된 요소의 배경색 변경 */
  background-color: #0069d9;
  color: white;
}

input[type="radio"] {
  display: none;
}

@keyframes wiggle {
  0%,
  7% {
    transform: rotateZ(0);
  }
  15% {
    transform: rotateZ(-5deg);
  }
  20% {
    transform: rotateZ(3deg);
  }
  25% {
    transform: rotateZ(-3deg);
  }
  30% {
    transform: rotateZ(1deg);
  }
  35% {
    transform: rotateZ(-1deg);
  }
  40%,
  100% {
    transform: rotateZ(0);
  }
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
