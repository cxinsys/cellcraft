<template>
  <div id="layout" @click="uncheckCheckbox">
    <div class="left-layout">
      <div class="setup-layout">
        <div class="setup-title">Recent Setting</div>
        <ul class="setup-list">
          <li class="setup-item">
            <!-- <div class="setup-date">Date</div> -->
            <div class="setup-filename">File Name</div>
          </li>
          <li
            class="setup-item"
            v-for="(file, idx) in sortedSetupList"
            :key="idx"
            @click="getOption(file)"
          >
            <!-- <div class="setup-date">{{ file | splitDateAndFileName(1) }}</div> -->
            <div class="setup-filename">
              {{ file | splitDateAndFileName(2) }}
            </div>
          </li>
        </ul>
      </div>
    </div>
    <div class="center-layout">
      <div class="algorithm-layout">
        <div class="algorithm-select">
          <!-- arrow left / right -->
          <!-- <div class="arrow__left" @click="algorithmSelect('left')"></div> -->
          <img
            v-if="current_algorithm == 'fasttenet'"
            class="algorithm-logo"
            :src="fasttenetLog"
            alt="fasttenet"
          />
          <img
            v-else-if="current_algorithm == 'tenet'"
            class="algorithm-logo"
            :src="tenetLogo"
            alt="tenet"
          />
          <!-- <div class="arrow__right" @click="algorithmSelect('right')"></div> -->
        </div>
        <div
          class="algorithm-select__tenet"
          v-if="current_algorithm === 'tenet'"
        >
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
            <span>TENET(TF)</span>
          </label>
        </div>
        <div class="algorithm-parts">
          <div>
            <div class="part-title">Select Input</div>
            <div class="parameters">
              <span class="parameter-id"
                >Option Name
                <span class="parameter-tooltip">Saved option name</span>
              </span>
              <input
                type="text"
                class="parameter__textInput"
                v-model="option_name"
              />
            </div>
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
                :placeholder="
                  this.$store.getters.getCurrentLinkedNodes[0].group
                "
                class="parameter__textInput"
                :disabled="true"
              />
            </div>
            <div class="parameters">
              <span class="parameter-id"
                >Clusters
                <span class="parameter-tooltip">Selected Clusters</span>
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
                <img
                  v-else
                  class="parameter__icon"
                  src="@/assets/checkbox.png"
                />
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
                <span class="parameter-tooltip">
                  Selected Pseudotime Column
                </span>
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
                  :step="parameter.step"
                  :disabled="parameter.disabled"
                  :min="parameter.min"
                  :max="parameter.max"
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
                <label
                  v-else-if="parameter.type === 'upload'"
                  class="parameter__button"
                  :style="buttonStyle"
                >
                  <img
                    v-if="parameter.value === ''"
                    class="parameter__button--icon"
                    src="@/assets/upload-file-black.png"
                  />
                  <div class="parameter__textInput" v-else>
                    {{ parameter.value }}
                  </div>
                  <input
                    class="form__input"
                    type="file"
                    name="file"
                    ref="selectFile"
                    @change.prevent="uploadFile"
                  />
                </label>
                <select
                  class="parameter__dropdown"
                  v-else-if="parameter.type === 'dropdown'"
                  v-model="parameter.value"
                >
                  <option class="parameter__menu" disabled value="">
                    Select Species
                  </option>
                  <option
                    class="parameter__menu"
                    v-for="(data, idx) in parameter.data"
                    :key="idx"
                    :value="data"
                  >
                    {{ data }}
                  </option>
                </select>
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
        <label class="form__button--setup" v-bind:class="{ apply: activate }">
          Set Up
          <input
            @click="setOption"
            class="form__input"
            type="submit"
            value="업로드"
          />
        </label>
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
        <div class="output-title">output result</div>
        <div>
          <label class="form__button--result">
            Check Result
            <input
              @click="refreshGenerated"
              class="form__input"
              type="submit"
              value="업로드"
            />
          </label>
        </div>
        <div
          v-for="outputDescription in outputDescriptions"
          :key="outputDescription.id"
          class="output-description"
          :class="{ generated: outputDescription.generated }"
          @click="downloadOutput(outputDescription.downloadLinkList)"
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
  </div>
</template>

<script>
import {
  getColumns,
  getClusters,
  setupAlgorithm,
  getOptions,
  checkOptions,
  uploadForm,
  getResultFile,
  getDownloadResult,
} from "@/api/index";

export default {
  data() {
    return {
      current_algorithm: "tenet",
      fasttenetLogo: require("@/assets/fasttenet.png"),
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
      ],
      //     id: "fasttenet parameters",
      //     parameters: [
      //       {
      //         id: "make_binary",
      //         type: "boolean",
      //         description: "if True, make binary expression and node name file",
      //         value: "False",
      //         required: false,
      //         disabled: false,
      //       },
      //     ],
      //   },
      //   {
      //     id: "worker run parameters",
      //     parameters: [
      //       {
      //         id: "Device",
      //         type: "string",
      //         description: "cpu or gpu",
      //         value: "gpu",
      //         required: false,
      //         disabled: false,
      //       },
      //       {
      //         id: "Device Ids",
      //         type: "list",
      //         description: "[0](cpu) or [list of whole gpu devices](gpu)",
      //         value: [0, 1, 2, 3, 4, 5, 6, 7],
      //         required: false,
      //         disabled: false,
      //       },
      //       {
      //         id: "Batch Size",
      //         type: "number",
      //         description: "batch size",
      //         value: 2 ** 16,
      //         required: true,
      //         disabled: false,
      //       },
      //       {
      //         id: "Kp",
      //         type: "number",
      //         description: "kernel percentail",
      //         value: 0.5,
      //         required: false,
      //         disabled: false,
      //       },
      //       {
      //         id: "Percentile",
      //         type: "number",
      //         description: "data crop percentile",
      //         value: 0,
      //         required: false,
      //         disabled: false,
      //       },
      //       {
      //         id: "Win length",
      //         type: "number",
      //         description: "smoothe func window length parameter",
      //         value: 10,
      //         required: false,
      //         disabled: false,
      //       },
      //       {
      //         id: "Polyorder",
      //         type: "number",
      //         description: "smoothe func polyorder parameter",
      //         value: 3,
      //         required: false,
      //         disabled: false,
      //       },
      //     ],
      //   },
      // ],
      tenetAlgorithm: [
        {
          id: "tenet parameters",
          parameters: [
            {
              id: "number of threads",
              type: "number",
              step: 1,
              description: "number of threads",
              value: 10,
              max: 10,
              min: 1,
              required: true,
              disabled: false,
            },
            {
              id: "history length",
              type: "number",
              step: 1,
              description: "history length",
              value: 1,
              max: 1,
              min: 1,
              required: true,
              disabled: false,
            },
            {
              id: "species",
              type: "dropdown",
              data: ["human", "mouse", "rat"],
              description: "human or mouse or rat",
              value: "",
              required: true,
              disabled: true,
            },
            {
              id: "gene list",
              type: "upload",
              description: "gene list",
              value: "",
              required: false,
              disabled: false,
            },
          ],
        },
        {
          id: "filter parameters",
          parameters: [
            {
              id: "cutoff for FDR",
              type: "number",
              step: 0.01,
              description: "A cutoff value for FDR by z-test",
              value: 0.01,
              min: 0,
              required: true,
              disabled: false,
            },
            {
              id: "number of links",
              type: "number",
              step: 100,
              description: "The number of links of the GRN",
              value: 1000,
              min: 0,
              required: true,
              disabled: false,
            },
            {
              id: "Trimming Indirect Edges",
              type: "number",
              step: 0.01,
              description: "The cutoff of trimming indirect edges",
              value: 0.01,
              min: 0,
              required: true,
              disabled: false,
            },
          ],
        },
      ],
      outputDescriptions: [
        {
          id: "Z-test filtered GRN",
          description: "Reconstructing GRN Output File",
          generated: false,
          requestFile: "fdr",
          downloadLinkList: [],
        },
        {
          id: "filtered GRN by TE score",
          description: "Reconstructing GRN Output File",
          generated: false,
          requestFile: "NumberOfLinks",
          downloadLinkList: [],
        },
        {
          id: "Indirect edges Trimmed GRN",
          description: "Trimming indirect edges Output File",
          generated: false,
          requestFile: "trimIndirect",
          downloadLinkList: [],
        },
        {
          id: "GRN Outdegree table",
          description: "Counting out-degree of a given GRN Output File",
          generated: false,
          requestFile: "outdegree",
          downloadLinkList: [],
        },
      ],
      annotations: [],
      pseudotime: [],
      clusters: [],
      current_option_file_name: "",
      option_name: "Untitled",
      current_file: null,
      commonOptions: {
        annotationColumn: "",
        pseudotimeColumn: "",
        clusterColumn: [],
      },
      activate: false,
      optionData: {},
      selectedTenetOption: "TENET",
      setupList: [],
      uploadPercentage: 0,
      uploadFailed: false,
      selectFile: null,
    };
  },
  async mounted() {
    // Vuex에서 algorithmOptions 가져오기
    const nodeIndex = this.$store.getters.getCurrentNodeIndex;
    console.log("현재 nodeIndex : " + nodeIndex + typeof nodeIndex);
    const algorithmOptions =
      this.$store.getters.getNodesAlgorithmOptions(nodeIndex);
    console.log("store algorithm option : " + algorithmOptions);
    if (
      algorithmOptions.algorithm !== null ||
      algorithmOptions.algorithm !== undefined
    )
      this.selectedTenetOption = algorithmOptions.algorithm;
    if (algorithmOptions.optionName !== "Untitled")
      this.option_name = algorithmOptions.optionName;
    if (algorithmOptions.optionFilePath !== null)
      this.current_option_file_name = algorithmOptions.optionFilePath;
    if (
      algorithmOptions.commonOptions !== null ||
      algorithmOptions.commonOptions !== undefined
    )
      this.commonOptions = algorithmOptions.commonOptions;
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

    this.current_file = this.filterAndAddSuffix(
      this.$store.getters.getCurrentFile.file
    );
    console.log("store file name is" + this.$store.getters.getCurrentFile.file);
    console.log("current file name is" + this.current_file);
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

    try {
      const result = await checkOptions();
      console.log(result.data);
      // this.current_file이 setupList 안에 있는 모든 아이템에 포함되도록 필터링
      this.setupList = result.data.option_files.filter((item) => {
        return item.includes(this.current_file);
      });
    } catch (error) {
      console.error(error);
    }
    await this.refreshGenerated();
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
    sortedSetupList() {
      return this.setupList.slice().sort((a, b) => {
        return b.localeCompare(a);
      });
    },
    buttonStyle() {
      // 업로드 실패 시 배경색을 빨간색으로 설정
      if (this.uploadFailed) {
        return {
          backgroundColor: "red",
        };
      }
      // 정상적인 업로드 진행률에 따라 배경색 설정
      return {
        backgroundColor: `rgba(76, 175, 80, ${this.uploadPercentage / 100})`,
      };
    },
  },
  watch: {
    selectedTenetOption: {
      handler(newVal) {
        this.$store.commit("setLinkedNodeAlgorithm", {
          nodeIndex: this.$store.getters.getCurrentNodeIndex,
          algorithm: newVal,
        });
      },
      deep: true,
    },
    commonOptions: {
      async handler(newVal, oldVal) {
        this.$store.commit("setLinkedNodeCommonOptions", {
          nodeIndex: this.$store.getters.getCurrentNodeIndex,
          commonOptions: newVal,
        });

        // commonOptions.annotationColumn 값이 변경되었는지 확인
        if (newVal.annotationColumn !== oldVal.annotationColumn) {
          // clusters 업데이트 함수 호출
          try {
            const result = await getClusters({
              file_name: this.current_file,
              anno_column: newVal.annotationColumn,
            });
            console.log(result.data);
            this.clusters = result.data.clusters;
          } catch (error) {
            console.error(error);
          }
        }
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
        this.$store.commit("setLinkedNodeTenetOptions", {
          nodeIndex: this.$store.getters.getCurrentNodeIndex,
          tenetOptions: newVal,
        });
      },
      deep: true,
    },
  },
  filters: {
    // 인자를 하나 더 받아서 1일 경우에 앞 문자열을 return, 2일 경우에 뒤 문자열을 return하도록 수정해줘
    // template에서 사용할 때는 {{ value | splitDateAndFileName(1) }} 이런 식으로 사용하면 됨
    splitDateAndFileName(value, num) {
      if (num === 1) {
        // 각 부분(년, 월, 일, 시, 분)을 추출합니다.
        const year = value.split("_")[0].substring(0, 4);
        const month = value.split("_")[0].substring(4, 6);
        const day = value.split("_")[0].substring(6, 8);
        const hour = value.split("_")[0].substring(8, 10);
        const minute = value.split("_")[0].substring(10, 12);
        // 추출한 부분을 '년/월/일 시:분' 형식으로 조합합니다.
        return `${year}/${month}/${day} ${hour}:${minute}`;
      } else if (num === 2) {
        return value.split("_")[1];
      }
    },
  },
  methods: {
    filterAndAddSuffix(inputString) {
      // Check if the inputString contains an underscore
      if (inputString.includes("_")) {
        // "_"로 구분된 문자열을 배열로 변환
        const segments = inputString.split("_");
        // 마지막 두 요소를 제외한 나머지를 합침
        const fileName = segments.slice(0, -2).join("_") + ".h5ad";
        return fileName;
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
    },
    async uploadFile() {
      if (this.$refs.selectFile[0].files.length > 0) {
        this.selectFile = new File(
          [this.$refs.selectFile[0].files[0]],
          this.$refs.selectFile[0].files[0].name
        );
        this.tenetAlgorithm[0].parameters[3].value =
          this.$refs.selectFile[0].files[0].name;
        const form = new FormData();
        form.append("files", this.selectFile);
        // 파일 업로드 진행률을 추적하기 위한 콜백
        const onUploadProgress = (progressEvent) => {
          this.uploadPercentage = parseInt(
            Math.round((progressEvent.loaded * 100) / progressEvent.total)
          );
        };
        try {
          const response = await uploadForm(form, onUploadProgress);
          console.log(response);
          this.uploadFailed = false; // 업로드 성공
        } catch (error) {
          console.error(error);
          this.uploadFailed = true; // 업로드 실패
        }
      }
    },
    async getOption(filename) {
      if (
        confirm(
          "Are you sure you want to import the selected options? Please make sure you have saved the currently selected options."
        )
      ) {
        console.log(filename);
        try {
          const result = await getOptions(filename);
          console.log(result.data);

          this.option_name = filename.split("_")[1];
          this.current_option_file_name = filename;
          this.$store.commit("setLinkedNodeOptionFilePath", {
            nodeIndex: this.$store.getters.getCurrentNodeIndex,
            optionFilePath: this.current_option_file_name,
          });

          this.selectedTenetOption = result.data.algorithm;

          // commonOptions 설정
          this.commonOptions = {
            annotationColumn: result.data.anno_of_interest,
            pseudotimeColumn: result.data.pseudo_of_interest,
            clusterColumn: result.data.clusters_of_interest,
          };

          this.tenetAlgorithm = this.tenetAlgorithm.map((group) => {
            return {
              ...group,
              parameters: group.parameters.map((param) => {
                switch (param.id) {
                  case "number of threads":
                    return { ...param, value: result.data.num_of_threads };
                  case "history length":
                    return { ...param, value: result.data.history_length };
                  case "species":
                    return { ...param, value: result.data.species };
                  case "gene list":
                    return { ...param, value: result.data.gene_list_file };
                  case "cutoff for FDR":
                    return { ...param, value: result.data.cutoff_for_fdr };
                  case "number of links":
                    return { ...param, value: result.data.num_of_links };
                  case "Trimming Indirect Edges":
                    return {
                      ...param,
                      value: result.data.trimming_indirect_edges,
                    };
                  default:
                    return param;
                }
              }),
            };
          });
        } catch (error) {
          console.error(error);
        }
      }
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
      console.log(this.commonOptions.pseudotimeColumn, this.isSelectedIndices);

      // 옵션들을 전부 선택했는지 확인 후, 다 선택되지 않았다면 alert
      if (
        this.commonOptions.pseudotimeColumn === "" ||
        (this.isSelectedIndices !== "selected" &&
          (this.commonOptions.annotationColumn === "" ||
            this.commonOptions.clusterColumn.length === 0))
      ) {
        alert("Please select all the options.");
        return;
      }

      // 선택한 옵션들을 저장하시겠습니까? 영어로 confirm
      if (
        !confirm(
          "Are you sure you want to save the selected options? Please make sure you have selected all the options."
        )
      ) {
        return;
      }

      this.optionData = {
        algorithm: this.selectedTenetOption,
        option_name: this.option_name,
        file_name: this.current_file,
        anno_of_interest: this.commonOptions.annotationColumn,
        pseudo_of_interest: this.commonOptions.pseudotimeColumn,
        clusters_of_interest: this.commonOptions.clusterColumn,
        num_of_threads: this.tenetAlgorithm[0].parameters[0].value,
        history_length: this.tenetAlgorithm[0].parameters[1].value,
        species: this.tenetAlgorithm[0].parameters[2].value,
        gene_list: this.tenetAlgorithm[0].parameters[3].value,
        cutoff_for_fdr: this.tenetAlgorithm[1].parameters[0].value,
        num_of_links: this.tenetAlgorithm[1].parameters[1].value,
        trimming_indirect_edges: this.tenetAlgorithm[1].parameters[2].value,
      };
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
        this.$store.commit("setLinkedNodeOptionFilePath", {
          nodeIndex: this.$store.getters.getCurrentNodeIndex,
          optionFilePath: result.data.file_name,
        });
        this.activate = true;
        try {
          const option_result = await checkOptions();
          console.log(option_result.data);
          // this.current_file이 setupList 안에 있는 모든 아이템에 포함되도록 필터링
          this.setupList = option_result.data.option_files.filter((item) => {
            return item.includes(this.current_file);
          });
        } catch (error) {
          console.error(error);
        }
      } catch (error) {
        console.error(error);
      }
    },
    async refreshGenerated() {
      // outputDescriptions를 순회하며 generated를 false로 설정 및 downloadLinkList 초기화
      this.outputDescriptions.forEach((desc) => {
        desc.generated = false;
        desc.downloadLinkList = [];
      });

      try {
        const option_file_name = this.$store.getters.getCurrentFile;
        const result = await getResultFile({
          file_name: this.current_file,
          option_file_name: option_file_name.algorithmOptions.optionFilePath,
        });
        const resultsList = result.data.result_files;

        if (resultsList.length === 0) {
          this.outputDescriptions.forEach((desc) => (desc.generated = false));
          // Recent Setting에서 실행한 알고리즘 세팅을 다시 선택해주세요. 혹은 분석을 실행한 알고리즘이 실행 중인지 Task Manager에서 확인해주세요. 영어로 경고
          alert(
            "Please select the algorithm setting you ran in Recent Setting again. Or check if the algorithm you ran in the analysis is running in Task Monitoring."
          );
          // alert("There is no generated output yet.");
          return;
        }

        // 이미 추가된 파일을 추적하기 위한 집합
        let addedFiles = new Set();

        // outputDescriptions를 역순으로 순회
        for (let i = this.outputDescriptions.length - 1; i >= 0; i--) {
          let desc = this.outputDescriptions[i];
          desc.generated = false; // 기본값으로 false 설정

          for (let file of resultsList) {
            // 이미 추가된 파일이 아니고, requestFile이 포함된 경우
            if (!addedFiles.has(file) && file.includes(desc.requestFile)) {
              desc.generated = true;
              desc.downloadLinkList.push(file);
              addedFiles.add(file); // 추가된 파일로 기록
            }
          }
        }
      } catch (error) {
        console.error("Error in refreshGenerated: ", error);
      }
    },
    async downloadOutput(downloadLinkList) {
      if (downloadLinkList.length === 0) {
        // 생성된 아웃풋이 없다고 영어로 경고
        alert("There is no generated output.");
        return;
      }

      // 해당 아웃풋을 다운 받으시겠습니까? 영어로 confirm
      if (
        !confirm(
          "Are you sure you want to download the selected output? Please make sure you have selected the output."
        )
      ) {
        return;
      }

      // 다운로드 링크를 모두 순회하며 다운로드
      downloadLinkList.forEach(async (downloadLink) => {
        try {
          const result = await getDownloadResult(downloadLink);
          console.log(result.data);
          // 다운로드 링크를 생성
          const url = window.URL.createObjectURL(new Blob([result.data]));
          // 링크를 클릭하도록 설정
          const linkElement = document.createElement("a");
          linkElement.href = url;
          linkElement.setAttribute("download", downloadLink);
          // 클릭 이벤트를 발생시킴
          linkElement.click();
        } catch (error) {
          console.error(error);
        }
      });
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
.setup-layout {
  width: 25%;
  height: 95%;
  margin-left: 1rem;
  display: flex;
  align-items: center;
  /* justify-content: center; */
  flex-direction: column;
  padding-left: 1rem;
  padding-right: 1rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
}
.output-layout {
  width: 25%;
  height: 95%;
  margin-left: 1rem;
  display: flex;
  align-items: center;
  /* justify-content: center; */
  flex-direction: column;
  padding-left: 1rem;
  padding-right: 1rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
}
.setup-title,
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
  color: #a1a1a1;
  background-color: rgb(224, 224, 224);
  border-radius: 1rem;
  margin: 0.5rem 0rem;
  border: 2px solid #e7e7e7;
  opacity: 0.9;
}
.output-description {
  cursor: not-allowed;
}
.input-description.linked,
.output-description.generated {
  opacity: 1;
  color: #fff;
  background-color: rgb(22, 152, 51);
  border: 2px solid #ecebff;
}
.output-description.generated {
  cursor: pointer;
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
.setup-list {
  width: 100%;
  height: 100%;
  overflow-y: auto;
  scrollbar-width: thin;
  scrollbar-color: #888 #f5f5f5;
}
.setup-item {
  width: 100%;
  display: flex;
  align-items: center;
  cursor: pointer;
  margin-bottom: 0.5rem;
}
.setup-item:first-child {
  pointer-events: none;
}
.setup-item:hover {
  background-color: rgb(202, 214, 255);
}
.setup-date,
.setup-filename {
  width: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
  font-size: 0.9rem;
  color: #353535;
  padding: 0.5rem;
}
.setup-item:first-child > .setup-date,
.setup-item:first-child > .setup-filename {
  font-weight: bold;
  font-size: 1rem;
  opacity: 1;
  background-color: rgb(224, 224, 224);
  border-radius: 1rem;
}
.center-layout {
  /* margin-top: 1rem;
  margin-bottom: 1rem; */
  width: 41%;
  height: 100%;
  margin-right: 1rem;
  /* background-color: blue; */
  display: flex;
  align-items: center;
  /* justify-content: center; */
  flex-direction: column;
  box-sizing: border-box;
  overflow-y: scroll;
  scrollbar-width: thin;
  scrollbar-color: #888 #f5f5f5;
}

/* 웹킷 브라우저(크롬, 사파리 등)에 대한 스크롤바 스타일 */
.center-layout::-webkit-scrollbar {
  width: 5px; /* 스크롤바 너비 */
}

.center-layout::-webkit-scrollbar-thumb {
  background-color: #cbcbcb; /* 스크롤바 색상 */
  border-radius: 1rem; /* 스크롤바 border-radius */
}

.center-layout::-webkit-scrollbar-thumb:hover {
  background-color: #8d8b8b; /* 스크롤바 호버 시 색상 */
}
.algorithm-layout {
  background-color: rgb(255, 255, 255);
  /* margin: 1rem; */
  padding: 1rem;
  border-radius: 1rem;
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
  text-transform: capitalize;
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

  text-transform: capitalize;
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
  font-size: 0.9rem;
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

.parameter__button {
  width: 9.3rem;
  height: 2rem;
  color: black;
  border: 0.5px solid black;
  border-radius: 3px;
  border-color: black;
  margin-bottom: 0px;
  display: flex;
  align-items: center;
  justify-content: center;
  position: relative;
  cursor: pointer;
}
.parameter__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
  opacity: 0.8;
  margin-right: 0.5rem;
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
  width: 8.5rem;
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
.left-layout {
  width: 23%;
  height: 100%;
  margin-left: -1rem;
  margin-right: 1rem;
  /* background-color: green; */
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
}
.right-layout {
  width: 27%;
  height: 95%;
  /* margin-right: 1rem; */
  /* background-color: green; */
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
}
.setup-layout,
.output-layout {
  width: 100%;
  margin-left: 1rem;
  margin-right: 0rem;
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
  display: none;
}

.form__button--setup {
  cursor: pointer;
  width: 8rem;
  height: 3rem;
  margin-left: 12.5rem;
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
  color: rgb(240, 240, 240);
}
.form__button--result {
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
  color: rgb(249, 249, 249);
  margin-bottom: 1rem;
}
.form__button--result:hover {
  background: rgb(40, 84, 197);
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

.reset-button {
  position: absolute;
  margin-top: -5px;
  width: 1.5rem;
  height: 1.5rem;
  opacity: 0.7;
}
.reset-button:hover {
  opacity: 0.5;
  cursor: pointer;
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
</style>
