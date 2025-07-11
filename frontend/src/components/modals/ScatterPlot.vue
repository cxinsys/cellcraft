<template>
  <div id="layout">
    <div class="plotly-layout">
      <div id="plotly__scatter"></div>
    </div>
    <div class="options-layout" v-if="!refreshOptionsLayout">
      <input type="text" placeholder="Title" class="options__textInput" @input="titleChangeFunc($event)" />
      <div class="options__item">
        <p class="options__title">X - axis</p>
        <select class="options__item--select" :value="selectedX" @change="setSelectX($event)">
          <option v-for="(item, index) in numList" :key="index" :value="item.value" :disabled="item.value === null">
            {{ item.name }}
          </option>
        </select>
      </div>
      <div class="options__item">
        <p class="options__title">Y - axis</p>
        <select class="options__item--select" :value="selectedY" @change="setSelectY($event)">
          <option v-for="(item, index) in numList" :key="index" :value="item.value" :disabled="item.value === null">
            {{ item.name }}
          </option>
        </select>
      </div>
      <div class="options__item">
        <p class="options__title">Group</p>
        <select class="options__item--select" :value="selectedCluster" @change="setSelectCluster($event)">
          <option v-for="(item, index) in clusterList" :key="index" :value="item.value" :disabled="item.value === null">
            {{ item.name }}
          </option>
        </select>
      </div>

      <div class="options__item">
        <p class="options__title">Marker Size</p>
        <div class="options__buttons">
          <img class="options__item__button__minus" src="@/assets/button_minus.png" alt="-"
            v-on:click="markerSizeMinus" />
          <label class="options__item__degree">{{ markerSize }}</label>
          <img class="options__item__button__plus" src="@/assets/button_plus.png" alt="-" v-on:click="markerSizePlus" />
        </div>
      </div>
      <div class="options__item">
        <p class="options__title">Download Plot Image</p>
        <img class="downloadPlot_button" src="@/assets/download.png" alt="Save Plot" @click="downloadPlot" />
      </div>
      <div class="options__item">
        <button id="reset-button" @click="resetSelect">reset</button>
        <button id="apply-button" @click="cellselect" :disabled="disableApplyButton">
          {{ disableApplyButton ? "Select Applied" : "Select Apply " }}
        </button>
      </div>
    </div>
  </div>
</template>

<script>
import { getFileData, saveWorkflowNodeFile, readWorkflowNodeFile, deleteWorkflowNodeFile } from "@/api/index";
import Plotly from "plotly.js-dist-min";

export default {
  // components: {
  //   Plotly,
  // },
  data() {
    return {
      plotData: null,
      current_file: null,
      file_num: null,
      keys: null, // key들의 리스트 - string[]
      areNum: null, // 각 column 타입의 num 여부 리스트 - bool[]
      lines: null, // data의 rows의 리스트(keys는 제거됨) - object[][]
      chartType: "scattergl",
      chartMode: "markers",
      chartTitle: "",
      selectedX: null, // 선택된 X축의 column index "string"
      selectedY: null, // 선택된 Y축의 column index "string"
      selectedName: null, // 선택된 Name의 column index
      selectedCluster: null, // 선택된 Cluster의 column index
      numClusterConstraint: 150, // Cluster가 너무 다양하면 성능에 저하가 생기므로 제한을 둠
      clusterQuantile: 0,
      clusterContrast: 1,
      markerSize: 3, // 점의 사이즈
      showGrid: false,
      showZeroLine: true,
      showLine: false,
      showLabel: true,
      numList: [{ name: "None", value: null }], // number type으로 x,y축에 들어가기 적합한 자료들의 option list
      keyList: [{ name: "None", value: null }], // key list
      clusterList: [{ name: "None", value: null, isTooVarious: null }], // cluster가 되기 적합한 list, unique한 자료수가 적은 column의 list
      refreshOptionsLayout: true,
      outputX: [],
      outputY: [],
      selectedIndices: [],
      appliedSelectedIndices: [],
      disableApplyButton: false,
      scatterResult: null,
      workflowId: this.$route.query.workflow_id,
      nodeId: this.$route.query.node,
    };
  },
  async mounted() {
    try {
      Plotly.newPlot("plotly__scatter", {
        data: [{ type: this.chartType }],
        layout: {
          autosize: true,
          automargin: true,
          width: 600, // !--조정필요
          height: 570, // !--조정필요
          legend: {
            itemsizing: "constant",
            font: {
              size: 15,
            },
          },
          paper_bgcolor: "#ffffff00",
          plot_bgcolor: "#ffffff00",
        },
        config: {
          responsive: true,
          scrollZoom: true,
          displayModeBar: true,
          displaylogo: false,
        },
      });

      // 저장된 selectedCluster 값을 가져오기
      const nodeInfo = this.$store.getters.getWorkflowNodeInfo(this.nodeId);
      if (nodeInfo.data && nodeInfo.data.selectedCluster !== undefined) {
        this.selectedCluster = nodeInfo.data.selectedCluster;
      }
    } catch (error) {
      console.error("Error initializing Plotly:", error);
    }

    try {
      this.current_file = this.$store.getters.getWorkflowNodeFileInfo(this.nodeId);
    } catch (error) {
      console.error("Error getting current file info:", error);
    }

    let savedIndices;
    try {
      savedIndices = await readWorkflowNodeFile({
        id: this.workflowId,
        node_id: String(this.nodeId),
        node_name: "ScatterPlot",
        file_extension: "json",
      });
    } catch (error) {
      console.error("Error reading workflow node file:", error);
    }

    const tmpIndices = savedIndices?.data?.file_content;
    console.log(tmpIndices);
    console.log(savedIndices);
    if (tmpIndices != null) {
      if (tmpIndices.length > 0) {
        this.appliedSelectedIndices = tmpIndices;
        this.disableApplyButton = true;
      }
    }

    if (this.current_file !== "") {
      let scatterResult;
      try {
        scatterResult = await getFileData(this.current_file);
      } catch (error) {
        console.error("Error getting file data:", error);
      }

      if (scatterResult) {
        this.scatterResult = scatterResult.data;
        console.log(this.scatterResult);

        this.processScatterResult(this.scatterResult);

        for (let i = 0; i < this.lines.length; i++) {
          this.lines[i].push(i);
        }

        if (this.appliedSelectedIndices != null) {
          if (this.appliedSelectedIndices.length > 0) {
            this.lines = this.lines.filter((_, index) =>
              this.appliedSelectedIndices.includes(index)
            );
            this.disableApplyButton = true;
          }
        }

        this.numList = [{ name: "None", value: null }];
        for (let i = 0; i < this.keys.length; i++) {
          if (this.areNum[i] == true) {
            this.numList.push({ name: this.keys[i], value: i });
          }
        }

        this.keyList = [{ name: "None", value: null }];
        for (let i = 0; i < this.keys.length; i++) {
          this.keyList.push({ name: this.keys[i], value: i });
        }

        this.clusterList = [{ name: "None", value: null, isTooVarious: null }];
        for (let i = 1; i < this.keys.length; i++) {
          const countUnique = new Set(this.lines.map((x) => x[i])).size;
          if (countUnique < this.numClusterConstraint) {
            this.clusterList.push({
              name: this.keys[i],
              value: i,
              isTooVarious: false,
            });
          } else {
            this.clusterList.push({
              name: this.keys[i],
              value: i,
              isTooVarious: true,
            });
          }
        }

        // Auto-select X and Y columns if they exist
        if (this.keys.indexOf("X") != -1) {
          this.selectedX = this.keys.indexOf("X");
        } else if (this.numList.length > 1) {
          // If no "X" column, select the first numeric column
          this.selectedX = this.numList[1].value;
        }
        
        if (this.keys.indexOf("Y") != -1) {
          this.selectedY = this.keys.indexOf("Y");
        } else if (this.numList.length > 2) {
          // If no "Y" column, select the second numeric column
          this.selectedY = this.numList[2].value;
        } else if (this.numList.length > 1) {
          // If only one numeric column exists, use it for Y as well
          this.selectedY = this.numList[1].value;
        }
        
        // if (this.keys.indexOf("leiden") != -1) {
        //   this.selectedCluster = this.keys.indexOf("leiden");
        // }
        if (this.keys.indexOf("INDEX") != -1) {
          this.selectedName = this.keys.indexOf("INDEX");
        }

        try {
          this.updateChart();
        } catch (error) {
          console.error("Error updating chart:", error);
        }
      }
    }

    this.refreshOptionsLayout = false;
  },
  methods: {
    // 차트 업데이트
    updateChart() {
      this.outputX = this.lines.map((x) => x[this.selectedX]);
      this.outputY = this.lines.map((x) => x[this.selectedY]);
      if (this.selectedCluster) {
        if (this.clusterList[this.selectedCluster].isTooVarious == false) {
          const clusterList = [
            ...new Set(this.lines.map((x) => x[this.selectedCluster])),
          ];
          var traces = [];
          // cluster 개수만큼 traces 생성
          for (let i = 0; i < clusterList.length; i++) {
            traces.push({
              x: [],
              y: [],
              text: [],
              name: clusterList[i] ?? "Undefined",
              type: this.chartType,
              mode: this.chartMode,
              marker: { size: this.markerSize },
            });
          }
          // lines를 순회하며 cluster에 맞는 traces에 x,y,text 값을 기입
          for (let i = 0; i < this.lines.length; i++) {
            traces[
              clusterList.indexOf(this.lines[i][this.selectedCluster])
            ].x.push(this.lines[i][this.selectedX]);
            traces[
              clusterList.indexOf(this.lines[i][this.selectedCluster])
            ].y.push(this.lines[i][this.selectedY]);
            traces[
              clusterList.indexOf(this.lines[i][this.selectedCluster])
            ].text.push(this.lines[i][this.selectedName]);
          }
          Plotly.newPlot("plotly__scatter", {
            data: traces,
            layout: {
              title: this.chartTitle,
              width: 600, // !--조정필요
              height: 570, // !--조정필요
              xaxis: {
                showgrid: this.showGrid,
                showticklabels: this.showLabel,
                zeroline: this.showZeroLine,
                showline: this.showLine,
              },
              yaxis: {
                showgrid: this.showGrid,
                showticklabels: this.showLabel,
                zeroline: this.showZeroLine,
                showline: this.showLine,
              },
              legend: {
                itemsizing: "constant",
                font: {
                  size: 15,
                },
              },
              paper_bgcolor: "#ffffff00",
              plot_bgcolor: "#ffffff00",
            },
          });
        } else {
          Plotly.newPlot("plotly__scatter", {
            data: [
              {
                x: this.lines.map((x) => x[this.selectedX]),
                y: this.lines.map((x) => x[this.selectedY]),
                text: this.lines.map((x) => x[this.selectedName]),
                type: this.chartType,
                mode: this.chartMode,
                marker: {
                  size: this.markerSize,
                  color: this.lines.map(
                    (x) =>
                      x[this.selectedCluster] * this.clusterContrast +
                      this.clusterQuantile
                  ),
                },
              },
            ],
            layout: {
              title: this.chartTitle,
              width: 600, // !--조정필요
              height: 570, // !--조정필요
              xaxis: {
                showgrid: this.showGrid,
                showticklabels: this.showLabel,
                zeroline: this.showZeroLine,
                showline: this.showLine,
              },
              yaxis: {
                showgrid: this.showGrid,
                showticklabels: this.showLabel,
                zeroline: this.showZeroLine,
                showline: this.showLine,
              },
              legend: {
                itemsizing: "constant",
                font: {
                  size: 15,
                },
              },
              paper_bgcolor: "#ffffff00",
              plot_bgcolor: "#ffffff00",
            },
          });
        }
      } else {
        Plotly.newPlot("plotly__scatter", {
          data: [
            {
              x: this.lines.map((x) => x[this.selectedX]),
              y: this.lines.map((x) => x[this.selectedY]),
              text: this.lines.map((x) => x[this.selectedName]),
              type: this.chartType,
              mode: this.chartMode,
              marker: { size: this.markerSize },
            },
          ],
          layout: {
            title: this.chartTitle,
            width: 600, // !--조정필요
            height: 570, // !--조정필요
            xaxis: {
              showgrid: this.showGrid,
              showticklabels: this.showLabel,
              zeroline: this.showZeroLine,
              showline: this.showLine,
            },
            yaxis: {
              showgrid: this.showGrid,
              showticklabels: this.showLabel,
              zeroline: this.showZeroLine,
              showline: this.showLine,
            },
            legend: {
              itemsizing: "constant",
              font: {
                size: 15,
              },
            },
            paper_bgcolor: "#ffffff00",
            plot_bgcolor: "#ffffff00",
          },
        });
      }
      var graphDiv = document.getElementById("plotly__scatter");
      graphDiv.on("plotly_selected", (eventData) => {
        this.selectedIndices = [];
        if (eventData.points.length > 0) {
          for (let i = 0; i < eventData.points.length; i++) {
            this.selectedIndices.push(parseInt(eventData.points[i].text));
          }
        }
      });
    },
    setSelectX(event) {
      this.selectedX = event.target.value;
      this.updateChart();
    },
    setSelectY(event) {
      this.selectedY = event.target.value;
      this.updateChart();
    },
    setSelectName(event) {
      this.selectedName = event.target.value;
      this.updateChart();
    },
    setSelectCluster(event) {
      this.selectedCluster = parseInt(event.target.value);

      // selectedCluster 값을 store에 저장
      const dataObject = {
        selectedCluster: this.selectedCluster
      };
      this.$store.commit("setWorkflowNodeDataObject", { nodeId: this.nodeId, dataObject });

      this.updateChart();
    },
    titleChangeFunc(event) {
      this.chartTitle = event.target.value;
      this.updateChart();
    },
    clusterQuantileMinus() {
      this.clusterQuantile--;
      this.updateChart();
    },

    clusterQuantilePlus() {
      this.clusterQuantile++;
      this.updateChart();
    },
    clusterQuantileReset() {
      this.clusterQuantile = 0;
      this.updateChart();
    },
    clusterContrastMinus() {
      if (this.clusterContrast > 0) {
        this.clusterContrast--;
        this.updateChart();
      }
    },
    clusterContrastPlus() {
      this.clusterContrast++;
      this.updateChart();
    },
    clusterContrastReset() {
      this.clusterContrast = 1;
      this.updateChart();
    },
    markerSizeMinus() {
      if (this.markerSize > 1) {
        this.markerSize--;
        this.updateChart();
      }
    },
    markerSizePlus() {
      if (this.markerSize < 21) {
        this.markerSize++;
        this.updateChart();
      }
    },
    markerSizeReset() {
      this.markerSize = 2;
      this.updateChart();
    },
    switchShowGrid() {
      this.showGrid = !this.showGrid;
      this.updateChart();
    },
    switchShowZeroLine() {
      this.showZeroLine = !this.showZeroLine;
      this.updateChart();
    },
    switchShowLine() {
      this.showLine = !this.showLine;
      this.updateChart();
    },
    switchShowLabel() {
      this.showLabel = !this.showLabel;
      this.updateChart();
    },
    downloadPlot() {
      Plotly.downloadImage("plotly__scatter", {
        format: "png",
        width: 600,
        height: 570,
        filename: "CELLCRAFT_Plot",
      });
    },
    async cellselect() {
      if (this.selectedIndices.length > 0) {
        // this.$store.commit("setSelectedIndices", [
        //   this.keys[this.selectedCluster],
        //   this.selectedIndices,
        // ]);
        const selectIndicesInfo = await saveWorkflowNodeFile({
          id: this.workflowId,
          node_id: String(this.nodeId),
          node_name: "ScatterPlot",
          file_content: this.selectedIndices,
          file_extension: "json",
        });
        console.log(selectIndicesInfo.data);
        const dataObject = {
          lasso_file_path: selectIndicesInfo.data.file_path
        }
        this.$store.commit("setWorkflowNodeDataObject", { nodeId: this.nodeId, dataObject });
        this.appliedSelectedIndices = this.selectedIndices;
        this.lines = this.lines.filter((_, index) =>
          this.appliedSelectedIndices.includes(index)
        );
        this.updateChart();
        this.disableApplyButton = true;
      }
    },
    async resetSelect() {
      try {
        // this.$store.commit("setSelectedIndices", ["", []]);
        const resetIndicesInfo = await deleteWorkflowNodeFile({
          id: this.workflowId,
          node_id: String(this.nodeId),
          node_name: "ScatterPlot",
          file_extension: "json",
        });
        console.log(resetIndicesInfo);
        this.appliedSelectedIndices = null;
        this.selectedIndices = [];
        this.disableApplyButton = false;
        this.markerSize = 3;

        this.processScatterResult(this.scatterResult);

        for (let i = 0; i < this.lines.length; i++) {
          this.lines[i].push(i);
        }

        this.numList = [{ name: "None", value: null }];
        for (let i = 0; i < this.keys.length; i++) {
          if (this.areNum[i] == true) {
            this.numList.push({ name: this.keys[i], value: i });
          }
        }

        this.keyList = [{ name: "None", value: null }];
        for (let i = 0; i < this.keys.length; i++) {
          this.keyList.push({ name: this.keys[i], value: i });
        }

        this.clusterList = [{ name: "None", value: null, isTooVarious: null }];
        for (let i = 1; i < this.keys.length; i++) {
          const countUnique = new Set(this.lines.map((x) => x[i])).size;
          if (countUnique < this.numClusterConstraint) {
            this.clusterList.push({
              name: this.keys[i],
              value: i,
              isTooVarious: false,
            });
          } else {
            this.clusterList.push({
              name: this.keys[i],
              value: i,
              isTooVarious: true,
            });
          }
        }

        if (this.keys.indexOf("X") != -1) {
          this.selectedX = this.keys.indexOf("X");
        }
        if (this.keys.indexOf("Y") != -1) {
          this.selectedY = this.keys.indexOf("Y");
        }
        // if (this.keys.indexOf("leiden") != -1) {
        //   this.selectedCluster = this.keys.indexOf("leiden");
        // }

        if (this.keys.indexOf("INDEX") != -1) {
          this.selectedName = this.keys.indexOf("INDEX");
        }
        this.updateChart();
      } catch (error) {
        console.error("Error resetting selection:", error);
      }
    },
    processScatterResult(scatterResult) {
      // Check if the data is a string (CSV) or an array (JSON)
      if (typeof scatterResult === 'string') {
        // Parse CSV data
        this.parseCSVData(scatterResult);
      } else if (Array.isArray(scatterResult)) {
        // Process JSON data as before
        // Extract keys from the first object
        this.keys = Object.keys(scatterResult[0]);
        this.keys.push("INDEX");

        // Map the scatterResult to an array of values and add the index
        this.lines = scatterResult.map((item, index) => {
          const values = Object.values(item);
          values.push(index);
          return values;
        });

        // Determine which columns are numeric
        this.areNum = this.lines[0].map((x) => !isNaN(x));
      } else {
        console.error("Unsupported data format");
      }
    },
    
    parseCSVData(csvString) {
      // Split CSV string into lines
      const lines = csvString.trim().split('\n');
      
      // Extract headers (first line)
      const headers = lines[0].split(',').map(h => h.trim());
      
      // Remove empty header if exists (for index column)
      if (headers[0] === '') {
        headers[0] = 'Row_Index';
      }
      
      this.keys = [...headers, "INDEX"];
      
      // Parse data rows
      this.lines = [];
      for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',').map(v => {
          const trimmed = v.trim();
          // Try to convert to number if possible
          const num = parseFloat(trimmed);
          return isNaN(num) ? trimmed : num;
        });
        
        // Skip the first value if it's a row index
        if (headers[0] === 'Row_Index') {
          values.shift();
        }
        
        values.push(i - 1); // Add INDEX
        this.lines.push(values);
      }
      
      // Determine which columns are numeric
      if (this.lines.length > 0) {
        this.areNum = this.lines[0].map((x) => !isNaN(x) && typeof x === 'number');
      }
    }
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

.plotly-layout {
  width: 70%;
  height: 95%;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
  border-radius: 1rem;
  margin: 1%;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
}

.options-layout {
  width: 30%;
  height: 95%;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
  gap: 0.5rem;
  z-index: 9997;
}

.options__textInput {
  width: 80%;
  padding: 0.5rem;
  color: black;
  border-radius: 0.5rem;
  border: 1px solid #e7eaff;
  font-size: medium;
  text-align: center;
}

.options__item {
  width: 90%;
  display: flex;
  align-items: center;
  justify-content: space-between;
}

.options__title {
  font-weight: 600;
  color: rgb(55, 55, 55);
}

.options__item--select {
  width: 8rem;
  padding: 0.5rem;
  border-radius: 8px;
  border-color: #e7eaff;
  color: #545454;
}

.options__item--select option:disabled {
  color: #ccc;
  background-color: #f5f5f5;
  cursor: not-allowed;
}

.options__buttons {
  width: 7rem;
  display: flex;
  align-items: center;
  justify-content: space-evenly;
}

.options__item__button__minus,
.options__item__button__plus,
.options__item__degree {
  width: 1rem;
  height: 1rem;
}

button {
  background-color: #ffffff;
  border: 1px solid #999999;
  border-radius: 0.3rem;
  color: #333333;
  /* font-size: 1rem;
  padding: 0.2rem 1rem; */
  text-align: center;
  text-decoration: none;
  transition: background-color 0.3s ease;
  margin: 0.2rem 0.1rem;
}

button:hover {
  background-color: #e7eaff;
  border-color: #b3b3b3;
}

/* The switch - the box around the slider_button */
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
  height: 15px;
  width: 15px;
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

.downloadPlot_button {
  width: 1.5rem;
  height: 1.5rem;
  opacity: 0.8;
}

.downloadPlot_button:hover {
  opacity: 1;
  cursor: pointer;
}

#apply-button {
  background-color: #2d2fbf;
  /* 버튼 배경색 */
  width: 9rem;
  color: white;
  /* 글자색 */
  padding: 10px 0px;
  /* 상하 10px, 좌우 20px의 여백 */
  border: none;
  /* 테두리 없앰 */
  border-radius: 4px;
  /* 테두리 모서리 둥글게 */
  cursor: pointer;
  /* 마우스 오버 시 커서 변경 */
  font-size: 16px;
  /* 글자 크기 */
  transition: background-color 0.3s;
  /* 배경색 변경시 트랜지션 효과 */
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.4);
}

#apply-button:hover {
  background-color: #4655ff;
  /* 마우스 오버시 버튼의 배경색 변경 */
}

#apply-button:disabled {
  background-color: #ccc;
  /* 비활성화 상태의 배경색 */
  color: #666;
  /* 비활성화 상태의 글자색 */
  cursor: not-allowed;
  /* 비활성화 상태에서의 커서 */
}

#reset-button {
  background-color: #616161;
  /* 버튼 배경색 */
  width: 4rem;
  color: white;
  /* 글자색 */
  padding: 10px 0px;
  /* 상하 10px, 좌우 20px의 여백 */
  border: none;
  /* 테두리 없앰 */
  border-radius: 4px;
  /* 테두리 모서리 둥글게 */
  cursor: pointer;
  /* 마우스 오버 시 커서 변경 */
  font-size: 16px;
  /* 글자 크기 */
  transition: background-color 0.3s;
  /* 배경색 변경시 트랜지션 효과 */
  margin-right: 0.5rem;
}

#reset-button:hover {
  background-color: #797979;
  /* 마우스 오버시 버튼의 배경색 변경 */
}
</style>
