<template>
  <div id="layout">
    <div class="plotly-layout">
      <div id="plotly-heatmap"></div>
    </div>
    <div class="options-layout">
      <div class="options__item">
        sif&nbsp;
        <select
          class="options__item__select"
          :value="selectedSif"
          @change="setSelectedSif($event)"
        >
          <option v-for="sif in resultSifs" :key="sif" :value="sif">
            {{ sif | splitUnderScore }}
          </option>
        </select>
      </div>
      <div class="options__item">
        Degree&nbsp;
        <select
          class="options__item__select"
          :value="selectedDegree"
          @change="setSelectedDegree($event)"
        >
          <option
            v-for="(item, index) in numList"
            :key="index"
            :value="item.value"
          >
            {{ item.name }}
          </option>
        </select>
      </div>
      <!-- <div class="options__item">
        p_min&nbsp;
        <img
          class="options__item__button__minus"
          src="@/assets/button_minus.png"
          alt="-"
          v-on:click="pMinMinus"
        />
        <label class="options__item__degree">{{ p_min / 10 }}</label>
        <img
          class="options__item__button__plus"
          src="@/assets/button_plus.png"
          alt="-"
          v-on:click="pMinPlus"
        />
      </div>
      <div class="options__item">
        p_max&nbsp;
        <img
          class="options__item__button__minus"
          src="@/assets/button_minus.png"
          alt="-"
          v-on:click="pMaxMinus"
        />
        <label class="options__item__degree">{{ p_max / 10 }}</label>
        <img
          class="options__item__button__plus"
          src="@/assets/button_plus.png"
          alt="-"
          v-on:click="pMaxPlus"
        />
      </div> -->
    </div>
  </div>
</template>

<script>
import { getResultFileOne, getResultFile } from "@/api/index";
// Importing the Plotly library
import Plotly from "plotly.js-dist-min";

export default {
  data() {
    return {
      selectedDegree: 0,
      numList: [
        { name: "Indegree", value: 0 },
        { name: "Outdegree", value: 2 },
      ],
      p_min: -7,
      p_max: 7,
      resultSifs: [],
      selectedSif: "",
      current_file: "",
      expression_file: "",
      pseudo_time_file: "",
      cell_select_file: "",
    };
  },
  name: "PlotlyHeatMap",

  async mounted() {
    // option_file_name.algorithmOptions.optionFilePath이 null이면 error 띄우기
    const option_file_name = this.$store.getters.getCurrentLinkedNodes;
    console.log(option_file_name[0].algorithmOptions.optionFilePath);
    if (option_file_name[0].algorithmOptions.optionFilePath === null) {
      alert("Please select an option file.");
    }

    try {
      this.current_file = this.filterAndAddSuffix(
        this.$store.getters.getCurrentFile.file
      );
      const result = await getResultFile({
        file_name: this.current_file,
        option_file_name: option_file_name[0].algorithmOptions.optionFilePath,
      });
      const resultsList = result.data.result_files;
      console.log(resultsList);
      // resultsList에서 .sif 파일만 추출
      // 추출할 때, "_"로 구분된 문자열을 배열로 변환하여 마지막 요소만 추출
      // 마지막 요소가 "sif"인지 확인
      // "sif"이면 배열을 문자열로 변환하여 resultSifs에 추가
      for (var i = 0; i < resultsList.length; i++) {
        const segments = resultsList[i].split("_");
        // segments 문자열 안에 sif가 포함되어 있으면 segments[segments.length - 1] push
        if (segments[segments.length - 1].includes("sif")) {
          this.resultSifs.push(resultsList[i]);
        }
      }

      // resultSifs에서 trimIndirect가 포함된 파일명 this.selectedSif에 할당
      // for (i = 0; i < this.resultSifs.length; i++) {
      //   if (this.resultSifs[i].includes("trimIndirect")) {
      //     this.selectedSif = this.resultSifs[i];
      //   }
      // }

      // expMatrix, pseudotime, cellSelect가 포함된 파일명 추출
      // 그냥 resultsList에서 순회하면서 해당 문자열이 포함된 파일명을 찾아서
      // this.expression_file, this.pseudo_time_file, this.cell_select_file에 할당
      for (i = 0; i < resultsList.length; i++) {
        if (resultsList[i].includes("expMatrix")) {
          this.expression_file = resultsList[i];
          console.log(this.expression_file);
        } else if (resultsList[i].includes("pseudotime")) {
          this.pseudo_time_file = resultsList[i];
          console.log(this.pseudo_time_file);
        } else if (resultsList[i].includes("cellSelect")) {
          this.cell_select_file = resultsList[i];
          console.log(this.cell_select_file);
        }
      }
    } catch (error) {
      console.log(error);
    }

    await this.updatePlot(this.selectedDegree); // indegree는 0 outdegree는 2
  },
  filters: {
    // 문자열 안에 "_"가 포함되어 있으면 "_"로 구분된 문자열을 중 마지막 요소를만 반환
    // "_"가 포함되어 있지 않으면 원래 문자열 반환
    splitUnderScore(inputString) {
      // Check if the inputString contains an underscore
      if (inputString.includes("_")) {
        // "_"로 구분된 문자열을 배열로 변환
        const segments = inputString.split("_");
        // 마지막 요소를 반환
        const fileName = segments[segments.length - 1];
        return fileName;
      }
      // If no underscore found, return the original string
      return inputString;
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
    async updatePlot(rowNum = 0) {
      var [
        algorithmResult,
        expressionResult,
        pseudoTimeResult,
        cellSelectResult,
      ] = ["", "", "", ""];
      try {
        // 비동기 처리 병렬화
        [
          algorithmResult,
          expressionResult,
          pseudoTimeResult,
          cellSelectResult,
        ] = await Promise.all([
          getResultFileOne(this.selectedSif),
          getResultFileOne(this.expression_file),
          getResultFileOne(this.pseudo_time_file),
          getResultFileOne(this.cell_select_file),
        ]);
      } catch (error) {
        console.log(error);
      }
      // y축 데이터 처리 최적화
      const frequencyMap = this.processAlgorithmResult(
        algorithmResult.data,
        rowNum
      );

      // z축 데이터 처리 최적화
      const { expressionTwoDimensionalList, expressionFirstRow } =
        this.processExpressionResult(expressionResult.data);
      this.appendPseudoTimes(
        expressionTwoDimensionalList,
        pseudoTimeResult.data
      );
      const filteredExpressionList = this.filterExpressionList(
        expressionTwoDimensionalList,
        cellSelectResult.data
      );

      var zMatrix = this.computeZMatrix(
        expressionFirstRow,
        filteredExpressionList,
        frequencyMap
      );

      // Plotly 데이터 준비
      const data = this.preparePlotData(frequencyMap, zMatrix);
      // Plotly 레이아웃 설정
      const layout = this.getPlotLayout();
      // Plotly 차트 렌더링
      Plotly.newPlot("plotly-heatmap", data, layout, {
        responsive: true,
        scrollZoom: true,
        displayModeBar: true,
        displaylogo: false,
      });
    },
    async setSelectedDegree(event) {
      this.selectedDegree = event.target.value;
      await this.updatePlot(this.selectedDegree);
    },
    async setSelectedSif(event) {
      this.selectedSif = event.target.value;
      console.log(this.selectedSif);
      await this.updatePlot(this.selectedDegree);
    },
    localRegression(data, bandwidth = 75) {
      const n = data.length;
      const smoothedData = [];

      for (let i = 0; i < n; i++) {
        let weightedSum = 0;
        let weightSum = 0;

        for (
          let j = Math.max(0, i - bandwidth);
          j < Math.min(n, i + bandwidth + 1);
          j++
        ) {
          const weight = 1; // 간단한 예시로 가중치를 1로 고정
          weightedSum += data[j] * weight;
          weightSum += weight;
        }

        const smoothedValue = weightedSum / weightSum;
        smoothedData.push(smoothedValue);
      }
      const nowMax = Math.max(...smoothedData);
      const nowMin = Math.min(...smoothedData);
      const nowDiff = nowMax - nowMin;
      let multipliedArray = smoothedData.map(
        (element) =>
          (element - (nowMax + nowMin) / 2) *
            ((this.p_max - this.p_min) / nowDiff / 10) +
          (this.p_max + this.p_min) / 20
      );
      return multipliedArray;
    },
    processAlgorithmResult(data, rowNum) {
      const lines = data.split("\n");
      let frequencyMap = {};

      lines.forEach((line) => {
        const elements = line.split("\t");
        const element = elements[rowNum];
        if (element) {
          frequencyMap[element] = (frequencyMap[element] || 0) + 1;
        }
      });

      return frequencyMap;
    },
    processExpressionResult(data) {
      const lines = data.split("\n");
      let twoDimensionalList = lines.map((line) => line.split(","));

      return {
        expressionTwoDimensionalList: twoDimensionalList,
        expressionFirstRow: twoDimensionalList[0],
      };
    },
    appendPseudoTimes(twoDimensionalList, pseudoTimeData) {
      const pseudoTimes = pseudoTimeData.split("\n");
      pseudoTimes.forEach((time, index) => {
        if (twoDimensionalList[index]) {
          twoDimensionalList[index].push(time);
        }
      });
    },
    filterExpressionList(twoDimensionalList, cellSelectData) {
      const cellSelect = cellSelectData.split("\r\n");
      let filteredList = [];

      twoDimensionalList.forEach((row, index) => {
        if (cellSelect[index] === "1") {
          filteredList.push(row);
        }
      });

      filteredList.sort((a, b) => a[a.length - 1] - b[b.length - 1]);
      return filteredList;
    },
    computeZMatrix(firstRow, expressionList, frequencyMap) {
      let zMatrix = [];

      Object.keys(frequencyMap).forEach((key) => {
        const cellIndex = firstRow.indexOf(key);
        let row = [];
        expressionList.forEach((expressionRow) => {
          row.push(expressionRow[cellIndex + 1]);
        });
        zMatrix.push(this.localRegression(row));
      });

      return zMatrix;
    },
    preparePlotData(frequencyMap, zMatrix) {
      return [
        {
          y: Object.keys(frequencyMap),
          z: zMatrix,
          type: "heatmap",
          orientation: "h",
        },
      ];
    },
    getPlotLayout() {
      return {
        width: 600,
        height: 570,
      };
    },
    // pMinMinus() {
    //   if (this.p_min <= -7) {
    //     return;
    //   }
    //   this.p_min -= 1;
    //   this.updatePlot(this.selectedDegree);
    // },
    // pMinPlus() {
    //   if (this.p_min >= -1) {
    //     return;
    //   }
    //   this.p_min += 1;
    //   this.updatePlot(this.selectedDegree);
    // },
    // pMaxMinus() {
    //   if (this.p_max <= 1) {
    //     return;
    //   }
    //   this.p_max -= 1;
    //   this.updatePlot(this.selectedDegree);
    // },
    // pMaxPlus() {
    //   if (this.p_max >= 7) {
    //     return;
    //   }
    //   this.p_max += 1;
    //   this.updatePlot(this.selectedDegree);
    // },
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
  /* background-color: blue; */
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
  padding: 1rem;
  margin: 2.5%;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
}
.options-layout {
  width: 25%;
  height: 95%;
  padding-right: 1%;
  /* background-color: red; */
  display: flex;
  align-items: flex-start;
  justify-content: center;
  flex-direction: column;
  z-index: 9997;
}
.options__textInput {
  color: black;
  padding: 5px;
  left: 10px;
  border-radius: 10px;
  border-color: #e7eaff;
  font-size: medium;
  text-align: center;
  margin-bottom: 5px;
}
.options__item {
  /* margin: auto; */
  font-weight: 600;
  color: rgb(55, 55, 55);
  /* padding: 2% 0; */
  left: 10px;
  margin-top: 15px;
  text-transform: capitalize;
}
.options__item__select {
  position: absolute;
  margin-top: -10px;
  right: 10px;
  padding: 5px;
  border-radius: 8px;
  width: 140px;
  border-color: #e7eaff;
  color: #545454;
}
.options__item__button__minus {
  position: absolute;
  margin-top: -3px;
  right: 90px;
  width: 15px;
  height: 15px;
}
.options__item__button__plus {
  position: absolute;
  margin-top: -3px;
  right: 19px;
  width: 15px;
  height: 15px;
}
.options__item__degree {
  position: absolute;
  right: 50px;
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

.downloadPlot_button {
  position: absolute;
  /* top: -0.2rem; */
  right: 1.5rem;
  margin-top: -0.2rem;
  width: 1.5rem;
  height: 1.5rem;
  opacity: 0.8;
}
.downloadPlot_button:hover {
  opacity: 1;
  cursor: pointer;
}

#apply-button {
  background-color: #2d2fbf; /* 버튼 배경색 */
  width: 8.5rem;
  color: white; /* 글자색 */
  padding: 10px 0px; /* 상하 10px, 좌우 20px의 여백 */
  border: none; /* 테두리 없앰 */
  border-radius: 4px; /* 테두리 모서리 둥글게 */
  cursor: pointer; /* 마우스 오버 시 커서 변경 */
  font-size: 16px; /* 글자 크기 */
  transition: background-color 0.3s; /* 배경색 변경시 트랜지션 효과 */
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.4);
}

#apply-button:hover {
  background-color: #4655ff; /* 마우스 오버시 버튼의 배경색 변경 */
}

#apply-button:disabled {
  background-color: #ccc; /* 비활성화 상태의 배경색 */
  color: #666; /* 비활성화 상태의 글자색 */
  cursor: not-allowed; /* 비활성화 상태에서의 커서 */
}

#reset-button {
  background-color: #616161; /* 버튼 배경색 */
  width: 3.5rem;
  color: white; /* 글자색 */
  padding: 10px 0px; /* 상하 10px, 좌우 20px의 여백 */
  border: none; /* 테두리 없앰 */
  border-radius: 4px; /* 테두리 모서리 둥글게 */
  cursor: pointer; /* 마우스 오버 시 커서 변경 */
  font-size: 16px; /* 글자 크기 */
  transition: background-color 0.3s; /* 배경색 변경시 트랜지션 효과 */
  /* margin-left: 10px; */
  margin-right: 0.5rem;
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.4);
}

#reset-button:hover {
  background-color: #797979; /* 마우스 오버시 버튼의 배경색 변경 */
}
</style>
