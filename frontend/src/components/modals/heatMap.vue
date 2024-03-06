<template>
  <div id="layout">
    <div class="plotly-layout">
      <div id="plotly-heatmap"></div>
      <div v-if="isLoading" class="loading-layout">
        <span> </span>
      </div>
      <div v-else-if="errorOccured">
        <span> SOME ERROR OCCURED</span>
      </div>
      <div v-else-if="!plotReady">
        <span> NO DATA FOR HEATMAP</span>
      </div>
    </div>
    <div class="options-layout">
      <div class="options__item">
        sif&nbsp;
        <select class="options__item__select" v-model="selectedSif">
          <option v-for="sif in resultSifs" :key="sif" :value="sif">
            {{ sif | splitUnderScore }}
          </option>
        </select>
      </div>
      <div class="options__item">
        Degree&nbsp;
        <select
          class="options__item__select"
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
      <div class="options__item">
        <button id="apply-button" @click="applySelect">
          {{ "Select Apply " }}
        </button>
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
      plotReady: false,
      isLoading: false,
      errorOccured: false,
    };
  },
  name: "PlotlyHeatMap",

  async mounted() {
    // option_file_name.algorithmOptions.optionFilePath이 null이면 error 띄우기
    const option_file_name = this.$store.getters.getCurrentLinkedNodes;
    console.log(option_file_name[0].algorithmOptions.optionFilePath);

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
      for (i = 0; i < this.resultSifs.length; i++) {
        if (this.resultSifs[i].includes("trimIndirect")) {
          this.selectedSif = this.resultSifs[i];
        }
      }

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
    async applySelect() {
      this.isLoading = true;
      await this.updatePlot();
      this.isLoading = false;
    },
    async updatePlot() {
      try {
        // ---- y축 ------
        const algorithmResult = await getResultFileOne(this.selectedSif);
        // 줄바꿈 문자("\n")를 기준으로 문자열을 분할하여 각 라인을 배열로 변환
        let lines = algorithmResult.data.split("\n");

        // 각 라인을 탭 문자("\t")를 기준으로 분할하여 2차원 배열 생성
        var twoDimensionalList = lines.map(function (line) {
          return line.split("\t");
        });

        // 결과 출력 (옵션)

        var frequencyMap = {};
        const tmpDegree = this.selectedDegree;
        twoDimensionalList.forEach(function (row) {
          var element = row[tmpDegree]; // 첫 번째 열의 요소
          // 해당 요소가 frequencyMap에 이미 존재하면 카운트 증가, 없으면 1로 초기화
          if (frequencyMap[element]) {
            frequencyMap[element]++;
          } else {
            frequencyMap[element] = 1;
          }
        });

        if (Object.prototype.hasOwnProperty.call(frequencyMap, "")) {
          delete frequencyMap[""];
        }

        // 결과 출력 (옵션)

        // 객체를 키-값 쌍의 배열로 변환
        var items = Object.keys(frequencyMap).map(function (key) {
          return [key, frequencyMap[key]];
        });

        // 배열을 값(value)에 따라 내림차순으로 정렬
        items.sort(function (first, second) {
          return first[1] - second[1];
        });

        // 정렬된 배열을 다시 객체로 변환 (선택적)
        var sortedFrequencyMap = {};
        items.forEach(function (item) {
          sortedFrequencyMap[item[0]] = item[1];
        });

        //sortedFrequencyMap에 key가 undefined인 요소가 있으면 삭제
        if (
          Object.prototype.hasOwnProperty.call(sortedFrequencyMap, "undefined")
        ) {
          delete sortedFrequencyMap["undefined"];
        }

        // 결과 출력
        // ---- y축 ------ 끝

        // ---- z ------
        const expressionResult = await getResultFileOne(this.expression_file);

        var expressionLines = expressionResult.data.split("\n");

        // 각 라인을 탭 문자("\t")를 기준으로 분할하여 2차원 배열 생성
        var expressionTwoDimensionalList = expressionLines.map(function (line) {
          return line.split(",");
        });

        // 결과 출력 (옵션)

        const pseudoTimeResult = await getResultFileOne(this.pseudo_time_file);
        const pseudoTimes = pseudoTimeResult.data.split("\n");

        for (var i = 0; i < pseudoTimes.length; i++) {
          expressionTwoDimensionalList[i].push(pseudoTimes[i]);
        }

        const cellSelectResult = await getResultFileOne(this.cell_select_file);
        const cellSelect = cellSelectResult.data.split("\r\n");

        // 새로운 2차원 배열 생성
        var filteredExpressionList = [];
        var expressionFirstRow = expressionTwoDimensionalList[0];
        // expressionTwoDimensionalList를 순회하면서
        // cellSelect에서 해당 인덱스가 "1"인 행만 새 배열에 추가
        for (i = 0; i < expressionTwoDimensionalList.length; i++) {
          if (cellSelect[i] === "1") {
            filteredExpressionList.push(expressionTwoDimensionalList[i + 1]);
          }
        }

        filteredExpressionList.sort(function (a, b) {
          // 각 행의 마지막 요소를 비교하여 정렬
          return a[a.length - 1] - b[b.length - 1];
        });

        var zMatrix = [];

        for (i = 0; i < Object.keys(sortedFrequencyMap).length; i++) {
          var cellIndex = expressionFirstRow.indexOf(
            Object.keys(sortedFrequencyMap)[i]
          );
          var row = [];
          for (var j = 0; j < filteredExpressionList.length; j++) {
            row.push(filteredExpressionList[j][cellIndex + 1]);
          }
          zMatrix.push(this.localRegression(row));
        }
        console.log(Object.keys(sortedFrequencyMap));
        console.log(zMatrix);
        // Data for the bar plot
        const data = [
          {
            // x: items.map((item) => item[1]),
            // y: items.map((item) => item[0]),
            // z: items.map((item) => item[1]),
            // x: [1, 2, 3],
            y: Object.keys(sortedFrequencyMap),
            z: zMatrix, // Heatmap values
            type: "heatmap",
            orientation: "h",
          },
        ];

        // Layout configuration
        const layout = {
          // title: "Simple Bar Plot",
          // xaxis: {
          //   title: "Categories",
          // },
          // yaxis: {
          //   title: "Values",
          // },
          width: 600, // !--조정필요
          height: 570, // !--조정필요
        };
        const config = {
          responsive: true,
          scrollZoom: true,
          displayModeBar: true,
          displaylogo: false,
        };
        // Rendering the bar plot
        Plotly.newPlot("plotly-heatmap", data, layout, config);
        this.plotReady = true;
      } catch (error) {
        this.plotReady = false;
        this.errorOccured = true;
        console.log(error);
      }
      console.log(
        111,
        this.selectedSif,
        222,
        this.current_file,
        333,
        this.expression_file,
        444,
        this.pseudo_time_file,
        555,
        this.cell_select_file
      );
    },
    async setSelectedDegree(event) {
      this.selectedDegree = event.target.value;
      if (typeof this.selectedDegree === "string") {
        this.selectedDegree = parseInt(this.selectedDegree);
      }
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
    pMinMinus() {
      if (this.p_min <= -7) {
        return;
      }
      this.p_min -= 1;
      this.updatePlot(this.selectedDegree);
    },
    pMinPlus() {
      if (this.p_min >= -1) {
        return;
      }
      this.p_min += 1;
      this.updatePlot(this.selectedDegree);
    },
    pMaxMinus() {
      if (this.p_max <= 1) {
        return;
      }
      this.p_max -= 1;
      this.updatePlot(this.selectedDegree);
    },
    pMaxPlus() {
      if (this.p_max >= 7) {
        return;
      }
      this.p_max += 1;
      this.updatePlot(this.selectedDegree);
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
  width: 13.5rem;
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

.loading-layout {
  display: flex;
  align-items: center;
  justify-content: center;
}

.loading-layout span {
  border: 4px solid #f3f3f3;
  border-top: 4px solid #3498db;
  border-radius: 50%;
  width: 40px;
  height: 40px;
  animation: spin 2s linear infinite;
}
</style>
