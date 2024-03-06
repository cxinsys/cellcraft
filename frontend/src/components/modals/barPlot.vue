<template>
  <div id="layout">
    <div class="plotly-layout">
      <div id="plotly-barplot"></div>
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
      resultSifs: [],
      selectedSif: "",
      current_file: "",
    };
  },
  name: "PlotlyBarPlot",

  async mounted() {
    // option_file_name.algorithmOptions.optionFilePath이 null이면 error 띄우기
    const option_file_name = this.$store.getters.getCurrentLinkedNodes;

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
      this.selectedSif = this.resultSifs[0];
    } catch (error) {
      console.log(error);
    }
    await this.updatePlot(this.selectedDegree);
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
      const algorithmResult = await getResultFileOne(this.selectedSif);
      // 줄바꿈 문자("\n")를 기준으로 문자열을 분할하여 각 라인을 배열로 변환
      var lines = algorithmResult.data.split("\n");

      // 각 라인을 탭 문자("\t")를 기준으로 분할하여 2차원 배열 생성
      var twoDimensionalList = lines.map(function (line) {
        return line.split("\t");
      });

      // 결과 출력 (옵션)

      var frequencyMap = {};

      twoDimensionalList.forEach(function (row) {
        var element = row[rowNum]; // 첫 번째 열의 요소

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

      // 결과 출력

      // Data for the bar plot
      const data = [
        {
          x: items.map((item) => item[1]),
          y: items.map((item) => item[0]),
          type: "bar",
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
      Plotly.newPlot("plotly-barplot", data, layout, config);
    },
    async setSelectedDegree(event) {
      this.selectedDegree = event.target.value;
      await this.updatePlot(this.selectedDegree);
    },
    async setSelectedSif(event) {
      this.selectedSif = event.target.value;
      await this.updatePlot(this.selectedDegree);
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
  right: 80px;
  width: 15px;
  height: 15px;
}
.options__item__button__plus {
  position: absolute;
  right: 10px;
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
