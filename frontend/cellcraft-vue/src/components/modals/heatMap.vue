<template>
  <div id="layout">
    <div class="plotly-layout">
      <div id="plotly__heatmap"></div>
    </div>
    <div class="options-layout">
      <input
        type="text"
        placeholder="Title"
        class="options__item"
        @input="searchChangeFunc($event)"
      />
      <div class="options__item">
        X - axis
        <select :value="selectedX" @change="setSelectX($event)">
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
        Y - axis
        <select :value="selectedY" @change="setSelectY($event)">
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
        Name
        <select :value="selectedName" @change="setSelectName($event)">
          <option
            v-for="(item, index) in keyList"
            :key="index"
            :value="item.value"
          >
            {{ item.name }}
          </option>
        </select>
      </div>
      <div class="options__item">
        Cluster
        <select :value="selectedCluster" @change="setSelectCluster($event)">
          <option
            v-for="(item, index) in clusterList"
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
import { getResult } from "@/api/index";
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
      selectedX: null, // 선택된 X축의 column index
      selectedY: null, // 선택된 Y축의 column index
      selectedName: null, // 선택된 Name의 column index
      selectedCluster: null, // 선택된 Cluster의 column index
      numClusterConstraint: 1000, // Cluster가 너무 다양하면 성능에 저하가 생기므로 제한을 둠
      markerSize: 1, // 점의 사이즈
      numList: [
        { name: "None", value: null },
        { name: "None", value: 0 }, // ! -- tab을 옮겨갔다와야 업데이트되는 문제를 해결하면 삭제할 라인
        { name: "None", value: 2 }, // ! -- tab을 옮겨갔다와야 업데이트되는 문제를 해결하면 삭제할 라인
        { name: "None", value: 3 }, // ! -- tab을 옮겨갔다와야 업데이트되는 문제를 해결하면 삭제할 라인
      ], // number type으로 x,y축에 들어가기 적합한 자료들의 option list
      keyList: [{ name: "None", value: null }], // key list
      clusterList: [{ name: "None", value: null }], // cluster가 되기 적합한 list, unique한 자료수가 적은 column의 list
    };
  },
  mounted() {
    Plotly.newPlot("plotly__heatmap", {
      data: [{ type: this.chartType }],
      layout: {
        autosize: true,
        automargin: true,
        width: 600, // !--조정필요
        height: 500, // !--조정필요
        // paper_bgcolor: rgb(255, 1, 255),
        // plot_bgcolor: rgb(1, 255, 255),
      },
    });
  },
  methods: {
    // 차트 업데이트
    updateChart() {
      if (this.selectedCluster) {
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
        Plotly.newPlot("plotly__heatmap", {
          data: traces,
          layout: {
            title: this.chartTitle,
          },
        });
      } else {
        Plotly.newPlot("plotly__heatmap", {
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
          },
        });
      }
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
      this.selectedCluster = event.target.value;
      this.updateChart();
    },
    searchChangeFunc(event) {
      this.chartTitle = event.target.value;
      this.updateChart();
    },
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
      // console.log(current_node);
      // console.log(this.current_file.file);
      if (current_node.name === "heatMap") {
        const filename = {
          filename: `${current_node.name}_${this.current_file.replace(
            ".csv",
            ""
          )}`,
        };
        console.log(filename);
        const dataTableResult = await getResult(filename);

        //백엔드에서 넘겨준 plot 데이터
        // dataTableResult.data;
        // lines, keys, areNum 업데이트
        this.lines = dataTableResult.data.split("\n").map((x) => x.split(","));
        this.keys = this.lines.splice(0, 1)[0];
        this.keys[0] = "INDEX"; // keys의 [0]을 ""로 받아오기 때문에 "INDEX로 변환"
        this.areNum = this.lines[0].map((x) => !isNaN(x));

        console.log(this.lines, this.keys);
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

        this.clusterList = [{ name: "None", value: null }];
        for (let i = 0; i < this.keys.length; i++) {
          const countUnique = new Set(this.lines.map((x) => x[i])).size;
          if (countUnique < this.numClusterConstraint) {
            this.clusterList.push({ name: this.keys[i], value: i });
          }
        }

        // 초기 x,y축 세팅하기
        if (this.numList.length == 3) {
          this.selectedX = this.numList[2].value;
          this.selectedY = this.numList[2].value;
        } else if (this.numList.length > 3) {
          this.selectedX = this.numList[2].value;
          this.selectedY = this.numList[3].value;
        }

        this.updateChart();
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
.plotly-layout {
  width: 75%;
  height: 95%;
  /* background-color: blue; */
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
  padding: 1rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
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
  .plotly-layout {
    background-color: rgb(41, 43, 48);
  }
  .options__item {
    /* color: rgb(255, 255, 255); */
  }
}
</style>
