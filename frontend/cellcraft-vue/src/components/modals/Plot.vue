<template>
  <div id="layout">
    <div class="plotly-layout">
      <div id="plotly__chart"></div>
    </div>
    <div class="options-layout">
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
      keys: null,
      areNum: null,
      lines: null,
      chartType: "scattergl",
      chartMode: "markers",
      selectedX: null,
      selectedY: null,
      selectedName: null,
      selectedCluster: null,
      numClusterConstraint: 1000,
      markerSize: 1,
      numList: [{ name: " ", value: null }],
      keyList: [{ name: " ", value: null }],
      clusterList: [{ name: " ", value: null }],
    };
  },
  mounted() {
    Plotly.newPlot("plotly__chart", {
      data: [{ type: this.chartType }],
      layout: {
        autosize: true,
        automargin: true,
        width: 600,
        height: 400,
      },
    });
  },
  methods: {
    updateChart() {
      if (this.selectedCluster) {
        const clusterList = [
          ...new Set(this.lines.map((x) => x[this.selectedCluster])),
        ];

        var traces = [];
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
        Plotly.newPlot("plotly__chart", {
          data: traces,
          layout: {},
        });
      } else {
        Plotly.newPlot("plotly__chart", {
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
          layout: {},
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
      if (current_node.name === "Plot") {
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
        this.lines = dataTableResult.data.split("\n").map((x) => x.split(","));
        this.keys = this.lines.splice(0, 1)[0];
        this.keys[0] = "INDEX";
        this.areNum = this.lines[0].map((x) => !isNaN(x));

        this.numList = [{ name: " ", value: null }];
        for (let i = 0; i < this.keys.length; i++) {
          if (this.areNum[i] == true) {
            this.numList.push({ name: this.keys[i], value: i });
          }
        }

        this.keyList = [{ name: " ", value: null }];
        for (let i = 0; i < this.keys.length; i++) {
          this.keyList.push({ name: this.keys[i], value: i });
        }

        this.clusterList = [{ name: " ", value: null }];
        for (let i = 0; i < this.keys.length; i++) {
          const countUnique = new Set(this.lines.map((x) => x[i])).size;
          if (countUnique < this.numClusterConstraint) {
            this.clusterList.push({ name: this.keys[i], value: i });
          }
        }

        this.updateChart();
      }
    },
  },
};
</script>

<style>
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
  background-color: blue;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
}
.options-layout {
  width: 25%;
  height: 95%;
  background-color: red;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
}
.options__item {
  /* margin: auto; */
}
</style>
