<template>
  <div id="layout">
    <div class="plotly-layout">
      <div id="plotly__chart"></div>
    </div>
    <div class="options-layout"></div>
  </div>
</template>

<script>
// import { Plotly } from "vue-plotly";
// import { changeData, changeLayout } from "./plotly/plotlypresenter.js";
import { getResult } from "@/api/index";
import Plotly from "plotly.js-dist-min";

export default {
  // components: {
  //   Plotly,
  // },
  data() {
    return {
      keys: null,
      lines: null,
      plotData: null,
      current_file: null,
      file_num: null,
    };
  },
  mounted() {
    Plotly.newPlot(
      "plotly__chart",
      /* JSON object */ {
        data: [{ type: "scattergl" }],
        layout: {
          autosize: true,
          automargin: true,
          width: 600,
          height: 400,
        },
        // config: {
        //   responsive: true,
        // },
      }
    );
  },
  methods: {},
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

        //얘를 여기 넣으면 안 될것 같은데 ㅜ
        Plotly.newPlot(
          "plotly__chart",
          /* JSON object */ {
            data: [
              {
                x: this.lines.map((x) => x[2]),
                y: this.lines.map((x) => x[3]),
                type: "scattergl",
                mode: "markers",
                marker: { size: 1 },
              },
            ],
            // layout: {
            //   autosize: true,
            //   automargin: true,
            // },
            // config: {
            //   responsive: true,
            // },
          }
        );
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
</style>
