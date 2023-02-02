<template>
  <div id="layout">
    <div class="plotly-layout">
      <Plotly :data="data" :layout="layout" :display-mode-bar="true"></Plotly>
    </div>
  </div>
</template>

<script>
import { Plotly } from "vue-plotly";
import { changeData, changeLayout } from "./plotly/plotlypresenter.js";
import { getResult } from "@/api/index";

export default {
  components: {
    Plotly,
  },
  data() {
    return {
      node_name: "Plot",
      type: "scatter",
      id: "",
      x: "",
      y: "",
      cluster: "",
      data: changeData("scatter", "", "", "", ""),
      layout: changeLayout("sample title"),
      // chartList: [
      //   { name: "None", properties: "" },
      //   { name: "Scatter", properties: "scatter" },
      //   { name: "Heatmap", properties: "heatmap" },
      // ],
      keys: null,
      plotData: null,
      current_file: null,
      file_num: null,
    };
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
        console.log(dataTableResult.data);

        //백엔드에서 넘겨준 plot 데이터
        this.plotData = dataTableResult.data;
      }
    },
    // checkCurrentNode(val) {
    //   const current_node = this.$store.getters.getNodeInfo(val);
    //   this.current_file = this.$store.getters.getCurrentFile.file;
    //   console.log(current_node, this.current_file);
    //   if (this.current_file === "sampledata.csv") {
    //     this.keys = ["NONE"].concat(keys);
    //     this.tag = this.keys[1];
    //     this.x = this.keys[2];
    //     this.y = this.keys[3];
    //     this.cluster = this.keys[4];
    //     this.file_num = 1;
    //     this.data = changeData(
    //       "scatter",
    //       this.tag,
    //       this.x,
    //       this.y,
    //       this.cluster,
    //       this.file_num
    //     );
    //   } else if (this.current_file === "sampledata2.csv") {
    //     this.keys = ["NONE"].concat(keys2);
    //     this.tag = this.keys[1];
    //     this.x = this.keys[2];
    //     this.y = this.keys[3];
    //     this.cluster = this.keys[4];
    //     this.file_num = 2;
    //     this.data = changeData(
    //       "scatter",
    //       this.tag,
    //       this.x,
    //       this.y,
    //       this.cluster,
    //       this.file_num
    //     );
    //   }
    // },
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
  flex-direction: column;
}
.plotly-layout {
  width: 95%;
  height: auto;
}
</style>
