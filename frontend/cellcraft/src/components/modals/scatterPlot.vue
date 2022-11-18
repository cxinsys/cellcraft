<template>
  <div class="layout">
    <div class="canvas">
      <svg :width="width" :height="height">
        <g class="graph" :width="width" :height="height" :d="scatter">
          <g class="xAxisG"></g>
          <g class="yAxisG"></g>
        </g>
        <text class="xLabel"></text>
        <text class="yLabel"></text>
      </svg>
    </div>
  </div>
</template>

<script>
/* eslint-disable */
import * as d3 from "d3";
import { getResult } from '@/api/index';
export default {
  props: {
    file_name: null
  },
  data() {
    return {
      csv: [],
      data: [],
      xLabel: "",
      yLabel: "",
      width: 400,
      height: 400,
      graphWidth: 0,
      graphHeight: 0,
      margin: {
        mt: 30,
        mb: 50,
        mr: 80,
        ml: 80
      }
    };
  },
  methods: {},
  async mounted() {
    const filename = { filename: `${this.node_name}_${this.file_name}` }
    console.log(filename)
    const PlotResult = await getResult(filename)
    console.log(PlotResult)
    this.data = PlotResult.data

    this.graphWidth = this.width - this.margin.ml - this.margin.mr;
    this.graphHeight = this.height - this.margin.mt - this.margin.mb;
    d3.select(".graph").attr(
      "transform",
      `translate(${this.margin.ml}, ${this.margin.mt})`
    );
    d3.select(".xAxisG").attr("transform", `translate(0, ${this.graphHeight})`);
    d3.select(".xLabel")
      .attr("text-anchor", "end")
      .attr("x", this.width / 2 + this.margin.ml / 2)
      .attr("y", this.height + this.margin.mt - 45);
    d3.select(".yLabel")
      .attr("text-anchor", "end")
      .attr("transform", "rotate(-90)")
      .attr("y", this.margin.ml / 2 - 5)
      .attr("x", -this.height / 2 + 20);
  },
  computed: {
    scatter() {
      console.log(this.csv[0]);
      d3.selectAll(".graph > circle").remove();
      this.data = this.csv.slice(1);
      // console.log(this.data);
      if (this.csv[0]) {
        this.xLabel = Object.values(this.csv[0])[0];
        this.yLabel = Object.values(this.csv[0])[1];
      }
      // console.log(d3.extent(this.data, (d) => d.AxisX));
      // console.log(d3.extent(this.data, (d) => d.AxisY));
      console.log(d3.min(this.data, d => d.AxisX));
      console.log(d3.max(this.data, d => d.AxisX));
      const x = d3
        .scaleLinear()
        .domain([
          parseInt(d3.min(this.data, d => d.AxisX)) - 10,
          parseInt(d3.max(this.data, d => d.AxisX)) + 10
        ])
        .range([0, this.graphWidth]);
      const y = d3
        .scaleLinear()
        .domain([
          parseInt(d3.min(this.data, d => d.AxisY)) - 10,
          parseInt(d3.max(this.data, d => d.AxisY)) + 10
        ])
        .range([this.graphHeight, 0]);
      const xAxis = d3.axisBottom(x).ticks(10);
      const yAxis = d3.axisLeft(y).ticks(10);
      d3.select(".xAxisG").call(xAxis);
      d3.select(".yAxisG").call(yAxis);
      d3.select(".xLabel").text(this.xLabel);
      d3.select(".yLabel").text(this.yLabel);
      return d3
        .select(".graph")
        .selectAll("circle")
        .data(this.data)
        .enter()
        .append("circle")
        .attr("cx", d => x(d.AxisX))
        .attr("cy", d => y(d.AxisY))
        .attr("r", 10);
    }
  }
};
</script>

<style scoped>
.layout {
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  justify-content: center;
}
.canvas {
}
.uploadCsv {
  width: 100%;
  height: 10%;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  margin: 10px 0;
}
</style>
