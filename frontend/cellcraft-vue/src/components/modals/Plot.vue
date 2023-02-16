<template>
  <div class="layout">
    <div class="options">
      <div class="options__column">
        <div class="options__name">Values:</div>
        <input
          type="text"
          class="options__button"
          name="value"
          v-model="value"
          @click="optionClick"
        />
        <ul class="options__menu">
          <li
            class="options__content"
            v-for="(header, idx) in headers"
            :key="idx"
            @click="menuClick"
          >
            {{ header }}
          </li>
        </ul>
      </div>
      <div class="options__column">
        <div class="options__name">Annotations:</div>
        <input
          type="text"
          class="options__button"
          name="annotation"
          v-model="annotation"
          @click="optionClick"
        />
        <ul class="options__menu">
          <li
            class="options__content"
            v-for="(header, idx) in headers"
            :key="idx"
            @click="menuClick"
          >
            {{ header }}
          </li>
        </ul>
      </div>
    </div>
    <div class="canvas">
      <svg :width="width" :height="height">
        <g class="graph" :width="width" :height="height">
          <g class="xAxisG"></g>
          <g class="yAxisG"></g>
        </g>
        <text class="xLabel"></text>
        <text class="yLabel"></text>
      </svg>
      <div class="tooltip"></div>
    </div>
  </div>
</template>

<script>
/* eslint-disable */
import * as d3 from "d3";
import { getResult } from "@/api/index";
export default {
  props: {
    file_name: null,
  },
  data() {
    return {
      node_name: "Plot",
      headers: [],
      num_headers: [],
      str_headers: [],
      bool_headers: [],
      value: "Values",
      annotation: "Annotations",
      csv: [],
      data: [],
      filtered_data: [],
      width: 480,
      height: 520,
      graphWidth: 0,
      graphHeight: 0,
      margin: {
        mt: 60,
        mb: 60,
        mr: 20,
        ml: 60,
      },
    };
  },
  methods: {
    optionClick(event) {
      event.target.style.opacity = "0.6";
      event.target.nextElementSibling.style.display = "block";
    },
    menuClick(event) {
      // console.dir(event.target.parentElement.previousElementSibling)
      const menuContent = event.target.innerText;
      event.target.parentElement.style.display = "none";
      event.target.parentElement.previousElementSibling.style.opacity = "1";
      if (
        event.target.parentElement.previousElementSibling.previousElementSibling.innerText.includes(
          "Annotation"
        )
      ) {
        this.annotation = menuContent;
      } else if (
        event.target.parentElement.previousElementSibling.previousElementSibling.innerText.includes(
          "Value"
        )
      ) {
        this.value = menuContent;
      }
      this.bar();
    },
    dotMouseover(ev, d, i) {
      console.log(d);
      d3.select(".tooltip")
        .html(
          `${this.value} : ${d.value} <br> ${this.annotation} : ${d.annotation}`
        )
        .style("display", "block")
        .style("left", ev.layerX + 10 + "px")
        .style("top", ev.layerY + 10 + "px");
    },
    dotMouseleave(ev, d, i) {
      console.log(ev);
      d3.select(".tooltip").style("display", "none");
    },
    dotMouseclick(ev, d, i) {
      console.log(ev.target.style);
      // ev.target.style.boxShadow = "rgba(0, 0, 0, 0.25) 0px 54px 55px, rgba(0, 0, 0, 0.12) 0px -12px 30px, rgba(0, 0, 0, 0.12) 0px 4px 6px, rgba(0, 0, 0, 0.17) 0px 12px 13px, rgba(0, 0, 0, 0.09) 0px -3px 5px";
    },
    bar() {
      d3.selectAll(".graph > rect").remove();
      // console.log(this.data[this.value][0]);
      this.filtered_data = [];

      console.log(this.value, this.annotation);
      for (let i = 0; i < Object.keys(this.data[this.value]).length; i++) {
        this.filtered_data.push({
          value: this.data[this.value][i],
          annotation: this.data[this.annotation][i],
        });
      }
      console.log(this.filtered_data);

      const x = d3
        .scaleBand()
        .domain(this.filtered_data.map((d) => d.annotation))
        .range([0, this.graphWidth])
        .padding(0.25);
      const y = d3
        .scaleLinear()
        .domain([0, d3.max(this.filtered_data, (d) => d.value)])
        .range([this.graphHeight, 0]);
      const xAxis = d3.axisBottom(x).ticks(5);
      const yAxis = d3.axisLeft(y).ticks(5);
      d3.select(".xAxisG").call(xAxis);
      d3.select(".yAxisG").call(yAxis);
      d3.select(".xLabel").text(this.annotation);
      d3.select(".yLabel").text(this.value);
      d3.select(".xAxisG")
        .selectAll("text")
        .attr("fill", "blue")
        .attr("transform", "rotate(-45)")
        .attr("text-anchor", "end");
      return d3
        .select(".graph")
        .selectAll("rect")
        .data(this.filtered_data)
        .enter()
        .append("rect")
        .attr("height", (d) => this.graphHeight - y(d.value))
        .attr("width", x.bandwidth)
        .attr("x", (d) => x(d.annotation))
        .attr("y", (d) => y(d.value))
        .on("mouseover", this.dotMouseover)
        .on("mouseleave", this.dotMouseleave)
        .on("click", this.dotMouseclick);
    },
  },
  async mounted() {
    // const filename = { filename: `${this.node_name}_${this.file_name}` };
    // // console.log(filename)
    // const PlotResult = await getResult(filename);
    // // console.log(PlotResult.data.club[0])
    // this.data = PlotResult.data;
    // console.log(this.data);
    // this.headers = Object.keys(PlotResult.data).slice(1);
    // this.headers.forEach((element) => {
    //   console.log(typeof this.data[element][0]);
    //   if (typeof this.data[element][0] === "number") {
    //     this.num_headers.push(element);
    //   } else if (typeof this.data[element][0] === "string") {
    //     this.str_headers.push(element);
    //   } else if (typeof this.data[element][0] === "boolean") {
    //     this.bool_headers.push(element);
    //   }
    // });
    // console.log(this.num_headers);
    // this.value = this.num_headers[0];
    // this.annotation = this.num_headers[1];
    // // console.log(typeof(this.value));

    // this.graphWidth = this.width - this.margin.ml - this.margin.mr;
    // this.graphHeight = this.height - this.margin.mt - this.margin.mb;
    // console.log(this.graphWidth, this.graphHeight);
    // d3.select(".graph").attr(
    //   "transform",
    //   `translate(${this.margin.ml}, ${this.margin.mt})`
    // );
    // d3.select(".xAxisG").attr("transform", `translate(0, ${this.graphHeight})`);
    // d3.select(".xLabel")
    //   .attr("text-anchor", "end")
    //   .attr("x", this.width / 2 + this.margin.ml / 2)
    //   .attr("y", this.height + this.margin.mt - 80);
    // d3.select(".yLabel")
    //   .attr("text-anchor", "end")
    //   .attr("transform", "rotate(-90)")
    //   .attr("y", this.margin.ml / 2 - 5)
    //   .attr("x", -this.height / 2 + 20);

    // for (let i = 0; i < Object.keys(this.data[this.value]).length; i++) {
    //   this.filtered_data.push({
    //     value: this.data[this.value][i],
    //     annotation: this.data[this.annotation][i],
    //   });
    // }

    // console.log(this.filtered_data);

    // const x = d3
    //   .scaleBand()
    //   .domain(this.filtered_data.map((d) => d.annotation))
    //   .range([0, this.graphWidth])
    //   .padding(0.25);
    // const y = d3
    //   .scaleLinear()
    //   .domain([0, d3.max(this.filtered_data, (d) => d.value)])
    //   .range([this.graphHeight, 0]);
    // const xAxis = d3.axisBottom(x).ticks(5);
    // const yAxis = d3.axisLeft(y).ticks(5);
    // d3.select(".xAxisG").call(xAxis);
    // d3.select(".yAxisG").call(yAxis);
    // d3.select(".xLabel").text(this.annotation);
    // d3.select(".yLabel").text(this.value);

    // d3.select(".xAxisG")
    //   .selectAll("text")
    //   .attr("fill", "blue")
    //   .attr("transform", "rotate(-45)")
    //   .attr("text-anchor", "end");
    // d3.select(".graph")
    //   .selectAll("rect")
    //   .data(this.filtered_data)
    //   .enter()
    //   .append("text")
    //   .attr("x", (d) => x(d.annotation))
    //   .attr("y", (d) => y(d.value) - 5)
    //   .text((d) => d.values)
    //   .attr("text-anchor", "start")
    //   .style("font-size", "12px");

    // d3.select(".graph")
    //   .selectAll("rect")
    //   .data(this.filtered_data)
    //   .enter()
    //   .append("rect")
    //   .attr("height", (d) => this.graphHeight - y(d.value))
    //   .attr("width", x.bandwidth)
    //   .attr("x", (d) => x(d.annotation))
    //   .attr("y", (d) => y(d.value))
    //   .on("mouseover", this.dotMouseover)
    //   .on("mouseleave", this.dotMouseleave)
    //   .on("click", this.dotMouseclick);
  },
};
</script>

<style scoped>
.options,
.options__column,
.options__name,
.options__button,
.options__menu,
.options__content {
  box-sizing: border-box;
  border-radius: 5px;
}
.layout {
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  justify-content: center;
  position: relative;
}
.options {
  width: 40%;
  height: 100%;
  padding: 10px;
  display: flex;
  flex-direction: column;
  align-items: center;
  background: rgb(235, 235, 235);
  position: absolute;
  left: 0;
}
.options__column {
  width: 100%;
  height: 5%;
  margin-bottom: 15px;
  display: flex;
  align-items: center;
  position: relative;
  background: rgb(210, 210, 210);
}
.options__name {
  width: 35%;
  height: 100%;
  padding: 5px;
  display: flex;
  align-items: center;
}
.options__button,
.options__menu {
  width: 65%;
  height: 100%;
  padding: 5px;
  display: flex;
  align-items: center;
  justify-content: flex-start;
  background: white;
}
.options__menu {
  display: none;
  height: auto;
  flex-direction: column;
  position: absolute;
  box-sizing: border-box;
  border-radius: 7px;
  top: 100%;
  right: 0;
  z-index: 1;
}
.options__content {
  margin-bottom: 4px;
}
.options__content:hover {
  background: rgb(125, 125, 255);
}
.layout {
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  justify-content: center;
  position: relative;
}
.options {
  width: 40%;
  height: 100%;
  padding: 10px;
  display: flex;
  flex-direction: column;
  align-items: center;
  background: rgb(235, 235, 235);
  position: absolute;
  left: 0;
}
.options__column {
  width: 100%;
  height: 5%;
  margin-bottom: 15px;
  display: flex;
  align-items: center;
  position: relative;
  background: rgb(210, 210, 210);
}
.options__name {
  width: 35%;
  height: 100%;
  padding: 5px;
  display: flex;
  align-items: center;
}
.options__button,
.options__menu {
  width: 65%;
  height: 100%;
  padding: 5px;
  display: flex;
  align-items: center;
  justify-content: flex-start;
  background: white;
}
.options__menu {
  display: none;
  height: auto;
  flex-direction: column;
  position: absolute;
  box-sizing: border-box;
  border-radius: 7px;
  top: 100%;
  right: 0;
  z-index: 1;
}
.options__content {
  margin-bottom: 4px;
}
.options__content:hover {
  background: rgb(125, 125, 255);
}
.canvas {
  width: 60%;
  height: 100%;
  position: absolute;
  right: 0;
  display: flex;
  align-items: center;
  justify-content: center;
}
.tooltip {
  position: absolute;
  display: none;
  width: 100px;
  height: 150px;
  padding: 10px;
  background: white;
  box-shadow: rgba(0, 0, 0, 0.24) 0px 3px 8px;
  border-radius: 20px;
}
</style>
