export const connectionRules = {
  InputFile: ["DataTable", "ScatterPlot", "Algorithm"],
  DataTable: ["ScatterPlot", "Algorithm"],
  ScatterPlot: ["DataTable", "Algorithm"],
  Algorithm: ["ResultFiles", "Visualization"],
  ResultFiles: ["Visualization"],
  Visualization: [],
};
