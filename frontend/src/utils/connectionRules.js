export const connectionRules = {
  InputFile: ["DataTable", "ScatterPlot", "Algorithm"],
  DataTable: ["ScatterPlot", "Algorithm"],
  ScatterPlot: ["DataTable", "Algorithm"],
  Algorithm: ["ResultFile", "ResultFiles", "Visualization"],
  ResultFile: ["Visualization"],
  ResultFiles: ["Visualization"],
  Visualization: [],
};
