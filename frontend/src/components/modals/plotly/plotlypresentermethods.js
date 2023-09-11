import { Trace } from "./trace";

// !!! need performance improvement
function newScatter(ids, x, y, markerSize, clusters = "") {
  // clusters, ids, x, y must have same length
  if (clusters == "") {
    return [new Trace("markers", "scattergl", markerSize, ids, x, y, "")];
  } else {
    let uniqueClusters = Array.from(new Set(clusters));
    var tmpTraces = [];
    for (let i = 0; i < uniqueClusters.length; i++) {
      tmpTraces.push(
        new Trace(
          "markers",
          "scattergl",
          markerSize,
          [],
          [],
          [],
          uniqueClusters[i]
        )
      );
    }
    for (let i = 0; i < clusters.length; i++) {
      tmpTraces[uniqueClusters.indexOf(clusters[i])].text.push(ids[i]);
      tmpTraces[uniqueClusters.indexOf(clusters[i])].x.push(x[i]);
      tmpTraces[uniqueClusters.indexOf(clusters[i])].y.push(y[i]);
    }
    return tmpTraces;
  }
}

function newHeatmap() {
  return [
    {
      z: [
        [1, 10, 30, 50, 1],
        [20, 30, 60, 80, 30],
        [30, 60, 1, -10, 20],
      ],
      x: ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday"],
      y: ["Morning", "Afternoon", "Evening"],
      type: "heatmap",
      hoverongaps: false,
    },
  ];
}

export { newScatter, newHeatmap };
