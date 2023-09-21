import { keys, datas, keys2, datas2 } from "./plotlymodel";
import { newScatter } from "./plotlypresentermethods";

// change type
function changeData(plotlyType, id, x, y, cluster = "", file) {
  if (plotlyType == "scatter") {
    // newScatter(clusters, ids, x, y, markerSize)
    // return newScatter(datas[keys[3]], datas[keys[0]], datas[keys[1]], datas[keys[2]], 2)
    if (file === 1) {
      if (cluster == "") {
        return newScatter(datas[id], datas[x], datas[y], 2);
      } else {
        return newScatter(datas[id], datas[x], datas[y], 2, datas[cluster]);
      }
    } else if (file === 2) {
      if (cluster == "") {
        return newScatter(datas2[id], datas2[x], datas2[y], 2);
      } else {
        return newScatter(datas2[id], datas2[x], datas2[y], 2, datas2[cluster]);
      }
    }
  } else {
    return [];
  }
}

function changeLayout(title) {
  return { title: title };
}

export { keys, keys2, changeData, changeLayout };
