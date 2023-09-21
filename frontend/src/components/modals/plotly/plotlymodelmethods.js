// json -> keys, datas
function getKeys(data) {
  return Object.keys(data[0]);
}

// function getDataArrays(keys, data) {
//     var datas = {}
//     keys.forEach(function(key) {
//          datas[key] = []
//          for (let i=0; i<data.length; i++) {
//              datas[key].push(data[i][key])
//          }
//     })
//     return datas
// }

function getDataArrays(keys, data) {
  var datas = {};
  keys.forEach(function (key) {
    datas[key] = data.map((x) => {
      return x[key];
    });
  });
  return datas;
}

// get max, min of array
function getMaxOfArray(numArray) {
  return Math.max.apply(null, numArray);
}

function getMinOfArray(numArray) {
  return Math.min.apply(null, numArray);
}

export { getKeys, getDataArrays, getMaxOfArray, getMinOfArray };
