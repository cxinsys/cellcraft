// import {getKeys, getDataArrays, getMaxOfArray, getMinOfArray} from "./plotlymethods"
import { getKeys, getDataArrays } from "./plotlymodelmethods";

// (tmp) input data.json --- will get from presenter
import { sampleData } from "./sampledata";
import { sampleData2 } from "./sampledata2";

let keys = getKeys(sampleData);
let datas = getDataArrays(keys, sampleData);

let keys2 = getKeys(sampleData2);
let datas2 = getDataArrays(keys2, sampleData2);

export { keys, datas, keys2, datas2 };
