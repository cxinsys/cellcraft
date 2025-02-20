<template>
    <div id="layout">
        <div class="plotly-layout">
            <div id="plotlyChart" ref="plotlyChart"></div>
        </div>
        <div class="options-layout">
            <div class="options__item">
                <p class="options__title">Method</p>
                <select class="options__item--select" v-model="selectedVisualizationTitle"
                    @change="setSelectVisualization($event)">
                    <option v-for="(item, index) in visualizationList" :key="index" :value="item.name">
                        {{ item.name }}
                    </option>
                </select>
            </div>
            <div class="options__item" v-for="(parameter, index) in selectedVisualizationParams" :key="index"
                v-show="parameter.type !== 'outputFile'">
                <p class="options__title" v-if="parameter.type === 'inputFile' && parameter.name.includes('target')">
                    Input File
                </p>
                <p class="options__title" v-else>
                    {{ parameter.name }}
                </p>
                <select class="options__item--select" v-model="parameter.defaultValue"
                    v-if="parameter.type === 'inputFile' && parameter.name.includes('target')">
                    <option v-for="(resultFile, index) in resultFileConnectionList" :key="index"
                        :value="resultFile.data.file.name">
                        {{ resultFile.data.file.name }}
                    </option>
                </select>
                <div class="checkbox-wrapper-18 options__item--input" v-else-if="parameter.type === 'inputFile'">
                    <div class="round">
                        <input type="checkbox" id="checkbox-18" v-model="checkStatuses[parameter.defaultValue]"
                            disabled />
                        <label for="checkbox-18"></label>
                    </div>
                </div>
                <input type="number" v-if="parameter.type === 'int' || parameter.type === 'float'"
                    class="options__item--input" v-model="parameter.defaultValue" :step="parameter.defaultValue"
                    :min="parameter.min" :max="parameter.max" />
                <input type="checkbox" v-else-if="parameter.type === 'boolean'" class="options__item--input"
                    v-model="parameter.defaultValue">
                <input type="text" v-else-if="parameter.type === 'string'" class="options__item--input"
                    v-model="parameter.defaultValue" />
            </div>
            <div class="options__item">
                <button id="reset-button" @click="resetSelect">Reset</button>
                <button id="apply-button" @click="runVisualization" :disabled="disableApplyButton || on_progress"
                    :class="{ 'failure': showFailure, 'success': showSuccess }">
                    <p v-if="on_progress">
                        <span class="button-loader"></span>
                        Processing...
                    </p>
                    <p v-else-if="showFailure">
                        Failure...
                    </p>
                    <p v-else-if="showSuccess">
                        Visualization
                    </p>
                    <p v-else>
                        {{ disableApplyButton ? "Select Applied" : "Select Apply " }}
                    </p>
                </button>
            </div>
        </div>
    </div>
</template>

<script>
import { getPluginInfo, getResults, runVisualization, getVisualizationResult, createTaskEventSource } from "@/api/index";
import Plotly from "plotly.js-dist-min";

export default {
    data() {
        return {
            workflowId: this.$route.query.workflow_id,
            nodeId: this.$route.query.node,
            algorithmId: null,
            plotlyData: null, // Plotly 데이터를 저장
            layout: {}, // Plotly 레이아웃 설정
            selectedPlugin: "", // 선택된 플러그인
            selectedVisualizationParams: [], // 선택된 각화 방법의 파라미터
            selectedVisualizationTitle: "", // 선택된 시각화 방법 제목
            visualizationList: [], // 시각화 방법 리스트
            resultFileConnectionList: [], // 결과 연결 리스트
            resultFileList: [], // 결과 파 리스트
            checkStatuses: {}, // 각 파일 체크 상태를 저장
            eventSources: {}, // EventSource 객체를 저장
            disableApplyButton: false, // Apply 버튼 비활성화 여부
            on_progress: false, // task 진행 중 여부
            taskStatus: '', // 추가: task 상태 저장
            showFailure: false, // 추가: 실패 메시지 표시 여부
            showSuccess: false, // 추가: 성공 메시지 표시 여부
        };
    },
    async mounted() {
        try {
            const current_node = this.$store.getters.getWorkflowNodeInfo(this.nodeId);
            const connectionList = current_node.inputs.input_1.connections

            const connectionAlgorithmNode = this.$store.getters.getAlgorithmNodeConnectedToInput(this.nodeId);
            console.log(connectionAlgorithmNode);
            this.algorithmId = connectionAlgorithmNode.id;

            this.resultFileList = await this.checkResultFiles();

            this.selectedPlugin = connectionAlgorithmNode.data.selectedPlugin.name;
            const pluginInfo = await getPluginInfo(this.selectedPlugin);

            // 1. pluginInfo에서 rules를 가져와서 isVisualization이 True인 것만 visualizationList에 저장
            const pluginRules = Object.values(pluginInfo.data.plugin.rules);

            const visualizationList = pluginRules.filter(rule => rule.isVisualization);
            this.visualizationList = visualizationList;
            console.log(this.visualizationList);

            // 2. 결과 연결 리스트를 가져와서 resultFileConnectionList에 저장
            // connectionList 순회하면서 getWorkflowNodeInfo(connectionItem.node) 실행한 결를 resultFileConnectionList에 저장
            this.resultFileConnectionList = connectionList.map(connectionItem => {
                return this.$store.getters.getWorkflowNodeInfo(connectionItem.node);
            });

            // 3. current_node의 data에서 title을 가져와서 visualizationList안에 존재한다면 selectedVisualization에 저장
            const node_title = current_node.data.title;
            const selectedVisualization = visualizationList.find(visualization => visualization.output[0] === node_title);


            // 4. current_node의 data에서 selectedVisualizationParams와 selectedVisualizationTitle을 가져와서 저장
            if (current_node.data["selectedVisualizationParams"] && current_node.data["selectedVisualizationTitle"]) {
                this.selectedVisualizationParams = current_node.data["selectedVisualizationParams"];
                this.selectedVisualizationTitle = current_node.data["selectedVisualizationTitle"];
            } else if (selectedVisualization) {
                this.selectedVisualizationParams = selectedVisualization.parameters
                this.selectedVisualizationTitle = selectedVisualization.name;
            } else {
                this.selectedVisualizationParams = this.visualizationList[0].parameters
                this.selectedVisualizationTitle = this.visualizationList[0].name;
            }

            // 5. checkStatuses 초기화
            this.initializeCheckStatuses();

            if (this.selectedVisualizationParams.length > 0) {
                const dataObject = {
                    "selectedVisualizationParams": this.selectedVisualizationParams,
                };
                const nodeId = this.nodeId;
                this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
            }

            if (this.selectedVisualizationTitle) {
                const dataObject = {
                    "selectedVisualizationTitle": this.selectedVisualizationTitle,
                };
                const nodeId = this.nodeId;
                this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
            }

            // 시각화 결과를 가져와서 plotlyData에 저장
            // const response = await getVisualizationResult();
            this.renderPlot();
        } catch (error) {
            console.error(error);
        }
    },
    watch: {
        selectedVisualizationParams: {
            handler(newVal) {
                if (newVal) {
                    const dataObject = {
                        "selectedVisualizationParams": newVal,
                    };
                    const nodeId = this.nodeId;
                    this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
                }
            },
            deep: true,
        },
        selectedVisualizationTitle: {
            handler(newVal) {
                if (newVal) {
                    const dataObject = {
                        "selectedVisualizationTitle": newVal,
                    };
                    const nodeId = this.nodeId;
                    this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
                }
            },
        },
    },
    beforeDestroy() {
        // Close all the event source connections before the component is destroyed
        for (let task_id in this.eventSources) {
            this.closeEventSource(task_id);
        }
    },
    methods: {
        async runVisualization() {
            try {
                const title = this.$store.getters.getTitle;
                const thumbnail = this.$store.getters.getThumbnail;
                const workflow_info = this.$store.getters.getWorkflowInfo;
                const workflow = {
                    id: this.workflowId,
                    current_node_id: this.nodeId,
                    title: title,
                    thumbnail: thumbnail,
                    workflow_info: workflow_info,
                };
                const workflow_data = await runVisualization(workflow);
                console.log(workflow_data);

                if (workflow_data.data.message === "Visualization result already exists") {
                    const workflowResult = {
                        id: this.workflowId,
                        algorithm_id: this.algorithmId,
                        filename: workflow_data.data.result_path,
                    }
                    const response = await getVisualizationResult(workflowResult);
                    this.plotlyData = response.data;
                    this.renderPlot();
                    this.taskStatus = 'SUCCESS';
                } else {
                    this.createEventSource(workflow_data.data.task_id);
                }
            } catch (error) {
                console.error(error);
            }
        },
        renderPlot() {
            // plotlyData가 존재할 때만 렌더링 수행
            if (this.plotlyData) {
                Plotly.newPlot(this.$refs.plotlyChart, this.plotlyData.data, this.plotlyData.layout);
            }
        },
        async resetSelect() {
            try {
                const pluginInfo = await getPluginInfo(this.selectedPlugin);
                const pluginRules = Object.values(pluginInfo.data.plugin.rules);
                const visualizationList = pluginRules.filter(rule => rule.isVisualization);
                this.visualizationList = visualizationList;
                const selectedVisualization = visualizationList.find(visualization => visualization.name === this.selectedVisualizationTitle);
                this.selectedVisualizationParams = selectedVisualization.parameters

                // plotly 랜더링 초기화
                this.plotlyData = null;
                this.renderPlot();
            } catch (error) {
                console.error(error);
            }
        },
        createEventSource(task_id) {
            this.on_progress = true;
            this.taskStatus = '';
            this.showFailure = false;

            this.eventSources[task_id] = createTaskEventSource(task_id, {
                onMessage: (event) => {
                    console.log("Received update: ", event.data);
                },
                onComplete: (status) => {
                    console.log("close", status);
                    if (status === 'FAILURE') {
                        this.taskStatus = 'FAILURE';
                        this.showFailure = true;
                        setTimeout(() => {
                            this.showFailure = false;
                        }, 3000);
                    } else if (status === 'SUCCESS') {
                        this.taskStatus = 'SUCCESS';
                        this.showSuccess = true;
                    }
                    this.on_progress = false;
                    this.closeEventSource(task_id);
                    clearInterval(this.timeInterval);
                },
                onError: (error) => {
                    console.error("SSE Error:", error);
                    this.taskStatus = 'FAILURE';
                    this.showFailure = true;
                    setTimeout(() => {
                        this.showFailure = false;
                    }, 3000);
                    this.on_progress = false;
                    this.closeEventSource(task_id);
                    clearInterval(this.timeInterval);
                }
            });
        },
        closeEventSource(task_id) {
            // If the event source exists, close it
            if (this.eventSources[task_id]) {
                this.eventSources[task_id].close();
                delete this.eventSources[task_id];
            }
        },
        setSelectVisualization(event) {
            const selectedVisualization = this.visualizationList.find(visualization => visualization.name === event.target.value);
            this.selectedVisualizationParams = selectedVisualization.parameters
            this.initializeCheckStatuses();
            this.renderPlot();
        },
        setSelectInputFile(event) {
            this.selectedInputFile = event.target.value;
            this.renderPlot();
        },
        initializeCheckStatuses() {
            if (Array.isArray(this.selectedVisualizationParams)) {
                this.selectedVisualizationParams.forEach((parameter) => {
                    if (parameter.type === 'inputFile' && !parameter.name.includes('target')) {
                        console.log(parameter.defaultValue);

                        this.checkStatuses[parameter.defaultValue] = this.resultFileList.some(resultFile => resultFile.name === parameter.defaultValue);
                        console.log(this.resultFileList);

                        this.disableApplyButton = !Object.values(this.checkStatuses).every(status => status);
                    }
                });
            }
        },
        async checkResultFiles() {
            const workflow_result = {
                id: this.workflowId,
                algorithm_id: this.algorithmId,
            };
            const response = await getResults(workflow_result);
            const resultFiles = response.data;
            console.log(resultFiles);

            return resultFiles;
        },
    }
}
</script>

<style scoped>
#layout {
    width: 100%;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: row;
    position: relative;
}

#plotlyChart {
    width: 100%;
    height: 100%;
}

.plotly-layout {
    width: 70%;
    height: 95%;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: column;
    border-radius: 1rem;
    margin: 1%;
    box-sizing: border-box;
    background-color: rgb(255, 255, 255);
}

.options-layout {
    width: 30%;
    height: 95%;
    display: flex;
    justify-content: center;
    flex-direction: column;
    gap: 0.5rem;
    z-index: 9997;
}

.options__item {
    width: 90%;
    display: flex;
    align-items: center;
    justify-content: space-between;
    border-radius: 3px;
    border-color: #f1f2fc;
}

.options__title {
    padding: 0.5rem 0;
    font-weight: 600;
    color: rgb(55, 55, 55);
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}

.options__item--select,
.options__item--input {
    width: 8rem;
    padding: 0.5rem;
    border-radius: 8px;
    border-color: #e7eaff;
    color: #545454;
}

.options__item--input {
    width: 7rem;
    border: none;
    display: flex;
    justify-content: center;
}

.options__item--select:focus,
.options__item--input:focus {
    outline: none;
}

.checkbox-wrapper-18 .round {
    position: relative;
}

.checkbox-wrapper-18 .round label {
    background-color: #fff;
    border: 1px solid #ccc;
    border-radius: 50%;
    cursor: pointer;
    height: 28px;
    width: 28px;
    display: block;
}

.checkbox-wrapper-18 .round label:after {
    border: 2px solid #fff;
    border-top: none;
    border-right: none;
    content: "";
    height: 6px;
    left: 8px;
    opacity: 0;
    position: absolute;
    top: 9px;
    transform: rotate(-45deg);
    width: 12px;
}

.checkbox-wrapper-18 .round input[type="checkbox"] {
    visibility: hidden;
    display: none;
    opacity: 0;
}

.checkbox-wrapper-18 .round input[type="checkbox"]:checked+label {
    background-color: #66bb6a;
    border-color: #66bb6a;
}

.checkbox-wrapper-18 .round input[type="checkbox"]:checked+label:after {
    opacity: 1;
}

#apply-button {
    background-color: #2d2fbf;
    /* 버튼 배경색 */
    width: 9rem;
    height: 2.5rem;
    color: white;
    /* 글자색 */
    padding: 10px 0px;
    /* 상하 10px, 좌우 20px의 여백 */
    border: none;
    /* 테두리 없앰 */
    border-radius: 4px;
    /* 테두리 모서리 둥글게 */
    cursor: pointer;
    /* 마우스 오버 시 커서 변경 */
    font-size: 16px;
    /* 글자 크기 */
    transition: background-color 0.3s;
    /* 배경색 변경시 트랜지션 효과 */
    display: flex;
    align-items: center;
    justify-content: center;
}

#apply-button:hover {
    background-color: #4655ff;
    /* 마우스 오버시 버튼의 배경색 변경 */
}

#apply-button:disabled {
    background-color: #ccc;
    /* 비활성화 상태의 배경색 */
    color: #666;
    /* 비활성화 상태의 글자색 */
    cursor: not-allowed;
    /* 비활성화 상태에서의 커서 */
}

#reset-button {
    background-color: #616161;
    /* 버튼 배경색 */
    width: 4rem;
    height: 2.5rem;
    color: white;
    /* 글자색 */
    padding: 10px 0px;
    /* 상하 10px, 좌우 20px의 여백 */
    border: none;
    /* 테두리 없앰 */
    border-radius: 4px;
    /* 테두리 모서리 둥글게 */
    cursor: pointer;
    /* 마우스 오버 시 커서 변경 */
    font-size: 16px;
    /* 글자 크기 */
    transition: background-color 0.3s;
    /* 배경색 변경시 트랜지션 효과 */
    margin-right: 0.5rem;
}

#reset-button:hover {
    background-color: #797979;
    /* 마우스 오버시 버튼의 배경색 변경 */
}

.loading-animation {
    position: fixed;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    background-color: rgba(0, 0, 0, 0.5);
    color: #fff;
    padding: 1rem;
    border-radius: 8px;
    z-index: 9999;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: column;
}

.loader {
    width: 48px;
    height: 48px;
    border: 5px solid #FFF;
    border-bottom-color: transparent;
    border-radius: 50%;
    margin-bottom: 0.5rem;
    display: inline-block;
    box-sizing: border-box;
    animation: rotation 1s linear infinite;
}

.loading-animation p {
    margin: 0;
    color: #fff;
}

.button-loader {
    width: 16px;
    height: 16px;
    border: 2px solid #FFF;
    border-bottom-color: transparent;
    border-radius: 50%;
    display: inline-block;
    margin-right: 8px;
    animation: rotation 1s linear infinite;
}

@keyframes rotation {
    0% {
        transform: rotate(0deg);
    }

    100% {
        transform: rotate(360deg);
    }
}

#apply-button.failure {
    background-color: #ff4444;
    transition: background-color 0.3s ease;
}

#apply-button.success {
    background-color: #66bb6a;
    transition: background-color 0.3s ease;
}
</style>