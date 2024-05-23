<template>
    <div class="rule-container">
        <button class="change-button" @click="toggleRuleView">
            <img src="@/assets/change_circle.png">
        </button>

        <button class="create-button" @click="showCreateModal">
            Create Node <i class="fas fa-plus"></i>
        </button>

        <!-- Drawflow 컴포넌트 -->
        <div v-show="!isRuleView" id="rule-drawflow" @drop="drop" @dragover="allowDrop"></div>

        <!-- Rule 보기 컴포넌트 -->
        <div v-if="isRuleView" class="rule-view-container">
            <div v-for="(rule, index) in rules" :key="index" class="rule-item">
                <h4>rule {{ rule.name }}:</h4>
                <div>
                    <strong>input:</strong>
                    <div v-for="(value, key) in rule.input" :key="'input-' + key">{{ key }}="{{ value }}"</div>
                </div>
                <div>
                    <strong>output:</strong>
                    <div v-for="(value, key) in rule.output" :key="'output-' + key">{{ key }}="{{ value }}"</div>
                </div>
                <div v-if="rule.script">
                    <strong>script:</strong>
                    <div>{{ rule.script }}</div>
                </div>
                <button @click="editRule(index)">Edit</button>
            </div>
        </div>

        <div v-if="isShowCreateModal" class="createModal-overlay">
            <div class="createModal-container">
                <form>
                    <div class="controller-group">
                        <label for="ruleTitle">Rule Title:</label>
                        <input type="text" id="ruleTitle" v-model="ruleTitle" />

                        <label for="scriptFile">Upload Script File:</label>
                        <input type="file" id="scriptFile" @change="handleFileUpload" />
                    </div>

                    <div class="script-code-container">
                        <pre>{{ completeRule }}</pre>
                    </div>

                    <div class="add-parameter-group">
                        <label>Parameters:</label>
                        <div class="parameter-inputs">
                            <input type="text" v-model="newParameter.name" placeholder="Parameter Name"
                                class="parameter-name" />
                            <select v-model="newParameter.type" class="parameter-type" @change="resetParameterOption">
                                <option value="inputFile">Input File</option>
                                <option value="outputFile">Output File</option>
                                <option value="string">String</option>
                                <option value="int">Integer</option>
                                <option value="float">Float</option>
                                <option value="boolean">Boolean</option>
                            </select>
                        </div>
                        <div v-if="newParameter.type === 'inputFile' || newParameter.type === 'outputFile'"
                            class="file-inputs">
                            <input type="text" v-model="newParameter.fileExtension" placeholder="File Type"
                                class="file-extension-input" />
                            <select v-model="selectedFileExtension" @change="updateFileExtension"
                                class="file-extension-select">
                                <option value="">Direct Input</option>
                                <option value=".txt">.txt</option>
                                <option value=".csv">.csv</option>
                                <option value=".json">.json</option>
                            </select>
                        </div>
                        <div v-else class="col-input-group">
                            <input type="text"
                                v-if="newParameter.type !== 'inputFile' && newParameter.type !== 'outputFile' && newParameter.type !== 'boolean'"
                                v-model="newParameter.defaultValue" placeholder="Default Value"
                                class="default-value-input" />
                            <select v-if="newParameter.type === 'boolean'" v-model="newParameter.defaultValue"
                                class="boolean-select">
                                <option value="">Select Value</option>
                                <option value="true">True</option>
                                <option value="false">False</option>
                            </select>
                            <div v-if="newParameter.type === 'int' || newParameter.type === 'float'"
                                class="range-inputs">
                                <input type="number" v-model.number="newParameter.min" placeholder="Min"
                                    class="range-min" />
                                <input type="number" v-model.number="newParameter.max" placeholder="Max"
                                    class="range-max" />
                            </div>
                        </div>
                        <button type="button" @click="addParameter">Add Parameter</button>
                    </div>


                    <div class="button-group">
                        <button type="button" @click="createNode">Create</button>
                        <button type="button" @click="isShowCreateModal = false">Cancel</button>
                    </div>
                </form>
            </div>
        </div>
    </div>
</template>

<script>
import Vue from "vue";
import ruleNode from "@/components/nodes/ruleNode.vue";

export default {
    data() {
        return {
            isRuleView: false,
            isShowCreateModal: false,
            ruleTitle: "",
            nodeInputs: 0,
            nodeOutputs: 0,
            inputFiles: [],
            outputFiles: [],
            scriptFile: null,
            newParameter: {
                name: '',
                type: 'string',
                defaultValue: '',
                min: null,
                max: null,
                fileExtension: ''
            },
            parameters: [],
            selectedFileExtension: '',
            completeRule: '',
            shellCommand: '',
            rules: [],
        };
    },
    mounted() {
        const id = document.getElementById("rule-drawflow");
        Vue.prototype.$df = new Drawflow(id, Vue, this);
        //this.$df == editor
        this.$df.start();
        this.$df.registerNode('ruleNode', ruleNode, {}, {});
        this.$df.on('nodeCreated', (id) => this.updateSnakefile(id));
        this.$df.on('nodeRemoved', (id) => this.updateSnakefile(id));
        this.$df.on('connectionCreated', (connection) => this.onConnectionCreated(connection));
        this.$df.on('connectionRemoved', (connection) => this.onConnectionRemoved(connection));
        this.$df.on('nodeDataChanged', (id, data) => this.onNodeDataChanged(id, data));
    },
    watch: {
        ruleTitle: 'updateCompleteRule',
        parameters: {
            handler: 'updateShellCommand',
            deep: true
        }
    },
    methods: {
        showCreateModal() {
            this.isShowCreateModal = true;
        },
        toggleRuleView() {
            this.isRuleView = !this.isRuleView;
        },
        allowDrop(event) {
            event.preventDefault();
        },
        drop(event) {
            event.preventDefault();
        },
        handleFileUpload(event) {
            this.scriptFile = event.target.files[0];
            if (this.scriptFile) {
                const reader = new FileReader();
                reader.onload = (e) => {
                    this.scriptContent = e.target.result;
                    this.updateShellCommand();
                };
                reader.readAsText(this.scriptFile);
            }
        },
        updateShellCommand() {
            if (!this.scriptFile) return;
            const fileName = this.scriptFile.name;
            let command = '';
            if (fileName.endsWith('.py')) {
                command = `/python ${fileName}`;
            } else if (fileName.endsWith('.R')) {
                command = `/Rscript ${fileName}`;
            } else {
                command = `/${fileName}`;
            }

            const paramStr = this.parameters.map(p => {
                if (p.type === 'inputFile' || p.type === 'outputFile') {
                    return `${p.defaultValue}(${p.type})`;
                } else {
                    return `${p.name}(${p.type}:${p.defaultValue})`;
                }
            }).join(' ');

            this.shellCommand = `${command} ${paramStr}`;
            this.updateCompleteRule();
        },
        updateCompleteRule() {
            const inputs = this.inputFiles.slice(1).join(', ');
            const outputs = this.outputFiles.slice(1).join(', ');
            this.completeRule = `rule ${this.ruleTitle}:\n  input: ${inputs}\n  output: ${outputs}\n  shell:\n    "${this.shellCommand}"`;
        },
        addParameter() {
            const newParam = { ...this.newParameter }; // 객체 복사
            if (newParam.type === 'inputFile' || newParam.type === 'outputFile') {
                newParam.defaultValue = newParam.name + newParam.fileExtension;
            }
            this.parameters.push(newParam);
            this.resetNewParameter();
            this.updateShellCommand();
        },
        resetNewParameter() {
            this.newParameter = {
                name: '',
                type: 'string',
                defaultValue: '',
                min: null,
                max: null,
                fileExtension: ''
            };
            this.selectedFileExtension = ''; // 선택된 파일 확장자 초기화
        },
        resetParameterOption() {
            this.newParameter.defaultValue = '';
        },
        updateFileExtension() {
            this.newParameter.fileExtension = this.selectedFileExtension;
        },
        createNode() {
            const nodeData = {
                title: this.ruleTitle,
                inputs: this.inputFiles.slice(1),
                outputs: this.outputFiles.slice(1),
                script: this.scriptFile,
            };
            const nodeId = this.$df.addNode('ruleNode', this.nodeInputs, this.nodeOutputs, 10, 10, 'ruleNode', nodeData, 'ruleNode', 'vue');

            console.log(nodeId, nodeData);
            this.onNodeCreated(nodeId, nodeData);
            this.closeModal();
        },
        onNodeCreated(id, nodeData) {
            const rule = {
                name: nodeData.title,
                input: nodeData.inputs,
                output: nodeData.outputs,
                script: nodeData.script ? nodeData.script.name : '',
                parameters: nodeData.parameters,
            };
            this.rules.push(rule);
            this.updateSnakefile();
        },
        closeModal() {
            this.isShowCreateModal = false;
            this.ruleTitle = "";
            this.nodeInputs = 0;
            this.nodeOutputs = 0;
            this.inputFiles = [];
            this.outputFiles = [];
        },
        onNodeDataChanged(id, data) {
            console.log("Node data changed", id, data);
        },
        editRule(index) {
            // Rule 편집 로직
            console.log("Editing rule", this.rules[index]);
        },
        updateSnakefile() {
            // Snakefile 업데이트 로직
        },
    },
};
</script>

<style>
.rule-container {
    width: 100%;
    height: 100%;
    position: relative;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
    background: rgb(0, 0, 0);
    background-image: url("@/assets/fantastic_background3.png");
    background-size: cover;
    background-position: center center;
    background-repeat: no-repeat;
    border-radius: 1rem;
    padding: 2rem 1rem;
    box-sizing: border-box;
}

.change-button {
    width: 2.5rem;
    height: 2.5rem;
    background-color: #007BFF;
    color: white;
    border: none;
    border-radius: 6px;
    cursor: pointer;
    position: absolute;
    top: 1rem;
    right: 1rem;
    z-index: 9999;
}

.create-button {
    width: 10rem;
    height: 2.5rem;
    background-color: #007BFF;
    color: white;
    border: none;
    border-radius: 6px;
    cursor: pointer;
    font-size: 1.1rem;
    position: absolute;
    top: 1rem;
    z-index: 9999;
}

.rule-view-container {
    width: calc(100% - 2rem);
    height: 36.5rem;
    border-radius: 1rem;
    box-shadow: 0px 0px 1px 1px rgba(0, 0, 0, 0.1);
    backdrop-filter: blur(10px);
    background: white;
    color: #333;
    overflow-y: auto;
}

.rule-item {
    padding: 1rem;
    margin-bottom: 1rem;
    border-bottom: 1px solid #eee;
    line-height: 1.4rem;
}

.rule-item h4 {
    color: #000;
    font-size: 1.1rem;
    margin-bottom: 5px;
}

.rule-item strong {
    margin-right: 5px;
}

.rule-edit-button {
    background-color: #007BFF;
    color: white;
    border: none;
    padding: 5px 10px;
    border-radius: 5px;
    cursor: pointer;
}

.rule-edit-button:hover {
    background-color: #0056b3;
}

img {
    width: 100%;
    height: 100%;
    object-fit: contain;
}

.createModal-overlay {
    position: fixed;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    display: flex;
    justify-content: center;
    align-items: center;
    z-index: 9999;
}

.createModal-container {
    background-color: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0px 0px 10px rgba(0, 0, 0, 0.1);
    width: 80%;
    max-width: 600px;
    height: 30rem;
    overflow-y: auto;
    display: flex;
    flex-direction: column;
    align-items: center;
    position: relative;
}

/* 스크롤 바 디자인 */
.createModal-container::-webkit-scrollbar {
    width: 8px;
}

.createModal-container::-webkit-scrollbar-track {
    background: #f1f1f1;
}

.createModal-container::-webkit-scrollbar-thumb {
    background: #888;
    border-radius: 4px;
}

.createModal-container::-webkit-scrollbar-thumb:hover {
    background: #555;
}

.dynamic-input-wrapper {
    display: flex;
    flex-wrap: wrap;
    gap: 1rem;
    width: 100%;
}

.dynamic-input-group {
    width: calc(50% - 1rem);
    display: flex;
    flex-direction: column;
}

.dynamic-input-group label {
    color: #000;
}

.dynamic-input-group input {
    padding: 8px;
    border: 1px solid #ccc;
    border-radius: 4px;
}

.createModal-container form {
    width: 100%;
    display: flex;
    flex-direction: column;
    align-items: center;
}

.controller-group {
    width: 100%;
    display: flex;
    align-items: center;
    margin-bottom: 0.9rem;
}

.controller-group input {
    padding: 8px;
    border: 1px solid #ccc;
    border-radius: 4px;
    width: calc(40% - 2rem);
    margin-right: 1rem;
}

.createModal-container form label {
    margin-bottom: 5px;
    color: #000;
    /* 검은색 폰트 색깔 */
    font-size: 0.9rem;
    line-height: 1rem;
}

.script-code-container {
    background: #f9f9f9;
    border: 1px solid #ccc;
    border-radius: 4px;
    padding: 10px;
    width: 100%;
    height: 10rem;
    /* 기본 높이 설정 */
    margin-bottom: 20px;
    font-size: 0.9rem;
    line-height: 1.4rem;
}

.add-parameter-group {
    width: 100%;
    display: flex;
    flex-direction: column;
    gap: 1rem;
}

.add-parameter-group .parameter-inputs,
.add-parameter-group .file-inputs,
.add-parameter-group .col-input-group {
    display: flex;
    gap: 1rem;
    width: 100%;
}

.parameter-name,
.parameter-type,
.file-extension-select,
.file-extension-input,
.default-value-input,
.boolean-select,
.range-min,
.range-max {
    flex: 1 1 48%;
    /* 가로 한 줄에 2개씩 배치 */
}

.range-inputs {
    display: flex;
    gap: 1rem;
    width: 100%;
}

.range-inputs input {
    flex: 1 1 calc(50% - 0.5rem);
}

.add-parameter-group .parameter-inputs select,
.add-parameter-group .parameter-inputs input[type="text"],
.add-parameter-group .file-inputs select,
.add-parameter-group .file-inputs input[type="text"] {
    width: 100%;
}

.add-parameter-group .range-inputs input {
    width: calc(50% - 0.5rem);
}

.button-group {
    display: flex;
    width: 5rem;
    justify-content: center;
    position: absolute;
    bottom: 1rem;
    gap: 1rem;
}

.button-group button {
    background-color: #007BFF;
    /* 파란색 버튼 */
    color: white;
    border: none;
    padding: 10px 20px;
    border-radius: 5px;
    cursor: pointer;
    margin-top: 10px;
}

.button-group button:hover {
    background-color: #0056b3;
}

.button-group button:nth-child(2) {
    background-color: #f82f2f;
    /* 파란색 버튼 */
    color: white;
    border: none;
    padding: 10px 20px;
    border-radius: 5px;
    cursor: pointer;
    margin-top: 10px;
}

.button-group button:nth-child(2):hover {
    background-color: rgb(204, 0, 0);
}

#rule-drawflow {
    width: calc(100% - 2rem);
    height: 36.5rem;
    border-radius: 1rem;
    box-shadow: 0px 0px 1px 1px rgba(255, 255, 255, 0.2);
    backdrop-filter: blur(10px);
    background: var(--dfBackgroundColor);
    background-size: var(--dfBackgroundSize) var(--dfBackgroundSize);
    background-image: var(--dfBackgroundImage);
    z-index: 9998;
}

.drawflow .drawflow-node {
    display: flex;
    background: #ffffff;
    color: #000000;
    border: 2px solid #000000;
    border-radius: 4px;
    min-height: 40px;
    width: auto;
    min-width: 160px;
    padding-top: 15px;
    padding-bottom: 15px;
    -webkit-box-shadow: 0px 2px 15px 2px #000000;
    box-shadow: 0px 2px 15px 2px #000000;
}

.drawflow .drawflow-node:hover {
    background: #ffffff;
    color: #000000;
    border: 2px solid #000000;
    border-radius: 4px;
    -webkit-box-shadow: 0px 2px 15px 2px rgba(255, 255, 255, 1);
    box-shadow: 0px 2px 15px 2px rgba(255, 255, 255, 1);
}

.drawflow .drawflow-node.selected {
    background: rgba(230, 230, 230, 0.75);
    color: rgba(0, 0, 0, 1);
    border: 2px solid #000000;
    border-radius: 4px;
    -webkit-box-shadow: 0px 2px 15px 2px rgba(0, 0, 0, 1);
    box-shadow: 0px 2px 15px 2px rgba(0, 0, 0, 1);
}

.drawflow .drawflow-node .input {
    left: -25px;
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
    height: 13px;
    width: 13px;
}

.drawflow .drawflow-node .input:hover {
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
}

.drawflow .drawflow-node .outputs {
    float: none;
}

.drawflow .drawflow-node .output {
    right: -8px;
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
    height: 13px;
    width: 13px;
}

.drawflow .drawflow-node .output:hover {
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
}

.drawflow .connection .main-path {
    stroke-width: 5px;
    stroke: #4682b4;
}

.drawflow .connection .main-path:hover {
    stroke: #4682b4;
}

.drawflow .connection .main-path.selected {
    stroke: #43b993;
}

.drawflow .connection .point {
    stroke: #000000;
    stroke-width: 2px;
    fill: #ffffff;
}

.drawflow .connection .point:hover {
    stroke: #000000;
    stroke-width: 2px;
    fill: #ffffff;
}

.drawflow-delete {
    display: block;
    color: #ffffff;
    background: #000000;
    border: 2px solid #ffffff;
    border-radius: 50px;
}

.parent-node .drawflow-delete {
    top: -15px;
}

.drawflow-delete:hover {
    color: #000000;
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
}
</style>