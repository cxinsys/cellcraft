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
                        <!-- rule title -->
                        <label for="ruleTitle">Rule Title:</label>
                        <input type="text" id="ruleTitle" v-model="ruleTitle" />

                        <label for="inputs">Number of Inputs:</label>
                        <input type="number" id="inputs" v-model.number="nodeInputs" min="0" max="6" />

                        <label for="outputs">Number of Outputs:</label>
                        <input type="number" id="outputs" v-model.number="nodeOutputs" min="0" max="6" />
                    </div>


                    <div class="dynamic-input-wrapper">
                        <div v-for="i in nodeInputs" :key="'input' + i" class="dynamic-input-group">
                            <label :for="'inputFile' + i">Input File Name {{ i }}:</label>
                            <input type="text" :id="'inputFile' + i" v-model="inputFiles[i]" />
                        </div>
                        <div v-for="i in nodeOutputs" :key="'output' + i" class="dynamic-input-group">
                            <label :for="'outputFile' + i">Output File Name {{ i }}:</label>
                            <input type="text" :id="'outputFile' + i" v-model="outputFiles[i]" />
                        </div>
                    </div>

                    <div class="script-upload">
                        <label for="scriptFile">Upload Script File:</label>
                        <input type="file" id="scriptFile" @change="handleFileUpload" />
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
            rules: [],
            editor: null,
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
                script: nodeData.script.name,
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
    width: calc(30% - 2rem);
    margin-right: 1rem;
}

.createModal-container form label {
    margin-bottom: 5px;
    color: #000;
    /* 검은색 폰트 색깔 */
    font-size: 0.9rem;
    line-height: 1rem;
}

.script-upload {
    width: 100%;
    display: flex;
    align-items: center;
    margin-top: 1rem;
}

.script-upload input {
    padding: 8px;
    border: 1px solid #ccc;
    border-radius: 4px;
    width: calc(80% - 2rem);
    margin-right: 1rem;
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

#drawflow {
  background: #ffffff;
  background-size: 0px 0px;
  background-image: none;
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