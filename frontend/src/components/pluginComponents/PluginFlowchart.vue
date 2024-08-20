<template>
    <div class="rule-container">
        <button class="change-button" @click="toggleRuleView">
            <img src="@/assets/change_circle.png">
        </button>

        <button class="create-button" @click="showCreateModal">
            Create Node <i class="fas fa-plus"></i>
        </button>

        <div v-if="!isRuleView" class="node-zoom-buttons">
            <button class="node-zoom-button" @click="zoomIn">
                <img src="@/assets/zoom_in.png">
            </button>
            <button class="node-zoom-button" @click="zoomOut">
                <img src="@/assets/zoom_out.png">
            </button>
        </div>

        <!-- Drawflow 컴포넌트 -->
        <div v-show="!isRuleView" ref="drawflow" id="rule-drawflow" @drop="drop" @dragover="allowDrop"></div>

        <!-- Rule 보기 컴포넌트 -->
        <div v-show="isRuleView" class="rule-view-container">
            <div class="warning-comment" v-if="rules.length === 0">
                <p>Rule content does not exist</p>
            </div>
            <div v-else>
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
                        <div class="script-container">
                            <div>{{ generateShellCommand(rule) }}</div>
                            <draggable class="script-drag-componenet" v-model="rule.parameters" @start="drag = true"
                                @end="drag = false">
                                <div class="script-drag-parameter" v-for="(param, paramIndex) in rule.parameters"
                                    :key="paramIndex">
                                    {{ param.name }} ({{ param.type }})
                                </div>
                            </draggable>
                        </div>
                    </div>
                    <button class="rule-remove-button" @click="removeRule(index)">remove</button>
                </div>
            </div>
        </div>

        <div v-if="isShowCreateModal" class="createModal-overlay">
            <div class="createModal-container">
                <div class="controller-group">
                    <label for="ruleTitle">Rule Title:</label>
                    <input type="text" id="ruleTitle" v-model="ruleTitle" />

                    <label for="scriptFile">Script File:</label>
                    <input type="file" id="scriptFile" @change="handleFileUpload" />
                </div>

                <div class="script-code-container">
                    <pre>{{ completeRule }}</pre>
                    <div class="checkbox-group" v-if="this.scriptFile && this.ruleTitle != ''">
                        <input type="checkbox" id="isVisualization" v-model="isVisualization"
                            @change="showAlertAndAddOutput" />
                        <label for="isVisualization">Visualization Node</label>
                    </div>
                </div>

                <div class="add-parameter-group">
                    <div class="parameter-inputs">
                        <input type="text" v-model="newParameter.name" placeholder="Parameter Name"
                            class="parameter-name" />
                        <select v-model="newParameter.type" class="parameter-type" @change="resetParameterOption">
                            <option value="inputFile">Input File</option>
                            <option value="outputFile">Output File</option>
                            <option v-if="isSelectedH5ad" value="h5adParameter">h5ad Parameter</option>
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
                            <option value="">Enter directly</option>
                            <option value=".h5ad">.h5ad</option>
                            <option value=".txt">.txt</option>
                            <option value=".csv">.csv</option>
                            <option value=".json">.json</option>
                        </select>
                    </div>
                    <div v-else-if="newParameter.type === 'h5adParameter'" class="col-input-group">
                        <select v-model="newParameter.name" class="file-extension-select">
                            <option value="">Please select h5ad parameter</option>
                            <option value="cellgroup">cellgroup</option>
                            <option value="clusters">clusters</option>
                            <option value="pseudotime">pseudotime column</option>
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
                        <div v-if="newParameter.type === 'int' || newParameter.type === 'float'" class="range-inputs">
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
                    <button type="button" @click="closeCreateModal">Cancel</button>
                </div>
            </div>
        </div>

        <alert-modal :show="alertContent.isShowAlertModal" :title="alertContent.title" :messages="alertContent.messages"
            @close="closeAlertModal" />
    </div>
</template>

<script>
/* eslint-disable */
import draggable from 'vuedraggable';
import Vue from "vue";
import ruleNode from "@/components/nodes/ruleNode.vue";
import AlertModal from '@/components/modals/AlertModal.vue';

export default {
    components: {
        draggable,
        AlertModal
    },
    props: {
        newRules: {
            type: Array,
            required: true
        },
        newDrawflow: {
            type: Object,
            required: true
        }
    },
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
            nodeData: {},
            rules: [...this.newRules],
            drawflow: { ...this.newDrawflow },
            allowRuleEdit: false,
            isVisualization: false,
            alertContent: {
                isShowAlertModal: false,
                title: "",
                messages: [],
            },
            drawflowInstance: null,
            isSelectedH5ad: false,
        };
    },
    mounted() {
        this.initDrawflow();
        const id = document.getElementById("rule-drawflow");
        Vue.prototype.$df = new Drawflow(id, Vue, this);
        //this.$df == editor
        this.$df.start();
        this.$df.registerNode('ruleNode', ruleNode, {}, {});
        this.$df.on('nodeCreated', (id) => this.onNodeCreated(id));
        this.$df.on('nodeRemoved', (id) => this.onNodeRemoved(id));
        this.$df.on('connectionCreated', (connection) => this.onConnectionCreated(connection));
        this.$df.on('connectionRemoved', (connection) => this.onConnectionRemoved(connection));
        this.$df.on('nodeDataChanged', (id) => this.onNodeDataChanged(id, data));
        this.importDrawflowData();
    },
    watch: {
        ruleTitle: 'updateCompleteRule',
        parameters: {
            handler: 'updateShellCommand',
            deep: true
        },
        newRules: {
            handler(newValue) {
                this.rules = [...newValue];
            },
            deep: true,
            immediate: true
        },
        newDrawflow: {
            handler(newValue) {
                this.drawflow = { ...newValue };
                console.log("newDrawflow", this.drawflow);
                this.importDrawflowData();
            },
            deep: true,
            immediate: true
        },
    },
    methods: {
        initDrawflow() {
            this.$nextTick(() => {
                const drawflowContainer = this.$refs.drawflow;
                if (drawflowContainer) {
                    this.drawflowInstance = new Drawflow(drawflowContainer);
                    this.drawflowInstance.start();
                    if (this.newDrawflow) {
                        this.importDrawflowData();
                    }
                } else {
                    console.error("Drawflow container not found.");
                }
            });
        },
        showCreateModal() {
            this.isShowCreateModal = true;
            // nodecreate 데이터 초기화
            this.resetNewParameter();
        },
        showAlertAndAddOutput() {
            if (this.isVisualization) {
                const title = "※ Please read the instructions for the visualization node"
                const messages = [
                    'Visualization nodes must always be configured to output a single .json file that can be uploaded to Plotly.',
                    'Visualization nodes can later select the input file from multiple files, so consider the current input file setting as a default value.'
                ]
                this.showAlertandFillContent(title, messages);

                // 기존 output 파라미터 제거
                this.parameters = this.parameters.filter(param => param.type !== 'outputFile');

                // 새로운 output 파라미터 추가
                const outputParameter = {
                    name: 'output',
                    type: 'outputFile',
                    defaultValue: `${this.ruleTitle}.json`,
                    fileExtension: '.json',
                };
                this.parameters.push(outputParameter);
            }
        },
        showAlertandFillContent(title, messages) {
            this.alertContent.title = title;
            this.alertContent.messages = messages;
            this.alertContent.isShowAlertModal = true;
        },
        closeAlertModal() {
            this.alertContent.isShowAlertModal = false;
        },
        toggleRuleView() {
            this.removeNodeIfNotExists();
            this.isRuleView = !this.isRuleView;
        },
        allowDrop(event) {
            event.preventDefault();
        },
        drop(event) {
            event.preventDefault();
        },
        zoomIn() {
            this.$df.zoom_in();
        },
        zoomOut() {
            this.$df.zoom_out();
        },
        removeNodeIfNotExists() {
            const nodeIds = this.$df.getNodesFromName('ruleNode');
            const rules = this.rules;
            nodeIds.forEach(nodeId => {
                const rule = rules.find(rule => rule.nodeId === nodeId);
                if (!rule) {
                    this.allowRuleEdit = true;
                    const removeFuncInput = "node-" + nodeId;
                    this.$df.removeNodeId(removeFuncInput);
                }
            });
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
                    return `${p.name}(${p.type})`;
                }
            }).join(' ');

            this.shellCommand = `${command} ${paramStr}`;
            this.updateCompleteRule();
        },
        updateParameters() {
            this.updateShellCommand();
        },
        updateCompleteRule() {
            const inputs = this.parameters.filter(p => p.type === 'inputFile').map(p => p.defaultValue).join(', ');
            const outputs = this.parameters.filter(p => p.type === 'outputFile').map(p => p.defaultValue).join(', ');

            this.completeRule = `rule ${this.ruleTitle}:\n  input: ${inputs}\n  output: ${outputs}\n  shell:\n    "${this.shellCommand}"`;
        },
        addParameter() {
            if (!this.scriptFile) {
                const title = "※ Please try again with the following in mind"
                const messages = [
                    'Please upload a script file first.',
                ]
                this.showAlertandFillContent(title, messages);
                return;
            }
            const newParam = { ...this.newParameter }; // 객체 복사
            // newParam.name, newParam.defaultValue 둘 중에 하나라도 비어있으면 추가하지 않음
            if (newParam.type === 'inputFile' || newParam.type === 'outputFile') {
                newParam.defaultValue = newParam.name + newParam.fileExtension;
            }
            if (!newParam.name) {
                const title = "※ Please try again with the following in mind"
                const messages = [
                    'Please fill in the parameter name.',
                ]
                this.showAlertandFillContent(title, messages);
                return;
            }
            this.parameters.push(newParam);
            if (newParam.type === 'inputFile' && newParam.fileExtension === '.h5ad') {
                this.isSelectedH5ad = true;
            }
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
            const isDuplicateTitle = this.rules.some(rule => rule.name === this.ruleTitle);
            if (isDuplicateTitle) {
                const title = "※ Please try again with the following in mind"
                const messages = [
                    'The rule title already exists. Please change the rule title.',
                ]
                this.showAlertandFillContent(title, messages);
                this.ruleTitle = "";
                return;
            }

            const inputFiles = this.parameters.filter(p => p.type === 'inputFile').map(p => p.defaultValue);
            const outputFiles = this.parameters.filter(p => p.type === 'outputFile').map(p => p.defaultValue);

            // Check if any required data is empty
            if (!this.ruleTitle || inputFiles.length === 0 || outputFiles.length === 0 || !this.scriptFile || this.parameters.length === 0) {
                const title = "※ Please try again with the following in mind"
                const messages = [
                    'All fields must be filled out.',
                ]
                this.showAlertandFillContent(title, messages);
                return;
            }

            const nodeData = {
                title: this.ruleTitle,
                inputs: inputFiles,
                outputs: outputFiles,
                script: this.scriptFile,
                parameters: this.parameters,
                isVisualization: this.isVisualization,
            };

            this.nodeData = nodeData;

            // 노드 x,y 좌표를 랜덤으로 생성
            let nodeX = Math.floor((Math.random() * 500) + 10)
            let nodeY = Math.floor((Math.random() * 500) + 10)

            if (this.isVisualization) {
                const nodeId = this.$df.addNode('ruleNode', inputFiles.length, outputFiles.length, nodeX, nodeY, 'visualizationNode', nodeData, 'ruleNode', 'vue');
                console.log(nodeId, nodeData);
            } else {
                const nodeId = this.$df.addNode('ruleNode', inputFiles.length, outputFiles.length, nodeX, nodeY, 'ruleNode', nodeData, 'ruleNode', 'vue');
                console.log(nodeId, nodeData);
            }
            this.isSelectedH5ad = false;
            this.closeCreateModal();
        },
        addNewRule(id, nodeData) {
            const rule = {
                nodeId: id,
                name: nodeData.title,
                input: nodeData.inputs,
                output: nodeData.outputs,
                // script: nodeData.script ? nodeData.script.name : '',
                script: nodeData.script,
                parameters: nodeData.parameters,
                isVisualization: nodeData.isVisualization,
            };
            this.rules.push(rule);

            // 노드 데이터 초기화
            this.nodeData = {};
        },
        closeCreateModal() {
            this.isShowCreateModal = false;
            this.isVisualization = false;
            this.ruleTitle = "";
            this.scriptFile = null;
            this.parameters = [];
            this.completeRule = "";
            this.shellCommand = "";
        },
        generateShellCommand(rule) {
            const paramStr = rule.parameters.map(p => {
                if (p.type === 'inputFile' || p.type === 'outputFile') {
                    return `${p.defaultValue}(${p.type})`;
                } else {
                    return `${p.name}(${p.type}:${p.defaultValue})`;
                }
            }).join(' ');

            const scriptName = rule.script.name || rule.script;
            let command = '';
            console.log("scriptName", scriptName);
            if (scriptName.endsWith('.py')) {
                command = `/python ${scriptName}`;
            } else if (scriptName.endsWith('.R', 'r')) {
                command = `/Rscript ${scriptName}`;
            } else {
                command = `/${scriptName}`;
            }

            // return `${command} ${paramStr}`;
            return `${command}`;
        },
        checkAndConnectNodes() {
            const nodeIds = this.$df.getNodesFromName('ruleNode');

            // 노드 ID와 규칙을 매핑
            const nodeToRuleMap = new Map();
            nodeIds.forEach(nodeId => {
                const rule = this.rules.find(rule => rule.nodeId === nodeId);
                if (rule) {
                    nodeToRuleMap.set(nodeId, rule);
                }
            });

            // 모든 노드 쌍을 비교하여 일치하는 output과 input을 연결
            nodeIds.forEach((nodeAId, i) => {
                const ruleA = nodeToRuleMap.get(nodeAId);
                if (!ruleA) return;

                nodeIds.forEach((nodeBId, j) => {
                    if (i === j) return;

                    const ruleB = nodeToRuleMap.get(nodeBId);
                    if (!ruleB) return;

                    // nodeA의 output과 nodeB의 input 비교
                    ruleA.output.forEach((output, outputIndex) => {
                        ruleB.input.forEach((input, inputIndex) => {
                            if (output === input) {
                                console.log("log connection check (A to B)", nodeAId, nodeBId, outputIndex + 1, inputIndex + 1);
                                this.allowRuleEdit = true;
                                this.$df.addConnection(nodeAId, nodeBId, `output_${outputIndex + 1}`, `input_${inputIndex + 1}`);
                            }
                        });
                    });

                    // nodeB의 output과 nodeA의 input 비교
                    ruleB.output.forEach((output, outputIndex) => {
                        ruleA.input.forEach((input, inputIndex) => {
                            if (output === input) {
                                console.log("log connection check (B to A)", nodeBId, nodeAId, outputIndex + 1, inputIndex + 1);
                                this.allowRuleEdit = true;
                                this.$df.addConnection(nodeBId, nodeAId, `output_${outputIndex + 1}`, `input_${inputIndex + 1}`);
                            }
                        });
                    });
                });
            });
        },
        onNodeCreated(id) {
            try {
                console.log("Node created", id);
                this.addNewRule(id, this.nodeData);
                this.checkAndConnectNodes();
                // this.emitRules();
                // this.emitDrawflow();
            } catch (error) {
                console.error(error);
            }
        },
        onConnectionCreated(connection) {
            if (this.allowRuleEdit) {
                try {
                    this.allowRuleEdit = false;
                    // this.emitRules();
                    // this.emitDrawflow();
                } catch (error) {
                    console.error(error);
                }
            } else {
                this.$df.removeSingleConnection(connection.output_id, connection.input_id, connection.output_class, connection.input_class);
                const title = "※ Please try again with the following in mind"
                const messages = [
                    'Connection is not allowed.',
                ]
                this.showAlertandFillContent(title, messages);
            }
        },
        onNodeRemoved(id) {
            if (this.allowRuleEdit) {
                try {
                    this.allowRuleEdit = false;
                    // 만약, drawflow 상태에서 직접 노드 제거가 허용되는 기능 제작되면 아래 주석 풀기
                    // nodeId 확인 후, 해당 노드에 대응되는 rule을 rules에서 제거
                    // const rule = this.rules.find(rule => rule.nodeId === id);
                    // const index = this.rules.indexOf(rule);
                    // this.rules.splice(index, 1);
                } catch (error) {
                    console.error(error);
                }
            } else {
                console.log("Node removed", id);
                // alert('Node removal is not allowed.');
                // // rules에서 제거된 노드에 해당하는 rule을 찾아서 다시 node를 추가
                // const rule = this.rules.find(rule => rule.nodeId === id);
                // let nodeX, nodeY = Math.floor(Math.random() * 100) + 10;
                // this.$df.addNode('ruleNode', rule.input.length, rule.output.length, nodeX, nodeY, 'ruleNode', rule, 'ruleNode', 'vue');
            }
        },
        removeRule(index) {
            const confirmed = confirm('Are you sure you want to delete this rule?');
            if (confirmed) {
                this.rules.splice(index, 1);

                // Drawflow에서 해당 노드를 찾아서 제거
                // const node = this.$df.getNodeFromId(this.rules[index].nodeId);
                // this.$df.removeNodeId(node.id);
                // console.log("node data:", node);
                this.emitRules();
            }
        },
        onConnectionRemoved(connection) {
            if (this.allowRuleEdit) {
                try {
                    this.allowRuleEdit = false;
                    // this.emitRules();
                    // this.emitDrawflow();
                } catch (error) {
                    console.error(error);
                }
            } else {
                console.log("Connection removed", connection);
                // alert('Connection removal is not allowed.');
                // // Drawflow에서 연결이 제거되지 않도록 처리
                // this.$df.addConnection(connection.output_id, connection.input_id, connection.output_class, connection.input_class);
            }
        },
        onNodeDataChanged(id) {
            // 노드 데이터가 변경되었을 때, 해당 노드에 대응되는 rule을 rules에서 찾아서 업데이트
            const node = this.$df.getNodeFromName(id);
            const rule = this.rules.find(rule => rule.nodeId === node.id);
            rule.name = node.data.title;
            // this.emitRules();
            // this.emitDrawflow();
        },
        emitRules() {
            this.$emit('update-rules', this.rules);
        },
        emitDrawflow() {
            this.drawflow = this.$df.export();
            this.$emit('update-drawflow', this.drawflow);
        },
        emitFlowchartData() {
            this.emitRules();
            this.emitDrawflow();
        },
        importDrawflowData() {
            if (this.drawflow && this.drawflow.drawflow) {
                try {
                    this.$df.clear();
                    console.log("importing drawflow data: ", this.drawflow);
                    this.$df.import(this.drawflow);
                } catch (error) {
                    console.error("Error importing drawflow data: ", error);
                }
            }
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

.node-zoom-buttons {
    position: absolute;
    bottom: 2.5rem;
    right: 2.5rem;
    display: flex;
    gap: 0.5rem;
    z-index: 9999;
}

.node-zoom-button {
    width: 2.5rem;
    height: 2.5rem;
    background-color: #007BFF;
    color: white;
    border: none;
    border-radius: 6px;
    cursor: pointer;
}

.node-zoom-button img {
    width: 100%;
    height: 100%;
    object-fit: contain;
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

.warning-comment {
    display: flex;
    justify-content: center;
    align-items: center;
    height: 100%;
}

.warning-comment p {
    font-size: 1.2rem;
    color: #333;
}

.visualizationNode {
    background-color: #4fc3f7 !important;
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

.rule-remove-button {
    /* red button */
    background-color: rgb(204, 0, 0);
    color: white;
    border: none;
    padding: 5px 10px;
    border-radius: 5px;
    cursor: pointer;
    margin-top: 0.5rem;
}

.rule-remove-button:hover {
    background-color: #ff5c5c;
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
    padding: 1.5rem;
    border-radius: 8px;
    box-shadow: 0px 0px 10px rgba(0, 0, 0, 0.1);
    width: 80%;
    max-width: 600px;
    overflow-y: auto;
    display: flex;
    flex-direction: column;
    align-items: center;
    position: relative;
    box-sizing: border-box;
    overflow-x: hidden;
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

.controller-group {
    width: 100%;
    display: flex;
    align-items: center;
    margin-bottom: 0.9rem;
    gap: 1rem;
}

.controller-group input {
    padding: 8px;
    border: 1px solid #ccc;
    border-radius: 4px;
    width: calc(40% - 4rem);
    height: 1.2rem;
    display: flex;
    align-items: center;
    justify-content: center;
}

.createModal-container label {
    color: #000;
    /* 검은색 폰트 색깔 */
    font-size: 0.9rem;
    line-height: 1rem;
}

.checkbox-group {
    position: absolute;
    top: 0.5rem;
    right: 0.5rem;
    display: flex;
    align-items: center;
}

.checkbox-group input[type="checkbox"] {
    margin-right: 0.5rem;
}

.script-code-container {
    background: #f9f9f9;
    border: 1px solid #ccc;
    border-radius: 4px;
    padding: 10px;
    width: 100%;
    height: 7.5rem;
    /* 기본 높이 설정 */
    margin-bottom: 1rem;
    font-size: 0.9rem;
    line-height: 1.4rem;
    overflow-x: auto;
    overflow-y: hidden;
    position: relative;
}

.script-container {
    display: flex;
    align-items: center;
    gap: 0.5rem
}

.script-drag-componenet {
    display: flex;
    align-items: center;
    gap: 0.5rem
}

.script-drag-parameter {
    padding: 4px;
    border: 1px solid #ccc;
    border-radius: 4px;
    background-color: #f9f9f9;
    cursor: grab;
    /* 드래그 가능함을 나타내는 커서 */
    transition: background-color 0.2s ease;
}

.script-drag-parameter:hover {
    background-color: #e0e0e0;
    /* 마우스를 올렸을 때 배경색 변경 */
    cursor: grabbing;
    /* 드래그 중일 때 커서 변경 */
}

.script-drag-parameter:active {
    background-color: #d0d0d0;
    /* 클릭 시 배경색 변경 */
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
    width: 5rem;
    display: flex;
    justify-content: center;
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
    z-index: 9997;
}

#rule-drawflow .drawflow .drawflow-node {
    display: flex;
    background: #ffffff;
    color: #000000;
    border: 2px solid #000000;
    border-radius: 1rem;
    min-height: 40px;
    width: auto;
    min-width: 160px;
    padding-top: 15px;
    padding-bottom: 15px;
    -webkit-box-shadow: 0px 2px 15px 2px #000000;
    box-shadow: 0px 2px 15px 2px #000000;
}

#rule-drawflow .drawflow .drawflow-node:hover {
    background: #ffffff;
    color: #000000;
    border: 2px solid #000000;
    border-radius: 1rem;
    -webkit-box-shadow: 0px 2px 15px 2px rgba(255, 255, 255, 1);
    box-shadow: 0px 2px 15px 2px rgba(255, 255, 255, 1);
}

#rule-drawflow .drawflow .drawflow-node.selected {
    background: rgba(230, 230, 230, 0.75);
    color: rgba(0, 0, 0, 1);
    border: 2px solid #000000;
    border-radius: 1rem;
    -webkit-box-shadow: 0px 2px 15px 2px rgba(0, 0, 0, 1);
    box-shadow: 0px 2px 15px 2px rgba(0, 0, 0, 1);
}

#rule-drawflow .drawflow .drawflow-node .input {
    left: -25px;
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
    height: 13px;
    width: 13px;
}

#rule-drawflow .drawflow .drawflow-node .input:hover {
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
}

#rule-drawflow .drawflow .drawflow-node .outputs {
    float: none;
}

#rule-drawflow .drawflow .drawflow-node .output {
    right: -8px;
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
    height: 13px;
    width: 13px;
}

#rule-drawflow .drawflow .drawflow-node .output:hover {
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
}

#rule-drawflow .drawflow .connection .main-path {
    stroke-width: 5px;
    stroke: #4682b4;
}

#rule-drawflow .drawflow .connection .main-path:hover {
    stroke: #4682b4;
}

#rule-drawflow .drawflow .connection .main-path.selected {
    stroke: #43b993;
}

#rule-drawflow .drawflow .connection .point {
    stroke: #000000;
    stroke-width: 2px;
    fill: #ffffff;
}

#rule-drawflow .drawflow .connection .point:hover {
    stroke: #000000;
    stroke-width: 2px;
    fill: #ffffff;
}

#rule-drawflow .drawflow-delete {
    /* display: block; */
    display: none;
    color: #ffffff;
    background: #000000;
    border: 2px solid #ffffff;
    border-radius: 50px;
}

#rule-drawflow .parent-node .drawflow-delete {
    top: -15px;
}

#rule-drawflow .drawflow-delete:hover {
    color: #000000;
    background: #ffffff;
    border: 2px solid #000000;
    border-radius: 50px;
}
</style>