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
            <div v-else class="rules-grid">
                <div v-for="(rule, index) in rules" :key="index" class="rule-card">
                    <!-- Rule Header -->
                    <div class="rule-header">
                        <div class="rule-title-section">
                            <div class="rule-icon">
                                <i class="fas fa-cogs"></i>
                            </div>
                            <div class="rule-title-content">
                                <h3 v-if="!rule.isEditing" class="rule-name">rule {{ rule.name }}</h3>
                                <div v-else class="rule-name-edit">
                                    <input type="text" v-model="rule.name" class="rule-name-input"
                                        @keyup.enter="toggleEditRule(index)" @keyup.escape="toggleEditRule(index)"
                                        placeholder="Rule name">
                                </div>
                                <div class="rule-meta">
                                    <span class="rule-node-id">Node ID: {{ rule.nodeId }}</span>
                                    <div v-if="rule.script" class="script-info">{{ getScriptName(rule.script) }}</div>
                                    <span v-if="rule.isVisualization" class="visualization-badge">
                                        <i class="fas fa-chart-bar"></i> Visualization
                                    </span>
                                </div>
                            </div>
                        </div>

                        <div class="rule-actions">
                            <button class="rule-action-btn move-btn" @click="moveRuleUp(index)" :disabled="index === 0"
                                title="Move Up">
                                <i class="fas fa-arrow-up"></i>
                            </button>
                            <button class="rule-action-btn move-btn" @click="moveRuleDown(index)"
                                :disabled="index === rules.length - 1" title="Move Down">
                                <i class="fas fa-arrow-down"></i>
                            </button>
                            <button class="rule-action-btn edit-btn" @click="toggleEditRule(index)"
                                :title="rule.isEditing ? 'Save' : 'Edit'">
                                <i :class="rule.isEditing ? 'fas fa-check' : 'fas fa-edit'"></i>
                            </button>
                            <button class="rule-action-btn delete-btn" @click="removeRule(index)" title="Delete Rule">
                                <i class="fas fa-trash"></i>
                            </button>
                        </div>
                    </div>

                    <!-- Rule Content -->
                    <div class="rule-content">
                        <!-- Input Section -->
                        <div class="rule-section">
                            <div class="section-header">
                                <div class="section-title">
                                    <i class="fas fa-arrow-right section-icon"></i>
                                    <strong>input:</strong>
                                </div>
                                <div class="section-count">{{ getInputFiles(rule).length }} files</div>
                            </div>
                            <div class="section-content">
                                <div v-if="getInputFiles(rule).length === 0" class="empty-state">
                                    No input files
                                </div>
                                <div v-else class="file-list">
                                    <div v-for="(inputFile, fileIndex) in getInputFiles(rule)"
                                        :key="'input-' + fileIndex" class="file-item input-file">
                                        <div class="file-icon">
                                            <i class="fas fa-file-import"></i>
                                        </div>
                                        <div class="file-info">
                                            <span class="file-name">{{ inputFile }}</span>
                                            <span class="file-type">{{ getFileExtension(inputFile) }}</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>

                        <!-- Output Section -->
                        <div class="rule-section">
                            <div class="section-header">
                                <div class="section-title">
                                    <i class="fas fa-arrow-left section-icon"></i>
                                    <strong>output:</strong>
                                </div>
                                <div class="section-count">{{ getOutputFiles(rule).length }} files</div>
                            </div>
                            <div class="section-content">
                                <div v-if="getOutputFiles(rule).length === 0" class="empty-state">
                                    No output files
                                </div>
                                <div v-else class="file-list">
                                    <div v-for="(outputFile, fileIndex) in getOutputFiles(rule)"
                                        :key="'output-' + fileIndex" class="file-item output-file">
                                        <div class="file-icon">
                                            <i class="fas fa-file-export"></i>
                                        </div>
                                        <div class="file-info">
                                            <span class="file-name">{{ outputFile }}</span>
                                            <span class="file-type">{{ getFileExtension(outputFile) }}</span>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>

                        <!-- Parameters Section -->
                        <div class="rule-section" v-if="rule.script">
                            <div class="section-header">
                                <div class="section-title">
                                    <i class="fas fa-sliders-h section-icon"></i>
                                    <strong>params:</strong>
                                    <button class="add-param-btn-inline" @click="addParameterToRule(rule, index)"
                                        title="Add Parameter">
                                        <i class="fas fa-plus"></i>
                                    </button>
                                </div>
                                <div class="section-count">{{ rule.parameters.length }} parameters</div>
                            </div>
                            <div class="section-content">
                                <div v-if="rule.parameters.length === 0" class="empty-state">
                                    No parameters configured
                                </div>
                                <div v-else class="parameters-container">
                                    <draggable class="script-drag-component" v-model="rule.parameters"
                                        @start="drag = true" @end="drag = false" handle=".drag-handle">
                                        <div class="script-drag-parameter"
                                            v-for="(param, paramIndex) in rule.parameters" :key="paramIndex"
                                            :class="{ 'editing': param.isEditing }">

                                            <!-- 일반 보기 모드 -->
                                            <div v-if="!param.isEditing" class="param-view-mode">
                                                <div class="drag-handle" title="Drag to reorder">
                                                    <i class="fas fa-grip-vertical"></i>
                                                </div>
                                                <span class="param-name">{{ param.name }}</span>
                                                <span class="param-type" :class="'type-' + param.type">{{ param.type
                                                    }}</span>
                                                <div class="param-actions-inline">
                                                    <button class="action-btn edit-btn"
                                                        @click="toggleParameterEdit(rule, paramIndex)" title="Edit">
                                                        <i class="fas fa-edit"></i>
                                                    </button>
                                                    <button class="action-btn duplicate-btn"
                                                        @click="duplicateParameter(rule, paramIndex)" title="Duplicate">
                                                        <i class="fas fa-copy"></i>
                                                    </button>
                                                    <button class="action-btn delete-btn"
                                                        @click="removeParameter(rule, paramIndex)" title="Delete">
                                                        <i class="fas fa-trash"></i>
                                                    </button>
                                                </div>
                                            </div>

                                            <!-- 편집 모드 -->
                                            <div v-else class="param-edit-mode">
                                                <div class="edit-form-compact">
                                                    <div class="edit-row">
                                                        <input type="text" v-model="param.name"
                                                            class="edit-input-compact" placeholder="Name"
                                                            @keyup.enter="toggleParameterEdit(rule, paramIndex)"
                                                            @keyup.escape="cancelParameterEdit(rule, paramIndex)">
                                                        <select v-model="param.type" class="edit-select-compact"
                                                            @change="resetParameterEditOption(param)">
                                                            <option value="inputFile">Input File</option>
                                                            <option value="optionalInputFile">Input File(Optional)
                                                            </option>
                                                            <option value="outputFile">Output File</option>
                                                            <option v-if="isSelectedH5adInEdit(rule)"
                                                                value="h5adParameter">h5ad Parameter</option>
                                                            <option value="string">String</option>
                                                            <option value="int">Integer</option>
                                                            <option value="float">Float</option>
                                                            <option value="boolean">Boolean</option>
                                                        </select>
                                                    </div>
                                                    <!-- File type parameters -->
                                                    <div v-if="param.type === 'inputFile' || param.type === 'outputFile' || param.type === 'optionalInputFile'"
                                                        class="edit-row">
                                                        <input type="text" v-model="param.fileExtension"
                                                            class="edit-input-compact" placeholder="File Type">
                                                        <select v-model="param.selectedFileExtension"
                                                            @change="updateEditFileExtension(param)"
                                                            class="edit-select-compact">
                                                            <option value="">Enter directly</option>
                                                            <option value=".h5ad">.h5ad</option>
                                                            <option value=".sif">.sif</option>
                                                            <option value=".txt">.txt</option>
                                                            <option value=".csv">.csv</option>
                                                            <option value=".json">.json</option>
                                                        </select>
                                                    </div>
                                                    <!-- h5ad parameter type -->
                                                    <div v-else-if="param.type === 'h5adParameter'" class="edit-row">
                                                        <select v-model="param.name" class="edit-select-compact">
                                                            <option value="">Please select h5ad parameter</option>
                                                            <option value="cell group">cell group</option>
                                                            <option value="pseudotime">pseudotime column</option>
                                                            <option value="clusters">clusters</option>
                                                            <option value="UMAP lasso">UMAP lasso</option>
                                                        </select>
                                                    </div>
                                                    <!-- Other types -->
                                                    <div v-else class="edit-row">
                                                        <input type="text"
                                                            v-if="param.type !== 'inputFile' && param.type !== 'outputFile' && param.type !== 'optionalInputFile' && param.type !== 'boolean'"
                                                            v-model="param.defaultValue" class="edit-input-compact"
                                                            placeholder="Default Value">
                                                        <select v-if="param.type === 'boolean'"
                                                            v-model="param.defaultValue" class="edit-select-compact">
                                                            <option value="">Select Value</option>
                                                            <option value="true">True</option>
                                                            <option value="false">False</option>
                                                        </select>
                                                        <div v-if="param.type === 'int' || param.type === 'float'"
                                                            class="edit-range-inputs">
                                                            <input type="number" v-model.number="param.min"
                                                                placeholder="Min" class="edit-input-compact" />
                                                            <input type="number" v-model.number="param.max"
                                                                placeholder="Max" class="edit-input-compact" />
                                                        </div>
                                                    </div>
                                                    <div class="edit-actions">
                                                        <button class="action-btn save-btn"
                                                            @click="toggleParameterEdit(rule, paramIndex)" title="Save">
                                                            <i class="fas fa-check"></i>
                                                        </button>
                                                        <button class="action-btn cancel-btn"
                                                            @click="cancelParameterEdit(rule, paramIndex)"
                                                            title="Cancel">
                                                            <i class="fas fa-times"></i>
                                                        </button>
                                                    </div>
                                                </div>
                                            </div>
                                        </div>
                                    </draggable>
                                </div>
                            </div>
                        </div>


                    </div>
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
                            @change="showAlertAndAddFile" />
                        <label for="isVisualization">Visualization Node</label>
                    </div>
                    <draggable class="script-drag-component" v-model="parameters" @start="drag = true"
                        @end="drag = false">
                        <div class="script-drag-parameter" v-for="(param, paramIndex) in parameters" :key="paramIndex">
                            {{ param.name }} ({{ param.type }})
                        </div>
                    </draggable>
                </div>

                <div class="add-parameter-group">
                    <div class="parameter-inputs">
                        <input type="text" v-model="newParameter.name" placeholder="Parameter Name"
                            class="parameter-name" />
                        <select v-model="newParameter.type" class="parameter-type" @change="resetParameterOption">
                            <option value="inputFile">Input File</option>
                            <option value="optionalInputFile">Input File(Optional)</option>
                            <option value="outputFile">Output File</option>
                            <option v-if="isSelectedH5ad" value="h5adParameter">h5ad Parameter</option>
                            <option value="string">String</option>
                            <option value="int">Integer</option>
                            <option value="float">Float</option>
                            <option value="boolean">Boolean</option>
                        </select>
                    </div>
                    <div v-if="newParameter.type === 'inputFile' || newParameter.type === 'outputFile' || newParameter.type === 'optionalInputFile'"
                        class="file-inputs">
                        <input type="text" v-model="newParameter.fileExtension" placeholder="File Type"
                            class="file-extension-input" />
                        <select v-model="selectedFileExtension" @change="updateFileExtension"
                            class="file-extension-select">
                            <option value="">Enter directly</option>
                            <option value=".h5ad">.h5ad</option>
                            <option value=".sif">.sif</option>
                            <option value=".txt">.txt</option>
                            <option value=".csv">.csv</option>
                            <option value=".json">.json</option>
                        </select>
                    </div>
                    <div v-else-if="newParameter.type === 'h5adParameter'" class="col-input-group">
                        <select v-model="newParameter.name" class="file-extension-select">
                            <option value="">Please select h5ad parameter</option>
                            <option value="cell group">cell group</option>
                            <option value="pseudotime">pseudotime column</option>
                            <option value="clusters">clusters</option>
                            <option value="UMAP lasso">UMAP lasso</option>
                        </select>
                    </div>
                    <div v-else class="col-input-group">
                        <input type="text"
                            v-if="newParameter.type !== 'inputFile' && newParameter.type !== 'outputFile' && newParameter.type !== 'optionalInputFile' && newParameter.type !== 'boolean'"
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
            rules: this.newRules.map(rule => ({
                ...rule,
                isEditing: false
            })),
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
            showParameterDetails: {},
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
        showAlertAndAddFile() {
            if (this.isVisualization) {
                const title = "※ Please read the instructions for the visualization node"
                const messages = [
                    'Visualization nodes must always be configured to output a single .json file that can be uploaded to Plotly.',
                    'Visualization nodes can later select the input file from multiple files, so consider the current input file setting as a default value.'
                ]
                this.showAlertandFillContent(title, messages);

                // 기존 output 파라미터 제거
                this.parameters = this.parameters.filter(param => param.type !== 'outputFile');

                // 새로운 파라미터 추가
                const outputParameter = {
                    name: `${this.ruleTitle}`,
                    type: 'outputFile',
                    defaultValue: `${this.ruleTitle}.json`,
                    fileExtension: '.json',
                };
                // const inputParameter = {
                //     name: 'target',
                //     type: 'inputFile', 
                //     defaultValue: 'target.sif',
                //     fileExtension: '.sif'
                // };
                this.parameters.push(outputParameter);
                // this.parameters.push(inputParameter);
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
                if (p.type === 'inputFile' || p.type === 'outputFile' || p.type === 'optionalInputFile') {
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
            // inputFile과 outputFile 타입의 요소를 필터링 후 defaultValue를 가져옴
            const inputs = this.parameters
                .filter(p => p.type === 'inputFile' || p.type === 'optionalInputFile')
                .map(p => p.defaultValue)
                .join(', ');

            const outputs = this.parameters
                .filter(p => p.type === 'outputFile')
                .map(p => p.defaultValue)
                .join(', ');

            // inputFile, outputFile 제외한 나머지 요소들의 defaultValue 가져오기
            const params = this.parameters
                .filter(p => p.type !== 'inputFile' && p.type !== 'outputFile' && p.type !== 'optionalInputFile')
                .map(p => `${p.name}(${p.type})`)
                .join(', ');

            // completeRule을 구성
            this.completeRule = `rule ${this.ruleTitle}:\n` +
                `  input: ${inputs}\n` +
                `  output: ${outputs}\n` +
                `  params: ${params}\n` +
                `  shell:\n`;
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
            if (newParam.type === 'inputFile' || newParam.type === 'outputFile' || newParam.type === 'optionalInputFile') {
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

            const inputFiles = this.parameters.filter(p => p.type === 'inputFile' || p.type === 'optionalInputFile').map(p => p.defaultValue);
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
                if (p.type === 'inputFile' || p.type === 'outputFile' || p.type === 'optionalInputFile') {
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
        toggleEditRule(index) {
            this.$set(this.rules, index, {
                ...this.rules[index],
                isEditing: !this.rules[index].isEditing
            });
        },
        moveRuleUp(index) {
            if (index > 0) {
                const temp = this.rules[index];
                this.$set(this.rules, index, this.rules[index - 1]);
                this.$set(this.rules, index - 1, temp);
                this.emitRules();
            }
        },
        moveRuleDown(index) {
            if (index < this.rules.length - 1) {
                const temp = this.rules[index];
                this.$set(this.rules, index, this.rules[index + 1]);
                this.$set(this.rules, index + 1, temp);
                this.emitRules();
            }
        },
        // 정리: 불필요한 메서드들 제거됨
        toggleParameterEdit(rule, paramIndex) {
            const param = rule.parameters[paramIndex];
            if (param.isEditing) {
                // 저장 로직
                this.saveParameterChanges(rule, paramIndex);
            } else {
                // 편집 모드 진입 전 백업
                this.backupParameterData(rule, paramIndex);
            }
            this.$set(param, 'isEditing', !param.isEditing);
        },
        saveParameterChanges(rule, paramIndex) {
            const param = rule.parameters[paramIndex];

            // 타입별 검증
            if (!this.validateParameter(param)) {
                return false;
            }

            // 파일 타입의 경우 defaultValue 업데이트
            if (param.type === 'inputFile' || param.type === 'outputFile' || param.type === 'optionalInputFile') {
                param.defaultValue = param.name + (param.fileExtension || '');
            }

            // h5adParameter 타입의 경우 defaultValue를 name과 동일하게 설정
            if (param.type === 'h5adParameter') {
                param.defaultValue = param.name;
            }

            // 백업 데이터 삭제
            delete param._backup;

            this.emitRules();
            return true;
        },
        backupParameterData(rule, paramIndex) {
            const param = rule.parameters[paramIndex];
            param._backup = {
                ...param,
                selectedFileExtension: param.selectedFileExtension || ''
            };
        },
        cancelParameterEdit(rule, paramIndex) {
            const param = rule.parameters[paramIndex];
            if (param._backup) {
                // 백업 데이터로 복원
                Object.keys(param._backup).forEach(key => {
                    if (key !== '_backup') {
                        param[key] = param._backup[key];
                    }
                });
                delete param._backup;
            }
            this.$set(param, 'isEditing', false);
        },
        validateParameter(param) {
            if (!param.name || param.name.trim() === '') {
                this.showAlertandFillContent(
                    "Validation Error",
                    ["Parameter name is required."]
                );
                return false;
            }

            if (param.type === 'int' || param.type === 'float') {
                if (param.min !== null && param.max !== null && param.min > param.max) {
                    this.showAlertandFillContent(
                        "Validation Error",
                        ["Minimum value cannot be greater than maximum value."]
                    );
                    return false;
                }
            }

            return true;
        },
        addParameterToRule(rule, ruleIndex) {
            const newParam = {
                name: 'new_parameter',
                type: 'string',
                defaultValue: '',
                isEditing: true,
                _isNew: true
            };

            rule.parameters.push(newParam);
            this.$set(this.showParameterDetails, `${ruleIndex}-${rule.parameters.length - 1}`, true);
        },
        removeParameter(rule, paramIndex) {
            const confirmed = confirm('이 파라미터를 삭제하시겠습니까?');
            if (confirmed) {
                rule.parameters.splice(paramIndex, 1);
                this.emitRules();
            }
        },
        duplicateParameter(rule, paramIndex) {
            const originalParam = rule.parameters[paramIndex];
            const duplicatedParam = {
                ...originalParam,
                name: originalParam.name + '_copy',
                isEditing: false,
                _backup: undefined
            };

            rule.parameters.splice(paramIndex + 1, 0, duplicatedParam);
            this.emitRules();
        },
        moveParameterUp(rule, paramIndex) {
            if (paramIndex > 0) {
                const temp = rule.parameters[paramIndex];
                this.$set(rule.parameters, paramIndex, rule.parameters[paramIndex - 1]);
                this.$set(rule.parameters, paramIndex - 1, temp);
                this.emitRules();
            }
        },
        moveParameterDown(rule, paramIndex) {
            if (paramIndex < rule.parameters.length - 1) {
                const temp = rule.parameters[paramIndex];
                this.$set(rule.parameters, paramIndex, rule.parameters[paramIndex + 1]);
                this.$set(rule.parameters, paramIndex + 1, temp);
                this.emitRules();
            }
        },
        // Helper methods for new rule display
        getInputFiles(rule) {
            if (!rule.parameters) return [];
            return rule.parameters
                .filter(param => param.type === 'inputFile' || param.type === 'optionalInputFile')
                .map(param => param.defaultValue || param.name)
                .filter(value => value);
        },
        getOutputFiles(rule) {
            if (!rule.parameters) return [];
            return rule.parameters
                .filter(param => param.type === 'outputFile')
                .map(param => param.defaultValue || param.name)
                .filter(value => value);
        },
        getFileExtension(filename) {
            if (!filename) return '';
            const parts = filename.split('.');
            return parts.length > 1 ? '.' + parts[parts.length - 1] : '';
        },
        getScriptName(script) {
            if (!script) return '';
            if (typeof script === 'string') return script;
            return script.name || 'Unknown script';
        },
        getParametersByType(rule, type) {
            if (!rule.parameters) return [];
            return rule.parameters.filter(param => param.type === type);
        },
        // 편집 모드용 헬퍼 메서드들
        resetParameterEditOption(param) {
            param.defaultValue = '';
            param.fileExtension = '';
            param.selectedFileExtension = '';
            param.min = null;
            param.max = null;
        },
        updateEditFileExtension(param) {
            param.fileExtension = param.selectedFileExtension;
        },
        isSelectedH5adInEdit(rule) {
            if (!rule.parameters) return false;
            return rule.parameters.some(param =>
                param.type === 'inputFile' && param.fileExtension === '.h5ad'
            );
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
    background: #191f26;
    /* background-image: url("@/assets/fantastic_background3.png"); */
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
    background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
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
    height: 9.5rem;
    /* 기본 높이 설정 */
    margin-bottom: 1rem;
    font-size: 0.9rem;
    line-height: 1.4rem;
    overflow-x: auto;
    overflow-y: auto;
    position: relative;
}

.script-container {
    display: flex;
    align-items: center;
    gap: 0.5rem
}

.script-drag-component {
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
    background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    color: #2c3e50;
    border: 1px solid rgba(0, 0, 0, 0.1);
    border-radius: 1rem;
    min-height: 40px;
    width: auto;
    min-width: 160px;
    padding-top: 15px;
    padding-bottom: 15px;
    backdrop-filter: blur(10px);
    -webkit-box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
}

#rule-drawflow .drawflow .drawflow-node:hover {
    background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%);
    color: #1565c0;
    border: 1px solid rgba(21, 101, 192, 0.3);
    border-radius: 1rem;
    transform: translateY(-2px);
    -webkit-box-shadow: 0 8px 24px rgba(0, 0, 0, 0.2);
    box-shadow: 0 8px 24px rgba(0, 0, 0, 0.2);
}

#rule-drawflow .drawflow .drawflow-node.selected {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    border: 1px solid rgba(102, 126, 234, 0.5);
    border-radius: 1rem;
    transform: translateY(-2px);
    -webkit-box-shadow: 0 12px 32px rgba(102, 126, 234, 0.4);
    box-shadow: 0 12px 32px rgba(102, 126, 234, 0.4);
}

#rule-drawflow .drawflow .drawflow-node .input {
    left: -25px;
    background: linear-gradient(135deg, #4fc3f7 0%, #29b6f6 100%);
    border: 2px solid rgba(79, 195, 247, 0.3);
    border-radius: 50px;
    height: 13px;
    width: 13px;
    box-shadow: 0 2px 8px rgba(79, 195, 247, 0.3);
    transition: all 0.2s ease;
}

#rule-drawflow .drawflow .drawflow-node .input:hover {
    background: linear-gradient(135deg, #29b6f6 0%, #0288d1 100%);
    border: 2px solid rgba(41, 182, 246, 0.6);
    border-radius: 50px;
    transform: scale(1.2);
    box-shadow: 0 4px 12px rgba(41, 182, 246, 0.5);
}

#rule-drawflow .drawflow .drawflow-node .outputs {
    float: none;
}

#rule-drawflow .drawflow .drawflow-node .output {
    right: -8px;
    background: linear-gradient(135deg, #ff7043 0%, #ff5722 100%);
    border: 2px solid rgba(255, 112, 67, 0.3);
    border-radius: 50px;
    height: 13px;
    width: 13px;
    box-shadow: 0 2px 8px rgba(255, 112, 67, 0.3);
    transition: all 0.2s ease;
}

#rule-drawflow .drawflow .drawflow-node .output:hover {
    background: linear-gradient(135deg, #ff5722 0%, #e64a19 100%);
    border: 2px solid rgba(255, 87, 34, 0.6);
    border-radius: 50px;
    transform: scale(1.2);
    box-shadow: 0 4px 12px rgba(255, 87, 34, 0.5);
}

#rule-drawflow .drawflow .connection .main-path {
    stroke-width: 3px;
    stroke: #667eea;
    filter: drop-shadow(0 2px 4px rgba(102, 126, 234, 0.3));
    transition: all 0.2s ease;
}

#rule-drawflow .drawflow .connection .main-path:hover {
    stroke: #5a6fd8;
    stroke-width: 4px;
    filter: drop-shadow(0 4px 8px rgba(102, 126, 234, 0.5));
}

#rule-drawflow .drawflow .connection .main-path.selected {
    stroke: #4fc3f7;
    stroke-width: 4px;
    filter: drop-shadow(0 4px 12px rgba(79, 195, 247, 0.6));
}

#rule-drawflow .drawflow .connection .point {
    stroke: #667eea;
    stroke-width: 2px;
    fill: #ffffff;
    filter: drop-shadow(0 2px 4px rgba(102, 126, 234, 0.3));
    transition: all 0.2s ease;
}

#rule-drawflow .drawflow .connection .point:hover {
    stroke: #4fc3f7;
    stroke-width: 3px;
    fill: #e3f2fd;
    filter: drop-shadow(0 4px 8px rgba(79, 195, 247, 0.5));
    transform: scale(1.2);
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

.rule-button-group {
    display: flex;
    gap: 0.5rem;
    margin-top: 0.5rem;
    align-items: center;
}

.rule-edit-button {
    background-color: #4CAF50;
    color: white;
    border: none;
    padding: 5px 10px;
    border-radius: 5px;
    cursor: pointer;
}

.rule-edit-button:hover {
    background-color: #45a049;
}

.rule-remove-button {
    /* red button */
    background-color: rgb(204, 0, 0);
    color: white;
    border: none;
    padding: 5px 10px;
    border-radius: 5px;
    cursor: pointer;
}

.rule-remove-button:hover {
    background-color: #ff5c5c;
}

.arrow-buttons {
    display: flex;
    gap: 0.25rem;
}

.rule-arrow-button {
    background-color: #6c757d;
    color: white;
    border: none;
    padding: 5px 10px;
    border-radius: 5px;
    cursor: pointer;
}

.rule-arrow-button:hover:not(:disabled) {
    background-color: #5a6268;
}

.rule-arrow-button:disabled {
    background-color: #ccc;
    cursor: not-allowed;
}

.edit-input {
    background: transparent;
    border: 1px solid #ccc;
    border-radius: 3px;
    padding: 2px 4px;
    font-size: inherit;
    width: auto;
    min-width: 50px;
}

.parameters-grid {
    display: flex;
    flex-direction: column;
    gap: 1rem;
    margin-top: 1rem;
}

.parameters-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 0.5rem;
    background-color: #f8f9fa;
    border-radius: 5px;
    border: 1px solid #dee2e6;
}

.parameters-header span {
    font-weight: bold;
    color: #495057;
}

.add-param-btn {
    background-color: #28a745;
    color: white;
    border: none;
    padding: 6px 12px;
    border-radius: 4px;
    cursor: pointer;
    font-size: 0.875rem;
    display: flex;
    align-items: center;
    gap: 0.25rem;
    transition: background-color 0.2s ease;
}

.add-param-btn:hover {
    background-color: #218838;
}

.parameters-list {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
    gap: 1rem;
}

.parameter-card {
    background: white;
    border: 1px solid #dee2e6;
    border-radius: 8px;
    padding: 1rem;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    transition: all 0.2s ease;
    position: relative;
}

.parameter-card:hover {
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.15);
    border-color: #007BFF;
}

.parameter-card.editing {
    border-color: #28a745;
    box-shadow: 0 0 0 2px rgba(40, 167, 69, 0.2);
}

.parameter-header {
    display: flex;
    align-items: flex-start;
    gap: 0.75rem;
    margin-bottom: 0.75rem;
}

.drag-handle {
    cursor: grab;
    color: #6c757d;
    padding: 0.25rem;
    border-radius: 3px;
    display: flex;
    align-items: center;
    justify-content: center;
    min-width: 20px;
}

.drag-handle:hover {
    background-color: #f8f9fa;
    color: #495057;
}

.drag-handle:active {
    cursor: grabbing;
}

.parameter-info {
    flex: 1;
    min-width: 0;
}

.param-name {
    font-weight: 600;
    color: #212529;
    font-size: 0.95rem;
    margin-bottom: 0.25rem;
    display: block;
}

.param-edit-input {
    width: 100%;
    padding: 0.375rem 0.5rem;
    border: 1px solid #ced4da;
    border-radius: 4px;
    font-size: 0.875rem;
    font-weight: 600;
}

.param-edit-input:focus {
    border-color: #007BFF;
    box-shadow: 0 0 0 2px rgba(0, 123, 255, 0.25);
    outline: none;
}

.param-type {
    font-size: 0.75rem;
    padding: 0.25rem 0.5rem;
    border-radius: 12px;
    font-weight: 500;
    text-transform: uppercase;
    letter-spacing: 0.025em;
}

.type-inputFile {
    background-color: #e3f2fd;
    color: #1565c0;
}

.type-outputFile {
    background-color: #f3e5f5;
    color: #7b1fa2;
}

.type-optionalInputFile {
    background-color: #e8f5e8;
    color: #2e7d32;
}

.type-string {
    background-color: #fff3e0;
    color: #f57c00;
}

.type-int,
.type-float {
    background-color: #fce4ec;
    color: #c2185b;
}

.type-boolean {
    background-color: #e0f2f1;
    color: #00695c;
}

.type-h5adParameter {
    background-color: #f1f8e9;
    color: #33691e;
}

.parameter-actions {
    display: flex;
    gap: 0.25rem;
    flex-wrap: wrap;
}

.parameter-actions button {
    padding: 0.375rem;
    border: none;
    border-radius: 4px;
    cursor: pointer;
    font-size: 0.75rem;
    width: 28px;
    height: 28px;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: all 0.2s ease;
}

.param-edit-btn {
    background-color: #007BFF;
    color: white;
}

.param-edit-btn:hover {
    background-color: #0056b3;
}

.param-cancel-btn {
    background-color: #6c757d;
    color: white;
}

.param-cancel-btn:hover {
    background-color: #545b62;
}

.param-duplicate-btn {
    background-color: #17a2b8;
    color: white;
}

.param-duplicate-btn:hover {
    background-color: #138496;
}

.param-move-btn {
    background-color: #ffc107;
    color: #212529;
}

.param-move-btn:hover:not(:disabled) {
    background-color: #e0a800;
}

.param-move-btn:disabled {
    background-color: #e9ecef;
    color: #6c757d;
    cursor: not-allowed;
}

.param-delete-btn {
    background-color: #dc3545;
    color: white;
}

.param-delete-btn:hover {
    background-color: #c82333;
}

.parameter-details {
    margin-top: 0.75rem;
    padding-top: 0.75rem;
    border-top: 1px solid #dee2e6;
}

.parameter-toggle {
    margin-top: 0.5rem;
    text-align: center;
}

.details-toggle-btn {
    background: none;
    border: none;
    padding: 0.25rem 0.5rem;
    cursor: pointer;
    color: #6c757d;
    font-size: 0.75rem;
    display: flex;
    align-items: center;
    gap: 0.25rem;
    margin: 0 auto;
    border-radius: 4px;
    transition: all 0.2s ease;
}

.details-toggle-btn:hover {
    background-color: #f8f9fa;
    color: #495057;
}

.param-field {
    margin-bottom: 0.75rem;
}

.param-field label {
    display: block;
    margin-bottom: 0.25rem;
    font-weight: 500;
    color: #495057;
    font-size: 0.875rem;
}

.param-field input,
.param-field select {
    width: 100%;
    padding: 0.375rem 0.5rem;
    border: 1px solid #ced4da;
    border-radius: 4px;
    font-size: 0.875rem;
    transition: border-color 0.15s ease-in-out, box-shadow 0.15s ease-in-out;
}

.param-field input:focus,
.param-field select:focus {
    border-color: #007BFF;
    box-shadow: 0 0 0 2px rgba(0, 123, 255, 0.25);
    outline: none;
}

.param-field input:read-only {
    background-color: #f8f9fa;
    color: #6c757d;
}

.range-inputs {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 0.5rem;
}

.range-input {
    width: 100% !important;
}

/* 반응형 디자인 */
@media (max-width: 768px) {
    .parameters-list {
        grid-template-columns: 1fr;
    }

    .parameter-actions {
        justify-content: center;
    }

    .parameter-header {
        flex-direction: column;
        align-items: stretch;
        gap: 0.5rem;
    }

    .drag-handle {
        align-self: flex-end;
    }
}

/* 드래그 중일 때 스타일 */
.sortable-ghost {
    opacity: 0.6;
    background-color: #e3f2fd;
}

.sortable-chosen {
    transform: scale(1.02);
}

/* 반응형 */
@media (max-width: 768px) {
    .script-drag-parameter {
        min-width: 200px;
    }

    .edit-row {
        flex-direction: column;
        align-items: stretch;
    }

    .edit-input-compact,
    .edit-select-compact {
        width: 100%;
    }
}

/* ============================================
   새로운 모던 Rule 카드 스타일
   ============================================ */

.rules-grid {
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
    padding: 1rem;
}

.rule-card {
    background: white;
    border-radius: 12px;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
    border: 1px solid #e1e5e9;
    overflow: hidden;
    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
}

.rule-card:hover {
    box-shadow: 0 8px 24px rgba(0, 0, 0, 0.15);
    transform: translateY(-2px);
}

/* Rule Header */
.rule-header {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 1.5rem;
    display: flex;
    justify-content: space-between;
    align-items: center;
}

.rule-title-section {
    display: flex;
    align-items: center;
    gap: 1rem;
    flex: 1;
}

.rule-icon {
    width: 3rem;
    height: 3rem;
    background: rgba(255, 255, 255, 0.2);
    border-radius: 10px;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 1.25rem;
}

.rule-title-content {
    flex: 1;
}

.rule-name {
    margin: 0;
    font-size: 1.5rem;
    font-weight: 600;
    margin-bottom: 0.5rem;
}

.rule-name-edit {
    margin-bottom: 0.5rem;
}

.rule-name-input {
    background: rgba(255, 255, 255, 0.9);
    border: 1px solid rgba(255, 255, 255, 0.3);
    border-radius: 6px;
    padding: 0.5rem;
    font-size: 1.25rem;
    font-weight: 600;
    color: #333;
    width: 100%;
    max-width: 300px;
}

.rule-meta {
    display: flex;
    align-items: center;
    gap: 1rem;
    opacity: 0.9;
}

.rule-node-id {
    font-size: 0.875rem;
    background: rgba(255, 255, 255, 0.2);
    padding: 0.25rem 0.5rem;
    border-radius: 12px;
    font-weight: 500;
}

.visualization-badge {
    background: #ff6b6b;
    color: white;
    padding: 0.25rem 0.75rem;
    border-radius: 12px;
    font-size: 0.75rem;
    font-weight: 600;
    display: flex;
    align-items: center;
    gap: 0.25rem;
}

.rule-actions {
    display: flex;
    gap: 0.5rem;
}

.rule-action-btn {
    width: 2.5rem;
    height: 2.5rem;
    background: rgba(255, 255, 255, 0.2);
    border: 1px solid rgba(255, 255, 255, 0.3);
    border-radius: 8px;
    color: white;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: all 0.2s ease;
    font-size: 0.875rem;
}

.rule-action-btn:hover:not(:disabled) {
    background: rgba(255, 255, 255, 0.3);
    transform: translateY(-1px);
}

.rule-action-btn:disabled {
    opacity: 0.5;
    cursor: not-allowed;
}

.rule-action-btn.delete-btn:hover:not(:disabled) {
    background: #ff4757;
    border-color: #ff4757;
}

/* Rule Content */
.rule-content {
    padding: 1.5rem;
}

.rule-section {
    margin-bottom: 2rem;
}

.rule-section:last-child {
    margin-bottom: 0;
}

.section-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 1rem;
    padding-bottom: 0.75rem;
    border-bottom: 2px solid #f1f3f4;
}

.section-title {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    font-size: 1.125rem;
    color: #2c3e50;
}

.section-icon {
    width: 1.5rem;
    height: 1.5rem;
    color: #667eea;
    font-size: 1rem;
}

.section-count {
    background: #667eea;
    color: white;
    padding: 0.25rem 0.75rem;
    border-radius: 12px;
    font-size: 0.75rem;
    font-weight: 600;
}

.script-info {
    background: #34495e;
    color: white;
    padding: 0.25rem 0.75rem;
    border-radius: 12px;
    font-size: 0.75rem;
    font-weight: 600;
}

.section-content {
    padding-left: 2.25rem;
}

/* File Lists */
.file-list {
    display: flex;
    flex-direction: column;
    gap: 0.75rem;
}

.file-item {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    padding: 0.75rem;
    background: #f8f9fa;
    border-radius: 8px;
    border-left: 4px solid #667eea;
    transition: all 0.2s ease;
}

.file-item:hover {
    background: #e9ecef;
    transform: translateX(4px);
}

.file-item.input-file {
    border-left-color: #28a745;
}

.file-item.output-file {
    border-left-color: #dc3545;
}

.file-icon {
    width: 2rem;
    height: 2rem;
    background: #667eea;
    color: white;
    border-radius: 6px;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 0.875rem;
}

.input-file .file-icon {
    background: #28a745;
}

.output-file .file-icon {
    background: #dc3545;
}

.file-info {
    flex: 1;
    display: flex;
    flex-direction: column;
    gap: 0.25rem;
}

.file-name {
    font-weight: 600;
    color: #2c3e50;
    /* color: white; */
    font-size: 0.95rem;
}

.file-type {
    font-size: 0.75rem;
    color: #6c757d;
    background: #e9ecef;
    padding: 0.125rem 0.5rem;
    border-radius: 8px;
    align-self: flex-start;
    font-weight: 500;
}

/* Empty States */
.empty-state {
    text-align: center;
    color: #6c757d;
    font-style: italic;
    padding: 2rem 1rem;
    background: #f8f9fa;
    border-radius: 8px;
    border: 2px dashed #dee2e6;
}

/* Shell Command */
.shell-command {
    background: #2c3e50;
    color: #ecf0f1;
    padding: 1rem;
    border-radius: 8px;
    overflow-x: auto;
    font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
}

.shell-command code {
    font-size: 0.875rem;
    line-height: 1.4;
    color: #ecf0f1;
    white-space: pre-wrap;
}

/* Parameters Container - 기존 스타일과 통합 */
.rule-section .parameters-container {
    margin-top: 0;
}

.rule-section .parameters-header-inline {
    justify-content: flex-end;
    margin-bottom: 1rem;
}

.rule-section .script-drag-component {
    background: #f8f9fa;
    border-radius: 8px;
    padding: 0.75rem;
    border: 1px solid #dee2e6;
}

/* 반응형 디자인 */
@media (max-width: 768px) {
    .rule-header {
        flex-direction: column;
        gap: 1rem;
        align-items: stretch;
    }

    .rule-title-section {
        justify-content: center;
        text-align: center;
    }

    .rule-actions {
        justify-content: center;
    }

    .section-header {
        flex-direction: column;
        gap: 0.5rem;
        align-items: stretch;
    }

    .section-content {
        padding-left: 0;
    }

    .file-list {
        gap: 0.5rem;
    }

    .file-item {
        padding: 0.5rem;
    }
}

/* 다크모드 지원 */
@media (prefers-color-scheme: dark) {
    .rule-card {
        background: #2c3e50;
        border-color: #34495e;
        color: #ecf0f1;
    }

    .section-title {
        color: #ecf0f1;
    }

    .file-item {
        background: #34495e;
        color: #ecf0f1;
    }

    .file-item:hover {
        background: #4a6741;
    }

    .empty-state {
        background: #34495e;
        color: #bdc3c7;
        border-color: #4a6741;
    }
}

/* 파라미터 컨테이너 - 컴팩트 디자인 스타일 */
.parameters-container {
    margin-top: 0.5rem;
}

.parameters-header-inline {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    margin-bottom: 0.5rem;
}

.add-param-btn-inline {
    background-color: #007BFF;
    color: white;
    border: none;
    width: 24px;
    height: 24px;
    border-radius: 4px;
    cursor: pointer;
    font-size: 0.75rem;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: background-color 0.2s ease;
}

.add-param-btn-inline:hover {
    background-color: #0056b3;
}

/* 기존 스크립트 드래그 컴포넌트 스타일 */
.script-drag-component {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    overflow-x: auto;
    overflow-y: hidden;
    padding: 0.25rem 0;
    min-height: 2.5rem;
}

.script-drag-parameter {
    display: flex;
    align-items: center;
    min-width: fit-content;
    background-color: #f9f9f9;
    border: 1px solid #ccc;
    border-radius: 6px;
    padding: 0.5rem;
    cursor: grab;
    transition: all 0.2s ease;
    white-space: nowrap;
    position: relative;
}

.script-drag-parameter:hover {
    background-color: #e9ecef;
    border-color: #007BFF;
    transform: translateY(-1px);
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

.script-drag-parameter.editing {
    background-color: #e3f2fd;
    border-color: #2196f3;
    min-width: 280px;
}

.script-drag-parameter:active {
    cursor: grabbing;
}

/* 일반 보기 모드 */
.param-view-mode {
    display: flex;
    align-items: center;
    gap: 0.375rem;
}

.drag-handle {
    color: #999;
    cursor: grab;
    display: flex;
    align-items: center;
    padding: 0.125rem;
    border-radius: 2px;
}

.drag-handle:hover {
    background-color: rgba(0, 0, 0, 0.1);
    color: #666;
}

.param-name {
    font-size: 0.85rem;
    font-weight: 500;
    color: #333;
    margin-right: 0.25rem;
}

.param-type {
    font-size: 0.7rem;
    padding: 0.125rem 0.375rem;
    border-radius: 8px;
    font-weight: 500;
    text-transform: uppercase;
    letter-spacing: 0.025em;
}

.param-actions-inline {
    display: flex;
    gap: 0.125rem;
    margin-left: 0.25rem;
}

.action-btn {
    width: 18px;
    height: 18px;
    border: none;
    border-radius: 3px;
    cursor: pointer;
    font-size: 0.7rem;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: all 0.2s ease;
    opacity: 0.7;
}

.script-drag-parameter:hover .action-btn {
    opacity: 1;
}

.edit-btn {
    background-color: #007BFF;
    color: white;
}

.edit-btn:hover {
    background-color: #0056b3;
}

.duplicate-btn {
    background-color: #17a2b8;
    color: white;
}

.duplicate-btn:hover {
    background-color: #138496;
}

.delete-btn {
    background-color: #dc3545;
    color: white;
}

.delete-btn:hover {
    background-color: #c82333;
}

.save-btn {
    background-color: #28a745;
    color: white;
}

.save-btn:hover {
    background-color: #218838;
}

.cancel-btn {
    background-color: #6c757d;
    color: white;
}

.cancel-btn:hover {
    background-color: #545b62;
}

/* 편집 모드 */
.param-edit-mode {
    width: 100%;
}

.edit-form-compact {
    display: flex;
    flex-direction: column;
    gap: 0.375rem;
    width: 100%;
}

.edit-row {
    display: flex;
    gap: 0.25rem;
    align-items: center;
}

.edit-input-compact,
.edit-select-compact {
    border: 1px solid #ced4da;
    border-radius: 3px;
    padding: 0.25rem 0.375rem;
    font-size: 0.75rem;
    background-color: white;
    transition: border-color 0.15s ease;
}

.edit-input-compact {
    flex: 1;
    min-width: 60px;
}

.edit-select-compact {
    min-width: 80px;
}

.edit-input-compact:focus,
.edit-select-compact:focus {
    border-color: #007BFF;
    outline: none;
    box-shadow: 0 0 0 1px rgba(0, 123, 255, 0.25);
}

.edit-actions {
    display: flex;
    gap: 0.25rem;
}

.edit-range-inputs {
    display: flex;
    gap: 0.25rem;
    width: 100%;
}

.edit-range-inputs .edit-input-compact {
    flex: 1;
}

/* 타입별 색상 유지 */
.type-inputFile {
    background-color: #e3f2fd;
    color: #1565c0;
}

.type-outputFile {
    background-color: #f3e5f5;
    color: #7b1fa2;
}

.type-optionalInputFile {
    background-color: #e8f5e8;
    color: #2e7d32;
}

.type-string {
    background-color: #fff3e0;
    color: #f57c00;
}

.type-int,
.type-float {
    background-color: #fce4ec;
    color: #c2185b;
}

.type-boolean {
    background-color: #e0f2f1;
    color: #00695c;
}

.type-h5adParameter {
    background-color: #f1f8e9;
    color: #33691e;
}

/* 스크롤바 개선 */
.script-drag-component::-webkit-scrollbar {
    height: 6px;
}

.script-drag-component::-webkit-scrollbar-track {
    background: #f1f1f1;
    border-radius: 3px;
}

.script-drag-component::-webkit-scrollbar-thumb {
    background: #c1c1c1;
    border-radius: 3px;
}

.script-drag-component::-webkit-scrollbar-thumb:hover {
    background: #a8a8a8;
}

/* 드래그 상태 스타일 */
.sortable-ghost {
    opacity: 0.6;
    background-color: #e3f2fd;
}

.sortable-chosen {
    transform: scale(1.02);
}
</style>