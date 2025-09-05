<template>
    <div id="layout">
        <!-- Plugin Selection Mode -->
        <div v-if="isPluginSelectionMode" class="plugin-selection-layout">
            <div class="plugin-list-section">
                <div class="section-header">
                    <h3 class="section-title">Visualization Plugins</h3>
                </div>
                <div class="plugin-list-container">
                    <div v-if="availableVisualizationPlugins.length === 0" class="empty-state">
                        <p>No visualization plugins available</p>
                    </div>
                    <div v-else class="plugin-list">
                        <button v-for="plugin in availableVisualizationPlugins" :key="plugin.name" class="plugin-item"
                            :class="{ 'selected': selectedVisualizationPlugin && selectedVisualizationPlugin.name === plugin.name }"
                            @click="selectVisualizationPlugin(plugin)" type="button">
                            <div class="plugin-name">{{ plugin.name }}</div>
                            <div class="plugin-description">{{ plugin.description || 'No description available' }}</div>
                        </button>
                    </div>
                </div>
            </div>

            <div class="visualization-list-section">
                <div class="section-header">
                    <h3 class="section-title">Available Visualizations</h3>
                </div>
                <div class="visualization-list-container">
                    <div v-if="!selectedVisualizationPlugin" class="empty-state">
                        <p>Please select a visualization plugin first</p>
                    </div>
                    <div v-else-if="visualizationScripts.length === 0" class="empty-state">
                        <p>No visualizations available in this plugin</p>
                    </div>
                    <div v-else class="visualization-list">
                        <button v-for="script in visualizationScripts" :key="script.name" class="visualization-item"
                            :class="{ 'selected': selectedScript && selectedScript.name === script.name }"
                            @click="selectVisualizationScript(script)" type="button">
                            <div class="visualization-name">{{ script.name }}</div>
                            <div class="visualization-info">
                                <span class="info-item">
                                    <svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24"
                                        fill="none" stroke="currentColor" stroke-width="2" aria-hidden="true">
                                        <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path>
                                        <polyline points="14 2 14 8 20 8"></polyline>
                                    </svg>
                                    {{ script.input ? script.input.length : 0 }} input(s)
                                </span>
                                <span class="info-item">
                                    <svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24"
                                        fill="none" stroke="currentColor" stroke-width="2" aria-hidden="true">
                                        <path d="M12 2v20M17 5H9.5a3.5 3.5 0 0 0 0 7h5a3.5 3.5 0 0 1 0 7H6"></path>
                                    </svg>
                                    {{script.parameters ? script.parameters.filter(p => p.type !== 'inputFile' &&
                                        p.type !== 'outputFile').length : 0 }} param(s)
                                </span>
                            </div>
                        </button>
                    </div>
                </div>
                <div class="section-footer" v-if="selectedScript">
                    <button class="continue-button" @click="proceedToConfiguration">
                        Continue to Configuration
                    </button>
                </div>
            </div>
        </div>

        <!-- Visualization Configuration Mode (existing layout with enhancements) -->
        <div v-else class="configuration-layout">
            <div class="plotly-layout">
                <div id="plotlyChart" ref="plotlyChart"></div>
            </div>
            <div class="options-layout">
                <!-- Back Button -->
                <div class="options__header">
                    <button class="back-button" @click="returnToPluginSelection">
                        <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none"
                            stroke="currentColor" stroke-width="2" aria-hidden="true">
                            <line x1="19" y1="12" x2="5" y2="12"></line>
                            <polyline points="12 19 5 12 12 5"></polyline>
                        </svg>
                        Change Visualization
                    </button>
                </div>

                <!-- Selected Visualization Info -->
                <div class="options__item selected-info">
                    <div class="selected-plugin">{{ selectedVisualizationPlugin?.name || '' }}</div>
                    <div class="selected-script">{{ selectedScript?.name || selectedVisualizationTitle }}</div>
                </div>

                <!-- Input Files Section -->
                <div class="options__section" v-if="inputFileParameters.length > 0">
                    <h4 class="section-subtitle">Input Files</h4>
                    <div class="options__item" v-for="(parameter, index) in inputFileParameters" :key="index">
                        <label class="options__title" :for="`input-file-${index}`">
                            {{ parameter.name }}
                            <span v-if="parameter.type === 'optionalInputFile'" class="optional-badge">(Optional)</span>
                        </label>
                        <select :id="`input-file-${index}`" class="options__item--select"
                            v-model="parameter.selectedFile" @change="updateInputFile(parameter)">
                            <option value="">{{ parameter.type === 'optionalInputFile' ? 'Not Selected' : 'Select File'
                                }}</option>
                            <option v-for="file in availableFiles" :key="file.name" :value="file.name">
                                {{ file.name }}
                            </option>
                        </select>
                    </div>
                </div>

                <!-- Parameters Section -->
                <div class="options__section" v-if="otherParameters.length > 0">
                    <h4 class="section-subtitle">Parameters</h4>
                    <div class="options__item" v-for="(parameter, index) in otherParameters" :key="index">
                        <label class="options__title" :for="`param-${index}`">{{ parameter.name }}</label>
                        <input :id="`param-${index}`" type="number"
                            v-if="parameter.type === 'int' || parameter.type === 'float'" class="options__item--input"
                            v-model="parameter.defaultValue" :step="parameter.type === 'float' ? '0.01' : '1'"
                            :min="parameter.min" :max="parameter.max" @change="onParameterChange" />
                        <input :id="`param-${index}`" type="checkbox" v-else-if="parameter.type === 'boolean'"
                            class="options__item--input" v-model="parameter.defaultValue" @change="onParameterChange">
                        <input :id="`param-${index}`" type="text" v-else-if="parameter.type === 'string'"
                            class="options__item--input" v-model="parameter.defaultValue" @change="onParameterChange" />
                    </div>
                </div>

                <!-- Action Buttons -->
                <div class="options__item options__actions">
                    <button id="reset-button" @click="resetParameters">Reset</button>
                    <button id="apply-button" @click="handleExecuteOrVisualize"
                        :disabled="(taskStatus !== 'SUCCESS' && !isConfigurationValid) || on_progress || loadingStates.execution"
                        :class="{ 'failure': showFailure, 'success': showSuccess }">
                        <p v-if="on_progress || loadingStates.execution">
                            <span class="button-loader"></span>
                            Processing...
                        </p>
                        <p v-else-if="showFailure">
                            Failed - Try Again
                        </p>
                        <p v-else-if="taskStatus === 'SUCCESS' && !hasParametersChanged">
                            Show Visualization
                        </p>
                        <p v-else>
                            Execute Visualization
                        </p>
                    </button>
                </div>
            </div>
        </div>
    </div>
</template>

<script>
import { getPlugins, getPluginInfo, getResults, runVisualization, getVisualizationResult, createTaskEventSource } from "@/api/index";
import { handleAPIError } from "@/api/common/interceptors";
import Plotly from "plotly.js-dist-min";

export default {
    data() {
        return {
            workflowId: this.$route.query.workflow_id,
            nodeId: this.$route.query.node,
            algorithmId: null,

            // Plugin Selection Mode
            isPluginSelectionMode: true,
            availableVisualizationPlugins: [],
            selectedVisualizationPlugin: null,
            visualizationScripts: [],
            selectedScript: null,

            // Configuration Mode
            selectedVisualizationParams: [],
            selectedVisualizationTitle: "",
            availableFiles: [],
            originalParameterValues: {},
            hasParametersChanged: false,

            // Visualization State
            plotlyData: null,
            layout: {},
            on_progress: false,
            taskStatus: '',
            showFailure: false,
            showSuccess: false,
            eventSources: {},

            // Error Handling
            errorState: {
                hasError: false,
                errorMessage: '',
                errorType: '',
                suggestedActions: [],
                showErrorDetails: false,
                errorDetails: ''
            },

            // Loading States
            loadingStates: {
                plugins: false,
                scripts: false,
                files: false,
                execution: false
            },

            // Legacy Support
            visualizationList: [],
            resultFileConnectionList: [],
            resultFileList: [],
            checkStatuses: {},
            disableApplyButton: false,
        };
    },
    computed: {
        inputFileParameters() {
            return this.selectedVisualizationParams.filter(p =>
                p.type === 'inputFile' || p.type === 'optionalInputFile'
            );
        },
        otherParameters() {
            return this.selectedVisualizationParams.filter(p =>
                p.type !== 'inputFile' && p.type !== 'optionalInputFile' && p.type !== 'outputFile'
            );
        },
        isConfigurationValid() {
            // Check if all required input files are selected
            const requiredInputs = this.inputFileParameters.filter(p => p.type === 'inputFile');
            const hasValidInputs = requiredInputs.every(p => p.selectedFile && p.selectedFile !== '');

            // Additional validation checks
            const hasValidPlugin = this.selectedVisualizationPlugin && this.selectedVisualizationPlugin.name;
            const hasValidScript = this.selectedScript && this.selectedScript.name;

            return hasValidInputs && hasValidPlugin && hasValidScript && !this.errorState.hasError;
        },

        canProceedToConfiguration() {
            return this.selectedVisualizationPlugin && this.selectedScript && !this.loadingStates.plugins && !this.loadingStates.scripts;
        },

        showErrorMessage() {
            return this.errorState.hasError && this.errorState.errorMessage;
        }
    },
    async mounted() {
        try {
            await this.initializeComponent();
        } catch (error) {
            console.error('Error initializing Visualization component:', error);
        }
    },
    watch: {
        selectedVisualizationParams: {
            handler(newVal) {
                if (newVal && this.selectedScript) {
                    this.saveNodeData();
                }
            },
            deep: true,
        }
    },
    beforeDestroy() {
        // Close all the event source connections before the component is destroyed
        for (let task_id in this.eventSources) {
            this.closeEventSource(task_id);
        }
    },
    methods: {
        async initializeComponent() {
            try {
                this.clearError();

                // Get current node info
                const current_node = this.$store.getters.getWorkflowNodeInfo(this.nodeId);

                if (!current_node) {
                    throw new Error(`Visualization node with ID ${this.nodeId} not found`);
                }

                // Check if node has previous selections
                if (current_node.data.selectedVisualizationPlugin && current_node.data.selectedScript) {
                    // Restore previous state
                    this.isPluginSelectionMode = false;
                    this.selectedVisualizationPlugin = current_node.data.selectedVisualizationPlugin;
                    this.selectedScript = current_node.data.selectedScript;
                    this.selectedVisualizationParams = current_node.data.selectedVisualizationParams || [];
                    this.selectedVisualizationTitle = current_node.data.selectedVisualizationTitle || this.selectedScript.name;

                    // Load available files and restore configuration
                    await this.loadAvailableFiles();
                    this.storeOriginalParameterValues();

                    // Check if visualization was already executed
                    if (current_node.data.taskStatus === 'SUCCESS') {
                        this.taskStatus = 'SUCCESS';
                        this.showSuccess = true;
                    }
                } else {
                    // New node - fetch visualization plugins
                    await this.fetchVisualizationPlugins();

                    // Check for legacy data
                    await this.checkLegacyConfiguration(current_node);
                }
            } catch (error) {
                this.handleError(error, {
                    operation: 'initialize_component',
                    nodeId: this.nodeId
                });
            }
        },

        async fetchVisualizationPlugins() {
            this.loadingStates.plugins = true;
            this.clearError();

            try {
                const response = await getPlugins();

                if (!response.data || !response.data.plugins) {
                    throw new Error('Invalid response format from server');
                }

                // Filter for visualization plugins
                this.availableVisualizationPlugins = response.data.plugins.filter(plugin =>
                    plugin.plugin_type === 'VISUALIZATION' ||
                    plugin.plugin_type === 'visualization'
                );

                if (this.availableVisualizationPlugins.length === 0) {
                    this.showWarning('No visualization plugins are currently available. Contact your administrator to install visualization plugins.');
                }

            } catch (error) {
                this.handleError(error, {
                    operation: 'fetch_visualization_plugins',
                    message: 'Failed to load visualization plugins'
                });
                this.availableVisualizationPlugins = [];
            } finally {
                this.loadingStates.plugins = false;
            }
        },

        async selectVisualizationPlugin(plugin) {
            if (!plugin || !plugin.name) {
                this.handleError(new Error('Invalid plugin selection'), {
                    operation: 'select_plugin'
                });
                return;
            }

            this.selectedVisualizationPlugin = plugin;
            this.visualizationScripts = [];
            this.selectedScript = null;
            this.loadingStates.scripts = true;
            this.clearError();

            try {
                const pluginInfo = await getPluginInfo(plugin.name);

                if (!pluginInfo.data || !pluginInfo.data.plugin) {
                    throw new Error('Invalid plugin information received from server');
                }

                if (pluginInfo.data.plugin.rules) {
                    // Extract visualization scripts
                    const rules = Object.values(pluginInfo.data.plugin.rules);

                    // For visualization plugins, all rules are considered visualization scripts
                    if (plugin.plugin_type === 'VISUALIZATION' || plugin.plugin_type === 'visualization') {
                        this.visualizationScripts = rules;
                    } else {
                        // For other plugins, only rules explicitly marked as visualization
                        this.visualizationScripts = rules.filter(rule =>
                            rule.isVisualization === true || rule.isVisualization === 'true'
                        );
                    }

                    if (this.visualizationScripts.length === 0) {
                        this.showWarning(`No visualization scripts found in plugin '${plugin.name}'. Please select a different plugin.`);
                    }

                    console.log('Plugin info:', pluginInfo.data.plugin);
                    console.log('Filtered visualization scripts:', this.visualizationScripts);
                } else {
                    this.showWarning(`Plugin '${plugin.name}' does not contain any visualization rules.`);
                }
            } catch (error) {
                this.handleError(error, {
                    operation: 'select_plugin',
                    plugin_name: plugin.name,
                    message: `Failed to load visualizations for ${plugin.name}`
                });
                this.visualizationScripts = [];
            } finally {
                this.loadingStates.scripts = false;
            }
        },

        selectVisualizationScript(script) {
            this.selectedScript = script;
        },

        async proceedToConfiguration() {
            if (!this.selectedScript) return;

            this.isPluginSelectionMode = false;
            this.selectedVisualizationParams = JSON.parse(JSON.stringify(this.selectedScript.parameters || []));
            this.selectedVisualizationTitle = this.selectedScript.name;

            // Initialize input file parameters
            this.selectedVisualizationParams.forEach(param => {
                if (param.type === 'inputFile' || param.type === 'optionalInputFile') {
                    param.selectedFile = param.defaultValue || '';
                }
            });

            await this.loadAvailableFiles();
            this.storeOriginalParameterValues();
            this.saveNodeData();
        },

        async returnToPluginSelection() {
            // Clear Plotly chart before switching modes
            if (this.$refs.plotlyChart) {
                Plotly.purge(this.$refs.plotlyChart);
            }

            // Reset UI mode
            this.isPluginSelectionMode = true;

            // Reset all selection data to initial state
            this.selectedVisualizationPlugin = null;
            this.selectedScript = null;
            this.visualizationScripts = [];
            this.selectedVisualizationParams = [];
            this.selectedVisualizationTitle = "";
            this.availableFiles = [];
            this.originalParameterValues = {};
            this.hasParametersChanged = false;

            // Reset visualization state
            this.plotlyData = null;
            this.layout = {};
            this.taskStatus = '';
            this.showSuccess = false;
            this.showFailure = false;
            this.on_progress = false;

            // Clear error states
            this.clearError();

            // Reset loading states
            this.loadingStates = {
                plugins: false,
                scripts: false,
                files: false,
                execution: false
            };

            // Clear Vuex store data by saving null/empty values
            const emptyDataObject = {
                selectedVisualizationPlugin: null,
                selectedScript: null,
                selectedVisualizationParams: [],
                selectedVisualizationTitle: "",
                taskStatus: "",
            };
            const nodeId = this.nodeId;
            this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject: emptyDataObject });

            // Reload visualization plugins to ensure fresh data
            await this.fetchVisualizationPlugins();
        },

        async loadAvailableFiles() {
            this.loadingStates.files = true;
            this.clearError();

            try {
                // Get connected ResultFiles nodes
                const current_node = this.$store.getters.getWorkflowNodeInfo(this.nodeId);

                if (!current_node) {
                    throw new Error(`Current node ${this.nodeId} not found`);
                }

                const connections = current_node.inputs.input_1?.connections || [];
                this.availableFiles = [];

                for (const connection of connections) {
                    if (!connection.node) {
                        console.warn('Invalid connection found:', connection);
                        continue;
                    }

                    const connectedNode = this.$store.getters.getWorkflowNodeInfo(connection.node);
                    if (connectedNode && connectedNode.name === 'ResultFiles') {
                        // Handle multi-file format
                        if (connectedNode.data.files && Array.isArray(connectedNode.data.files)) {
                            const selectedFiles = connectedNode.data.files.filter(f => f.selected);
                            this.availableFiles.push(...selectedFiles);
                        }
                    }
                }

                // Also check for algorithm node results (legacy support)
                if (this.availableFiles.length === 0) {
                    await this.loadLegacyFiles();
                }

                if (this.availableFiles.length === 0) {
                    this.showWarning('No input files are available. Please ensure the workflow is properly connected and executed.');
                }

            } catch (error) {
                this.handleError(error, {
                    operation: 'load_available_files',
                    nodeId: this.nodeId
                });
            } finally {
                this.loadingStates.files = false;
            }
        },

        async loadLegacyFiles() {
            try {
                // const current_node = this.$store.getters.getWorkflowNodeInfo(this.nodeId);
                const connectionAlgorithmNode = this.$store.getters.getAlgorithmNodeConnectedToInput(this.nodeId);

                if (!connectionAlgorithmNode) {
                    console.debug('No algorithm node connected for legacy file loading');
                    return;
                }

                this.algorithmId = connectionAlgorithmNode.id;

                if (!this.workflowId) {
                    throw new Error('Workflow ID is required for loading files');
                }

                const workflow_result = {
                    id: this.workflowId,
                    algorithm_id: this.algorithmId,
                };

                const response = await getResults(workflow_result);

                if (!response.data) {
                    console.warn('No data received from getResults API');
                    return;
                }

                this.availableFiles = Array.isArray(response.data) ? response.data : [];

            } catch (error) {
                // Don't show error to user for legacy file loading, just log
                console.warn('Error loading legacy files:', error);
            }
        },

        updateInputFile(parameter) {
            // Don't overwrite defaultValue as it's the placeholder name
            // The selectedFile property already contains the user's selection
            console.log('Updating input file:', parameter);
            this.onParameterChange();
        },

        onParameterChange() {
            this.hasParametersChanged = this.checkIfParametersChanged();
            if (this.hasParametersChanged && this.taskStatus === 'SUCCESS') {
                this.showSuccess = false;
            }
        },

        checkIfParametersChanged() {
            return this.selectedVisualizationParams.some(param => {
                const originalValue = this.originalParameterValues[param.name];
                return param.defaultValue !== originalValue || param.selectedFile !== originalValue;
            });
        },

        storeOriginalParameterValues() {
            this.originalParameterValues = {};
            this.selectedVisualizationParams.forEach(param => {
                this.originalParameterValues[param.name] = param.selectedFile || param.defaultValue;
            });
        },

        resetParameters() {
            // Reset to original values
            this.selectedVisualizationParams.forEach(param => {
                if (this.selectedScript) {
                    const originalParam = this.selectedScript.parameters.find(p => p.name === param.name);
                    if (originalParam) {
                        param.defaultValue = originalParam.defaultValue;
                        if (param.type === 'inputFile' || param.type === 'optionalInputFile') {
                            param.selectedFile = originalParam.defaultValue || '';
                        }
                    }
                }
            });

            this.storeOriginalParameterValues();
            this.hasParametersChanged = false;
            this.plotlyData = null;
            this.renderPlot();
        },

        async handleExecuteOrVisualize() {
            if (this.taskStatus === 'SUCCESS' && !this.hasParametersChanged) {
                // Show visualization
                this.renderPlot();
            } else {
                // Execute visualization
                await this.runVisualization();
            }
        },

        async runVisualization() {
            // Validate inputs before sending request
            if (!this.validateVisualizationInputs()) {
                return;
            }

            this.loadingStates.execution = true;
            this.on_progress = true;
            this.clearError();

            try {
                const title = this.$store.getters.getTitle;
                const thumbnail = this.$store.getters.getThumbnail;
                const workflow_info = this.$store.getters.getWorkflowInfo;

                // Resolve algorithm_id using the existing getter
                const connectionAlgorithmNode = this.$store.getters.getAlgorithmNodeConnectedToInput(this.nodeId);
                const resolvedAlgorithmId = connectionAlgorithmNode ? connectionAlgorithmNode.id : null;

                if (!resolvedAlgorithmId && this.availableFiles.length === 0) {
                    throw new Error('No algorithm connection found and no files available. Please ensure the workflow is properly connected.');
                }

                // Prepare the new WorkflowVisualizationRequest data structure
                const workflowVisualizationRequest = {
                    id: this.workflowId,
                    current_node_id: this.nodeId,
                    algorithm_id: resolvedAlgorithmId,
                    selectedVisualizationPlugin: this.selectedVisualizationPlugin,
                    selectedScript: this.selectedScript,
                    selectedVisualizationParams: this.selectedVisualizationParams,
                    availableFiles: this.availableFiles,
                    title: title,
                    thumbnail: thumbnail,
                    workflow_info: workflow_info,
                };

                console.log('Sending visualization request:', workflowVisualizationRequest);
                const workflow_data = await runVisualization(workflowVisualizationRequest);
                console.log('Visualization response:', workflow_data);

                if (!workflow_data.data) {
                    throw new Error('Invalid response from visualization service');
                }

                if (workflow_data.data.message === "Visualization result already exists" || workflow_data.data.cached) {
                    const workflowResult = {
                        id: this.workflowId,
                        algorithm_id: this.nodeId,  // Use visualization node ID for cached results
                        filename: workflow_data.data.result_path,
                    };

                    const response = await getVisualizationResult(workflowResult);
                    if (!response.data) {
                        throw new Error('Failed to retrieve cached visualization result');
                    }

                    this.plotlyData = response.data;
                    this.renderPlot();
                    this.taskStatus = 'SUCCESS';
                    this.showSuccess = true;
                    this.hasParametersChanged = false;
                    this.storeOriginalParameterValues();
                } else if (workflow_data.data.task_id) {
                    this.createEventSource(workflow_data.data.task_id);
                } else {
                    throw new Error('No task ID received from visualization service');
                }

            } catch (error) {
                this.handleError(error, {
                    operation: 'run_visualization',
                    plugin: this.selectedVisualizationPlugin?.name,
                    script: this.selectedScript?.name
                });
                this.handleTaskFailure();
            } finally {
                if (!this.eventSources || Object.keys(this.eventSources).length === 0) {
                    this.on_progress = false;
                    this.loadingStates.execution = false;
                }
            }
        },

        renderPlot() {
            try {
                if (!this.plotlyData) {
                    console.warn('No plotly data available for rendering');
                    return;
                }

                if (!this.$refs.plotlyChart) {
                    console.error('Plotly chart container not found');
                    return;
                }

                // Validate plotly data structure
                if (!this.plotlyData.data || !Array.isArray(this.plotlyData.data)) {
                    throw new Error('Invalid plotly data format: data array is missing');
                }

                Plotly.newPlot(
                    this.$refs.plotlyChart,
                    this.plotlyData.data,
                    this.plotlyData.layout || {},
                    { responsive: true }
                );

                console.log('Plot rendered successfully');

            } catch (error) {
                this.handleError(error, {
                    operation: 'render_plot',
                    message: 'Failed to render visualization plot'
                });
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
                onComplete: async (status) => {
                    console.log("Task complete:", status);
                    if (status === 'FAILURE') {
                        this.handleTaskFailure();
                    } else if (status === 'SUCCESS') {
                        await this.handleTaskSuccess(task_id);
                    }
                    this.on_progress = false;
                    this.closeEventSource(task_id);
                },
                onError: (error) => {
                    console.error("SSE Error:", error);
                    this.handleTaskError(error);
                    this.on_progress = false;
                    this.closeEventSource(task_id);
                }
            });
        },

        closeEventSource(task_id) {
            if (this.eventSources[task_id]) {
                this.eventSources[task_id].close();
                delete this.eventSources[task_id];
            }
        },

        async handleTaskSuccess(task_id) {
            try {
                this.taskStatus = 'SUCCESS';
                this.showSuccess = true;
                this.hasParametersChanged = false;
                this.storeOriginalParameterValues();
                this.clearError(); // Clear any previous errors

                // Fetch and render visualization result
                const connectionAlgorithmNode = this.$store.getters.getAlgorithmNodeConnectedToInput(this.nodeId);
                const resolvedAlgorithmId = connectionAlgorithmNode ? connectionAlgorithmNode.id : null;

                if (!resolvedAlgorithmId) {
                    throw new Error('Unable to resolve algorithm ID for result retrieval');
                }

                const workflowResult = {
                    id: this.workflowId,
                    algorithm_id: this.nodeId,  // Use visualization node ID for result retrieval
                    filename: `${this.selectedVisualizationTitle}.json`,
                };

                const response = await getVisualizationResult(workflowResult);

                if (!response.data) {
                    throw new Error('No visualization data received from server');
                }

                this.plotlyData = response.data;
                this.renderPlot();

                console.log('Visualization completed successfully');

            } catch (error) {
                this.handleError(error, {
                    operation: 'handle_task_success',
                    task_id: task_id,
                    message: 'Failed to load visualization result after successful task completion'
                });

                // Keep success state but show error for result loading
                this.taskStatus = 'SUCCESS';
                this.showSuccess = true;
            } finally {
                this.on_progress = false;
                this.loadingStates.execution = false;
                // Critical: Save the task status to persist it
                this.saveNodeData();
            }
        },

        handleTaskFailure() {
            this.taskStatus = 'FAILURE';
            this.showFailure = true;
            this.on_progress = false;
            this.loadingStates.execution = false;

            // Show failure state temporarily, but keep error details if available
            setTimeout(() => {
                this.showFailure = false;
            }, 3000);

            this.saveNodeData();
        },

        handleTaskError(error) {
            this.handleError(error, {
                operation: 'task_execution',
                message: 'Task execution failed'
            });

            this.taskStatus = 'FAILURE';
            this.showFailure = true;
            this.on_progress = false;
            this.loadingStates.execution = false;

            setTimeout(() => {
                this.showFailure = false;
            }, 3000);

            this.saveNodeData();
        },

        // Enhanced Error Handling Methods
        handleError(error, context = {}) {
            console.error('Visualization Error:', error, context);

            const errorInfo = handleAPIError(error, context);

            this.errorState = {
                hasError: true,
                errorMessage: errorInfo.message,
                errorType: errorInfo.type,
                suggestedActions: errorInfo.suggestedActions,
                showErrorDetails: false,
                errorDetails: errorInfo.details || error.stack || JSON.stringify(context)
            };

            // Clear loading states on error
            this.loadingStates = {
                plugins: false,
                scripts: false,
                files: false,
                execution: false
            };
        },

        clearError() {
            this.errorState = {
                hasError: false,
                errorMessage: '',
                errorType: '',
                suggestedActions: [],
                showErrorDetails: false,
                errorDetails: ''
            };
        },

        showWarning(message) {
            console.warn('Visualization Warning:', message);
            // Could be extended to show warning notifications
        },

        toggleErrorDetails() {
            this.errorState.showErrorDetails = !this.errorState.showErrorDetails;
        },

        retryLastOperation() {
            this.clearError();
            // Try to determine what operation to retry based on current state
            if (this.isPluginSelectionMode) {
                if (this.availableVisualizationPlugins.length === 0) {
                    this.fetchVisualizationPlugins();
                } else if (this.selectedVisualizationPlugin && this.visualizationScripts.length === 0) {
                    this.selectVisualizationPlugin(this.selectedVisualizationPlugin);
                }
            } else {
                if (this.availableFiles.length === 0) {
                    this.loadAvailableFiles();
                } else if (this.taskStatus === 'FAILURE') {
                    this.runVisualization();
                }
            }
        },

        validateVisualizationInputs() {
            if (!this.selectedVisualizationPlugin) {
                this.handleError(new Error('Please select a visualization plugin'), {
                    operation: 'validation',
                    field: 'plugin'
                });
                return false;
            }

            if (!this.selectedScript) {
                this.handleError(new Error('Please select a visualization script'), {
                    operation: 'validation',
                    field: 'script'
                });
                return false;
            }

            const requiredInputs = this.inputFileParameters.filter(p => p.type === 'inputFile');
            const missingInputs = requiredInputs.filter(p => !p.selectedFile || p.selectedFile === '');

            if (missingInputs.length > 0) {
                const fieldNames = missingInputs.map(p => p.name).join(', ');
                this.handleError(new Error(`Please select files for required inputs: ${fieldNames}`), {
                    operation: 'validation',
                    field: 'input_files',
                    missing_inputs: fieldNames
                });
                return false;
            }

            return true;
        },

        saveNodeData() {
            const dataObject = {
                selectedVisualizationPlugin: this.selectedVisualizationPlugin ? {
                    name: this.selectedVisualizationPlugin.name,
                    source: this.selectedVisualizationPlugin.source || null,
                    plugin_type: this.selectedVisualizationPlugin.plugin_type || null
                } : null,
                selectedScript: this.selectedScript,
                selectedVisualizationParams: this.selectedVisualizationParams,
                selectedVisualizationTitle: this.selectedVisualizationTitle,
                taskStatus: this.taskStatus,
            };
            const nodeId = this.nodeId;
            this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
        },

        // Legacy support methods
        async checkLegacyConfiguration(current_node) {
            // Check if this is a legacy node with pre-selected visualization
            const connectionAlgorithmNode = this.$store.getters.getAlgorithmNodeConnectedToInput(this.nodeId);

            if (connectionAlgorithmNode && connectionAlgorithmNode.data.selectedPlugin) {
                const pluginName = connectionAlgorithmNode.data.selectedPlugin.name;
                const pluginInfo = await getPluginInfo(pluginName);

                if (pluginInfo.data.plugin.rules) {
                    const rules = Object.values(pluginInfo.data.plugin.rules);
                    const visualizations = rules.filter(rule => rule.isVisualization);

                    if (visualizations.length > 0 && current_node.data.title) {
                        // Try to match by title
                        const matchingViz = visualizations.find(v => v.output[0] === current_node.data.title);
                        if (matchingViz) {
                            // Auto-select the matching visualization
                            const plugin = this.availableVisualizationPlugins.find(p => p.name === pluginName);
                            if (plugin) {
                                await this.selectVisualizationPlugin(plugin);
                                this.selectVisualizationScript(matchingViz);
                                await this.proceedToConfiguration();
                            }
                        }
                    }
                }
            }
        }
    }
}
</script>

<style scoped>
#layout {
    width: 100%;
    height: 100%;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
}

/* Error Handling Styles */
.error-container {
    width: 100%;
    padding: 1rem;
    margin-bottom: 1rem;
    box-sizing: border-box;
}

.error-message {
    background-color: #fee;
    border: 1px solid #fcc;
    border-radius: 0.5rem;
    padding: 1rem;
    color: #a00;
}

.error-header {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    margin-bottom: 0.5rem;
}

.error-icon {
    color: #d32f2f;
    flex-shrink: 0;
}

.error-title {
    font-weight: 600;
    font-size: 1rem;
    flex: 1;
    text-transform: capitalize;
}

.error-close {
    background: none;
    border: none;
    color: #a00;
    cursor: pointer;
    padding: 0.25rem;
    border-radius: 0.25rem;
    transition: background-color 0.2s;
}

.error-close:hover {
    background-color: rgba(0, 0, 0, 0.1);
}

.error-text {
    margin: 0.5rem 0;
    line-height: 1.4;
}

.error-suggestions {
    margin: 0.75rem 0;
}

.suggestions-title {
    font-weight: 500;
    margin: 0 0 0.25rem 0;
    font-size: 0.9rem;
}

.suggestions-list {
    margin: 0;
    padding-left: 1.25rem;
    font-size: 0.875rem;
}

.suggestions-list li {
    margin-bottom: 0.25rem;
}

.error-actions {
    display: flex;
    gap: 0.5rem;
    margin-top: 0.75rem;
}

.retry-button,
.details-button {
    display: flex;
    align-items: center;
    gap: 0.25rem;
    padding: 0.5rem 0.75rem;
    border: 1px solid #a00;
    background-color: #fff;
    color: #a00;
    border-radius: 0.25rem;
    font-size: 0.875rem;
    cursor: pointer;
    transition: all 0.2s;
}

.retry-button:hover {
    background-color: #a00;
    color: #fff;
}

.details-button:hover {
    background-color: rgba(170, 0, 0, 0.1);
}

.error-details {
    margin-top: 0.75rem;
    padding: 0.75rem;
    background-color: #f5f5f5;
    border-radius: 0.25rem;
    border: 1px solid #ddd;
}

.error-details pre {
    margin: 0;
    font-size: 0.75rem;
    line-height: 1.3;
    white-space: pre-wrap;
    word-wrap: break-word;
    color: #333;
}

/* Loading States */
.loading-state {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    padding: 2rem;
    color: #666;
}

.loader {
    width: 24px;
    height: 24px;
    border: 2px solid #f3f3f3;
    border-top: 2px solid #2d2fbf;
    border-radius: 50%;
    animation: spin 1s linear infinite;
    margin-bottom: 0.5rem;
}

@keyframes spin {
    0% {
        transform: rotate(0deg);
    }

    100% {
        transform: rotate(360deg);
    }
}

/* Plugin Selection Mode Styles */
.plugin-selection-layout {
    width: 100%;
    height: 100%;
    display: flex;
    gap: 1rem;
    padding: 1rem;
}

.plugin-list-section,
.visualization-list-section {
    flex: 1;
    background-color: #ffffff;
    border-radius: 1rem;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    display: flex;
    flex-direction: column;
    overflow: hidden;
}

.section-header {
    padding: 1.5rem;
    border-bottom: 1px solid #e5e5e5;
}

.section-title {
    margin: 0;
    font-size: 1.2rem;
    font-weight: 600;
    color: #333;
}

.plugin-list-container,
.visualization-list-container {
    flex: 1;
    overflow-y: auto;
    padding: 1rem;
}

.empty-state {
    display: flex;
    align-items: center;
    justify-content: center;
    height: 100%;
    color: #999;
    font-style: italic;
}

.plugin-list,
.visualization-list {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
}

.plugin-item,
.visualization-item {
    width: 100%;
    padding: 1rem;
    background-color: #f8f9fa;
    border: 2px solid transparent;
    border-radius: 0.5rem;
    cursor: pointer;
    transition: all 0.2s ease;
    text-align: left;
    font-family: inherit;
    font-size: inherit;
}

.plugin-item:hover,
.visualization-item:hover {
    background-color: #e9ecef;
    border-color: #dee2e6;
}

.plugin-item.selected,
.visualization-item.selected {
    background-color: #e7eaff;
    border-color: #2d2fbf;
}

.plugin-name,
.visualization-name {
    font-weight: 600;
    color: #333;
    margin-bottom: 0.25rem;
}

.plugin-description {
    font-size: 0.875rem;
    color: #666;
    overflow: hidden;
    text-overflow: ellipsis;
    display: -webkit-box;
    -webkit-line-clamp: 2;
    -webkit-box-orient: vertical;
}

.visualization-info {
    display: flex;
    gap: 1rem;
    margin-top: 0.5rem;
}

.info-item {
    display: flex;
    align-items: center;
    gap: 0.25rem;
    font-size: 0.875rem;
    color: #666;
}

.info-item svg {
    flex-shrink: 0;
}

.section-footer {
    padding: 1rem;
    border-top: 1px solid #e5e5e5;
}

.continue-button {
    width: 100%;
    padding: 0.75rem;
    background-color: #2d2fbf;
    color: white;
    border: none;
    border-radius: 0.5rem;
    font-size: 1rem;
    font-weight: 500;
    cursor: pointer;
    transition: background-color 0.2s ease;
}

.continue-button:hover:not(:disabled) {
    background-color: #4655ff;
}

.continue-button:disabled {
    background-color: #ccc;
    color: #666;
    cursor: not-allowed;
}

.continue-button:disabled .button-loader {
    border-top-color: #666;
}

/* Configuration Mode Styles */
.configuration-layout {
    width: 100%;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: row;
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

#plotlyChart {
    width: 100%;
    height: 100%;
}

.options-layout {
    width: 30%;
    height: 95%;
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
    padding: 1rem;
    background-color: #f8f9fa;
    border-radius: 1rem;
}

.options__header {
    margin-bottom: 1rem;
}

.back-button {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.5rem 1rem;
    background-color: #fff;
    border: 1px solid #dee2e6;
    border-radius: 0.5rem;
    color: #333;
    font-size: 0.875rem;
    cursor: pointer;
    transition: all 0.2s ease;
}

.back-button:hover {
    background-color: #e9ecef;
    border-color: #adb5bd;
}

.selected-info {
    background-color: #fff;
    padding: 1rem;
    border-radius: 0.5rem;
    margin-bottom: 1rem;
}

.selected-plugin {
    font-size: 0.875rem;
    color: #666;
    margin-bottom: 0.25rem;
}

.selected-script {
    font-size: 1.1rem;
    font-weight: 600;
    color: #333;
}

.options__section {
    margin-bottom: 1.5rem;
}

.section-subtitle {
    margin: 0 0 1rem 0;
    font-size: 1rem;
    font-weight: 600;
    color: #333;
}

.options__item {
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin-bottom: 0.75rem;
    gap: 1rem;
}

.options__title {
    flex: 1;
    font-weight: 500;
    color: #333;
    margin: 0;
}

.optional-badge {
    font-size: 0.75rem;
    color: #666;
    font-weight: normal;
}

.options__item--select,
.options__item--input {
    width: 50%;
    padding: 0.5rem;
    border: 1px solid #dee2e6;
    border-radius: 0.375rem;
    font-size: 0.875rem;
    background-color: #fff;
}

.options__item--select:focus,
.options__item--input:focus {
    outline: none;
    border-color: #2d2fbf;
    box-shadow: 0 0 0 2px rgba(45, 47, 191, 0.1);
}

.options__actions {
    margin-top: auto;
    gap: 0.5rem;
}

#reset-button {
    background-color: #616161;
    color: white;
    padding: 0.5rem 1rem;
    border: none;
    border-radius: 0.375rem;
    cursor: pointer;
    font-size: 0.875rem;
    transition: background-color 0.2s;
}

#reset-button:hover {
    background-color: #797979;
}

#apply-button {
    flex: 1;
    background-color: #2d2fbf;
    color: white;
    padding: 0.75rem;
    border: none;
    border-radius: 0.375rem;
    cursor: pointer;
    font-size: 1rem;
    font-weight: 500;
    transition: all 0.2s ease;
    display: flex;
    align-items: center;
    justify-content: center;
}

#apply-button:hover:not(:disabled) {
    background-color: #4655ff;
}

#apply-button:disabled {
    background-color: #ccc;
    color: #666;
    cursor: not-allowed;
}

#apply-button.failure {
    background-color: #ff4444;
}

#apply-button.success {
    background-color: #66bb6a;
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
    vertical-align: middle;
}

@keyframes rotation {
    0% {
        transform: rotate(0deg);
    }

    100% {
        transform: rotate(360deg);
    }
}

/* Scrollbar styling */
.plugin-list-container::-webkit-scrollbar,
.visualization-list-container::-webkit-scrollbar {
    width: 6px;
}

.plugin-list-container::-webkit-scrollbar-track,
.visualization-list-container::-webkit-scrollbar-track {
    background: #f1f1f1;
    border-radius: 3px;
}

.plugin-list-container::-webkit-scrollbar-thumb,
.visualization-list-container::-webkit-scrollbar-thumb {
    background: #c1c1c1;
    border-radius: 3px;
}

.plugin-list-container::-webkit-scrollbar-thumb:hover,
.visualization-list-container::-webkit-scrollbar-thumb:hover {
    background: #a8a8a8;
}
</style>