<template>
  <div id="layout">
    <div class="side-layout">
      <div class="setup-layout">
        <div class="setup-title">Plugin List</div>

        <!-- Enhanced Plugin Type Tabs -->
        <div class="plugin-type-tabs">
          <button class="plugin-type-tab" :class="{ active: activePluginType === 'official' }"
            @click="setActivePluginType('official')">
            <span class="tab-text">Official</span>
          </button>
          <button class="plugin-type-tab" :class="{ active: activePluginType === 'local' }"
            @click="setActivePluginType('local')">
            <span class="tab-text">Local</span>
          </button>
        </div>

        <ul class="setup-list">
          <li class="setup-item" v-for="(plugin, idx) in filteredPlugins" :key="idx" @click="selectPlugin(plugin)"
            :class="{ 'selected-plugin': selectedPlugin && selectedPlugin.name === plugin.name }">
            <div class="setup-filename">
              {{ plugin.name }}
              <span class="plugin-source-badge"
                :class="plugin.source === 'official' ? 'badge-official' : 'badge-local'">
                {{ plugin.source.toUpperCase() }}
              </span>
            </div>
          </li>
        </ul>
      </div>
    </div>
    <div class="center-layout">
      <div class="algorithm-layout">
        <div class="algorithm-select">
          <div class="algorithm-logo">{{ selectedPlugin ? selectedPlugin.name : 'Select a Plugin' }}</div>
        </div>
        <div class="algorithm-parts">
          <div v-for="(rule, ruleIndex) in selectedPluginRules" :key="`rule-${ruleIndex}-${rule.name}`">
            <div class="part-title" v-show="rule.parameters.length != 0">{{ rule.name }}</div>
            <div v-for="(parameter, paramIndex) in rule.parameters"
              :key="`param-${ruleIndex}-${paramIndex}-${parameter.name}`" v-show="rule.parameters.length != 0">
              <div class="parameters">
                <span class="parameter-id">
                  {{ parameter.name }}
                </span>
                <input type="number" v-if="parameter.type === 'int' || parameter.type === 'float'"
                  class="parameter__input" v-model="parameter.defaultValue" :step="parameter.defaultValue"
                  :min="parameter.min" :max="parameter.max" />
                <input type="text" v-else-if="parameter.type === 'string'" class="parameter__input"
                  v-model="parameter.defaultValue" :class="{ 'red-text': !parameter.defaultValue }" />
                <input type="checkbox" v-else-if="parameter.type === 'boolean'" class="parameter__input"
                  v-model="parameter.defaultValue">
                <input v-else-if="parameter.type === 'h5adParameter' && parameter.name.includes('UMAP')" type="text"
                  class="parameter__input umap-status-input"
                  :value="parameter.defaultValue && parameter.defaultValue !== '' ? 'Activated' : 'Not Activated'"
                  readonly
                  :class="{ 'umap-activated': parameter.defaultValue && parameter.defaultValue !== '', 'umap-deactivated': !parameter.defaultValue || parameter.defaultValue === '' }" />

                <div v-else-if="parameter.type === 'h5adParameter'" class="parameter__input">
                  <select v-if="parameter.name === 'cell group'" class="parameter__dropdown"
                    v-model="parameter.defaultValue" @change="selectColumns($event)">
                    <option class="parameter__menu" disabled value="">
                      Select Cell Group
                    </option>
                    <option v-for="(column, index) in annotations" :key="index" :value="column" class="parameter__menu">
                      {{ column }}
                    </option>
                  </select>
                  <select v-else-if="parameter.name === 'pseudotime'" class="parameter__dropdown"
                    v-model="parameter.defaultValue">
                    <option class="parameter__menu" disabled value="">
                      Select Pseudotime
                    </option>
                    <option class="parameter__menu" v-for="(column, index) in pseudotime" :key="index" :value="column">
                      {{ column }}
                    </option>
                  </select>
                  <div v-else-if="parameter.name === 'clusters'" class="parameter__dropdown--checkbox"
                    @click="activateClusters" :class="{ isactive: dropdownIsActive }">
                    Select Clusters
                    <ul class="parameter__dropdown--menu">
                      <li v-for="(column, index) in clusters" :key="index">
                        <label>
                          <input type="checkbox" :name="column" @change="clusterToggle($event, parameter, column)"
                            :checked="parameter.defaultValue && parameter.defaultValue.includes(column)" />{{
                              column }}</label>
                      </li>
                    </ul>
                  </div>
                </div>
              </div>
            </div>
          </div>
          <div class="algorithm-alert" v-show="allParametersEmpty">
            No parameters in this plugin.
          </div>
        </div>
      </div>
    </div>
    <div class="side-layout">
      <div class="setup-layout">
        <div class="plugin-io-section">
          <div class="setup-title">Plugin Inputs</div>
          <div class="plugin-inputs-container">
            <div class="plugin-input-item" v-for="(item, idx) in selectedPluginInputOutput" :key="idx"
              v-show="item.type === 'inputFile' || item.type === 'optionalInputFile'">
              <label class="plugin-input-label">
                {{ item.name }}
                <span v-if="item.type === 'optionalInputFile'" class="optional-badge">(Optional)</span>
              </label>
              <select class="plugin-input-dropdown" v-model="item.selected_file" @change="updateInputFile(item)">
                <option value="">{{ item.type === 'optionalInputFile' ? 'Not Selected' : 'Select File' }}</option>
                <option v-for="file in getAvailableFilesForInput(item)" :key="file.node_id + '_' + file.file_name"
                  :value="file">
                  {{ file.display_name }}
                </option>
              </select>
            </div>
          </div>
        </div>
        <div class="plugin-io-section">
          <div class="setup-title">Plugin Outputs</div>
          <div class="plugin-outputs-container">
            <div class="output-file-item" v-for="(item, idx) in selectedPluginInputOutput" :key="idx"
              v-show="item.type === 'outputFile'">
              <div class="output-file-name">{{ item.defaultValue || item.name }}</div>
              <div class="output-file-bottom">
                <div class="output-file-icon">
                  <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none"
                    stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                    <path d="M13 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V9z"></path>
                    <polyline points="13 2 13 9 20 9"></polyline>
                  </svg>
                </div>
                <div class="output-file-type">{{ getFileTypeDescription(item.fileExtension) }}</div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import {
  getPlugins,
  getColumns,
  getClusters,
} from "@/api/index";

export default {
  data() {
    return {
      annotations: [],
      pseudotime: [],
      clusters: [],
      plugins: [],
      workflowId: this.$route.query.workflow_id,
      nodeId: this.$route.query.node,
      selectedPlugin: null,
      selectedPluginRules: [],
      selectedPluginInputOutput: [],
      currentNodeConnection: [],
      dropdownIsActive: false,
      nodeInfo: {},
      currentCellGroup: '',
      activePluginType: 'official', // 기본적으로 Official 플러그인 표시
    };
  },
  async mounted() {
    try {
      this.checkCurrentNodeConnection();

      const plugins = await getPlugins();
      // Filter to only show analysis algorithms (plugin_type = 'ANALYSIS')
      this.plugins = plugins.data.plugins.filter(plugin => {
        return plugin.plugin_type === 'ANALYSIS' || plugin.plugin_type === 'analysis' ||
          (!plugin.plugin_type && this.isAnalysisPlugin(plugin));
      });

      const nodeInfo = this.$store.getters.getWorkflowNodeInfo(this.nodeId);

      if (nodeInfo && nodeInfo.data && nodeInfo.data["selectedPlugin"]) {
        this.selectedPlugin = this.plugins.find(plugin => plugin.name === nodeInfo.data.selectedPlugin.name);
      } else if (this.nodeInfo && this.nodeInfo.data && this.nodeInfo.data.title) {
        this.plugins.forEach((plugin) => {
          if (plugin.name === this.nodeInfo.data.title) {
            this.selectedPlugin = plugin;
          }
        });
      }

      // Ensure we have a selected plugin before proceeding
      if (!this.selectedPlugin && this.plugins.length > 0) {
        this.selectedPlugin = this.plugins[0]; // Default to first analysis plugin
      }

      if (this.selectedPlugin) {
        const selectedPlugin = this.plugins.find((plugin) => plugin.name === this.selectedPlugin.name);
        if (selectedPlugin && selectedPlugin.rules) {
          const result = this.filterRules(selectedPlugin.rules);

          if (nodeInfo && nodeInfo.data && nodeInfo.data["selectedPluginRules"]) {
            this.selectedPluginRules = nodeInfo.data.selectedPluginRules;
          } else {
            this.selectedPluginRules = result.filteredRules;
          }

          // 저장된 상태가 있으면 복원, 없으면 새로 초기화
          if (nodeInfo && nodeInfo.data && nodeInfo.data["selectedPluginInputOutput"]) {
            // 저장된 상태를 복원하면서 현재 연결 상태와 동기화
            this.selectedPluginInputOutput = nodeInfo.data.selectedPluginInputOutput.map(savedItem => {
              const updatedItem = this.activatePlugin([savedItem], this.currentNodeConnection)[0];
              // 저장된 selected_file이 있고 여전히 유효한 경우 유지
              if (savedItem.selected_file && this.currentNodeConnection.some(conn =>
                conn.id === savedItem.selected_file.node_id)) {
                updatedItem.selected_file = savedItem.selected_file;
                updatedItem.activate = true;
                updatedItem.file_name = savedItem.selected_file.file_name;
              }
              return updatedItem;
            });
          } else {
            this.selectedPluginInputOutput = this.activatePlugin(result.filteredInputOutput, this.currentNodeConnection);
          }

          if (this.selectedPluginInputOutput.some((item) => item.type === "inputFile" && item.fileExtension === ".h5ad")) {
            await this.loadColumns();

            // cell group 파라미터를 찾아서 defaultValue가 있으면 clusters를 할당
            const cellGroupParam = this.selectedPluginRules.find(rule =>
              rule.parameters.some(param => param.name === 'cell group' && param.defaultValue)
            )?.parameters.find(param => param.name === 'cell group');

            if (cellGroupParam?.defaultValue) {
              const clusters = await this.getCurrnetClusters(cellGroupParam.defaultValue);
              this.clusters = clusters;
            }
          }
        }
      }
      // Update store with selected plugin information
      if (this.selectedPlugin) {
        const dataObject = {
          "selectedPlugin": {
            name: this.selectedPlugin.name,
            source: this.selectedPlugin.source,
          },
        };
        const nodeId = this.nodeId;
        this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
      }

      if (this.selectedPluginRules) {
        const dataObject = {
          "selectedPluginRules": this.selectedPluginRules,
        };
        const nodeId = this.nodeId;
        this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
      }

      if (this.selectedPluginInputOutput) {
        const dataObject = {
          "selectedPluginInputOutput": this.selectedPluginInputOutput,
        };
        const nodeId = this.nodeId;
        this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
      }
    } catch (error) {
      console.error(error);
    }
  },
  watch: {
    selectedPlugin: {
      handler(newVal) {
        if (newVal) {
          const dataObject = {
            "selectedPlugin": {
              name: newVal.name,
              source: newVal.source,
            },
          };
          const nodeId = this.nodeId;
          this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
        }
      },
      deep: true,
    },
    selectedPluginRules: {
      handler(newVal) {
        if (newVal) {
          // ScatterPlot 파라미터 처리
          newVal.forEach(rule => {
            rule.parameters.forEach(param => {
              if ((param.type === 'h5adParameter' || param.type === 'string') && (param.name === 'ScatterPlot' || param.name.includes('UMAP'))) {
                const scatterPlotNode = this.currentNodeConnection.find(
                  node => node.class === 'ScatterPlot' && node.data.lasso_file_path
                );
                if (scatterPlotNode) {
                  param.defaultValue = scatterPlotNode.data.lasso_file_path;
                }
              }
            });
          });

          // 기존 store commit 로직
          const dataObject = {
            "selectedPluginRules": newVal,
          };
          const nodeId = this.nodeId;
          this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
        }
      },
      deep: true,
    },
    selectedPluginInputOutput: {
      handler(newVal) {
        if (newVal) {
          const dataObject = {
            "selectedPluginInputOutput": newVal,
          };
          const nodeId = this.nodeId;
          this.$store.commit("setWorkflowNodeDataObject", { nodeId, dataObject });
        }
      },
      deep: true,
    },
    annotations: {
      handler() {
        // annotations가 변경되었을 때 clusters 초기화
        this.clusters = [];
      },
      deep: true
    },
    pseudotime: {
      handler() {
        // pseudotime이 변경되었을 때 clusters 초기화
        this.clusters = [];
      },
      deep: true
    },
  },
  computed: {
    allParametersEmpty() {
      // selectedPluginRules를 순회하면서 모든 parameters가 비어 있는지 확인
      return this.selectedPluginRules.every(rule => rule.parameters.length === 0);
    },
    filteredPlugins() {
      return this.plugins.filter(plugin => {
        return plugin.source === this.activePluginType;
      });
    },
    availableFilesByExtension() {
      // 파일 확장자별로 사용 가능한 파일들을 미리 매핑
      const fileMap = {};

      this.currentNodeConnection.forEach(connection => {
        if (connection.data && connection.data.file) {
          // 여러 확장자를 처리할 수 있도록 함
          const extensions = ['.h5ad', '.csv', '.tsv', '.txt', '.json'];

          extensions.forEach(ext => {
            if (connection.data.file.includes(ext)) {
              if (!fileMap[ext]) {
                fileMap[ext] = [];
              }

              const fileInfo = {
                node_id: connection.id,
                file_name: connection.data.file,
                display_name: connection.data.title || connection.data.file,
                node_class: connection.class
              };

              // 중복 방지
              if (!fileMap[ext].some(f => f.node_id === connection.id)) {
                fileMap[ext].push(fileInfo);
              }
            }
          });
        }
      });

      return fileMap;
    }
  },
  methods: {
    getAvailableFilesForInput(inputItem) {
      // 성능 최적화: computed property에서 미리 계산된 파일 맵을 사용
      if (!inputItem.fileExtension) {
        return [];
      }

      // 해당 확장자에 맞는 파일들을 반환
      const files = this.availableFilesByExtension[inputItem.fileExtension] || [];

      // InputFile, DataTable, ScatterPlot 노드만 필터링
      return files.filter(file =>
        file.node_class === 'InputFile' ||
        file.node_class === 'DataTable' ||
        file.node_class === 'ScatterPlot'
      );
    },
    updateInputFile(inputItem) {
      // 드롭다운에서 파일이 선택되었을 때 호출
      if (inputItem.selected_file) {
        inputItem.activate = true;
        inputItem.file_name = inputItem.selected_file.file_name;
      } else {
        // 선택 해제된 경우
        inputItem.activate = false;
        inputItem.file_name = null;
      }

      // Vue의 반응성을 위해 배열 항목을 교체
      const index = this.selectedPluginInputOutput.findIndex(item =>
        item.name === inputItem.name && item.type === inputItem.type
      );
      if (index !== -1) {
        this.selectedPluginInputOutput.splice(index, 1, inputItem);
      }
    },
    getFileTypeDescription(extension) {
      // 파일 확장자에 따른 설명 반환
      const fileTypes = {
        '.h5ad': 'H5AD Data File',
        '.csv': 'CSV Data File',
        '.tsv': 'TSV Data File',
        '.txt': 'Text File',
        '.json': 'JSON File',
        '.pdf': 'PDF Document',
        '.png': 'PNG Image',
        '.jpg': 'JPEG Image',
        '.svg': 'SVG Image'
      };

      return fileTypes[extension] || `${extension} File`;
    },
    isAnalysisPlugin(plugin) {
      // Helper method to determine if a plugin is an analysis algorithm
      // Fallback logic for plugins without explicit plugin_type
      const name = plugin.name.toLowerCase();
      const visualizationKeywords = ['plot', 'chart', 'graph', 'visual', 'heatmap', 'scatter'];
      return !visualizationKeywords.some(keyword => name.includes(keyword));
    },
    setActivePluginType(type) {
      this.activePluginType = type;
      // 탭이 변경되면 선택된 플러그인 초기화
      this.selectedPlugin = null;
      this.selectedPluginRules = [];
      this.selectedPluginInputOutput = [];

      // 변경된 탭의 첫 번째 플러그인을 자동 선택
      if (this.filteredPlugins.length > 0) {
        this.selectPlugin(this.filteredPlugins[0]);
      }
    },
    activateClusters() {
      if (this.clusters.length === 0) {
        alert("There is no cluster column. Please select the cellgroup first.");
        return;
      }
      this.dropdownIsActive = !this.dropdownIsActive;
    },
    async selectColumns(event) {
      const anno_column = event.target.value;
      this.currentCellGroup = anno_column;

      const selectedInputFile = this.selectedPluginInputOutput.find(
        (item) => item.type === "inputFile" && item.fileExtension === ".h5ad" && item.activate
      );
      if (selectedInputFile) {
        const clusters = await this.getCurrnetClusters(anno_column);
        this.clusters = clusters;
      }
    },
    async getCurrnetClusters(anno_column) {
      // H5AD 파일을 loadColumns()와 동일한 방식으로 찾기
      const selectedInputFile = this.selectedPluginInputOutput.find(
        (item) => item.type === "inputFile" && item.fileExtension === ".h5ad" && (item.activate || item.selected_file)
      );
      
      let h5adFileName;
      if (selectedInputFile) {
        // 1. selected_file에서 파일명 가져오기
        if (selectedInputFile.selected_file && selectedInputFile.selected_file.file_name) {
          h5adFileName = selectedInputFile.selected_file.file_name;
        } else if (selectedInputFile.file_name) {
          h5adFileName = selectedInputFile.file_name;
        }
      }
      
      // 2. currentNodeConnection에서 h5ad 파일 찾기
      if (!h5adFileName) {
        const h5adConnection = this.currentNodeConnection.find(conn => 
          conn.data && conn.data.file && conn.data.file.includes('.h5ad')
        );
        if (h5adConnection) {
          h5adFileName = h5adConnection.data.file;
        }
      }
      
      // 3. Algorithm 노드의 파일 정보에서 찾기 (fallback)
      if (!h5adFileName) {
        const nodeFilesInfo = this.$store.getters.getWorkflowNodeFilesInfo(this.nodeId);
        if (nodeFilesInfo) {
          if (Array.isArray(nodeFilesInfo)) {
            const h5adFile = nodeFilesInfo.find(file => file.name && file.name.includes('.h5ad'));
            h5adFileName = h5adFile?.name;
          } else if (typeof nodeFilesInfo === 'object') {
            h5adFileName = Object.values(nodeFilesInfo).find(fileName => 
              fileName && fileName.includes('.h5ad')
            );
          }
        }
      }

      if (!h5adFileName) {
        console.warn('No H5AD file found for getCurrnetClusters');
        return [];
      }

      try {
        const result = await getClusters({
          file_name: h5adFileName,
          anno_column: anno_column,
        });
        console.log('Clusters loaded:', result.data);
        return result.data.clusters;
      } catch (error) {
        console.error('Error loading clusters:', error);
        return [];
      }
    },
    activatePlugin(selectedPluginInputOutput, currentNodeConnection) {
      // ResultFile인 항목의 개수를 계산
      const resultFileCount = currentNodeConnection.filter(item => item.name === 'ResultFile').length;

      // ResultFile 항목을 처리하기 위한 인덱스
      let resultFileIndex = 0;

      return selectedPluginInputOutput.map(item => {
        let activate = false;
        let file_name = null; // file_name 초기값 설정
        let selected_file = null; // 드롭다운에서 선택된 파일

        if (item.type === 'inputFile' || item.type === 'optionalInputFile') {
          // inputFile 타입의 경우
          const matchingConnection = currentNodeConnection.find(connection => {
            if (!connection.data || !connection.data.file) return false;
            
            // 1. 확장자 기반 매칭 (우선순위)
            if (item.fileExtension && connection.data.file.includes(item.fileExtension)) {
              return true;
            }
            
            // 2. InputFile의 title 기반 매칭 (기존 로직)
            if (connection.class === 'InputFile' && connection.data.title && item.defaultValue) {
              return connection.data.title.includes(item.defaultValue) || 
                     item.defaultValue.includes(connection.data.title) ||
                     connection.data.file.includes(item.defaultValue);
            }
            
            // 3. DataTable, ScatterPlot의 확장자 매칭
            if ((connection.class === 'DataTable' || connection.class === 'ScatterPlot') && 
                item.fileExtension && connection.data.file.includes(item.fileExtension)) {
              return true;
            }
            
            return false;
          });

          if (matchingConnection) {
            activate = true; // 파일이 일치하면 활성화
            file_name = matchingConnection.data.file; // 해당 connection의 data.file을 file_name으로 할당
            // 드롭다운 선택을 위한 객체 생성
            selected_file = {
              node_id: matchingConnection.id,
              file_name: matchingConnection.data.file,
              display_name: matchingConnection.data.title || matchingConnection.data.file,
              node_class: matchingConnection.class
            };
          }
        } else if (item.type === 'outputFile') {
          // outputFile 타입의 경우
          if (resultFileIndex < resultFileCount) {
            activate = true;
            resultFileIndex++;
          }
        }

        return {
          ...item,
          activate,
          file_name,
          selected_file
        };
      });
    },
    checkCurrentNodeConnection() {
      this.nodeInfo = this.$store.getters.getWorkflowNodeInfo(this.nodeId);
      console.log(this.nodeInfo);

      // nodeInfo.inputs.input_1.connections, nodeInfo.outputs.output_1.connections 배열을 순회
      // 순회하며 각 아이템의 node를 기반으로 노드 정보 조회
      // 조회한 노드 정보를 currentNodeConnection 배열에 추가
      this.nodeInfo.inputs.input_1.connections.forEach((item) => {
        const node = this.$store.getters.getWorkflowNodeInfo(item.node);
        this.currentNodeConnection.push(node);
      });
      this.nodeInfo.outputs.output_1.connections.forEach((item) => {
        const node = this.$store.getters.getWorkflowNodeInfo(item.node);
        this.currentNodeConnection.push(node);
      });

      console.log(this.currentNodeConnection);
    },
    filterUniqueNames(inputOutputs) {
      const nameCounts = inputOutputs.reduce((acc, param) => {
        acc[param.name] = (acc[param.name] || 0) + 1;
        return acc;
      }, {});

      return inputOutputs.filter(param => nameCounts[param.name] === 1);
    },
    selectPlugin(plugin) {
      if (confirm("Are you sure you want to select this plugin?")) {
        this.selectedPlugin = plugin;
        const result = this.filterRules(this.selectedPlugin.rules);

        this.$nextTick(async () => {
          this.selectedPluginRules = result.filteredRules;
          this.selectedPluginInputOutput = this.activatePlugin(
            result.filteredInputOutput,
            this.currentNodeConnection
          );

          // h5ad 파일이 있으면 columns 로드
          if (this.selectedPluginInputOutput.some((item) => item.type === "inputFile" && item.fileExtension === ".h5ad")) {
            await this.loadColumns();
          }
        });
      }
    },
    filterRules(rules) {
      const filteredInputOutput = new Set();

      const filteredRules = rules
        .filter(rule => !rule.isVisualization)
        .map(rule => {
          const filteredParameters = rule.parameters.filter(param => {
            if (param.type === 'inputFile' || param.type === 'optionalInputFile' || param.type === 'outputFile') {
              filteredInputOutput.add(JSON.stringify(param));
              return false;
            }
            return true;
          });

          return {
            name: rule.name,
            parameters: filteredParameters
          };
        });

      return {
        filteredRules: filteredRules,
        filteredInputOutput: this.filterUniqueNames(Array.from(filteredInputOutput).map(item => JSON.parse(item)))
      };
    },
    clusterToggle(event, parameter, column) {
      if (typeof parameter.defaultValue !== "object") {
        parameter.defaultValue = [];
      }

      if (event.target.checked) {
        // 중복 체크 후 추가
        if (!parameter.defaultValue.includes(column)) {
          parameter.defaultValue.push(column);
        }
      } else {
        parameter.defaultValue = parameter.defaultValue.filter((item) => item !== column);
      }

      console.log(parameter.defaultValue, typeof parameter.defaultValue);
    },
    selectOption(selectedKey, values) {
      console.log(selectedKey, values);
      Object.keys(values).forEach((key) => {
        values[key] = key === selectedKey;
      });
    },
    async loadColumns() {
      // selectedPluginInputOutput 배열을 순회하며 type이 'inputFile'이고 h5ad 파일을 찾음
      const selectedInputFile = this.selectedPluginInputOutput.find(
        (item) => item.type === "inputFile" && item.fileExtension === ".h5ad" && (item.activate || item.selected_file)
      );
      
      if (selectedInputFile) {
        console.log('Selected H5AD file:', selectedInputFile);
        
        // 1. 먼저 selected_file에서 파일명 가져오기 시도
        let h5adFileName;
        if (selectedInputFile.selected_file && selectedInputFile.selected_file.file_name) {
          h5adFileName = selectedInputFile.selected_file.file_name;
        } else if (selectedInputFile.file_name) {
          h5adFileName = selectedInputFile.file_name;
        } else {
          // 2. currentNodeConnection에서 h5ad 파일 찾기
          const h5adConnection = this.currentNodeConnection.find(conn => 
            conn.data && conn.data.file && conn.data.file.includes('.h5ad')
          );
          if (h5adConnection) {
            h5adFileName = h5adConnection.data.file;
          }
        }
        
        // 3. 마지막으로 Algorithm 노드의 파일 정보에서 찾기 (fallback)
        if (!h5adFileName) {
          const nodeFilesInfo = this.$store.getters.getWorkflowNodeFilesInfo(this.nodeId);
          if (nodeFilesInfo) {
            if (Array.isArray(nodeFilesInfo)) {
              const h5adFile = nodeFilesInfo.find(file => file.name && file.name.includes('.h5ad'));
              h5adFileName = h5adFile?.name;
            } else if (typeof nodeFilesInfo === 'object') {
              h5adFileName = Object.values(nodeFilesInfo).find(fileName => 
                fileName && fileName.includes('.h5ad')
              );
            }
          }
        }

        console.log('H5AD filename resolved:', h5adFileName);

        if (!h5adFileName) {
          console.warn('No H5AD file found in any source');
          alert("No H5AD file found. Please ensure an H5AD file is connected to this algorithm.");
          return;
        }

        try {
          const result = await getColumns({
            file_name: h5adFileName,
          });
          console.log('Columns loaded:', result.data);
          this.annotations = result.data.anno_columns;
          this.pseudotime = result.data.pseudo_columns;
        } catch (error) {
          console.error('Error loading columns:', error);
          alert(`Error loading H5AD file columns: ${error.message || error}`);
        }
      } else {
        // H5AD 파일이 필요하지 않은 플러그인일 수도 있으므로 경고만 출력
        console.warn("No H5AD input file found in plugin configuration.");
      }
    },
  },
};
</script>

<style scoped>
#layout {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: row;
}

/* Plugin Type Tabs */
.plugin-type-tabs {
  display: flex;
  margin-bottom: 1rem;
  border-radius: 8px;
  overflow: hidden;
  /* box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1); */
  gap: 10px;
  padding: 0 10px;
}

.plugin-type-tab {
  flex: 1;
  padding: 10px 20px;
  border: none;
  background-color: #34495e;
  color: #ecf0f1;
  font-size: 14px;
  font-weight: 500;
  cursor: pointer;
  transition: all 0.3s ease;
  border-radius: 4px;
  display: flex;
  align-items: center;
  justify-content: center;
}


.plugin-type-tab:hover {
  background-color: #3498db;
  color: white;
}

.plugin-type-tab.active {
  background-color: #2980b9;
  color: white;
  font-weight: 600;
}


/* Plugin Source Badge */
.plugin-source-badge {
  font-size: 10px;
  padding: 2px 6px;
  border-radius: 12px;
  margin-left: 8px;
  font-weight: 600;
  letter-spacing: 0.5px;
}

.badge-official {
  background-color: #3498db;
  color: white;
}

.badge-local {
  background-color: #e67e22;
  color: white;
}

.setup-layout,
.setup-layout__node,
.algorithm-layout {
  width: 95%;
  min-height: 100%;
  margin: auto;
  display: flex;
  align-items: center;
  /* justify-content: center; */
  flex-direction: column;
  padding: 1rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
  /* overflow-y: auto; */
  scrollbar-width: thin;
  scrollbar-color: #888 #f5f5f5;
}

.setup-layout__node {
  min-height: 49%;
  margin: 0;
}

/* Plugin I/O Section Styles */
.plugin-io-section {
  flex: 1;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
  min-height: 48%;
  margin-bottom: 1rem;
  overflow: hidden;
  width: 100%;
}

.plugin-io-section:first-child {
  border-bottom: 1px solid #e7e7e7;
}

.plugin-io-section:last-child {
  margin-bottom: 0;
}

.output-layout {
  width: 25%;
  height: 95%;
  margin-left: 1rem;
  display: flex;
  align-items: center;
  /* justify-content: center; */
  flex-direction: column;
  padding-left: 1rem;
  padding-right: 1rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
}

.setup-title {
  text-transform: capitalize;
  font-weight: bold;
  font-size: larger;
  color: #494949;
  margin: 0.5rem 0;
}

.input-description,
.output-description {
  width: 100%;
  height: 5rem;
  text-align: center;
  justify-content: center;
  display: flex;
  flex-direction: column;
  color: #a1a1a1;
  background-color: rgb(224, 224, 224);
  border-radius: 1rem;
  margin: 0.5rem 0rem;
  border: 2px solid #e7e7e7;
  opacity: 0.9;
}

.output-description {
  cursor: not-allowed;
}

.input-description.linked,
.output-description.generated {
  opacity: 1;
  color: #fff;
  background-color: rgb(22, 152, 51);
  border: 2px solid #ecebff;
}

.output-description.generated {
  cursor: pointer;
}

.description-id {
  position: relative;
  display: inline-block;
  margin: 10px;
}

.description-id .description-tooltip {
  visibility: hidden;
  width: 200px;
  background-color: #555;
  color: #fff;
  text-align: center;
  padding: 15px 10px;
  border-radius: 10px;

  position: absolute;
  z-index: 1;
  bottom: 125%;
  left: 50%;
  margin-left: -60px;

  opacity: 0;
  transition: opacity 0.3s;
}

.description-id .description-tooltip::after {
  content: "";
  position: absolute;
  top: 100%;
  left: 30%;
  margin-left: -5px;
  border-width: 5px;
  border-style: solid;
  border-color: #555 transparent transparent transparent;
}

/* Show the tooltip text when you mouse over the tooltip container */
.description-id:hover .description-tooltip {
  visibility: visible;
  opacity: 1;
}

.setup-list {
  width: 100%;
  height: 100%;
  overflow-y: auto;
  scrollbar-width: thin;
  scrollbar-color: #888 #f5f5f5;
}

.setup-list__node {
  width: 100%;
  display: flex;
  justify-content: center;
  flex-wrap: wrap;
  /* 너비 초과 시 다음 줄로 넘김 */
  gap: 10px;
  /* 항목 간의 간격을 설정할 수 있습니다 */
}

.setup-item {
  width: 100%;
  display: flex;
  align-items: center;
  cursor: pointer;
  margin-bottom: 0.5rem;
}

.setup-item__node {
  width: 5rem;
  height: 5rem;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  text-align: center;
  background-color: rgb(224, 224, 224);
  border-radius: 1rem;
  position: relative;
}

.setup-item:hover {
  background-color: rgb(202, 214, 255);
}

.setup-item.selected-plugin>.setup-filename {
  font-weight: bold;
  font-size: 1rem;
  opacity: 1;
  background-color: rgb(224, 224, 224);
  border-radius: 1rem;
}

.setup-date,
.setup-filename {
  width: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
  font-size: 0.9rem;
  color: #353535;
  padding: 0.5rem;
}

.setup-filename__node {
  font-size: 0.8rem;
  color: #353535;
  max-width: 95%;
  height: 0.9rem;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.algorithm-select {
  width: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}

.algorithm-logo {
  font-size: 4rem;
  font-weight: bold;
  color: #353535;
}

.algorithm-parts {
  position: relative;
  align-content: start;
  width: 100%;
  margin-bottom: 2rem;
}

.algorithm-select__tenet {
  width: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 1rem;
  box-sizing: border-box;
}

.option {
  background: #fff;
  height: 100%;
  width: 100%;
  display: flex;
  align-items: center;
  justify-content: space-evenly;
  margin: 0 10px;
  border-radius: 5px;
  cursor: pointer;
  padding: 0 10px;
  border: 2px solid lightgrey;
  transition: all 0.3s ease;
}

.dot {
  height: 20px;
  width: 20px;
  background: #d9d9d9;
  border-radius: 50%;
  position: relative;
}

.dot::before {
  position: absolute;
  content: "";
  top: 4px;
  left: 4px;
  width: 12px;
  height: 12px;
  background: #0069d9;
  border-radius: 50%;
  opacity: 0;
  transform: scale(1.5);
  transition: all 0.3s ease;
}

#option-1:checked:checked~.option-1,
#option-2:checked:checked~.option-2 {
  border-color: #0069d9;
  background: #0069d9;
}

#option-1:checked:checked~.option-1 .dot,
#option-2:checked:checked~.option-2 .dot {
  background: #fff;
}

#option-1:checked:checked~.option-1 .dot::before,
#option-2:checked:checked~.option-2 .dot::before {
  opacity: 1;
  transform: scale(1);
}

.wrapper .option span {
  font-size: 20px;
  color: #808080;
}

#option-1:checked:checked~.option-1 span,
#option-2:checked:checked~.option-2 span {
  color: #fff;
}

.part-title {
  text-transform: capitalize;
  font-weight: bold;
  font-size: large;
  color: #353535;

  margin: 20px 0px 10px 0px;
}

.parameters {
  display: flex;
  direction: row;
  justify-content: space-between;
  width: 100%;
  height: 2rem;
  color: #353535;
  position: relative;
  margin: 2px;
  margin-bottom: 0.5rem;
  text-transform: capitalize;
}

.parameter-id {
  position: relative;
  display: flex;
  align-items: center;
}

.parameter-id .parameter-tooltip {
  visibility: hidden;
  width: 200px;
  background-color: #555;
  color: #fff;
  text-align: center;
  padding: 15px 10px;
  border-radius: 10px;

  position: absolute;
  z-index: 1;
  bottom: 125%;
  left: 50%;
  margin-left: -60px;

  opacity: 0;
  transition: opacity 0.3s;

  text-transform: capitalize;
}

.parameter-id .parameter-tooltip::after {
  content: "";
  position: absolute;
  top: 100%;
  left: 30%;
  margin-left: -5px;
  border-width: 5px;
  border-style: solid;
  border-color: #555 transparent;
}

/* Show the tooltip text when you mouse over the tooltip container */
.parameter-id:hover .parameter-tooltip {
  visibility: visible;
  opacity: 1;
}

.parameter__dropdown {
  width: 100%;
  height: 100%;
  color: black;
  border: 1px solid #aaa;
  border-radius: 3px;
  font-size: small;
  text-align: left;
  margin-bottom: 0px;
  display: flex;
  align-items: center;
  justify-content: left;
  position: relative;
}

.parameter__dropdown:focus {
  outline: none;
  border-color: #aaa;
}

.parameter__dropdown span {
  color: black;
  font-size: 13px;
  font-weight: 400;
  cursor: pointer;
  display: block;
  font-family: Arial, Helvetica, sans-serif;
  margin-left: 0.8rem;
}

.parameter__dropdown--checkbox {
  width: 100%;
  height: 100%;
  border: 1px solid #aaa;
  border-radius: 3px;
  position: relative;
  margin: 0;
  display: flex;
  align-items: center;
  padding-left: 5px;
  box-sizing: border-box;

  user-select: none;
}

.parameter__dropdown--menu {
  list-style: none;
  margin: 0;
  padding: 0;
  position: absolute;
  top: 100%;
  /* align the dropdown right below the dropdown text */
  border: inherit;
  border-top: none;
  left: -1px;
  /* align the dropdown to the left */
  right: -1px;
  /* align the dropdown to the right */
  opacity: 0;
  /* hide the dropdown */

  transition: opacity 0.4s ease-in-out;
  /* height: 100px; */
  overflow: scroll;
  overflow-x: hidden;
  pointer-events: none;
  max-height: 15rem;

  z-index: 9999;
}

.parameter__dropdown--menu li label {
  display: block;
  border-bottom: 1px solid silver;
  padding: 10px;
  background-color: white;

  transition: all 0.1s ease-out;
}

.parameter__dropdown--menu li label:first-child {
  border-top-left-radius: 3px;
  border-top-right-radius: 3px;
  border-top: 1px solid silver;
}

.parameter__dropdown--menu li label:hover {
  background-color: #555;
  color: white;
}

.isactive .parameter__dropdown--menu {
  opacity: 1;
  pointer-events: auto;
}

.parameter__dropdown--checkbox:after {
  content: '';
  height: 0;
  position: absolute;
  width: 0;
  border: 6px solid transparent;
  border-top-color: #000;
  top: 50%;
  right: 3px;
  margin-top: -3px;
}

.parameter__dropdown--checkbox .isactive:after {
  border-bottom-color: #000;
  border-top-color: #fff;
  margin-top: -9px;
}

.parameter__drop {
  position: absolute;
  opacity: 0;
  height: 0px;
}

.parameter__menu:checked::before {
  content: "✔";
}

.parameter__radioBox {
  width: 9.3rem;
  display: flex;
  align-items: center;
  justify-content: center;
}

.parameter__radio {
  width: 50%;
  border: 1px solid #8b8887;
  padding: 0.5rem 0;
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: pointer;
  font-size: 0.9rem;
}

.parameter__checkbox {
  display: none;
  position: absolute;
  clear: both;
  width: 9.3rem;
  max-height: 300px;
  overflow-y: auto;
  overflow-x: hidden;
  text-align: center;
  background: #6c6867;
  transition: height 0.4s ease;
  border-radius: 0.5rem;
}

/* 웹킷 브라우저(크롬, 사파리 등)에 대한 스크롤바 스타일 */
.parameter__checkbox::-webkit-scrollbar {
  width: 5px;
  /* 스크롤바 너비 */
}

.parameter__checkbox::-webkit-scrollbar-thumb {
  background-color: #cbcbcb;
  /* 스크롤바 색상 */
  border-radius: 1rem;
  /* 스크롤바 border-radius */
}

.parameter__checkbox::-webkit-scrollbar-thumb:hover {
  background-color: #8d8b8b;
  /* 스크롤바 호버 시 색상 */
}

.parameter__checkbox div {
  height: 26px;
  line-height: 1.2rem;
  color: white;
  display: flex;
  align-items: center;
  justify-content: center;
}

#touch:checked+.parameter__checkbox {
  border: 1px solid #8b8887;
  display: block;
  top: calc(100%);
  right: 0;
  z-index: 9999;
}

.parameter__check {
  transform: scale(1.5);
  margin-bottom: 1rem;
  margin-right: 0.5rem;
}

.parameter__check:nth-child(2) {
  margin-left: 0rem;
}

.parameter__icon {
  position: absolute;
  right: 0.6rem;
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  padding-left: 2rem;
}

.parameter__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
  opacity: 0.8;
  margin-right: 0.5rem;
}

.parameter__menu {
  width: 85%;
  height: 100%;
  color: black;
  padding: 3px 7.5%;
  /* right: 10px; */
  border-radius: 3px;
  border-color: #e7eaff;
  font-size: small;
  text-align: center;
  margin-bottom: 0px;
  cursor: pointer;
  position: relative;
  word-wrap: break-word;
  /* 긴 단어가 div의 너비를 넘어갈 때 줄바꿈 */
  overflow-wrap: break-word;
  /* CSS3에서 word-wrap 대신 사용 */
}

.parameter__input {
  width: 8.5rem;
  color: black;
  margin: 0;
  /* right: 10px; */
  border-radius: 3px;
  border-color: #f1f2fc;
  font-size: small;
  text-align: center;
  margin-bottom: 0px;
  box-sizing: border-box;
}

.parameter__input:disabled {
  background-color: lightgray;
}

.parameter__input:disabled::placeholder {
  color: black;
}

.side-layout,
.center-layout {
  height: 95%;
  display: flex;
  align-items: center;
  flex-direction: column;
  justify-content: space-between;
  box-sizing: border-box;
  overflow-y: hidden;
}

.side-layout {
  width: 25%;
}

.center-layout {
  width: 50%;
}

.setup-layout::-webkit-scrollbar,
.setup-layout__node::-webkit-scrollbar,
.algorithm-layout::-webkit-scrollbar {
  width: 5px;
  /* 스크롤바 너비 */
}

.setup-layout::-webkit-scrollbar-thumb,
.setup-layout__node::-webkit-scrollbar-thumb,
.algorithm-layout::-webkit-scrollbar-thumb {
  background-color: #cbcbcb;
  /* 스크롤바 색상 */
  border-radius: 1rem;
  /* 스크롤바 border-radius */
}

.setup-layout::-webkit-scrollbar-thumb:hover,
.setup-layout__node::-webkit-scrollbar-thumb:hover,
.algorithm-layout::-webkit-scrollbar-thumb:hover {
  background-color: #8d8b8b;
  /* 스크롤바 호버 시 색상 */
}

.options-layout {
  width: 25%;
  height: 95%;
  /* background-color: red; */
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
}

.options__item {
  /* margin: auto; */
  color: black;
  /* padding: 2% 0; */
  margin: 2%;
}

.form__input {
  position: absolute;
  width: 0;
  height: 0;
  padding: 0;
  overflow: hidden;
  border: 0;
  display: none;
}

.form__button--setup {
  cursor: pointer;
  width: 8rem;
  height: 3rem;
  margin-left: 12.5rem;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.5);
  border-radius: 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  /* background: rgb(163, 163, 163); */
  background: rgb(59, 108, 221);
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.2rem;
  line-height: 1rem;
  color: rgb(240, 240, 240);
}

.form__button--result {
  cursor: pointer;
  width: 12rem;
  height: 3rem;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.5);
  border-radius: 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  /* background: rgb(163, 163, 163); */
  background: rgb(59, 108, 221);
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.2rem;
  line-height: 1rem;
  color: rgb(249, 249, 249);
  margin-bottom: 1rem;
}

.form__button--result:hover {
  background: rgb(40, 84, 197);
}

.apply {
  cursor: pointer;
  background: rgb(40, 84, 197);
  color: white;
}

.wiggle {
  animation: wiggle 2s linear infinite;
}

.selected {
  /* 선택된 요소의 배경색 변경 */
  background-color: #0069d9;
  color: white;
}

input[type="radio"] {
  display: none;
}

.reset-button {
  position: absolute;
  margin-top: -5px;
  width: 1.5rem;
  height: 1.5rem;
  opacity: 0.7;
}

.reset-button:hover {
  opacity: 0.5;
  cursor: pointer;
}

.algorithm-alert {
  margin-top: 2rem;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-direction: column;
  font-size: 1.5rem;
  color: #353535;
}

.checkbox-wrapper-9 {
  margin-bottom: 0.5rem;
}

.checkbox-wrapper-9 .tgl {
  display: none;
}

.checkbox-wrapper-9 .tgl,
.checkbox-wrapper-9 .tgl:after,
.checkbox-wrapper-9 .tgl:before,
.checkbox-wrapper-9 .tgl *,
.checkbox-wrapper-9 .tgl *:after,
.checkbox-wrapper-9 .tgl *:before,
.checkbox-wrapper-9 .tgl+.tgl-btn {
  box-sizing: border-box;
}

.checkbox-wrapper-9 .tgl::-moz-selection,
.checkbox-wrapper-9 .tgl:after::-moz-selection,
.checkbox-wrapper-9 .tgl:before::-moz-selection,
.checkbox-wrapper-9 .tgl *::-moz-selection,
.checkbox-wrapper-9 .tgl *:after::-moz-selection,
.checkbox-wrapper-9 .tgl *:before::-moz-selection,
.checkbox-wrapper-9 .tgl+.tgl-btn::-moz-selection,
.checkbox-wrapper-9 .tgl::selection,
.checkbox-wrapper-9 .tgl:after::selection,
.checkbox-wrapper-9 .tgl:before::selection,
.checkbox-wrapper-9 .tgl *::selection,
.checkbox-wrapper-9 .tgl *:after::selection,
.checkbox-wrapper-9 .tgl *:before::selection,
.checkbox-wrapper-9 .tgl+.tgl-btn::selection {
  background: none;
}

.checkbox-wrapper-9 .tgl+.tgl-btn {
  outline: 0;
  display: block;
  width: 3.5em;
  height: 1.5em;
  position: relative;
  cursor: pointer;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}

.checkbox-wrapper-9 .tgl+.tgl-btn:after,
.checkbox-wrapper-9 .tgl+.tgl-btn:before {
  position: relative;
  display: block;
  content: "";
  width: 50%;
  height: 100%;
}

.checkbox-wrapper-9 .tgl+.tgl-btn:after {
  left: 0;
}

.checkbox-wrapper-9 .tgl+.tgl-btn:before {
  display: none;
}

.checkbox-wrapper-9 .tgl:checked+.tgl-btn:after {
  left: 50%;
}

.checkbox-wrapper-9 .tgl-flat+.tgl-btn {
  padding: 2px;
  transition: all 0.2s ease;
  background: #fff;
  border: 4px solid #f2f2f2;
  border-radius: 2em;
}

.checkbox-wrapper-9 .tgl-flat+.tgl-btn:after {
  transition: all 0.2s ease;
  background: #f2f2f2;
  content: "";
  border-radius: 1em;
}

.checkbox-wrapper-9 .tgl-flat:checked+.tgl-btn {
  border: 4px solid #7FC6A6;
}

.checkbox-wrapper-9 .tgl-flat:checked+.tgl-btn:after {
  left: 50%;
  background: #7FC6A6;
}

.checked {
  background-color: #a0eac9;
  box-shadow: 0 0 0 2px #a0eac9 inset;
}

/* Plugin Inputs Dropdown Styles */
.plugin-inputs-container {
  width: 100%;
  flex: 1;
  display: flex;
  flex-direction: column;
  gap: 1rem;
  padding: 0.5rem;
  padding-right: 0.75rem;
  /* 스크롤바 공간 확보 */
  overflow-y: auto;
  overflow-x: hidden;
  box-sizing: border-box;
}

.plugin-input-item {
  display: flex;
  flex-direction: column;
  gap: 0.5rem;
  flex-shrink: 0;
  width: 100%;
  box-sizing: border-box;
}

.plugin-input-label {
  font-size: 0.9rem;
  font-weight: 500;
  color: #353535;
  display: flex;
  align-items: center;
  gap: 0.5rem;
}

.optional-badge {
  font-size: 0.75rem;
  color: #888;
  font-weight: normal;
}

.plugin-input-dropdown {
  width: 100%;
  padding: 0.5rem;
  border: 1px solid #ddd;
  border-radius: 0.25rem;
  background-color: white;
  font-size: 0.85rem;
  color: #353535;
  cursor: pointer;
  transition: all 0.2s ease;
  box-sizing: border-box;
  min-width: 0;
}

.plugin-input-dropdown:hover {
  border-color: #3498db;
}

.plugin-input-dropdown:focus {
  outline: none;
  border-color: #2980b9;
  box-shadow: 0 0 0 2px rgba(41, 128, 185, 0.1);
}

/* Plugin Outputs File List Styles */
.plugin-outputs-container {
  width: 100%;
  flex: 1;
  display: flex;
  flex-direction: column;
  gap: 0.75rem;
  padding: 0.5rem;
  padding-right: 0.75rem;
  /* 스크롤바 공간 확보 */
  overflow-y: auto;
  overflow-x: hidden;
  box-sizing: border-box;
}

.output-file-item {
  display: flex;
  flex-direction: column;
  gap: 0.5rem;
  padding: 0.75rem;
  background-color: #f8f9fa;
  border-radius: 0.5rem;
  border: 1px solid #e9ecef;
  transition: all 0.2s ease;
  width: 100%;
  box-sizing: border-box;
  min-width: 0;
}

.output-file-item:hover {
  background-color: #e9ecef;
  border-color: #dee2e6;
}

.output-file-name {
  font-size: 0.9rem;
  font-weight: 500;
  color: #353535;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  width: 100%;
  line-height: 1.3;
  cursor: pointer;
  transition: all 0.2s ease;
  position: relative;
}

.output-file-name:hover {
  white-space: normal;
  word-wrap: break-word;
  overflow-wrap: break-word;
  overflow: visible;
  background-color: rgba(255, 255, 255, 0.95);
  padding: 0.25rem;
  margin: -0.25rem;
  border-radius: 0.25rem;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.15);
  z-index: 10;
  max-height: none;
}

.output-file-bottom {
  display: flex;
  align-items: center;
  gap: 0.5rem;
  width: 100%;
}

.output-file-icon {
  flex-shrink: 0;
  color: #6c757d;
  width: 16px;
  height: 16px;
}

.output-file-type {
  font-size: 0.75rem;
  color: #6c757d;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  flex: 1;
  min-width: 0;
}

/* Scrollbar styling for containers */
.plugin-inputs-container::-webkit-scrollbar,
.plugin-outputs-container::-webkit-scrollbar {
  width: 10px;
  /* 5px에서 10px로 증가 */
}

.plugin-inputs-container::-webkit-scrollbar-track {
  background: transparent;
  /* 투명한 트랙 배경 */
  border-radius: 1rem;
}

.plugin-inputs-container::-webkit-scrollbar-thumb,
.plugin-outputs-container::-webkit-scrollbar-thumb {
  background-color: rgba(203, 203, 203, 0.6);
  /* 반투명 처리 */
  border-radius: 1rem;
  border: 2px solid transparent;
  /* 투명한 테두리 추가 */
  background-clip: content-box;
  /* 테두리 안쪽만 색상 적용 */
}

.plugin-inputs-container::-webkit-scrollbar-thumb:hover,
.plugin-outputs-container::-webkit-scrollbar-thumb:hover {
  background-color: rgba(141, 139, 139, 0.8);
  /* 호버 시 더 진한 색 */
}

/* Firefox 지원 */
.plugin-inputs-container,
.plugin-outputs-container {
  scrollbar-width: thin;
  scrollbar-color: rgba(203, 203, 203, 0.6) transparent;
}

/* UMAP 파라미터 상태 스타일 */
.umap-status-input {
  width: 8.5rem;
  height: auto;
  padding: 0;
  margin: 0;
  border-radius: 3px;
  font-size: small;
  text-align: center;
  box-sizing: border-box;
}

.umap-activated {
  background-color: #d4edda;
  border-color: #c3e6cb;
  color: #155724;
  font-weight: 500;
}

.umap-deactivated {
  background-color: #f8d7da;
  border-color: #f5c6cb;
  color: #721c24;
  font-weight: 500;
}
</style>