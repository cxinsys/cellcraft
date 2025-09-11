<template>
  <div class="dag-modal" v-if="visible" @click.self="closeModal">
    <div class="dag-content" @click.stop>
      <!-- Header -->
      <div class="dag-header">
        <div class="dag-title">
          <h3>Workflow Progress: {{ taskName }}</h3>
          <div class="dag-progress" v-if="progressInfo">
            <div class="progress-bar">
              <div class="progress-fill" :style="{ width: progressInfo.percentage + '%' }"
                :class="{ 'has-failed': progressInfo.failed_rules > 0 }"></div>
            </div>
            <div class="progress-text">
              {{ progressInfo.completed_rules }}/{{ progressInfo.total_rules }} completed
              <span v-if="progressInfo.failed_rules > 0" class="failed-count">
                ({{ progressInfo.failed_rules }} failed)
              </span>
              <span v-if="progressInfo.running_rules > 0" class="running-count">
                ({{ progressInfo.running_rules }} running)
              </span>
            </div>
          </div>
        </div>

        <!-- Controls -->
        <div class="dag-controls">
          <button @click="refreshStatus" class="control-btn" title="Refresh Status">
            <span>Refresh</span>
          </button>
          <button @click="exportPNG" class="control-btn" title="Export as PNG">
            <span>PNG</span>
          </button>
          <button @click="exportSVG" class="control-btn" title="Export as SVG">
            <span>SVG</span>
          </button>
          <button @click="closeModal" class="close-btn" title="Close">
            ×
          </button>
        </div>
      </div>

      <!-- Main Content -->
      <div class="dag-body">
        <!-- Plot Container -->
        <div ref="dagPlot" class="dag-plot" :class="{ 'loading': isLoading }">
          <div v-if="isLoading" class="loading-overlay">
            <div class="spinner"></div>
            <p>Loading DAG...</p>
          </div>
          <div v-if="errorMessage" class="error-message">
            <p>{{ errorMessage }}</p>
            <button @click="retryLoad" class="retry-btn">Retry</button>
          </div>
        </div>

        <!-- Side Panel -->
        <div class="dag-sidebar" v-if="selectedNode">
          <div class="node-details">
            <h4>{{ selectedNode.label }} Details</h4>

            <div class="detail-section">
              <div class="detail-label">Status:</div>
              <div class="detail-content">
                <span class="status-badge" :class="'status-' + (nodeStatuses[selectedNode.id] || 'pending')">
                  {{ nodeStatuses[selectedNode.id] || 'pending' }}
                </span>
              </div>
            </div>

            <div class="detail-section" v-if="selectedNode.description">
              <div class="detail-label">Description:</div>
              <div class="detail-content">{{ selectedNode.description }}</div>
            </div>

            <div class="detail-section" v-if="selectedNode.inputs && selectedNode.inputs.length">
              <div class="detail-label">Input Files:</div>
              <div class="detail-content">
                <div v-for="input in selectedNode.inputs" :key="input" class="file-item">
                  • {{ getFileName(input) }}
                </div>
              </div>
            </div>

            <div class="detail-section" v-if="selectedNode.outputs && selectedNode.outputs.length">
              <div class="detail-label">Output Files:</div>
              <div class="detail-content">
                <div v-for="output in selectedNode.outputs" :key="output" class="file-item">
                  • {{ getFileName(output) }}
                </div>
              </div>
            </div>

            <div class="detail-section" v-if="selectedNode.params && selectedNode.params.length">
              <div class="detail-label">Parameters:</div>
              <div class="detail-content">
                <div v-for="param in selectedNode.params" :key="param" class="param-item">
                  • {{ param }}
                </div>
              </div>
            </div>

            <div class="detail-section" v-if="selectedNode.script">
              <div class="detail-label">Script:</div>
              <div class="detail-content script-content">
                {{ selectedNode.script }}
              </div>
            </div>

            <div class="detail-actions">
              <button @click="showRuleLogs" class="action-btn">
                📋 View Logs
              </button>
              <button @click="selectedNode = null" class="action-btn secondary">
                ✖ Close Details
              </button>
            </div>
          </div>
        </div>
      </div>

      <!-- Status Legend -->
      <div class="dag-footer">
        <div class="status-legend">
          <div class="legend-item">
            <span class="status-dot pending"></span> Pending
          </div>
          <div class="legend-item">
            <span class="status-dot running"></span> Running
          </div>
          <div class="legend-item">
            <span class="status-dot success"></span> Completed
          </div>
          <div class="legend-item">
            <span class="status-dot failed"></span> Failed
          </div>
        </div>
      </div>
    </div>

  </div>
</template>

<script>
import Plotly from 'plotly.js-dist-min';
import HierarchicalLayout from '@/utils/HierarchicalLayout';
import { getDAGStructure, getRuleStatuses, getRuleLogs } from '@/api';

export default {
  name: 'DAGVisualization',
  props: {
    visible: {
      type: Boolean,
      default: false
    },
    taskId: {
      type: String,
      required: false,
      default: null
    },
    taskName: {
      type: String,
      default: 'Unknown Task'
    },
    taskStatus: {
      type: String,
      required: false,
      default: 'UNKNOWN'
    }
  },
  data() {
    return {
      // DAG 데이터
      dagData: null,
      nodeStatuses: {},
      progressInfo: null,

      // 시각화 상태
      isLoading: false,
      errorMessage: '',
      selectedNode: null,

      // 레이아웃 설정 (고정값)
      currentDirection: 'LR',  // Left to Right로 고정
      textPosition: 'outside', // Outside로 고정


      // 색상 매핑
      statusColors: {
        'pending': '#9E9E9E',    // 회색
        'running': '#f39c12',    // 주황색
        'success': '#4CAF50',    // 초록색
        'failed': '#F44336',     // 빨간색
        'skipped': '#FFC107'     // 노란색
      }
    };
  },
  watch: {
    visible(newValue) {
      console.log('DAGVisualization visible changed:', newValue, 'taskId:', this.taskId);
      if (newValue && this.taskId) {
        this.loadDAGData();
      } else {
        this.selectedNode = null;
      }
    },
    taskId(newValue) {
      console.log('DAGVisualization taskId changed:', newValue, 'visible:', this.visible);
      if (this.visible && newValue) {
        this.loadDAGData();
      }
    },
    taskStatus(newValue) {
      console.log('DAGVisualization taskStatus changed:', newValue);
      // Refresh rule statuses when task status changes
      if (this.visible && this.taskId) {
        this.updateRuleStatuses();
      }
    }
  },
  mounted() {
    console.log('DAGVisualization mounted - visible:', this.visible, 'taskId:', this.taskId, 'taskStatus:', this.taskStatus);
    // 컴포넌트가 마운트된 시점에 visible=true이고 taskId가 있으면 데이터 로드
    if (this.visible && this.taskId) {
      this.loadDAGData();
      // 추가로 refreshStatus 함수를 한 번 더 실행
      this.refreshStatus();
    }
  },
  methods: {
    /**
     * DAG 데이터 로드
     */
    async loadDAGData() {
      console.log('loadDAGData called with taskId:', this.taskId);

      if (!this.taskId) {
        console.warn('loadDAGData: taskId is null or undefined');
        this.errorMessage = 'Task ID is required';
        return;
      }

      this.isLoading = true;
      this.errorMessage = '';

      try {
        console.log('Calling getDAGStructure with taskId:', this.taskId);
        const response = await getDAGStructure(this.taskId);
        console.log('getDAGStructure response:', response);

        this.dagData = response.data.dag_structure;
        console.log('DAG structure loaded:', this.dagData);
        console.log('DAG nodes:', this.dagData?.nodes?.map(n => ({ id: n.id, label: n.label })));

        // 상태 정보도 함께 로드
        await this.updateRuleStatuses();

        // DAG 시각화 렌더링
        await this.renderDAG();

        this.isLoading = false;
      } catch (error) {
        console.error('Failed to load DAG data:', error);
        if (error.response) {
          // 서버 응답이 있는 경우
          const statusCode = error.response.status;
          if (statusCode === 404) {
            this.errorMessage = 'Workflow not found. Please check if the task exists.';
          } else if (statusCode === 403) {
            this.errorMessage = 'Access denied. You do not have permission to view this workflow.';
          } else {
            this.errorMessage = `Server error (${statusCode}). Please try again later.`;
          }
        } else {
          // 네트워크 오류 등
          this.errorMessage = 'Failed to load DAG data. Please check your connection and try again.';
        }
        this.isLoading = false;
      }
    },

    /**
     * Rule 상태 업데이트
     */
    async updateRuleStatuses() {
      console.log('updateRuleStatuses called with taskId:', this.taskId);

      if (!this.taskId) {
        console.warn('Cannot update rule statuses: taskId is not available');
        return;
      }

      try {
        console.log('Calling getRuleStatuses with taskId:', this.taskId, 'taskStatus:', this.taskStatus);
        const response = await getRuleStatuses(this.taskId, this.taskStatus);
        console.log('getRuleStatuses response:', response);
        console.log('Raw rule_statuses from API:', response.data.rule_statuses);
        console.log('Raw rule_statuses keys:', Object.keys(response.data.rule_statuses || {}));

        this.nodeStatuses = response.data.rule_statuses;
        this.progressInfo = response.data.progress;

        console.log('nodeStatuses after assignment:', this.nodeStatuses);
        console.log('nodeStatuses keys:', Object.keys(this.nodeStatuses || {}));
        console.log('nodeStatuses type:', typeof this.nodeStatuses);
        console.log('nodeStatuses is array:', Array.isArray(this.nodeStatuses));

        // 이미 렌더링된 경우 색상 업데이트
        if (this.dagData && this.$refs.dagPlot && this.$refs.dagPlot.data) {
          this.updateNodeColors();
        }
      } catch (error) {
        console.error('Failed to update rule statuses:', error);
        // 사용자에게 알리지 않고 조용히 실패 - 폴링 중에는 너무 많은 오류 메시지가 방해가 될 수 있음
      }
    },

    /**
     * DAG 시각화 렌더링
     */
    async renderDAG() {
      if (!this.dagData || !this.$refs.dagPlot) return;

      const { nodes, edges } = this.dagData;

      // 레이아웃 계산
      const layout = new HierarchicalLayout(nodes, edges, {
        direction: this.currentDirection,
        rankSep: 180,
        nodeSep: 250
      });

      const positions = layout.calculatePositions();

      // Plotly 데이터 생성
      const plotData = this.createPlotlyData(nodes, edges, positions);

      // Plotly 레이아웃 설정
      const plotLayout = {
        showlegend: false,
        hovermode: 'closest',
        margin: { t: 30, b: 30, l: 30, r: 30 },
        xaxis: {
          showgrid: false,
          zeroline: false,
          showticklabels: false,
          fixedrange: false
        },
        yaxis: {
          showgrid: false,
          zeroline: false,
          showticklabels: false,
          fixedrange: false
        },
        dragmode: 'pan',
        plot_bgcolor: '#34495e',
        paper_bgcolor: '#2c3e50',
        font: {
          family: 'Arial, sans-serif',
          color: '#ecf0f1'
        }
      };

      const config = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['select2d', 'lasso2d', 'toggleSpikelines'],
        displaylogo: false
      };

      // Plotly 렌더링
      await Plotly.newPlot(this.$refs.dagPlot, plotData, plotLayout, config);

      // 이벤트 리스너 설정
      this.setupPlotlyEvents();
    },

    /**
     * Plotly 데이터 생성
     */
    createPlotlyData(nodes, edges, positions) {
      const traces = [];

      // 1. Edge Trace (연결선)
      const edgeTrace = {
        x: [],
        y: [],
        mode: 'lines',
        line: {
          width: 3,
          color: 'rgba(150, 150, 150, 0.4)'
        },
        hoverinfo: 'none',
        type: 'scatter',
        name: 'edges',
        showlegend: false
      };

      edges.forEach(edge => {
        const sourcePos = positions[edge.source];
        const targetPos = positions[edge.target];

        if (sourcePos && targetPos) {
          edgeTrace.x.push(sourcePos.x, targetPos.x, null);
          edgeTrace.y.push(sourcePos.y, targetPos.y, null);
        }
      });

      // 2. Arrow Markers (방향 표시)
      const arrowTrace = this.createArrowMarkers(edges, positions);

      // 3. Node Trace (노드)
      const nodeTrace = {
        x: [],
        y: [],
        text: [],
        mode: 'markers+text',
        textposition: this.textPosition === 'center' ? 'middle center' : 'top center',
        textfont: {
          size: 20,
          color: this.textPosition === 'center' ? 'white' : '#ecf0f1',
          family: 'Arial, sans-serif',
          weight: 'bold'
        },
        marker: {
          size: 80,
          symbol: 'circle',
          color: [],
          line: {
            color: 'white',
            width: 3
          },
          opacity: 0.95
        },
        hoverinfo: 'text',
        hovertext: [],
        type: 'scatter',
        name: 'nodes',
        showlegend: false
      };

      // 노드 데이터 채우기
      nodes.forEach(node => {
        const pos = positions[node.id];
        if (pos) {
          nodeTrace.x.push(pos.x);
          nodeTrace.y.push(pos.y);

          // 텍스트 위치에 따른 라벨 선택
          const label = this.textPosition === 'center' ?
            (node.shortLabel || node.label) : node.label;
          nodeTrace.text.push(label);

          // 노드 색상 (상태 기반)
          const status = this.nodeStatuses[node.id] || 'pending';
          console.log(`Node ${node.id} status lookup:`, status, 'available statuses:', Object.keys(this.nodeStatuses || {}));
          nodeTrace.marker.color.push(this.statusColors[status]);

          // 호버 텍스트
          nodeTrace.hovertext.push(this.createHoverText(node));
        }
      });

      traces.push(edgeTrace, arrowTrace, nodeTrace);

      return traces;
    },

    /**
     * 화살표 마커 생성
     */
    createArrowMarkers(edges, positions) {
      const arrowTrace = {
        x: [],
        y: [],
        mode: 'markers',
        marker: {
          size: 12,
          symbol: 'triangle-up',
          color: '#2196F3',
          angle: []
        },
        hoverinfo: 'none',
        type: 'scatter',
        name: 'arrows',
        showlegend: false
      };

      edges.forEach(edge => {
        const sourcePos = positions[edge.source];
        const targetPos = positions[edge.target];

        if (sourcePos && targetPos) {
          const dx = targetPos.x - sourcePos.x;
          const dy = targetPos.y - sourcePos.y;
          const angle = Math.atan2(dy, dx) * 180 / Math.PI - 90;
          const dist = Math.sqrt(dx * dx + dy * dy);

          // 노드 경계에서 떨어진 위치
          const arrowOffset = 50;
          const ratio = (dist - arrowOffset) / dist;

          if (ratio > 0 && ratio < 1) {
            arrowTrace.x.push(sourcePos.x + dx * ratio);
            arrowTrace.y.push(sourcePos.y + dy * ratio);
            arrowTrace.marker.angle.push(angle);
          }
        }
      });

      return arrowTrace;
    },

    /**
     * 호버 텍스트 생성
     */
    createHoverText(node) {
      const status = this.nodeStatuses[node.id] || 'pending';
      let hoverText = `<b>${node.label}</b><br>`;
      hoverText += `Status: <b>${status.toUpperCase()}</b><br>`;

      if (node.description) {
        hoverText += `${node.description}<br>`;
      }

      if (node.inputs && node.inputs.length > 0) {
        hoverText += `Inputs: ${node.inputs.length} files<br>`;
      }

      if (node.outputs && node.outputs.length > 0) {
        hoverText += `Outputs: ${node.outputs.length} files<br>`;
      }

      hoverText += '<br><i>Click for details</i>';

      return hoverText;
    },

    /**
     * Plotly 이벤트 설정
     */
    setupPlotlyEvents() {
      const plotElement = this.$refs.dagPlot;

      // 노드 클릭 이벤트
      plotElement.on('plotly_click', (data) => {
        if (data.points[0] && data.points[0].data.name === 'nodes') {
          const nodeIndex = data.points[0].pointIndex;
          this.handleNodeClick(nodeIndex);
        }
      });

      // 더블클릭으로 뷰 리셋
      plotElement.on('plotly_doubleclick', () => {
        this.resetView();
      });
    },

    /**
     * 노드 클릭 처리
     */
    handleNodeClick(nodeIndex) {
      if (this.dagData && this.dagData.nodes) {
        this.selectedNode = this.dagData.nodes[nodeIndex];
        this.highlightConnectedNodes(nodeIndex);
      }
    },

    /**
     * 연결된 노드 하이라이트
     */
    highlightConnectedNodes(nodeIndex) {
      if (!this.dagData) return;

      const selectedNodeId = this.dagData.nodes[nodeIndex].id;
      const connected = new Set([nodeIndex]);

      // 연결된 노드 찾기
      this.dagData.edges.forEach(edge => {
        if (edge.source === selectedNodeId) {
          const targetIndex = this.dagData.nodes.findIndex(n => n.id === edge.target);
          if (targetIndex !== -1) connected.add(targetIndex);
        }
        if (edge.target === selectedNodeId) {
          const sourceIndex = this.dagData.nodes.findIndex(n => n.id === edge.source);
          if (sourceIndex !== -1) connected.add(sourceIndex);
        }
      });

      // 불투명도 조정
      const opacity = this.dagData.nodes.map((_, idx) =>
        connected.has(idx) ? 1 : 0.3
      );

      Plotly.restyle(this.$refs.dagPlot, {
        'marker.opacity': [opacity]
      }, [2]); // nodes trace index
    },

    /**
     * 노드 색상 업데이트
     */
    updateNodeColors() {
      if (!this.dagData || !this.$refs.dagPlot) return;

      const colors = this.dagData.nodes.map(node => {
        const status = this.nodeStatuses[node.id] || 'pending';
        return this.statusColors[status];
      });

      Plotly.restyle(this.$refs.dagPlot, {
        'marker.color': [colors]
      }, [2]); // nodes trace index
    },



    /**
     * 상태 수동 새로고침
     */
    async refreshStatus() {
      console.log('refreshStatus called');
      await this.updateRuleStatuses();
    },

    /**
     * PNG 내보내기
     */
    exportPNG() {
      if (this.$refs.dagPlot) {
        Plotly.downloadImage(this.$refs.dagPlot, {
          format: 'png',
          width: 1920,
          height: 1080,
          filename: `dag-${this.taskName}-${Date.now()}`
        });
      }
    },

    /**
     * SVG 내보내기
     */
    exportSVG() {
      if (this.$refs.dagPlot) {
        Plotly.downloadImage(this.$refs.dagPlot, {
          format: 'svg',
          width: 1920,
          height: 1080,
          filename: `dag-${this.taskName}-${Date.now()}`
        });
      }
    },


    /**
     * 뷰 리셋
     */
    resetView() {
      this.selectedNode = null;

      if (this.dagData && this.$refs.dagPlot) {
        const opacity = new Array(this.dagData.nodes.length).fill(1);
        Plotly.restyle(this.$refs.dagPlot, {
          'marker.opacity': [opacity]
        }, [2]);
      }
    },

    /**
     * Rule 로그 보기
     */
    async showRuleLogs() {
      if (!this.selectedNode || !this.taskId) return;

      try {
        const response = await getRuleLogs(this.taskId, this.selectedNode.id);

        // 로그 데이터를 별도 모달이나 탭에서 표시
        console.log('Rule logs:', response.data);
        // TODO: 로그 뷰어 모달 구현
      } catch (error) {
        console.error('Failed to load rule logs:', error);
        // 사용자 친화적인 오류 메시지 표시
        alert('Failed to load rule logs. Please try again.');
      }
    },

    /**
     * 파일명 추출
     */
    getFileName(filepath) {
      return filepath.split('/').pop();
    },

    /**
     * 재시도
     */
    async retryLoad() {
      console.log('retryLoad called');
      await this.loadDAGData();
    },


    /**
     * 모달 닫기
     */
    closeModal() {
      this.$emit('close');
    }
  }
};
</script>

<style scoped>
/* 모달 스타일 */
.dag-modal {
  position: fixed;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  background-color: rgba(0, 0, 0, 0.7);
  display: flex;
  justify-content: center;
  align-items: center;
  z-index: 10000;
}

.dag-content {
  width: 80%;
  height: 70%;
  background-color: #2c3e50;
  border-radius: 16px;
  display: flex;
  flex-direction: column;
  box-shadow: 0 4px 10px rgba(0, 0, 0, 0.3);
}

/* 헤더 스타일 */
.dag-header {
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
  padding: 16px;
  border-bottom: 1px solid #576574;
  background-color: #34495e;
  border-radius: 16px 16px 0 0;
}

.dag-title h3 {
  margin: 0 0 10px 0;
  color: #ecf0f1;
  font-size: 1.3em;
}

.dag-progress {
  width: 400px;
}

.progress-bar {
  width: 100%;
  height: 8px;
  background-color: #576574;
  border-radius: 4px;
  overflow: hidden;
  margin-bottom: 5px;
}

.progress-fill {
  height: 100%;
  background: linear-gradient(90deg, #4CAF50 0%, #2196F3 100%);
  transition: width 0.3s ease;
}

.progress-fill.has-failed {
  background: linear-gradient(90deg, #4CAF50 0%, #F44336 50%, #2196F3 100%);
}

.progress-text {
  font-size: 0.9em;
  color: #ecf0f1;
}

.failed-count {
  color: #F44336;
  font-weight: bold;
}

.running-count {
  color: #2196F3;
  font-weight: bold;
}

/* 컨트롤 스타일 */
.dag-controls {
  display: flex;
  gap: 8px;
  align-items: flex-start;
}

.control-btn {
  display: flex;
  align-items: center;
  gap: 4px;
  padding: 6px 10px;
  background: #2196F3;
  color: white;
  border: none;
  border-radius: 8px;
  cursor: pointer;
  font-size: 0.8em;
  transition: background-color 0.2s;
}

.control-btn:hover {
  background: #1976D2;
}

.close-btn {
  padding: 6px 10px;
  background: #f44336;
  color: white;
  border: none;
  border-radius: 8px;
  cursor: pointer;
  font-size: 1.5rem;
}

.close-btn:hover {
  background: #d32f2f;
}


/* 바디 스타일 */
.dag-body {
  display: flex;
  flex: 1;
  overflow: hidden;
}

.dag-plot {
  flex: 1;
  position: relative;
  background-color: #34495e;
}

.dag-plot.loading {
  display: flex;
  align-items: center;
  justify-content: center;
}

.loading-overlay {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 20px;
  color: #666;
}

.spinner {
  width: 40px;
  height: 40px;
  border: 4px solid #e3f2fd;
  border-top: 4px solid #2196f3;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  0% {
    transform: rotate(0deg);
  }

  100% {
    transform: rotate(360deg);
  }
}

.error-message {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 10px;
  color: #f44336;
}

.retry-btn {
  padding: 8px 16px;
  background: #2196f3;
  color: white;
  border: none;
  border-radius: 4px;
  cursor: pointer;
}

/* 사이드바 스타일 */
.dag-sidebar {
  width: 280px;
  border-left: 1px solid #576574;
  background-color: #2c3e50;
  overflow-y: auto;
}

.node-details {
  padding: 16px;
}

.node-details h4 {
  margin: 0 0 16px 0;
  color: #ecf0f1;
  border-bottom: 2px solid #576574;
  padding-bottom: 8px;
}

.detail-section {
  margin-bottom: 15px;
}

.detail-label {
  font-weight: bold;
  color: #ecf0f1;
  margin-bottom: 5px;
}

.detail-content {
  color: #bdc3c7;
  line-height: 1.4;
}

.status-badge {
  padding: 4px 8px;
  border-radius: 4px;
  font-size: 0.8em;
  font-weight: bold;
  text-transform: uppercase;
}

.status-pending {
  background-color: #9E9E9E;
  color: white;
}

.status-running {
  background-color: #2196F3;
  color: white;
}

.status-success {
  background-color: #4CAF50;
  color: white;
}

.status-failed {
  background-color: #F44336;
  color: white;
}

.file-item,
.param-item {
  padding: 2px 0;
  font-family: monospace;
  font-size: 0.9em;
}

.script-content {
  background-color: #1e2833;
  padding: 10px;
  border-radius: 8px;
  font-family: monospace;
  font-size: 0.8em;
  max-height: 100px;
  overflow-y: auto;
  color: #ecf0f1;
}

.detail-actions {
  margin-top: 20px;
  display: flex;
  gap: 10px;
}

.action-btn {
  padding: 8px 12px;
  border: none;
  border-radius: 4px;
  cursor: pointer;
  font-size: 0.9em;
  transition: background-color 0.2s;
}

.action-btn:not(.secondary) {
  background: #2196F3;
  color: white;
}

.action-btn:not(.secondary):hover {
  background: #1976D2;
}

.action-btn.secondary {
  background: #576574;
  color: #ecf0f1;
}

.action-btn.secondary:hover {
  background: #495A68;
}

/* 푸터 스타일 */
.dag-footer {
  padding: 12px 16px;
  border-top: 1px solid #576574;
  background-color: #34495e;
  border-radius: 0 0 16px 16px;
}

.status-legend {
  display: flex;
  justify-content: center;
  gap: 25px;
}

.legend-item {
  display: flex;
  align-items: center;
  gap: 8px;
  font-size: 0.9em;
  color: #ecf0f1;
}

.status-dot {
  width: 12px;
  height: 12px;
  border-radius: 50%;
}

.status-dot.pending {
  background-color: #9E9E9E;
}

.status-dot.running {
  background-color: #f39c12;
}

.status-dot.success {
  background-color: #4CAF50;
}

.status-dot.failed {
  background-color: #F44336;
}


/* 반응형 디자인 */
@media (max-width: 768px) {
  .dag-content {
    width: 95%;
    height: 85%;
  }

  .dag-header {
    flex-direction: column;
    gap: 12px;
    padding: 12px;
  }

  .dag-progress {
    width: 100%;
  }

  .dag-body {
    flex-direction: column;
  }

  .dag-sidebar {
    width: 100%;
    max-height: 200px;
  }

  .dag-controls {
    flex-wrap: wrap;
    gap: 4px;
  }

  .control-btn {
    padding: 4px 8px;
    font-size: 0.7em;
  }
}
</style>