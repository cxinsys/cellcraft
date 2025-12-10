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
import { getDAGStructure, getRuleStatuses } from '@/api';

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
        // v-if로 인해 DOM이 렌더링되기 전에 watcher가 실행될 수 있음
        // $nextTick으로 DOM 렌더링 완료 후 데이터 로드
        this.$nextTick(() => {
          this.loadDAGData();
        });
      }
    },
    taskId(newValue) {
      console.log('DAGVisualization taskId changed:', newValue, 'visible:', this.visible);
      if (this.visible && newValue) {
        this.$nextTick(() => {
          this.loadDAGData();
        });
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
      // $nextTick으로 DOM이 완전히 준비된 후 데이터 로드
      this.$nextTick(() => {
        this.loadDAGData();
      });
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

        // 로딩 상태 해제 (렌더링 전에 해제하여 CSS가 플롯에 영향을 주지 않도록)
        this.isLoading = false;

        // DOM 업데이트 대기 후 DAG 시각화 렌더링
        await this.$nextTick();
        await this.renderDAG();
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
      if (!this.dagData) {
        console.warn('renderDAG: dagData is not available');
        return;
      }

      // DOM이 준비될 때까지 대기
      await this.$nextTick();

      if (!this.$refs.dagPlot) {
        console.warn('renderDAG: dagPlot ref is not available, retrying...');
        // 한 번 더 시도
        await new Promise(resolve => setTimeout(resolve, 100));
        await this.$nextTick();
        if (!this.$refs.dagPlot) {
          console.error('renderDAG: dagPlot ref is still not available after retry');
          return;
        }
      }

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

      // 컨테이너 크기에 맞게 플롯 리사이즈 (모달이 완전히 표시된 후)
      await this.$nextTick();
      if (this.$refs.dagPlot) {
        Plotly.Plots.resize(this.$refs.dagPlot);
      }
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
     * 노드 타입에 따른 동적 설명 생성
     */
    getTypeDescription(nodeType) {
      const typeDescriptions = {
        'input_processing': 'Data input and preprocessing',
        'output': 'Result output generation',
        'analysis': 'Data analysis step',
        'network_analysis': 'Network reconstruction and analysis',
        'preprocessing': 'Data filtering and cleaning',
        'visualization': 'Result visualization',
        'process': 'Processing step'
      };
      return typeDescriptions[nodeType] || 'Processing step';
    },

    /**
     * 호버 텍스트 생성 (동적 정보 표시)
     */
    createHoverText(node) {
      const status = this.nodeStatuses[node.id] || 'pending';
      let hoverText = `<b>${node.label}</b><br>`;
      hoverText += `Status: <b>${status.toUpperCase()}</b><br>`;

      // 노드 타입 기반 동적 설명
      if (node.type) {
        hoverText += `<i>${this.getTypeDescription(node.type)}</i><br>`;
      }

      hoverText += '<br>';

      // Input 파일명 표시 (최대 3개)
      if (node.inputs && node.inputs.length > 0) {
        hoverText += `<b>Inputs:</b><br>`;
        const displayInputs = node.inputs.slice(0, 3);
        displayInputs.forEach(input => {
          hoverText += `  • ${this.getFileName(input)}<br>`;
        });
        if (node.inputs.length > 3) {
          hoverText += `  <i>+${node.inputs.length - 3} more</i><br>`;
        }
      }

      // Output 파일명 표시 (최대 3개)
      if (node.outputs && node.outputs.length > 0) {
        hoverText += `<b>Outputs:</b><br>`;
        const displayOutputs = node.outputs.slice(0, 3);
        displayOutputs.forEach(output => {
          hoverText += `  • ${this.getFileName(output)}<br>`;
        });
        if (node.outputs.length > 3) {
          hoverText += `  <i>+${node.outputs.length - 3} more</i><br>`;
        }
      }

      return hoverText;
    },

    /**
     * Plotly 이벤트 설정
     */
    setupPlotlyEvents() {
      // 현재는 hover만 사용하므로 추가 이벤트 설정 없음
      // Plotly의 기본 hover 동작 사용
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