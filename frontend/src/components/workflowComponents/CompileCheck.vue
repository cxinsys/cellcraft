<template>
    <div class="modal-content">
        <div class="modal-container" :class="{ 'loading-active': isLoadingResources }">
            <!-- Full modal loading overlay -->
            <div v-if="isLoadingResources" class="full-loading-overlay">
                <div class="loading-spinner"></div>
                <p class="loading-text">Loading resource information...</p>
            </div>
            
            <h2 class="modal-title">Confirm Task Execution</h2>

            <div class="task-info">
                <h3 class="task-info__title">Task Information</h3>
                <div v-for="(task, index) in taskInfoList" :key="index" class="task-info__item">
                    <div class="task-plugin">{{ task.pluginName }}</div>
                    <div class="task-container">
                        <!-- Input Section -->
                        <div class="task-inputs">
                            <div v-for="(input, inputIndex) in task.inputs" :key="inputIndex" class="task-input">
                                {{ input }}
                            </div>
                        </div>

                        <!-- Arrow -->
                        <div class="task-arrow">→</div>

                        <!-- Output Section -->
                        <div class="task-outputs">
                            <div v-for="(output, outputIndex) in task.outputs" :key="outputIndex" class="task-output">
                                {{ output }}
                            </div>
                        </div>
                    </div>
                </div>
            </div>

            <div class="resource-info">
                <h3 class="resource-info__title">Server Resource Status</h3>

                <!-- CPU 정보 표시 -->
                <div class="resource-bar">
                    <label>CPU Usage: {{ Number(serverResources.total_cpu_usage_percent).toFixed(2) }}%</label>
                    <div class="bar">
                        <div class="fill"
                            :style="{ width: Math.min(serverResources.total_cpu_usage_percent, 100) + '%' }">
                        </div>
                    </div>
                    <div class="resource-details">
                        <p><strong>Total Memory Usage:</strong> {{ formatBytes(serverResources.total_memory_usage_bytes)
                            }}</p>
                        <p><strong>Total Memory Limit:</strong> {{ formatBytes(serverResources.total_memory_limit_bytes)
                            }}</p>
                    </div>
                </div>

                <!-- 메모리 정보 표시 -->
                <div class="resource-bar">
                    <label>Memory Usage: {{ Number(serverResources.total_memory_usage_percent).toFixed(2) }}%</label>
                    <div class="bar">
                        <div class="fill"
                            :style="{ width: Math.min(serverResources.total_memory_usage_percent, 100) + '%' }"></div>
                    </div>
                    <div class="resource-details">
                        <p>
                            <strong>Total Memory:</strong> {{ formatBytes(serverResources.total_memory_limit_bytes) }}
                        </p>
                        <p>
                            <strong>Used Memory:</strong>
                            <span :class="getMemoryUsageClass(serverResources.total_memory_usage_percent)">
                                {{ formatBytes(serverResources.total_memory_usage_bytes) }}
                            </span>
                        </p>
                        <p>
                            <strong>Available Memory:</strong> {{ formatBytes(serverResources.total_memory_limit_bytes -
                                serverResources.total_memory_usage_bytes) }}
                        </p>
                    </div>
                </div>

                <!-- 컨테이너별 정보 표시 -->
                <div v-for="(container, index) in serverResources.containers" :key="index" class="container-info">
                    <h4>{{ container.container_info.name }}</h4>

                    <!-- 컨테이너 CPU 정보 -->
                    <div class="resource-sub-bar">
                        <label>CPU Usage: {{ Number(container.cpu.usage_percent).toFixed(2) }}%</label>
                        <div class="bar">
                            <div class="fill" :style="{ width: Math.min(container.cpu.usage_percent, 100) + '%' }">
                            </div>
                        </div>
                    </div>

                    <!-- 컨테이너 메모리 정보 -->
                    <div class="resource-sub-bar">
                        <label>Memory Usage: {{ Number(container.memory.percent).toFixed(2) }}%</label>
                        <div class="bar">
                            <div class="fill" :style="{ width: Math.min(container.memory.percent, 100) + '%' }"></div>
                        </div>
                    </div>

                    <!-- 컨테이너 GPU 정보 -->
                    <div v-if="container.gpu" class="gpu-info">
                        <div v-for="(gpu, gpuIndex) in container.gpu" :key="gpuIndex">
                            <label>{{ gpu.name }} (GPU {{ gpu.id }})</label>
                            <div class="resource-sub-bar">
                                <label>GPU Usage: {{ gpu.utilization_percent }}%</label>
                                <div class="bar">
                                    <div class="fill" :style="{ width: Math.min(gpu.utilization_percent, 100) + '%' }">
                                    </div>
                                </div>
                            </div>
                            <div class="resource-sub-bar">
                                <label>Memory Usage: {{ gpu.memory.utilization_percent.toFixed(2) }}%</label>
                                <div class="bar">
                                    <div class="fill"
                                        :style="{ width: Math.min(gpu.memory.utilization_percent, 100) + '%' }"></div>
                                </div>
                            </div>
                            <div class="resource-details">
                                <p>
                                    <strong>Memory:</strong>
                                    {{ formatBytes(gpu.memory.used_bytes) }} / {{ formatBytes(gpu.memory.total_bytes) }}
                                </p>
                                <p>
                                    <strong>Temperature:</strong> {{ gpu.temperature_c }}°C
                                </p>
                                <p>
                                    <strong>Power:</strong> {{ gpu.power.draw_watts }}W / {{ gpu.power.limit_watts }}W
                                </p>
                            </div>
                        </div>
                    </div>
                </div>
            </div>

            <!-- 실행/취소 버튼 -->
            <div class="modal-actions">
                <button class="btn confirm" @click="confirmTask">Execute</button>
                <button class="btn cancel" @click="closeModal">Cancel</button>
            </div>
        </div>
    </div>
</template>

<script>
import { getSystemResources } from '@/api/index';
import { formatBytes } from '@/utils/formatters';

export default {
    data() {
        return {
            taskInfoList: [],
            serverResources: {
                total_containers: 0,
                total_cpu_usage_percent: 0,
                total_memory_usage_bytes: 0,
                total_memory_limit_bytes: 0,
                total_memory_usage_percent: 0,
                containers: []
            },
            isLoadingResources: true,
            isPolling: false,
            resourceTimeoutId: null,
        };
    },
    async mounted() {
        try {
            const response = await getSystemResources();
            this.serverResources = response.data;
            this.isLoadingResources = false;
        } catch (error) {
            console.error('Failed to load initial resources:', error);
            this.isLoadingResources = false;
        }
        
        this.startPolling();

        // workflow 정보를 통해 Algorithm 노드들의 정보 가져오기
        try {
            const workflow_info = this.$store.getters.getWorkflowInfo;
            console.log(workflow_info);

            const nodes_list = Object.values(workflow_info.drawflow.Home.data);
            const algorithm_nodes = nodes_list.filter(node => node.class === 'Algorithm');
            console.log(algorithm_nodes);

            // 각 노드의 정보를 taskInfoList에 추가
            algorithm_nodes.forEach(node => {
                const taskInfo = {
                    pluginName: node.data.selectedPlugin.name,
                    // node.data.selectedPluginInputOutput에서 activate가 true고 type이 inputFile인 것만 고르기
                    inputs: node.data.selectedPluginInputOutput.filter(input => input.activate && input.type === 'inputFile').map(input => input.defaultValue),
                    outputs: node.data.selectedPluginInputOutput.filter(output => output.activate && output.type === 'outputFile').map(output => output.defaultValue)
                };
                this.taskInfoList.push(taskInfo);
            });
        } catch (error) {
            console.error(error);
        }
    },
    beforeDestroy() {
        this.stopPolling();
    },
    methods: {
        startPolling() {
            this.isPolling = true;
            this.pollResources();
        },
        stopPolling() {
            this.isPolling = false;
            if (this.resourceTimeoutId) {
                clearTimeout(this.resourceTimeoutId);
                this.resourceTimeoutId = null;
            }
        },
        async pollResources() {
            if (!this.isPolling) return;

            try {
                const response = await getSystemResources();
                this.serverResources = response.data;
                console.log(this.serverResources);
            } catch (error) {
                console.error('Failed to refresh resources:', error);
            } finally {
                if (this.isPolling) {
                    this.resourceTimeoutId = setTimeout(this.pollResources, 10000);
                }
            }
        },
        formatBytes,
        confirmTask() {
            alert("Task is being executed...");
            this.closeModal();
            this.$emit('run-workflow');
        },
        closeModal() {
            this.stopPolling();
            this.$emit('deactivate-compile-check');
        },
        formatCPUUsage(usage) {
            return usage ? Number(usage).toLocaleString() : '0';
        },
        getCPUUsageClass(usage) {
            if (usage >= 90) return 'usage-critical';
            if (usage >= 70) return 'usage-warning';
            return 'usage-normal';
        },
        getMemoryUsageClass(percent) {
            if (percent >= 90) return 'usage-critical';
            if (percent >= 70) return 'usage-warning';
            return 'usage-normal';
        },
        calculateCPUPercentage(usage) {
            if (!this.serverResources.cpu.system_usage) return 0;
            return (usage / (this.serverResources.cpu.system_usage / this.serverResources.cpu.num_cpus)) * 100;
        }
    },
    computed: {
        getTotalMemory() {
            return this.serverResources.available_memory_bytes / (1 - this.serverResources.memory_usage_percent / 100);
        },
        getUsedMemory() {
            return this.getTotalMemory - this.serverResources.available_memory_bytes;
        }
    }
};
</script>

<style scoped>
/* 모달 콘텐츠 */
.modal-content {
    position: fixed;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    background-color: rgba(0, 0, 0, 0.5);
    display: flex;
    justify-content: center;
    align-items: center;
    z-index: 10000; /* 높은 z-index 설정 */
}

.modal-container {
    background-color: #2c3e50;
    color: #ecf0f1;
    padding: 1rem;
    border-radius: 1rem;
    width: 480px;
    max-width: 90%;
    text-align: center;
    position: relative;
    z-index: 10001; /* modal-content보다 높은 z-index */
}

.modal-container.loading-active > *:not(.full-loading-overlay) {
    filter: blur(2px);
    opacity: 0.6;
    pointer-events: none;
}

.modal-title {
    font-size: 1.5em;
    margin-bottom: 20px;
}

.task-info,
.resource-info {
    height: 240px;
    background-color: #34495e;
    padding: 0.5rem;
    border-radius: 1rem;
    overflow-y: auto;
    text-align: left;
    margin-bottom: 1.5rem;

    /* 스크롤바의 색상 설정 (Firefox) */
    scrollbar-color: #2c3e50 #34495e;
    /* 스크롤바 핸들, 트랙 색상 */
    scrollbar-width: thin;
    /* 스크롤바 두께 설정 (Firefox) */
}

/* WebKit 기반 브라우저용 스크롤바 스타일링 (Chrome, Edge, Safari 등) */
.task-info::-webkit-scrollbar,
.resource-info::-webkit-scrollbar {
    width: 8px;
    /* 스크롤바의 너비 */
}

.task-info::-webkit-scrollbar-track,
.resource-info::-webkit-scrollbar-track {
    background: #34495e;
    /* 스크롤바 트랙(배경) 색상 */
    border-radius: 1rem;
}

.task-info::-webkit-scrollbar-thumb,
.resource-info::-webkit-scrollbar-thumb {
    background-color: #2c3e50;
    /* 스크롤바 핸들 색상 */
    border-radius: 1rem;
}

.task-info::-webkit-scrollbar-thumb:hover,
.resource-info::-webkit-scrollbar-thumb:hover {
    background-color: #1f2a38;
    /* 스크롤바 핸들 hover 색상 */
}


.task-info__title,
.resource-info__title {
    font-size: 1.2rem;
    margin-bottom: 1rem;
}

.task-info__item {
    margin-bottom: 10px;
    background-color: #2c3e50;
    padding: 0.5rem;
    border-radius: 1rem;
}

.task-plugin {
    font-size: 18px;
    margin-bottom: 10px;
    display: flex;
    justify-content: center;
    align-items: center;
}

.task-container {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 10px;
    border-radius: 8px;
}

.task-inputs,
.task-outputs {
    display: flex;
    flex-direction: column;
    gap: 10px;
}

.task-input,
.task-output {
    padding: 10px;
    background-color: #1abc9c;
    color: white;
    border-radius: 5px;
    text-align: center;
    font-size: 0.8rem;
    width: 144px;
    text-overflow: ellipsis;
    white-space: nowrap;
    overflow: hidden;
}

.task-arrow {
    font-size: 24px;
    color: white;
    margin: 0 10px;
}

.resource-bar {
    margin-bottom: 10px;
}

.bar {
    width: 100%;
    height: 12px;
    background-color: #46627e;
    border-radius: 5px;
    overflow: hidden;
    margin-top: 5px;
}

.fill {
    height: 100%;
    background-color: #3498db;
}

.resource-details {
    margin-top: 15px;
    color: #ecf0f1;
}

.resource-details p {
    margin: 5px 0;
}

.modal-actions {
    display: flex;
    justify-content: space-between;
}

.btn {
    padding: 12px 24px;
    border: none;
    border-radius: 5px;
    cursor: pointer;
    color: white;
}

.confirm {
    background-color: #27ae60;
}

.confirm:hover {
    background-color: #2ecc71;
}

.cancel {
    background-color: #e74c3c;
}

.cancel:hover {
    background-color: #c0392b;
}

.gpu-info {
    margin-bottom: 15px;
    padding: 10px;
    background-color: #2c3e50;
    border-radius: 8px;
}

.resource-sub-bar {
    margin: 8px 0;
}

.gpu-info label {
    font-weight: bold;
    margin-bottom: 5px;
    display: block;
}

.usage-critical {
    color: #e74c3c;
}

.usage-warning {
    color: #f39c12;
}

.usage-normal {
    color: #2ecc71;
}

.cpu-cores {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
    gap: 8px;
    margin-top: 10px;
}

.usage-normal {
    color: #2ecc71;
}

.usage-warning {
    color: #f1c40f;
}

.usage-critical {
    color: #e74c3c;
}

.resource-details p {
    margin: 8px 0;
    display: flex;
    justify-content: space-between;
    align-items: center;
}

.container-info {
    margin-bottom: 15px;
    padding: 10px;
    background-color: #2c3e50;
    border-radius: 8px;
}

.container-info h4 {
    margin-bottom: 10px;
    color: #ecf0f1;
}

/* Full modal loading overlay */
.full-loading-overlay {
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    z-index: 10002; /* 가장 높은 z-index로 설정 */
    background-color: rgba(44, 62, 80, 0.3);
    border-radius: 1rem;
}

.loading-spinner {
    width: 40px;
    height: 40px;
    border: 4px solid #34495e;
    border-top: 4px solid #3498db;
    border-radius: 50%;
    animation: spin 1s linear infinite;
    margin-bottom: 15px;
}

.loading-text {
    color: #ecf0f1;
    font-size: 1rem;
    margin: 0;
}

@keyframes spin {
    0% { transform: rotate(0deg); }
    100% { transform: rotate(360deg); }
}
</style>