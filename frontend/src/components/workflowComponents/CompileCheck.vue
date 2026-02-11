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
                        <div class="task-arrow">&rarr;</div>

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
                <h3 class="resource-info__title">Resource Allocation Status</h3>

                <!-- Error state -->
                <div v-if="hasError" class="error-state">
                    <p>Resource monitoring temporarily unavailable</p>
                </div>

                <!-- Worker not initialized -->
                <div v-else-if="resources && resources.cpu.total === 0 && resources.gpu.total === 0" class="empty-state">
                    <p>Waiting for worker initialization...</p>
                </div>

                <!-- Resource overview -->
                <div v-else-if="resources" class="resource-overview">
                    <!-- CPU Slots -->
                    <div class="resource-bar">
                        <label>CPU Slots: {{ resources.cpu.available }} / {{ resources.cpu.total }} available</label>
                        <div class="bar">
                            <div class="fill" :class="getUsageClass(cpuPercent)"
                                :style="{ width: cpuPercent + '%' }">
                            </div>
                        </div>
                        <span class="bar-percent" :class="getUsageClass(cpuPercent)">{{ cpuPercent.toFixed(0) }}% used</span>
                    </div>

                    <!-- Memory -->
                    <div v-if="resources.memory" class="resource-bar">
                        <label>Memory: {{ formatBytes(resources.memory.used_bytes) }} / {{ formatBytes(resources.memory.total_bytes) }}</label>
                        <div class="bar">
                            <div class="fill" :class="getUsageClass(resources.memory.percent)"
                                :style="{ width: Math.min(resources.memory.percent, 100) + '%' }">
                            </div>
                        </div>
                        <span class="bar-percent" :class="getUsageClass(resources.memory.percent)">{{ resources.memory.percent.toFixed(1) }}%</span>
                        <div class="resource-detail-line">
                            Available: {{ formatBytes(resources.memory.available_bytes) }}
                        </div>
                    </div>

                    <!-- GPU Slots (only when gpu.total > 0) -->
                    <div v-if="resources.gpu.total > 0" class="resource-bar">
                        <label>GPU Slots: {{ resources.gpu.available }} / {{ resources.gpu.total }} available</label>
                        <div class="bar">
                            <div class="fill" :class="getUsageClass(gpuPercent)"
                                :style="{ width: gpuPercent + '%' }">
                            </div>
                        </div>
                        <span class="bar-percent" :class="getUsageClass(gpuPercent)">{{ gpuPercent.toFixed(0) }}% used</span>
                    </div>

                    <!-- GPU Device details (only when gpu_devices has items) -->
                    <div v-if="resources.gpu_devices && resources.gpu_devices.length > 0" class="gpu-devices">
                        <div v-for="gpu in resources.gpu_devices" :key="gpu.id" class="gpu-device">
                            <label>[GPU {{ gpu.id }}] {{ gpu.name }} - {{ gpu.load_percent }}%</label>
                            <div class="resource-detail-line">
                                VRAM: {{ formatBytes(gpu.memory_used_bytes) }} / {{ formatBytes(gpu.memory_total_bytes) }}
                            </div>
                            <div v-if="gpu.temperature_c !== null" class="resource-detail-line">
                                Temp: {{ gpu.temperature_c }}&deg;C
                            </div>
                        </div>
                    </div>

                    <!-- Running Tasks -->
                    <div class="task-list-section">
                        <label class="task-list-label">Running Tasks ({{ resources.tasks.length }})</label>
                        <div v-if="resources.tasks.length > 0" class="task-list">
                            <div v-for="task in resources.tasks" :key="task.task_id" class="task-item">
                                <span class="task-dot">&#9679;</span>
                                <span class="task-name">{{ task.plugin_name || 'Unknown' }}</span>
                                <span class="resource-badge" :class="'badge-' + task.resource_type">
                                    {{ task.resource_type.toUpperCase() }} x{{ task.resource_slots }}
                                </span>
                                <span class="task-elapsed">{{ getElapsedTime(task.started_at) }}</span>
                            </div>
                        </div>
                        <div v-else class="empty-state">
                            <p>No active tasks</p>
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
import { getTaskResources } from '@/api/index';
import { formatBytes, getRunningTime } from '@/utils/formatters';

export default {
    data() {
        return {
            taskInfoList: [],
            resources: null,
            isLoadingResources: true,
            isPolling: false,
            resourceTimeoutId: null,
            currentTime: new Date(),
            currentTimeIntervalId: null,
            hasError: false,
        };
    },
    computed: {
        cpuPercent() {
            if (!this.resources || this.resources.cpu.total === 0) return 0;
            return (this.resources.cpu.used / this.resources.cpu.total) * 100;
        },
        gpuPercent() {
            if (!this.resources || this.resources.gpu.total === 0) return 0;
            return (this.resources.gpu.used / this.resources.gpu.total) * 100;
        },
    },
    async mounted() {
        try {
            const response = await getTaskResources();
            this.resources = response.data;
            this.hasError = false;
        } catch (error) {
            console.error('Failed to load initial resources:', error);
            if (error.response && error.response.status === 503) {
                this.hasError = true;
            }
        } finally {
            this.isLoadingResources = false;
        }

        this.startPolling();

        // 1초마다 currentTime 갱신 (경과시간 실시간 표시)
        this.currentTimeIntervalId = setInterval(() => {
            this.currentTime = new Date();
        }, 1000);

        // workflow 정보를 통해 Algorithm 노드들의 정보 가져오기
        try {
            const workflow_info = this.$store.getters.getWorkflowInfo;

            const nodes_list = Object.values(workflow_info.drawflow.Home.data);
            const algorithm_nodes = nodes_list.filter(node => node.class === 'Algorithm');

            // 각 노드의 정보를 taskInfoList에 추가
            algorithm_nodes.forEach(node => {
                const taskInfo = {
                    pluginName: node.data.selectedPlugin.name,
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
        if (this.currentTimeIntervalId) {
            clearInterval(this.currentTimeIntervalId);
            this.currentTimeIntervalId = null;
        }
    },
    methods: {
        formatBytes,
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
                const response = await getTaskResources();
                this.resources = response.data;
                this.hasError = false;
            } catch (error) {
                if (error.response && error.response.status === 503) {
                    this.hasError = true;
                }
                console.error('Failed to refresh resources:', error);
            } finally {
                if (this.isPolling) {
                    this.resourceTimeoutId = setTimeout(this.pollResources, 5000);
                }
            }
        },
        getElapsedTime(startedAt) {
            if (!startedAt) return '--:--:--';
            return getRunningTime(startedAt, this.currentTime);
        },
        getUsageClass(percent) {
            if (percent >= 90) return 'usage-critical';
            if (percent >= 70) return 'usage-warning';
            return 'usage-normal';
        },
        confirmTask() {
            alert("Task is being executed...");
            this.closeModal();
            this.$emit('run-workflow');
        },
        closeModal() {
            this.stopPolling();
            this.$emit('deactivate-compile-check');
        },
    },
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
    z-index: 10000;
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
    z-index: 10001;
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

    scrollbar-color: #2c3e50 #34495e;
    scrollbar-width: thin;
}

.task-info::-webkit-scrollbar,
.resource-info::-webkit-scrollbar {
    width: 8px;
}

.task-info::-webkit-scrollbar-track,
.resource-info::-webkit-scrollbar-track {
    background: #34495e;
    border-radius: 1rem;
}

.task-info::-webkit-scrollbar-thumb,
.resource-info::-webkit-scrollbar-thumb {
    background-color: #2c3e50;
    border-radius: 1rem;
}

.task-info::-webkit-scrollbar-thumb:hover,
.resource-info::-webkit-scrollbar-thumb:hover {
    background-color: #1f2a38;
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

/* Resource bars */
.resource-overview {
    display: flex;
    flex-direction: column;
    gap: 12px;
}

.resource-bar {
    margin-bottom: 2px;
}

.resource-bar label {
    font-size: 0.85rem;
    color: #bdc3c7;
    display: block;
    margin-bottom: 4px;
}

.bar {
    width: 100%;
    height: 12px;
    background-color: #46627e;
    border-radius: 5px;
    overflow: hidden;
}

.fill {
    height: 100%;
    background-color: #3498db;
    border-radius: 5px;
    transition: width 0.3s ease;
}

.fill.usage-warning {
    background-color: #f39c12;
}

.fill.usage-critical {
    background-color: #e74c3c;
}

.bar-percent {
    font-size: 0.75rem;
    margin-top: 2px;
    display: inline-block;
}

.resource-detail-line {
    font-size: 0.8rem;
    color: #95a5a6;
    margin-top: 2px;
}

/* GPU devices */
.gpu-devices {
    margin-top: 4px;
}

.gpu-device {
    background-color: #2c3e50;
    padding: 8px 10px;
    border-radius: 6px;
    margin-bottom: 6px;
}

.gpu-device label {
    font-size: 0.85rem;
    font-weight: bold;
    color: #ecf0f1;
    display: block;
    margin-bottom: 4px;
}

/* Running tasks list */
.task-list-section {
    margin-top: 8px;
}

.task-list-label {
    font-size: 0.9rem;
    color: #bdc3c7;
    font-weight: bold;
    display: block;
    margin-bottom: 6px;
}

.task-list {
    display: flex;
    flex-direction: column;
    gap: 6px;
}

.task-item {
    display: flex;
    align-items: center;
    gap: 8px;
    background-color: #2c3e50;
    padding: 6px 10px;
    border-radius: 6px;
    font-size: 0.85rem;
}

.task-dot {
    color: #2ecc71;
    font-size: 0.6rem;
}

.task-name {
    flex: 1;
    font-weight: 500;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}

.resource-badge {
    padding: 2px 8px;
    border-radius: 4px;
    font-size: 0.7rem;
    font-weight: bold;
    white-space: nowrap;
}

.badge-cpu {
    background-color: #2980b9;
    color: #fff;
}

.badge-gpu {
    background-color: #8e44ad;
    color: #fff;
}

.task-elapsed {
    font-family: monospace;
    font-size: 0.8rem;
    color: #95a5a6;
    white-space: nowrap;
}

/* States */
.empty-state {
    text-align: center;
    color: #7f8c8d;
    padding: 12px 0;
    font-size: 0.9rem;
}

.error-state {
    text-align: center;
    color: #e67e22;
    padding: 12px 0;
    font-size: 0.9rem;
}

/* Usage color classes for text */
.usage-normal {
    color: #2ecc71;
}

.usage-warning {
    color: #f39c12;
}

.usage-critical {
    color: #e74c3c;
}

/* Modal actions */
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
    z-index: 10002;
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
