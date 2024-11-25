<template>
    <div class="modal-content">
        <div class="modal-container">
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
                    <label>CPU Usage: {{ Number(serverResources.cpu.usage_percent).toFixed(2) }}% ({{
                        serverResources.cpu.num_cpus }} Cores)</label>
                    <div class="bar">
                        <div class="fill" :style="{ width: Math.min(serverResources.cpu.usage_percent, 100) + '%' }">
                        </div>
                    </div>
                    <div class="resource-details">
                        <p><strong>Total CPU Usage:</strong> {{ formatCPUUsage(serverResources.cpu.total_usage) }}
                            cycles</p>
                        <p><strong>System CPU Usage:</strong> {{ formatCPUUsage(serverResources.cpu.system_usage) }}
                            cycles</p>
                        <div class="cpu-cores">
                            <p v-for="(usage, index) in serverResources.cpu.per_cpu_usage" :key="index">
                                <strong>Core {{ index }}:</strong>
                                <span :class="getCPUUsageClass(calculateCPUPercentage(usage))">
                                    {{ calculateCPUPercentage(usage).toFixed(2) }}%
                                </span>
                            </p>
                        </div>
                    </div>
                </div>

                <!-- 메모리 정보 표시 -->
                <div class="resource-bar">
                    <label>Memory Usage: {{ Number(serverResources.memory.percent).toFixed(2) }}%</label>
                    <div class="bar">
                        <div class="fill" :style="{ width: Math.min(serverResources.memory.percent, 100) + '%' }"></div>
                    </div>
                    <div class="resource-details">
                        <p>
                            <strong>Total Memory:</strong> {{ formatBytes(serverResources.memory.total_bytes) }}
                        </p>
                        <p>
                            <strong>Used Memory:</strong>
                            <span :class="getMemoryUsageClass(serverResources.memory.percent)">
                                {{ formatBytes(serverResources.memory.used_bytes) }}
                            </span>
                        </p>
                        <p>
                            <strong>Available Memory:</strong> {{ formatBytes(serverResources.memory.available_bytes) }}
                        </p>
                    </div>
                </div>

                <!-- GPU 정보 표시 -->
                <div class="resource-bar" v-if="serverResources.gpu">
                    <div v-for="(gpu, index) in serverResources.gpu" :key="index" class="gpu-info">
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

export default {
    data() {
        return {
            taskInfoList: [],
            serverResources: {
                container_info: {
                    id: '',
                    name: '',
                    status: '',
                    created: ''
                },
                cpu: {
                    usage_percent: 0,
                    num_cpus: 0,
                    total_usage: 0,
                    system_usage: 0,
                    per_cpu_usage: []
                },
                memory: {
                    total_bytes: 0,
                    used_bytes: 0,
                    available_bytes: 0,
                    percent: 0
                },
                gpu: null,
                network: {}
            },
            intervalId: null,
        };
    },
    async mounted() {
        const response = await getSystemResources();
        this.serverResources = response.data;
        this.intervalId = setInterval(async () => {
            const response = await getSystemResources();
            this.serverResources = response.data;
            console.log(this.serverResources);
        }, 1000);

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
        if (this.intervalId) {
            clearInterval(this.intervalId);
        }
    },
    methods: {
        formatBytes(bytes, decimals = 2) {
            if (bytes === 0) return '0 Bytes';
            const k = 1024;
            const dm = decimals < 0 ? 0 : decimals;
            const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
            const i = Math.floor(Math.log(bytes) / Math.log(k));
            return parseFloat((bytes / Math.pow(k, i)).toFixed(dm)) + ' ' + sizes[i];
        },
        confirmTask() {
            alert("Task is being executed...");
            if (this.intervalId) {
                clearInterval(this.intervalId);
            }
            this.closeModal();
            this.$emit('run-workflow');
        },
        closeModal() {
            if (this.intervalId) {
                clearInterval(this.intervalId);
            }
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
}

.modal-container {
    background-color: #2c3e50;
    color: #ecf0f1;
    padding: 1rem;
    border-radius: 1rem;
    width: 480px;
    max-width: 90%;
    text-align: center;
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
</style>