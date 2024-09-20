<template>
    <div class="modal-content">
        <h2 class="modal-title">Confirm Task Execution</h2>

        <div class="task-info">
            <h3 class="task-info__title">Task Information</h3>
            <div v-for="(task, index) in taskInfoList" :key="index" class="task-info__item">
                <div><strong>Selected Plugin:</strong> {{ task.pluginName }}</div>
                <div v-for="(input, inputIndex) in task.inputs" :key="inputIndex">{{ input }}</div>
                <div v-for="(output, outputIndex) in task.outputs" :key="outputIndex">{{ output }}</div>
            </div>
        </div>

        <div class="resource-info">
            <h3 class="resource-info__title">Server Resource Status</h3>

            <!-- CPU 사용량 표시 -->
            <div class="resource-bar">
                <label>CPU Usage: {{ serverResources.cpu_usage_percent }}%</label>
                <div class="bar">
                    <div class="fill" :style="{ width: serverResources.cpu_usage_percent + '%' }"></div>
                </div>
            </div>

            <!-- 메모리 사용량 표시 -->
            <div class="resource-bar">
                <label>Memory Usage: {{ serverResources.memory_usage_percent }}%</label>
                <div class="bar">
                    <div class="fill" :style="{ width: serverResources.memory_usage_percent + '%' }"></div>
                </div>
            </div>

            <!-- 총 메모리와 사용된 메모리, 가용 메모리 표시 -->
            <div class="resource-details">
                <p><strong>Total Memory:</strong> {{ formatBytes(serverResources.total_memory_bytes) }}</p>
                <p><strong>Used Memory:</strong> {{ formatBytes(serverResources.used_memory_bytes) }}</p>
                <p><strong>Available Memory:</strong> {{ formatBytes(serverResources.available_memory_bytes) }}</p>
            </div>
        </div>

        <!-- 실행/취소 버튼 -->
        <div class="modal-actions">
            <button class="btn confirm" @click="confirmTask">Execute</button>
            <button class="btn cancel" @click="closeModal">Cancel</button>
        </div>
    </div>
</template>

<script>
import { getSystemResources } from '@/api/index';

export default {
    data() {
        return {
            taskInfoList: [
                {
                    pluginName: "Plugin A",
                    inputs: ["Input 1", "Input 2"],
                    outputs: ["Output 1", "Output 2"]
                },
                {
                    pluginName: "Plugin B",
                    inputs: ["Input 3", "Input 4"],
                    outputs: ["Output 3", "Output 4"]
                }
            ],
            serverResources: {
                cpu_usage_percent: 0,
                memory_usage_percent: 0,
                total_memory_bytes: 0,
                used_memory_bytes: 0,
                available_memory_bytes: 0
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
        }, 1000);
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
            // 작업 실행 로직
            alert("Task is being executed...");
            this.closeModal();
            if (this.intervalId) {
                clearInterval(this.intervalId);
            }
        },
        closeModal() {
            this.$emit('deactivate-compile-check');
            if (this.intervalId) {
                clearInterval(this.intervalId);
            }
        }
    }
};
</script>

<style scoped>
/* 모달 콘텐츠 */
.modal-content {
    background-color: #2c3e50;
    color: #ecf0f1;
    padding: 1rem;
    border-radius: 1rem;
    width: 400px;
    max-width: 90%;
    text-align: center;

    position: absolute;
    top: calc(50% - 300px);
    left: calc(50% - 200px);
}

.modal-title {
    font-size: 1.5em;
    margin-bottom: 20px;
}

.task-info,
.resource-info {
    height: 200px;
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

.resource-bar {
    margin-bottom: 10px;
}

.bar {
    width: 100%;
    height: 10px;
    background-color: #34495e;
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
    padding: 10px 20px;
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
</style>