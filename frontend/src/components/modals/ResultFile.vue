<template>
    <div class="result__layout">
        <div class="result__download-container">
            <select v-model="selectedFile" class="result__file-dropdown" @change="setupFile = ''">
                <option value="">Please Select Result File</option>
                <option v-for="file in fileList" :key="file.name" :value="file">
                    {{ file.name }}
                </option>
            </select>

            <div class="set__button" v-if="setupFile == ''" @click="setFile">
                <div class="set__button--wrapper">
                    <div class="set__text">Set up File</div>
                </div>
            </div>

            <div class="result__button" v-else @click="downloadFile"
                :data-tooltip="selectedFile ? formatFileSize(selectedFile.size) : 'Select a file'"
                :disabled="!selectedFile">
                <div class="result__button--wrapper">
                    <div class="result__text">Download</div>
                    <span class="result__icon">
                        <svg xmlns="http://www.w3.org/2000/svg" aria-hidden="true" role="img" width="2em" height="2em"
                            preserveAspectRatio="xMidYMid meet" viewBox="0 0 24 24">
                            <path fill="none" stroke="currentColor" stroke-linecap="round" stroke-linejoin="round"
                                stroke-width="2"
                                d="M12 15V3m0 12l-4-4m4 4l4-4M2 17l.621 2.485A2 2 0 0 0 4.561 21h14.878a2 2 0 0 0 1.94-1.515L22 17">
                            </path>
                        </svg>
                    </span>
                </div>
            </div>
        </div>
    </div>
</template>

<script>
import { getResults, getResult } from "@/api/index";

export default {
    data() {
        return {
            workflowId: this.$route.query.workflow_id,
            nodeId: this.$route.query.node,
            algorithmId: null,
            fileList: [],
            selectedFile: "",
            setupFile: "",
        };
    },
    async mounted() {
        const current_node = this.$store.getters.getWorkflowNodeInfo(this.nodeId);
        console.log(current_node.inputs.input_1.connections[0].node);
        this.algorithmId = current_node.inputs.input_1.connections[0].node;
        const workflow_result = {
            id: this.workflowId,
            algorithm_id: this.algorithmId,
        }
        try {
            const response = await getResults(workflow_result);
            console.log(response);
            this.fileList = response.data;
        } catch (error) {
            console.log(error);
        }

        // if (this.fileList.length !== 0) {
        //     const algorithmNodeInfo = this.$store.getters.getWorkflowNodeInfo(this.algorithmId);
        //     // algorithmNodeInfo의 data.selectedPluginInputOutput들에서 type이 outputFile이고, activate가 true인 것들만 List로 추출
        //     const outputFiles = algorithmNodeInfo.data.selectedPluginInputOutput.filter(output => output.type === 'outputFile' && output.activate);

        //     console.log(outputFiles);

        //     // outputFiles 안에 요소들의 name이 fileList의 요소들의 name과 일치하는 것들만 fileList에 남기기
        //     this.fileList = this.fileList.filter(file => outputFiles.some(outputFile => outputFile.defaultValue === file.name));
        // }

        const setFileName = this.$store.getters.getWorkflowNodeFileInfo(this.nodeId);
        // this.fileList 안에 setFileName이 존재하면 selectedFile에 setFileName을 할당
        if (this.fileList.some(file => file.name === setFileName)) {
            this.selectedFile = setFileName;
            this.setupFile = setFileName;
        } else {
            const node_title = current_node.data['title'];
            this.selectedFile = this.fileList.find(file => file.name.includes(node_title));
            this.setupFile = setFileName;
            this.setFile();
        }
    },
    methods: {
        async downloadFile() {
            try {
                const workflow_result = {
                    id: this.workflowId,
                    algorithm_id: this.algorithmId,
                    filename: this.selectedFile.name,
                }
                const response = await getResult(workflow_result);
                console.log(response);
                // 파일을 다운로드하도록 설정
                const url = window.URL.createObjectURL(new Blob([response.data]));
                const link = document.createElement('a');
                link.href = url;
                link.setAttribute('download', this.selectedFile.name); // 파일 이름 설정
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);

                // 메모리 해제
                window.URL.revokeObjectURL(url);
            } catch (error) {
                console.error('Error downloading the file:', error);
            }
        },
        setFile() {
            this.setupFile = this.selectedFile;
            this.$store.commit('setWorkflowFile', {
                id: this.nodeId,
                file_name: this.selectedFile
            });
        },
        formatFileSize(bytes) {
            if (bytes === 0) return '0 Bytes';
            const k = 1024;
            const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB'];
            const i = Math.floor(Math.log(bytes) / Math.log(k));
            return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
        },
    }
}
</script>

<style>
.result__layout {
    width: 100%;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
}

.result__download-container {
    display: flex;
    align-items: center;
    justify-content: center;
    background-color: rgb(255, 255, 255);
    padding: 1rem;
    border-radius: 1rem;
}

.result__file-dropdown {
    margin-right: 10px;
    padding: 5px;
}

.set__button {
    --width: 100px;
    --height: 35px;
    --tooltip-height: 35px;
    --tooltip-width: 90px;
    --gap-between-tooltip-to-button: 18px;
    --button-color: #1ac951;
    --tooltip-color: #fff;
    width: var(--width);
    height: var(--height);
    background: var(--button-color);
    position: relative;
    text-align: center;
    border-radius: 0.45em;
    font-family: "Arial";
    transition: background 0.3s;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
}

.set__button:hover {
    background: #0cb843;
}

.set__button--wrapper {
    display: flex;
    align-items: center;
    justify-content: center;
}

.set__text {
    color: #fff;
}

/* 기존 button 스타일을 대체하는 커스텀 버튼 스타일 */
.result__button {
    --width: 100px;
    --height: 35px;
    --tooltip-height: 35px;
    --tooltip-width: 90px;
    --gap-between-tooltip-to-button: 18px;
    --button-color: #007BFF;
    --tooltip-color: #fff;
    width: var(--width);
    height: var(--height);
    background: var(--button-color);
    position: relative;
    text-align: center;
    border-radius: 0.45em;
    font-family: "Arial";
    transition: background 0.3s;
    cursor: pointer;
}

.result__button::before {
    position: absolute;
    content: attr(data-tooltip);
    width: var(--tooltip-width);
    height: var(--tooltip-height);
    background-color: var(--tooltip-color);
    font-size: 0.9rem;
    color: #111;
    border-radius: .25em;
    line-height: var(--tooltip-height);
    bottom: calc(var(--height) + var(--gap-between-tooltip-to-button) + 10px);
    left: calc(50% - var(--tooltip-width) / 2);
}

.result__button::after {
    position: absolute;
    content: '';
    width: 0;
    height: 0;
    border: 10px solid transparent;
    border-top-color: var(--tooltip-color);
    left: calc(50% - 10px);
    bottom: calc(100% + var(--gap-between-tooltip-to-button) - 10px);
}

.result__button::after,
.result__button::before {
    opacity: 0;
    visibility: hidden;
    transition: all 0.5s;
}

.result__text,
.set__text {
    display: flex;
    align-items: center;
    justify-content: center;
}

.result__button-wrapper,
.result__text,
.result__icon {
    overflow: hidden;
    position: absolute;
    width: 100%;
    height: 100%;
    left: 0;
    color: #fff;
}

.result__text {
    top: 0;
    opacity: 1;
    visibility: visible;
}

.result__text,
.result__icon {
    transition: all 0.5s;
}

.result__icon {
    color: #fff;
    top: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    opacity: 0;
    visibility: hidden;
}

.result__icon svg {
    width: 24px;
    height: 24px;
}

.result__button:hover {
    background: #6200ff;
}

.result__button:hover .result__text {
    top: -100%;
    opacity: 0;
    visibility: hidden;
}

.result__button:hover .result__icon {
    top: 0;
    opacity: 1;
    visibility: visible;
}

.result__button:hover:before,
.result__button:hover:after {
    opacity: 1;
    visibility: visible;
}

.result__button:hover:after {
    bottom: calc(var(--height) + var(--gap-between-tooltip-to-button) - 20px);
}

.result__button:hover:before {
    bottom: calc(var(--height) + var(--gap-between-tooltip-to-button));
}
</style>