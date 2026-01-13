<template>
    <div class="result__layout">
        <div class="result__files-container">
            <div class="files__header">
                <h3 class="files__title">Select Outputs</h3>
                <div class="files__actions">
                    <button class="action__button action__button--secondary" @click="toggleSelectAll"
                        :disabled="allFiles.length === 0">
                        {{ allSelected ? 'Deselect All' : 'Select All' }}
                    </button>
                    <button
                        :class="['action__button', setupSuccess ? 'action__button--success' : 'action__button--primary']"
                        @click="setFiles"
                        :disabled="selectedFiles.length === 0">
                        {{ setupSuccess ? '✓ Files Set Up' : 'Set up Files' }}
                    </button>
                    <button class="action__button action__button--download" @click="downloadSelectedFiles"
                        :disabled="selectedFiles.length === 0" v-if="selectedFiles.length > 1">
                        Download Selected ({{ selectedFiles.length }})
                    </button>
                </div>
            </div>

            <div class="files__content">
                <!-- Final Files Section -->
                <div class="files__section files__section--final" v-if="finalFileList.length > 0">
                <div class="section__header">
                    <div class="section__title">
                        <svg class="section__icon" xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <polyline points="9,11 12,14 22,4"></polyline>
                            <path d="m21,12v7a2,2 0 0,1 -2,2H5a2,2 0 0,1 -2,-2V5a2,2 0 0,1 2,-2h11"></path>
                        </svg>
                        <h4>Primary Outputs</h4>
                        <span class="section__count">({{ finalFileList.length }})</span>
                    </div>
                    <button class="action__button action__button--small" @click="toggleFinalSelectAll"
                        :disabled="finalFileList.length === 0">
                        {{ allFinalSelected ? 'Deselect All' : 'Select All' }}
                    </button>
                </div>

                <div class="files__list">
                    <div class="file__item file__item--final" v-for="file in finalFileList" :key="'final-' + file.name"
                        :class="{ 'file__item--selected': isFinalFileSelected(file.name) }">
                        <label class="file__checkbox-label">
                            <input type="checkbox" class="file__checkbox" :checked="isFinalFileSelected(file.name)"
                                @change="toggleFinalFileSelection(file.name)" />
                            <span class="file__checkbox-custom"></span>
                            <div class="file__info">
                                <span class="file__name">{{ file.name }}</span>
                                <span class="file__size">{{ formatFileSize(file.size) }}</span>
                            </div>
                        </label>
                        <button class="file__download-btn" @click="downloadIndividualFile(file)"
                            :data-tooltip="formatFileSize(file.size)" title="Download this file">
                            <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none"
                                stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M12 15V3m0 12l-4-4m4 4l4-4M2 17l.621 2.485A2 2 0 0 0 4.561 21h14.878a2 2 0 0 0 1.94-1.515L22 17"></path>
                            </svg>
                        </button>
                    </div>
                </div>
            </div>

                <!-- Intermediate Files Section -->
                <div class="files__section files__section--intermediate" v-if="intermediateFileList.length > 0">
                <div class="section__header section__header--expandable" @click="showIntermediateFiles = !showIntermediateFiles">
                    <div class="section__title">
                        <svg class="section__icon section__icon--expandable" :class="{ 'section__icon--expanded': showIntermediateFiles }" 
                             xmlns="http://www.w3.org/2000/svg" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <polyline points="6,9 12,15 18,9"></polyline>
                        </svg>
                        <h4>Intermediate Outputs</h4>
                        <span class="section__count">({{ intermediateFileList.length }})</span>
                    </div>
                    <button class="action__button action__button--small" @click.stop="toggleIntermediateSelectAll"
                        :disabled="intermediateFileList.length === 0" v-if="showIntermediateFiles">
                        {{ allIntermediateSelected ? 'Deselect All' : 'Select All' }}
                    </button>
                </div>

                <div class="files__list files__list--intermediate" v-show="showIntermediateFiles">
                    <div class="file__item file__item--intermediate" v-for="file in intermediateFileList" :key="'intermediate-' + file.name"
                        :class="{ 'file__item--selected': isIntermediateFileSelected(file.name) }">
                        <label class="file__checkbox-label">
                            <input type="checkbox" class="file__checkbox" :checked="isIntermediateFileSelected(file.name)"
                                @change="toggleIntermediateFileSelection(file.name)" />
                            <span class="file__checkbox-custom"></span>
                            <div class="file__info">
                                <span class="file__name">{{ file.name }}</span>
                                <span class="file__size">{{ formatFileSize(file.size) }}</span>
                            </div>
                        </label>
                        <button class="file__download-btn" @click="downloadIndividualFile(file)"
                            :data-tooltip="formatFileSize(file.size)" title="Download this file">
                            <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none"
                                stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                <path d="M12 15V3m0 12l-4-4m4 4l4-4M2 17l.621 2.485A2 2 0 0 0 4.561 21h14.878a2 2 0 0 0 1.94-1.515L22 17"></path>
                            </svg>
                        </button>
                    </div>
                </div>
                </div>

                <!-- Empty State -->
                <div class="files__empty" v-if="finalFileList.length === 0 && intermediateFileList.length === 0">
                    <p>No result files available</p>
                </div>
            </div>

            <!-- Footer -->
            <div class="files__footer">
                <div class="files__selection-info">
                    <span v-if="selectedFiles.length > 0">
                        {{ selectedFiles.length }} of {{ allFiles.length }} files selected
                        <span v-if="selectedFinalFiles.length > 0" class="selection__breakdown">
                            ({{ selectedFinalFiles.length }} final
                            <span v-if="selectedIntermediateFiles.length > 0">
                                + {{ selectedIntermediateFiles.length }} intermediate
                            </span>)
                        </span>
                        <span v-else-if="selectedIntermediateFiles.length > 0" class="selection__breakdown">
                            ({{ selectedIntermediateFiles.length }} intermediate)
                        </span>
                        <span v-if="hasBeenSetup" class="setup-success-message" style="margin-left: 0.75rem;">
                            ✓ {{ configuredFilesCount }} file{{ configuredFilesCount !== 1 ? 's' : '' }} configured for workflow
                        </span>
                    </span>
                    <span v-else-if="allFiles.length > 0">
                        No files selected
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
            finalFileList: [],
            intermediateFileList: [],
            selectedFinalFiles: [],
            selectedIntermediateFiles: [],
            isSetup: false,
            showIntermediateFiles: false,
            setupSuccess: false,
            setupSuccessTimeout: null,
            hasBeenSetup: false,
        };
    },
    async mounted() {
        await this.initializeComponent();
    },
    computed: {
        // 전체 파일 리스트 (최종 + 중간)
        allFiles() {
            return [...this.finalFileList, ...this.intermediateFileList];
        },
        // 전체 선택된 파일 리스트
        selectedFiles() {
            return [...this.selectedFinalFiles, ...this.selectedIntermediateFiles];
        },
        // 최종 파일 전체 선택 여부
        allFinalSelected() {
            return this.finalFileList.length > 0 && this.selectedFinalFiles.length === this.finalFileList.length;
        },
        // 중간 파일 전체 선택 여부
        allIntermediateSelected() {
            return this.intermediateFileList.length > 0 && this.selectedIntermediateFiles.length === this.intermediateFileList.length;
        },
        // 전체 파일 선택 여부
        allSelected() {
            return this.allFiles.length > 0 && this.selectedFiles.length === this.allFiles.length;
        },
        hasMultipleSelected() {
            return this.selectedFiles.length > 1;
        },
        // Vuex store에 저장된 실제 설정된 파일 개수
        configuredFilesCount() {
            const storedFiles = this.$store.getters.getSelectedWorkflowFiles(this.nodeId);
            return storedFiles ? storedFiles.length : 0;
        }
    },
    methods: {
        async initializeComponent() {
            try {
                await this.loadFileList();
                this.loadExistingSelection();
            } catch (error) {
                console.error('Error initializing ResultFiles component:', error);
            }
        },

        async loadFileList() {
            const current_node = this.$store.getters.getWorkflowNodeInfo(this.nodeId);

            if (!current_node.inputs.input_1.connections[0]) {
                throw new Error('No algorithm node connected');
            }

            this.algorithmId = current_node.inputs.input_1.connections[0].node;

            const workflow_result = {
                id: this.workflowId,
                algorithm_id: this.algorithmId,
            };

            const response = await getResults(workflow_result);
            const allFiles = response.data;

            // Get algorithm output configuration
            const algorithmNodeInfo = this.$store.getters.getWorkflowNodeInfo(this.algorithmId);
            const outputFiles = algorithmNodeInfo.data.selectedPluginInputOutput
                .filter(output => output.type === 'outputFile' && output.activate);

            // 최종 생성 파일 (기존 필터링 로직)
            this.finalFileList = allFiles.filter(file =>
                outputFiles.some(outputFile => outputFile.defaultValue === file.name)
            );

            // 중간 생성 파일 (전체 파일에서 최종 파일 제외)
            this.intermediateFileList = allFiles.filter(file =>
                !outputFiles.some(outputFile => outputFile.defaultValue === file.name)
            );
        },

        loadExistingSelection() {
            // Algorithm의 현재 executionId와 저장된 executionId 비교
            const algorithmNodeInfo = this.$store.getters.getWorkflowNodeInfo(this.algorithmId);
            const currentExecutionId = algorithmNodeInfo?.data?.lastExecutionId;
            const current_node = this.$store.getters.getWorkflowNodeInfo(this.nodeId);
            const configuredExecutionId = current_node?.data?.configuredForExecutionId;

            // executionId가 다르면 상태 초기화 (Algorithm 재실행 감지)
            if (currentExecutionId && configuredExecutionId && currentExecutionId !== configuredExecutionId) {
                console.info('ResultFiles: Algorithm re-executed, resetting state');
                this.$store.commit('setWorkflowFiles', { id: this.nodeId, files: [] });
                this.$store.commit('setWorkflowNodeDataObject', {
                    nodeId: this.nodeId,
                    dataObject: { configuredForExecutionId: null }
                });
                this.selectedFinalFiles = [];
                this.selectedIntermediateFiles = [];
                this.isSetup = false;
                this.hasBeenSetup = false;
                return;
            }

            // Check if this is already a multi-file node
            const existingFiles = this.$store.getters.getWorkflowNodeFilesInfo(this.nodeId);

            if (existingFiles && Array.isArray(existingFiles)) {
                // Multi-file format: restore selections
                const selectedFileNames = existingFiles
                    .filter(f => f.selected)
                    .map(f => f.name);

                // 현재 파일 목록에 존재하는 선택된 파일들만 필터링
                const validFinalFiles = selectedFileNames.filter(fileName =>
                    this.finalFileList.some(file => file.name === fileName)
                );
                const validIntermediateFiles = selectedFileNames.filter(fileName =>
                    this.intermediateFileList.some(file => file.name === fileName)
                );

                // Algorithm 재실행으로 인해 저장된 파일이 없어진 경우 상태 초기화
                const validFilesCount = validFinalFiles.length + validIntermediateFiles.length;
                if (selectedFileNames.length > 0 && validFilesCount === 0) {
                    // 모든 저장된 파일이 더 이상 존재하지 않음 - 상태 초기화
                    console.info('ResultFiles: Stored files no longer exist, resetting state');
                    this.$store.commit('setWorkflowFiles', { id: this.nodeId, files: [] });
                    this.selectedFinalFiles = [];
                    this.selectedIntermediateFiles = [];
                    this.isSetup = false;
                    this.hasBeenSetup = false;
                    return;
                }

                this.selectedFinalFiles = validFinalFiles;
                this.selectedIntermediateFiles = validIntermediateFiles;

                this.isSetup = validFilesCount > 0;
                this.hasBeenSetup = validFilesCount > 0;
            } else {
                // Handle legacy single file format
                const existingFile = this.$store.getters.getWorkflowNodeFileInfo(this.nodeId);
                if (existingFile && this.allFiles.some(file => file.name === existingFile)) {
                    // 기존 파일이 최종 파일인지 중간 파일인지 확인
                    if (this.finalFileList.some(file => file.name === existingFile)) {
                        this.selectedFinalFiles = [existingFile];
                    } else {
                        this.selectedIntermediateFiles = [existingFile];
                    }
                    this.isSetup = true;
                    this.hasBeenSetup = true;
                } else {
                    // Auto-select based on node title if available
                    const current_node = this.$store.getters.getWorkflowNodeInfo(this.nodeId);
                    const node_title = current_node.data['title'];
                    const matchingFile = this.allFiles.find(file => file.name === node_title);

                    if (matchingFile) {
                        // 매칭된 파일이 최종 파일인지 중간 파일인지 확인
                        if (this.finalFileList.some(file => file.name === matchingFile.name)) {
                            this.selectedFinalFiles = [matchingFile.name];
                        } else {
                            this.selectedIntermediateFiles = [matchingFile.name];
                        }
                        this.isSetup = true;
                        this.hasBeenSetup = true;
                    }
                }
            }
        },

        // Final files selection methods
        isFinalFileSelected(fileName) {
            return this.selectedFinalFiles.includes(fileName);
        },

        toggleFinalFileSelection(fileName) {
            const index = this.selectedFinalFiles.indexOf(fileName);
            if (index > -1) {
                this.selectedFinalFiles.splice(index, 1);
            } else {
                this.selectedFinalFiles.push(fileName);
            }
        },

        toggleFinalSelectAll() {
            if (this.allFinalSelected) {
                this.selectedFinalFiles = [];
            } else {
                this.selectedFinalFiles = this.finalFileList.map(file => file.name);
            }
        },

        // Intermediate files selection methods
        isIntermediateFileSelected(fileName) {
            return this.selectedIntermediateFiles.includes(fileName);
        },

        toggleIntermediateFileSelection(fileName) {
            const index = this.selectedIntermediateFiles.indexOf(fileName);
            if (index > -1) {
                this.selectedIntermediateFiles.splice(index, 1);
            } else {
                this.selectedIntermediateFiles.push(fileName);
            }
        },

        toggleIntermediateSelectAll() {
            if (this.allIntermediateSelected) {
                this.selectedIntermediateFiles = [];
            } else {
                this.selectedIntermediateFiles = this.intermediateFileList.map(file => file.name);
            }
        },

        // Global selection methods
        toggleSelectAll() {
            if (this.allSelected) {
                this.selectedFinalFiles = [];
                this.selectedIntermediateFiles = [];
            } else {
                this.selectedFinalFiles = this.finalFileList.map(file => file.name);
                this.selectedIntermediateFiles = this.intermediateFileList.map(file => file.name);
            }
        },

        async downloadIndividualFile(file) {
            try {
                const workflow_result = {
                    id: this.workflowId,
                    algorithm_id: this.algorithmId,
                    filename: file.name,
                };

                const response = await getResult(workflow_result);
                this.downloadBlob(response.data, file.name);
            } catch (error) {
                console.error('Error downloading individual file:', error);
            }
        },

        async downloadSelectedFiles() {
            if (this.selectedFiles.length === 0) return;

            try {
                // Cache computed properties to prevent reactivity issues during async operations
                const allFilesList = this.allFiles || [];
                const selectedFilesList = this.selectedFiles || [];

                if (selectedFilesList.length === 1) {
                    // Single file download
                    const file = allFilesList.find(f => f.name === selectedFilesList[0]);
                    if (file) {
                        await this.downloadIndividualFile(file);
                    }
                } else {
                    // Multiple files download (batch)
                    for (const fileName of selectedFilesList) {
                        const file = allFilesList.find(f => f.name === fileName);
                        if (file) {
                            await this.downloadIndividualFile(file);
                            // Add small delay between downloads to prevent overwhelming the server
                            await new Promise(resolve => setTimeout(resolve, 100));
                        }
                    }
                }
            } catch (error) {
                console.error('Error downloading selected files:', error);
            }
        },

        downloadBlob(data, filename) {
            const url = window.URL.createObjectURL(new Blob([data]));
            const link = document.createElement('a');
            link.href = url;
            link.setAttribute('download', filename);
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            window.URL.revokeObjectURL(url);
        },

        setFiles() {
            if (this.selectedFiles.length === 0) return;

            // Algorithm의 lastExecutionId 가져오기
            const algorithmNodeInfo = this.$store.getters.getWorkflowNodeInfo(this.algorithmId);
            const executionId = algorithmNodeInfo?.data?.lastExecutionId || null;

            // Convert to multi-file format for Vuex store
            const filesData = this.allFiles.map(file => ({
                name: file.name,
                size: file.size,
                selected: this.selectedFiles.includes(file.name)
            }));

            // Store in multi-file format
            this.$store.commit('setWorkflowFiles', {
                id: this.nodeId,
                files: filesData
            });

            // configuredForExecutionId 저장 (Algorithm 재실행 감지용)
            this.$store.commit('setWorkflowNodeDataObject', {
                nodeId: this.nodeId,
                dataObject: { configuredForExecutionId: executionId }
            });

            // Share files with connected nodes
            this.$store.commit('shareWorkflowFiles', this.nodeId);

            this.isSetup = true;
            this.hasBeenSetup = true;

            // Show success state
            if (this.setupSuccessTimeout) {
                clearTimeout(this.setupSuccessTimeout);
            }
            this.setupSuccess = true;
            this.setupSuccessTimeout = setTimeout(() => {
                this.setupSuccess = false;
            }, 3000);
        },

        formatFileSize(bytes) {
            if (bytes === 0) return '0 Bytes';
            const k = 1024;
            const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB'];
            const i = Math.floor(Math.log(bytes) / Math.log(k));
            return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
        },
    },
    beforeDestroy() {
        if (this.setupSuccessTimeout) {
            clearTimeout(this.setupSuccessTimeout);
        }
    }
};
</script>

<style scoped>
.result__layout {
    width: 100%;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    box-sizing: border-box;
    position: relative;
}

.result__files-container {
    background-color: rgb(255, 255, 255);
    padding: 1rem;
    border-radius: 0.5rem;
    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    width: 750px;
    max-width: calc(100% - 2rem);
    height: 600px;
    max-height: calc(100% - 2rem);
    display: flex;
    flex-direction: column;
    position: relative;
    box-sizing: border-box;
}

.files__header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 1rem;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid #e5e5e5;
    flex-shrink: 0;
}

.files__content {
    flex: 1;
    overflow-y: auto;
    margin-bottom: 1rem;
    padding-right: 4px;
}

/* Section Styles */
.files__section {
    margin-bottom: 1rem;
    border-radius: 8px;
    overflow: hidden;
}

.files__section--final {
    border: 1px solid #007bff8a;
    background: linear-gradient(135deg, #f8f9ff 0%, #ffffff 100%);
}

.files__section--intermediate {
    border: 1px solid #dee2e6;
    background: #f8f9fa;
}

.section__header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 12px 16px;
    background: rgba(255, 255, 255, 0.8);
    border-bottom: 1px solid rgba(0, 0, 0, 0.1);
}

.files__section--final .section__header {
    background: rgba(0, 123, 255, 0.05);
    border-bottom-color: rgba(0, 123, 255, 0.2);
}

.files__section--intermediate .section__header {
    background: rgba(108, 117, 125, 0.05);
    border-bottom-color: rgba(108, 117, 125, 0.2);
}

.section__header--expandable {
    cursor: pointer;
    user-select: none;
    transition: background-color 0.2s ease;
}

.section__header--expandable:hover {
    background: rgba(0, 0, 0, 0.05);
}

.section__title {
    display: flex;
    align-items: center;
    gap: 8px;
    font-weight: 600;
    color: #333;
}

.files__section--final .section__title {
    color: #0056b3;
}

.files__section--intermediate .section__title {
    color: #6c757d;
}

.section__title h4 {
    margin: 0;
    font-size: 1rem;
    font-weight: 600;
}

.section__count {
    font-size: 0.9rem;
    opacity: 0.7;
    margin-left: 4px;
}

.section__icon {
    flex-shrink: 0;
}

.section__icon--expandable {
    transition: transform 0.3s ease;
}

.section__icon--expanded {
    transform: rotate(180deg);
}

.files__title {
    margin: 0;
    color: #333;
    font-size: 1.2rem;
    font-weight: 600;
}

.files__actions {
    display: flex;
    gap: 0.5rem;
    flex-wrap: wrap;
    align-items: center;
}

.files__list {
    overflow-y: auto;
    max-height: 300px;
    padding: 0.5rem;
}

.file__item {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 0.75rem;
    border: 1px solid #e5e5e5;
    border-radius: 0.5rem;
    margin-bottom: 0.5rem;
    transition: all 0.2s ease;
    background: white;
    width: 100%;
    box-sizing: border-box;
    min-width: 0;
}

.file__item:last-child {
    margin-bottom: 0;
}

.file__item:hover {
    border-color: #007BFF;
    box-shadow: 0 2px 4px rgba(0, 123, 255, 0.1);
}

.file__item--selected {
    border-color: #007BFF;
    background: #f8f9ff;
}

.file__item--final {
    border-left: 4px solid #007BFF;
}

.file__item--final:hover {
    box-shadow: 0 3px 6px rgba(0, 123, 255, 0.15);
    transform: translateY(-1px);
}

.file__item--intermediate {
    border-left: 3px solid #6c757d;
    opacity: 0.9;
}

.file__item--intermediate:hover {
    opacity: 1;
    box-shadow: 0 2px 4px rgba(108, 117, 125, 0.15);
}

.files__list--intermediate {
    max-height: 200px;
    overflow-y: auto;
}

.file__checkbox-label {
    display: flex;
    align-items: center;
    flex: 1;
    cursor: pointer;
    gap: 0.75rem;
    min-width: 0;
    overflow: hidden;
}

.file__checkbox {
    display: none;
}

.file__checkbox-custom {
    width: 18px;
    height: 18px;
    border: 2px solid #ddd;
    border-radius: 3px;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: all 0.2s ease;
    flex-shrink: 0;
}

.file__checkbox:checked+.file__checkbox-custom {
    background: #007BFF;
    border-color: #007BFF;
}

.file__checkbox:checked+.file__checkbox-custom::after {
    content: '✓';
    color: white;
    font-size: 12px;
    font-weight: bold;
}

.file__info {
    display: flex;
    flex-direction: column;
    gap: 0.25rem;
    flex: 1;
    min-width: 0;
    overflow: hidden;
}

.file__name {
    font-weight: 500;
    color: #333;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    max-width: 100%;
}

.file__size {
    font-size: 0.85rem;
    color: #666;
}

.file__download-btn {
    background: #28a745;
    color: white;
    border: none;
    border-radius: 4px;
    padding: 0.5rem;
    cursor: pointer;
    transition: all 0.2s ease;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-shrink: 0;
}

.file__download-btn:hover {
    background: #218838;
    transform: translateY(-1px);
}

.file__download-btn:disabled {
    background: #ccc;
    cursor: not-allowed;
    transform: none;
}

.files__empty {
    text-align: center;
    padding: 2rem;
    color: #666;
    font-style: italic;
}

.files__footer {
    border-top: 1px solid #e5e5e5;
    padding-top: 1rem;
    display: flex;
    align-items: center;
    flex-shrink: 0;
    margin-top: auto;
}

.files__selection-info {
    font-size: 0.9rem;
    color: #666;
}

.action__button {
    padding: 0.5rem 1rem;
    border: none;
    border-radius: 0.375rem;
    font-size: 0.9rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s ease;
    display: flex;
    align-items: center;
    gap: 0.5rem;
}

.action__button:disabled {
    opacity: 0.5;
    cursor: not-allowed;
}

.action__button--primary {
    background: #007BFF;
    color: white;
}

.action__button--primary:hover:not(:disabled) {
    background: #0056b3;
}

.action__button--secondary {
    background: #f8f9fa;
    color: #333;
    border: 1px solid #ddd;
}

.action__button--secondary:hover:not(:disabled) {
    background: #e9ecef;
}

.action__button--download {
    background: #1ac951;
    color: white;
}

.action__button--download:hover:not(:disabled) {
    background: #0cb843;
}

.action__button--small {
    padding: 0.375rem 0.75rem;
    font-size: 0.8rem;
    border-radius: 0.25rem;
}

.selection__breakdown {
    font-size: 0.85rem;
    color: #6c757d;
    margin-left: 0.5rem;
}

/* Scrollbar styling for webkit browsers */
.files__list::-webkit-scrollbar {
    width: 6px;
}

.files__list::-webkit-scrollbar-track {
    background: #f1f1f1;
    border-radius: 3px;
}

.files__list::-webkit-scrollbar-thumb {
    background: #c1c1c1;
    border-radius: 3px;
}

.files__list::-webkit-scrollbar-thumb:hover {
    background: #a8a8a8;
}

/* Success state styles */
.action__button--success {
    background: #28a745 !important;
    color: white;
    transition: all 0.3s ease;
}

.action__button--success:hover:not(:disabled) {
    background: #218838 !important;
}

.setup-success-message {
    color: #28a745;
    font-weight: 600;
    animation: fadeIn 0.3s ease-in;
}

@keyframes fadeIn {
    from {
        opacity: 0;
        transform: translateY(-5px);
    }
    to {
        opacity: 1;
        transform: translateY(0);
    }
}
</style>