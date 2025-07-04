<template>
    <div class="validation-container" v-on:click.stop="validationLoading">
        <!-- Header Section -->
        <div class="validation-header">
            <div class="header-content">
                <div class="header-icon">
                    <i class="fas fa-check-circle"></i>
                </div>
                <div class="header-text">
                    <h2>Plugin Validation & Upload</h2>
                    <p>Please review and upload your plugin</p>
                </div>
            </div>
        </div>

        <!-- Plugin Information Card -->
        <div class="info-card plugin-info-card">
            <div class="card-header">
                <div class="card-title">
                    <i class="fas fa-info-circle card-icon"></i>
                    <h3>Plugin Information</h3>
                </div>
                <div class="card-badge">Basic Information</div>
            </div>
            <div class="card-content">
                <div class="info-grid">
                    <div class="info-item">
                        <label>Plugin Name</label>
                        <div class="info-value plugin-name">{{ plugin.name }}</div>
                    </div>
                    <div class="info-item description">
                        <label>Description</label>
                        <div class="info-value">{{ plugin.description }}</div>
                    </div>
                </div>

                <!-- Dependency Files -->
                <div v-if="plugin.dependencyFiles && plugin.dependencyFiles.length > 0" class="subsection">
                    <div class="subsection-header">
                        <i class="fas fa-file-code"></i>
                        <span>Dependency Files</span>
                        <div class="count-badge">{{ plugin.dependencyFiles.length }}</div>
                    </div>
                    <div class="file-list">
                        <div v-for="(dependency, index) in plugin.dependencyFiles" :key="index"
                            class="file-item dependency-file">
                            <div class="file-icon">
                                <i class="fas fa-file-alt"></i>
                            </div>
                            <span class="file-name">{{ dependency.type }}</span>
                        </div>
                    </div>
                </div>

                <!-- Reference Folders -->
                <div v-if="plugin.referenceFolders && plugin.referenceFolders.length > 0" class="subsection">
                    <div class="subsection-header">
                        <i class="fas fa-folder-open"></i>
                        <span>Reference Folders</span>
                        <div class="count-badge">{{ plugin.referenceFolders.length }}</div>
                    </div>
                    <div class="folder-tree">
                        <div v-for="(folder, index) in plugin.referenceFolders" :key="index" class="folder-item">
                            <div class="folder-header">
                                <i class="fas fa-folder"></i>
                                <span class="folder-name">{{ folder.folderName }}</span>
                                <div class="file-count">{{ folder.files.length }} files</div>
                            </div>
                            <div v-if="folder.files.length > 0" class="folder-files">
                                <div v-for="(file, i) in folder.files" :key="i" class="file-item script-file">
                                    <div class="file-icon">
                                        <i class="fas fa-file-code"></i>
                                    </div>
                                    <span class="file-name">{{ file.name }}</span>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Rules Section -->
        <div class="info-card rules-card">
            <div class="card-header">
                <div class="card-title">
                    <i class="fas fa-cogs card-icon"></i>
                    <h3>Workflow Rules</h3>
                </div>
                <div class="card-badge">{{ rules.length }} rules</div>
            </div>
            <div class="card-content">
                <div v-if="rules.length === 0" class="empty-state">
                    <i class="fas fa-exclamation-triangle"></i>
                    <p>No rules defined</p>
                </div>
                <div v-else class="rules-grid">
                    <div v-for="(rule, index) in rules" :key="index" class="rule-card">
                        <!-- Rule Header -->
                        <div class="rule-header">
                            <div class="rule-title-section">
                                <div class="rule-icon">
                                    <i class="fas fa-play"></i>
                                </div>
                                <div class="rule-info">
                                    <h4>{{ rule.name }}</h4>
                                    <div class="rule-meta">
                                        <span class="rule-id">ID: {{ rule.nodeId }}</span>
                                        <span v-if="rule.isVisualization" class="visualization-badge">
                                            <i class="fas fa-chart-bar"></i> Visualization
                                        </span>
                                    </div>
                                </div>
                            </div>
                        </div>

                        <!-- Rule Content -->
                        <div class="rule-content">
                            <!-- Input Files -->
                            <div class="rule-section">
                                <div class="section-header">
                                    <i class="fas fa-arrow-right"></i>
                                    <span>Input Files</span>
                                    <div class="count-badge">{{ rule.input.length }}</div>
                                </div>
                                <div class="file-list">
                                    <div v-for="(input, iIndex) in rule.input" :key="iIndex"
                                        class="file-item input-file">
                                        <div class="file-icon">
                                            <i class="fas fa-file-import"></i>
                                        </div>
                                        <span class="file-name">{{ input }}</span>
                                        <span class="file-extension">{{ getFileExtension(input) }}</span>
                                    </div>
                                </div>
                            </div>

                            <!-- Output Files -->
                            <div class="rule-section">
                                <div class="section-header">
                                    <i class="fas fa-arrow-left"></i>
                                    <span>Output Files</span>
                                    <div class="count-badge">{{ rule.output.length }}</div>
                                </div>
                                <div class="file-list">
                                    <div v-for="(output, oIndex) in rule.output" :key="oIndex"
                                        class="file-item output-file">
                                        <div class="file-icon">
                                            <i class="fas fa-file-export"></i>
                                        </div>
                                        <span class="file-name">{{ output }}</span>
                                        <span class="file-extension">{{ getFileExtension(output) }}</span>
                                    </div>
                                </div>
                            </div>

                            <!-- Script Information -->
                            <div v-if="rule.script" class="rule-section">
                                <div class="section-header">
                                    <i class="fas fa-code"></i>
                                    <span>Script</span>
                                </div>
                                <div class="script-info">
                                    <div class="script-item">
                                        <div class="script-icon">
                                            <i class="fas fa-file-code"></i>
                                        </div>
                                        <span class="script-name">{{ getScriptName(rule.script) }}</span>
                                    </div>
                                    <div class="command-preview">
                                        <code>{{ generateShellCommand(rule) }}</code>
                                    </div>
                                </div>
                            </div>

                            <!-- Parameters -->
                            <div v-if="rule.parameters && rule.parameters.length > 0" class="rule-section">
                                <div class="section-header">
                                    <i class="fas fa-sliders-h"></i>
                                    <span>Parameters</span>
                                    <div class="count-badge">{{ rule.parameters.length }}</div>
                                </div>
                                <div class="parameters-grid">
                                    <div v-for="(param, pIndex) in rule.parameters" :key="pIndex" class="param-item">
                                        <div class="param-header">
                                            <span class="param-name">{{ param.name }}</span>
                                            <span class="param-type" :class="'type-' + param.type">{{ param.type
                                            }}</span>
                                        </div>
                                        <div v-if="param.defaultValue" class="param-value">{{ param.defaultValue }}
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Action Section -->
        <div class="action-section">
            <div class="action-content">
                <div class="confirmation-checklist">
                    <h4>Pre-upload Checklist</h4>
                    <ul>
                        <li><i class="fas fa-check"></i> Verified that all scripts meet dependencies</li>
                        <li><i class="fas fa-check"></i> Verified that shell commands execute correctly</li>
                        <li><i class="fas fa-check"></i> Verified that scripts function properly</li>
                        <li><i class="fas fa-check"></i> Verified that final output is correctly visualized</li>
                    </ul>
                </div>
                <div class="upload-section">
                    <!-- GPU Configuration -->
                    <div class="gpu-config-option">
                        <label class="gpu-toggle-label">
                            <input type="checkbox" class="gpu-toggle-input" v-model="useGpu">
                            <span class="gpu-toggle-slider"></span>
                            <div class="gpu-toggle-text">
                                <span class="gpu-toggle-title">Enable GPU Support</span>
                                <span class="gpu-toggle-desc">Use GPU acceleration for this plugin</span>
                            </div>
                        </label>
                    </div>
                    
                    <button @click="uploadPluginData" class="upload-button">
                        <i class="fas fa-cloud-upload-alt"></i>
                        <span>Upload Plugin</span>
                    </button>
                </div>
            </div>
        </div>

        <!-- Loading Animation -->
        <div v-if="validationLoading" class="loading-overlay">
            <div class="loading-content">
                <div class="loader"></div>
                <div class="loading-text">
                    <h3 v-if="uploadingStep == 1">Validating plugin...</h3>
                    <h3 v-else-if="uploadingStep == 2">Uploading plugin data...</h3>
                    <h3 v-else-if="uploadingStep == 3">Building plugin environment...</h3>
                    <h3 v-else-if="uploadingStep == 4">Uploading plugin scripts...</h3>
                    <p>Please wait a moment</p>
                </div>
            </div>
        </div>
    </div>
</template>

<script>
import { validationPlugin, uploadPluginMetadata, uploadPluginScripts, uploadPluginPackage, buildPluginDocker, syncPluginData } from '@/api/index';

export default {
    props: {
        plugin: {
            type: Object,
            required: true
        },
        rules: {
            type: Array,
            required: true
        },
        drawflow: {
            type: Object,
            required: true
        }
    },
    data() {
        return {
            processedPlugin: { ...this.plugin },
            validationLoading: false,
            uploadingStep: 0,
            buildTimer: null,
            useGpu: false,
        };
    },
    async mounted() {
        await this.processDependencyFiles();
        await this.processReferenceFolderFiles();
    },
    beforeDestroy() {
        // 컴포넌트 제거 시 타이머 정리
        if (this.buildTimer) {
            clearTimeout(this.buildTimer);
        }
    },
    methods: {
        async processDependencyFiles() {
            const promises = this.processedPlugin.dependencyFiles.map(async (dependency) => {
                if (dependency.file instanceof File) {
                    dependency.file = await this.readFileContent(dependency.file);
                }
            });
            await Promise.all(promises);
        },
        async processReferenceFolderFiles() {
            const processFolder = async (folder) => {
                const filePromises = folder.files.map(async (file, idx) => {
                    if (file.file instanceof File) {
                        folder.files[idx] = {
                            file: await this.readFileContent(file.file),
                            fileName: file.file.name,
                            type: file.file.type,
                        };
                    }
                });
                await Promise.all(filePromises);

                const subFolderPromises = folder.subFolders.map((subFolder) => processFolder(subFolder));
                await Promise.all(subFolderPromises);
            };

            const folderPromises = this.plugin.referenceFolders.map((folder) => processFolder(folder));
            await Promise.all(folderPromises);
        },
        readFileContent(file) {
            return new Promise((resolve, reject) => {
                const reader = new FileReader();
                reader.onload = () => resolve(reader.result);
                reader.onerror = reject;
                reader.readAsText(file);
            });
        },
        async uploadPluginData() {
            const confirmationMessage = `
                Have you confirmed the following:
                - All uploaded scripts meet the dependencies?
                - The complete shell command can correctly execute the specified script?
                - The script works correctly?
                - The final output is visualized correctly?
                `;

            const userConfirmed = confirm(confirmationMessage);

            if (userConfirmed) {
                try {
                    this.validationLoading = true;
                    this.uploadingStep = 1;

                    // 1. 규칙 데이터 준비
                    const rules = this.rules.map(rule => ({
                        name: rule.name,
                        input: rule.input,
                        output: rule.output,
                        script: rule.script.name,
                        parameters: rule.parameters,
                        nodeId: rule.nodeId,
                        isVisualization: rule.isVisualization
                    }));

                    // 2. 플러그인 데이터 준비
                    const plugin = {
                        name: this.plugin.name,
                        description: this.plugin.description,
                        dependencyFiles: this.plugin.dependencyFiles,
                        referenceFolders: this.plugin.referenceFolders,
                        useGpu: this.useGpu
                    };

                    // 3. 플러그인 검증
                    const validation_response = await validationPlugin(plugin, rules, this.drawflow);
                    console.log('Validation Response:', validation_response);

                    // 4. 메타데이터 업로드
                    this.uploadingStep = 2;
                    const pluginCreate = validation_response.data.plugin;
                    const scriptNameList = validation_response.data.scripts;

                    const db_plugin = await uploadPluginMetadata(pluginCreate);
                    console.log("Plugin Metadata Uploaded:", db_plugin.data);

                    // 5. 파일 업로드 준비
                    const prepareFiles = async () => {
                        try {
                            // 스크립트 파일 업로드 준비
                            const scriptFormData = new FormData();
                            scriptFormData.append('plugin_name', pluginCreate.name);

                            // 스크립트 파일 처리
                            await Promise.all(this.rules
                                .filter(rule => rule.script && scriptNameList.includes(rule.script.name))
                                .map(async (rule) => {
                                    try {
                                        const fileContent = await this.readFileAsText(rule.script);
                                        const file = new Blob([fileContent], { type: rule.script.type });
                                        scriptFormData.append('files', file, rule.script.name);
                                    } catch (error) {
                                        console.error(`Error processing script ${rule.script.name}:`, error);
                                        throw error;
                                    }
                                })
                            );

                            // 패키지 파일 업로드 준비
                            const packageFormData = this.plugin.packageFiles?.length > 0 ? new FormData() : null;
                            if (packageFormData) {
                                packageFormData.append("plugin_name", pluginCreate.name);
                                await Promise.all(this.plugin.packageFiles.map(async (packageFile) => {
                                    try {
                                        if (!packageFile.file) return;
                                        const fileBlob = new Blob([await packageFile.file.arrayBuffer()], {
                                            type: packageFile.file.type
                                        });
                                        packageFormData.append("files", fileBlob, packageFile.fileName);
                                    } catch (error) {
                                        console.error(`Error processing package file ${packageFile.fileName}:`, error);
                                        throw error;
                                    }
                                }));
                            }

                            return { scriptFormData, packageFormData };
                        } catch (error) {
                            console.error('Error preparing files:', error);
                            throw error;
                        }
                    };

                    const preparedFiles = await prepareFiles();

                    // 6. 스크립트 및 패키지 파일 업로드
                    this.uploadingStep = 4;
                    const uploadPromises = [];

                    if (preparedFiles.scriptFormData) {
                        uploadPromises.push(
                            uploadPluginScripts(preparedFiles.scriptFormData)
                                .then(() => console.log('Scripts uploaded successfully'))
                                .catch(error => {
                                    console.error('Error uploading scripts:', error);
                                    throw error;
                                })
                        );
                    }

                    if (preparedFiles.packageFormData) {
                        uploadPromises.push(
                            uploadPluginPackage(preparedFiles.packageFormData)
                                .then(() => console.log('Packages uploaded successfully'))
                                .catch(error => {
                                    console.error('Error uploading packages:', error);
                                    throw error;
                                })
                        );
                    }

                    // 모든 파일 업로드 완료 대기
                    await Promise.all(uploadPromises);
                    console.log('All files uploaded successfully');

                    // 7. DB 동기화 (파일 업로드 후 DB 정보 최신화)
                    // 스크립트나 텍스트 의존성 파일 변경사항을 DB에 반영
                    try {
                        console.log('Starting DB synchronization...');
                        const syncResponse = await syncPluginData(pluginCreate.name);
                        console.log('DB synchronization completed:', syncResponse.data);

                        if (syncResponse.data.updated_fields) {
                            const fields = syncResponse.data.updated_fields;
                            console.log(`Synchronized - Dependencies: ${fields.dependencies}, Rules: ${fields.rules}, Drawflow nodes: ${fields.drawflow_nodes}`);
                        }
                    } catch (error) {
                        console.warn('DB synchronization failed (but continuing):', error);
                        // 동기화 실패는 경고만 출력하고 계속 진행
                        // 수동으로 동기화할 수 있도록 안내 메시지 추가 가능
                    }

                    // 8. Docker 이미지 빌드 (마지막 단계)
                    this.uploadingStep = 3;
                    console.log('Starting Docker build...');

                    const buildResponse = await buildPluginDocker(pluginCreate.name, this.useGpu);
                    console.log('Docker build completed:', buildResponse.data);

                    // 업로드 완료 처리
                    this.validationLoading = false;
                    alert('Plugin uploaded and built successfully!');
                    this.$emit('close');

                } catch (error) {
                    console.error('Error while uploading plugin:', error);
                    this.validationLoading = false;
                    this.uploadingStep = 0;
                    if (this.buildTimer) {
                        clearTimeout(this.buildTimer);
                    }
                    alert(error.response?.data?.detail?.message || error.response?.data?.detail || 'Error while uploading plugin. Please try again.');
                }
            } else {
                alert('Please verify the plugin details and try again.');
            }
        },
        generateShellCommand(rule) {
            const paramStr = rule.parameters.map(p => {
                if (p.type === 'inputFile' || p.type === 'outputFile') {
                    return `${p.defaultValue}(${p.type})`;
                } else {
                    return `${p.name}(${p.type}:${p.defaultValue})`;
                }
            }).join(' ');

            const scriptName = rule.script.name || rule.script;
            let command = '';
            if (scriptName.endsWith('.py')) {
                command = `/python ${scriptName}`;
            } else if (scriptName.endsWith('.R')) {
                command = `/Rscript ${scriptName}`;
            } else {
                command = `/${scriptName}`;
            }

            return `${command} ${paramStr}`;
        },
        getFileExtension(filename) {
            if (!filename) return '';
            const parts = filename.split('.');
            return parts.length > 1 ? '.' + parts[parts.length - 1] : '';
        },
        getScriptName(script) {
            if (!script) return '';
            if (typeof script === 'string') return script;
            return script.name || 'Unknown script';
        },
        // 파일을 텍스트로 읽는 헬퍼 함수
        readFileAsText(file) {
            return new Promise((resolve, reject) => {
                const reader = new FileReader();
                reader.onload = (event) => resolve(event.target.result);
                reader.onerror = (error) => reject(error);
                reader.readAsText(file);
            });
        },
    },
};
</script>

<style scoped>
.validation-container {
    width: 100%;
    max-width: 1200px;
    margin: 0 auto;
    padding: 2rem;
    background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    min-height: 100vh;
    box-sizing: border-box;
}

/* Header Section */
.validation-header {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 2rem;
    border-radius: 16px;
    margin-bottom: 2rem;
    box-shadow: 0 8px 32px rgba(102, 126, 234, 0.3);
}

.header-content {
    display: flex;
    align-items: center;
    gap: 1.5rem;
}

.header-icon {
    width: 4rem;
    height: 4rem;
    background: rgba(255, 255, 255, 0.2);
    border-radius: 12px;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 1.5rem;
}

.header-text h2 {
    margin: 0 0 0.5rem 0;
    font-size: 2rem;
    font-weight: 600;
}

.header-text p {
    margin: 0;
    opacity: 0.9;
    font-size: 1.125rem;
}

/* Card Styles */
.info-card {
    background: white;
    border-radius: 16px;
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1);
    margin-bottom: 2rem;
    overflow: hidden;
    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
}

.info-card:hover {
    box-shadow: 0 8px 32px rgba(0, 0, 0, 0.15);
    transform: translateY(-2px);
}

.card-header {
    background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    padding: 1.5rem 2rem;
    display: flex;
    justify-content: space-between;
    align-items: center;
    border-bottom: 1px solid #e1e5e9;
}

.card-title {
    display: flex;
    align-items: center;
    gap: 0.75rem;
}

.card-icon {
    color: #667eea;
    font-size: 1.25rem;
}

.card-title h3 {
    margin: 0;
    font-size: 1.5rem;
    font-weight: 600;
    color: #2c3e50;
}

.card-badge {
    background: #667eea;
    color: white;
    padding: 0.5rem 1rem;
    border-radius: 20px;
    font-size: 0.875rem;
    font-weight: 500;
}

.card-content {
    padding: 2rem;
}

/* Plugin Information Styles */
.info-grid {
    display: grid;
    gap: 1.5rem;
    margin-bottom: 2rem;
}

.info-item {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
}

.info-item.description {
    grid-column: 1 / -1;
}

.info-item label {
    font-size: 0.875rem;
    font-weight: 600;
    color: #6c757d;
    text-transform: uppercase;
    letter-spacing: 0.025em;
}

.info-value {
    font-size: 1rem;
    color: #2c3e50;
    padding: 0.75rem;
    background: #f8f9fa;
    border-radius: 8px;
    border-left: 4px solid #667eea;
}

.info-value.plugin-name {
    font-weight: 600;
    font-size: 1.25rem;
    background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%);
}

/* Subsection Styles */
.subsection {
    margin-top: 2rem;
    padding-top: 2rem;
    border-top: 1px solid #e1e5e9;
}

.subsection-header {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    margin-bottom: 1rem;
    font-weight: 600;
    color: #2c3e50;
}

.subsection-header i {
    color: #667eea;
    font-size: 1.125rem;
}

.count-badge {
    background: #e3f2fd;
    color: #1565c0;
    padding: 0.25rem 0.75rem;
    border-radius: 12px;
    font-size: 0.75rem;
    font-weight: 600;
}

/* File List Styles */
.file-list {
    display: flex;
    flex-direction: column;
    gap: 0.75rem;
}

.file-item {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    padding: 0.75rem;
    background: #f8f9fa;
    border-radius: 8px;
    border-left: 4px solid #28a745;
    transition: all 0.2s ease;
}

.file-item:hover {
    background: #e9ecef;
    transform: translateX(4px);
}

.file-item.dependency-file {
    border-left-color: #ffc107;
}

.file-item.script-file {
    border-left-color: #17a2b8;
}

.file-item.input-file {
    border-left-color: #28a745;
}

.file-item.output-file {
    border-left-color: #dc3545;
}

.file-icon {
    width: 2rem;
    height: 2rem;
    background: #28a745;
    color: white;
    border-radius: 6px;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 0.875rem;
}

.dependency-file .file-icon {
    background: #ffc107;
}

.script-file .file-icon {
    background: #17a2b8;
}

.input-file .file-icon {
    background: #28a745;
}

.output-file .file-icon {
    background: #dc3545;
}

.file-name {
    font-weight: 500;
    color: #2c3e50;
    flex: 1;
}

.file-extension {
    font-size: 0.75rem;
    color: #6c757d;
    background: #e9ecef;
    padding: 0.25rem 0.5rem;
    border-radius: 8px;
    font-weight: 500;
}

/* Folder Tree Styles */
.folder-tree {
    display: flex;
    flex-direction: column;
    gap: 1rem;
}

.folder-item {
    background: #f8f9fa;
    border-radius: 8px;
    padding: 1rem;
    border-left: 4px solid #6f42c1;
}

.folder-header {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    margin-bottom: 0.75rem;
    font-weight: 600;
    color: #2c3e50;
}

.folder-header i {
    color: #6f42c1;
    font-size: 1.125rem;
}

.folder-name {
    flex: 1;
}

.file-count {
    background: #e3f2fd;
    color: #1565c0;
    padding: 0.25rem 0.5rem;
    border-radius: 8px;
    font-size: 0.75rem;
    font-weight: 500;
}

.folder-files {
    padding-left: 1.5rem;
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
}

/* Rules Grid */
.rules-grid {
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
}

.rule-card {
    background: white;
    border-radius: 12px;
    box-shadow: 0 2px 12px rgba(0, 0, 0, 0.08);
    border: 1px solid #e1e5e9;
    overflow: hidden;
    transition: all 0.3s ease;
}

.rule-card:hover {
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.12);
    transform: translateY(-2px);
}

.rule-header {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 1.5rem;
}

.rule-title-section {
    display: flex;
    align-items: center;
    gap: 1rem;
}

.rule-icon {
    width: 2.5rem;
    height: 2.5rem;
    background: rgba(255, 255, 255, 0.2);
    border-radius: 8px;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 1rem;
}

.rule-info h4 {
    margin: 0 0 0.5rem 0;
    font-size: 1.25rem;
    font-weight: 600;
}

.rule-meta {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    opacity: 0.9;
}

.rule-id {
    background: rgba(255, 255, 255, 0.2);
    padding: 0.25rem 0.5rem;
    border-radius: 8px;
    font-size: 0.75rem;
    font-weight: 500;
}

.visualization-badge {
    background: #ff6b6b;
    color: white;
    padding: 0.25rem 0.75rem;
    border-radius: 12px;
    font-size: 0.75rem;
    font-weight: 600;
    display: flex;
    align-items: center;
    gap: 0.25rem;
}

.rule-content {
    padding: 1.5rem;
}

.rule-section {
    margin-bottom: 1.5rem;
}

.rule-section:last-child {
    margin-bottom: 0;
}

.section-header {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    margin-bottom: 1rem;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid #e1e5e9;
    font-weight: 600;
    color: #2c3e50;
}

.section-header i {
    color: #667eea;
    font-size: 1rem;
}

/* Script Information */
.script-info {
    display: flex;
    flex-direction: column;
    gap: 1rem;
}

.script-item {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    padding: 0.75rem;
    background: #f8f9fa;
    border-radius: 8px;
    border-left: 4px solid #17a2b8;
}

.script-icon {
    width: 2rem;
    height: 2rem;
    background: #17a2b8;
    color: white;
    border-radius: 6px;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 0.875rem;
}

.script-name {
    font-weight: 500;
    color: #2c3e50;
}

.command-preview {
    padding: 0.75rem;
    background: #2c3e50;
    border-radius: 8px;
    overflow-x: auto;
}

.command-preview code {
    color: #ecf0f1;
    font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
    font-size: 0.875rem;
    white-space: nowrap;
}

/* Parameters Grid */
.parameters-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
    gap: 0.75rem;
}

.param-item {
    padding: 0.75rem;
    background: #f8f9fa;
    border-radius: 8px;
    border-left: 4px solid #6f42c1;
}

.param-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 0.5rem;
}

.param-name {
    font-weight: 600;
    color: #2c3e50;
    font-size: 0.875rem;
}

.param-type {
    font-size: 0.75rem;
    padding: 0.25rem 0.5rem;
    border-radius: 12px;
    font-weight: 500;
    text-transform: uppercase;
    letter-spacing: 0.025em;
}

.type-inputFile {
    background-color: #e3f2fd;
    color: #1565c0;
}

.type-outputFile {
    background-color: #f3e5f5;
    color: #7b1fa2;
}

.type-optionalInputFile {
    background-color: #e8f5e8;
    color: #2e7d32;
}

.type-string {
    background-color: #fff3e0;
    color: #f57c00;
}

.type-int,
.type-float {
    background-color: #fce4ec;
    color: #c2185b;
}

.type-boolean {
    background-color: #e0f2f1;
    color: #00695c;
}

.type-h5adParameter {
    background-color: #f1f8e9;
    color: #33691e;
}

.param-value {
    font-size: 0.875rem;
    color: #6c757d;
    font-weight: 500;
}

/* Empty State */
.empty-state {
    text-align: center;
    padding: 3rem 2rem;
    color: #6c757d;
}

.empty-state i {
    font-size: 3rem;
    margin-bottom: 1rem;
    color: #dee2e6;
}

.empty-state p {
    margin: 0;
    font-size: 1.125rem;
    font-weight: 500;
}

/* Action Section */
.action-section {
    background: white;
    border-radius: 16px;
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1);
    padding: 2rem;
}

.action-content {
    display: flex;
    justify-content: space-between;
    align-items: flex-start;
    gap: 2rem;
}

.confirmation-checklist {
    flex: 1;
}

.confirmation-checklist h4 {
    margin: 0 0 1rem 0;
    font-size: 1.125rem;
    color: #2c3e50;
}

.confirmation-checklist ul {
    list-style: none;
    padding: 0;
    margin: 0;
}

.confirmation-checklist li {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    padding: 0.5rem 0;
    color: #495057;
    font-size: 0.875rem;
}

.confirmation-checklist li i {
    color: #28a745;
    font-size: 1rem;
}

.upload-button {
    background: linear-gradient(135deg, #28a745 0%, #20c997 100%);
    color: white;
    border: none;
    padding: 1rem 2rem;
    border-radius: 12px;
    cursor: pointer;
    font-size: 1.125rem;
    font-weight: 600;
    display: flex;
    align-items: center;
    gap: 0.75rem;
    transition: all 0.3s ease;
    box-shadow: 0 4px 16px rgba(40, 167, 69, 0.3);
    min-width: 200px;
    justify-content: center;
}

.upload-button:hover {
    background: linear-gradient(135deg, #218838 0%, #1ea085 100%);
    transform: translateY(-2px);
    box-shadow: 0 6px 24px rgba(40, 167, 69, 0.4);
}

.upload-button:active {
    transform: translateY(0);
}

/* Loading Overlay */
.loading-overlay {
    position: fixed;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    background: rgba(0, 0, 0, 0.8);
    display: flex;
    align-items: center;
    justify-content: center;
    z-index: 9999;
    backdrop-filter: blur(10px);
}

.loading-content {
    background: white;
    padding: 3rem;
    border-radius: 16px;
    box-shadow: 0 8px 32px rgba(0, 0, 0, 0.3);
    text-align: center;
    max-width: 400px;
    width: 90%;
}

.loader {
    width: 64px;
    height: 64px;
    border: 4px solid #e1e5e9;
    border-top: 4px solid #667eea;
    border-radius: 50%;
    margin: 0 auto 2rem auto;
    animation: rotation 1s linear infinite;
}

@keyframes rotation {
    0% {
        transform: rotate(0deg);
    }

    100% {
        transform: rotate(360deg);
    }
}

.loading-text h3 {
    margin: 0 0 0.5rem 0;
    font-size: 1.25rem;
    color: #2c3e50;
}

.loading-text p {
    margin: 0;
    color: #6c757d;
    font-size: 0.875rem;
}

/* Responsive Design */
@media (max-width: 768px) {
    .validation-container {
        padding: 1rem;
    }

    .header-content {
        flex-direction: column;
        text-align: center;
        gap: 1rem;
    }

    .header-text h2 {
        font-size: 1.5rem;
    }

    .card-header {
        padding: 1rem;
        flex-direction: column;
        gap: 1rem;
        align-items: flex-start;
    }

    .card-content {
        padding: 1rem;
    }

    .action-content {
        flex-direction: column;
        gap: 2rem;
    }

    .upload-button {
        width: 100%;
    }

    .parameters-grid {
        grid-template-columns: 1fr;
    }

    .info-grid {
        grid-template-columns: 1fr;
    }
}

/* Upload Section Styles */
.upload-section {
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
    align-items: center;
}

/* GPU Configuration Styles */
.gpu-config-option {
    display: flex;
    align-items: center;
    justify-content: center;
    margin-bottom: 0.5rem;
}

.gpu-toggle-label {
    display: flex;
    align-items: center;
    cursor: pointer;
    user-select: none;
    gap: 1rem;
    padding: 0.75rem 1.5rem;
    background: #f8f9fa;
    border-radius: 12px;
    border: 1px solid #e1e5e9;
    transition: all 0.3s ease;
}

.gpu-toggle-label:hover {
    background: #e9ecef;
    border-color: #667eea;
}

.gpu-toggle-input {
    position: absolute;
    opacity: 0;
    cursor: pointer;
    height: 0;
    width: 0;
}

.gpu-toggle-slider {
    position: relative;
    width: 48px;
    height: 24px;
    background-color: #ccc;
    border-radius: 24px;
    transition: all 0.3s ease;
}

.gpu-toggle-slider::before {
    content: "";
    position: absolute;
    height: 20px;
    width: 20px;
    left: 2px;
    top: 2px;
    background-color: white;
    border-radius: 50%;
    transition: all 0.3s ease;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
}

.gpu-toggle-input:checked + .gpu-toggle-slider {
    background-color: #667eea;
}

.gpu-toggle-input:checked + .gpu-toggle-slider::before {
    transform: translateX(24px);
}

.gpu-toggle-text {
    display: flex;
    flex-direction: column;
    gap: 0.25rem;
}

.gpu-toggle-title {
    font-weight: 600;
    color: #2c3e50;
    font-size: 1rem;
}

.gpu-toggle-desc {
    color: #6c757d;
    font-size: 0.875rem;
}

@media (max-width: 480px) {
    .validation-container {
        padding: 0.5rem;
    }

    .header-text h2 {
        font-size: 1.25rem;
    }

    .rule-header {
        padding: 1rem;
    }

    .rule-content {
        padding: 1rem;
    }
    
    .upload-section {
        gap: 1rem;
    }
    
    .gpu-toggle-label {
        padding: 0.5rem 1rem;
        gap: 0.75rem;
    }
}
</style>