<template>
    <div class="validation-container" v-on:click.stop="validationLoading">
        <h2>Validate and Upload Plugin</h2>
        <div class="section plugin-info">
            <h3>Plugin Information</h3>
            <p><strong>Name:</strong> {{ plugin.name }}</p>
            <p><strong>Description:</strong> {{ plugin.description }}</p>
            <p>
                <strong>Dependency Files:</strong>
            <ul>
                <li v-for="(dependency, index) in plugin.dependencyFiles" :key="index">
                    {{ dependency.type }}
                </li>
            </ul>
            </p>
            <p><strong>Reference Folders:</strong>
            <ul>
                <li v-for="(folder, index) in plugin.referenceFolders" :key="index">
                    {{ folder.folderName }}
                    <ul>
                        <li v-for="(file, i) in folder.files" :key="i">{{ file.name }}</li>
                    </ul>
                </li>
            </ul>
            </p>
        </div>
        <div class="section rules-info">
            <h3>Rules</h3>
            <div v-for="(rule, index) in rules" :key="index" class="rule-item">
                <h4>Rule {{ rule.name }}</h4>
                <div>
                    <strong>Input:</strong>
                    <div v-for="(input, iIndex) in rule.input" :key="iIndex">{{ input }}</div>
                </div>
                <div>
                    <strong>Output:</strong>
                    <div v-for="(output, oIndex) in rule.output" :key="oIndex">{{ output }}</div>
                </div>
                <div v-if="rule.script">
                    <strong>Script:</strong>
                    <div>{{ rule.script.name }}</div>
                </div>
                <div>
                    <strong>Parameters:</strong>
                    <div v-for="(param, pIndex) in rule.parameters" :key="pIndex">
                        {{ param.name }} ({{ param.type }}): {{ param.defaultValue }}
                    </div>
                </div>
                <div v-if="rule.script">
                    <strong>script execution command:</strong>
                    <div>{{ generateShellCommand(rule) }}</div>
                </div>
            </div>
        </div>
        <div class="validation-actions">
            <button @click="uploadPluginData">Upload Plugin</button>
        </div>
        <div v-if="validationLoading" class="loading-animation">
            <div class="loader"></div>
            <p v-if="uploadingStep == 1">Validating Plugin...</p>
            <p v-else-if="uploadingStep == 2">Uploading Plugin Data...</p>
            <p v-else-if="uploadingStep == 3">Building Plugin Environment...</p>
            <p v-else-if="uploadingStep == 4">Uploading Plugin Scripts...</p>
        </div>
    </div>
</template>

<script>
import { validationPlugin, uploadPluginMetadata, uploadPluginScripts, uploadPluginPackage, buildPluginDocker } from '@/api/index';

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
                        referenceFolders: this.plugin.referenceFolders
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

                    // 7. Docker 이미지 빌드 (마지막 단계)
                    this.uploadingStep = 3;
                    console.log('Starting Docker build...');

                    const buildResponse = await buildPluginDocker(pluginCreate.name);
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
    padding: 2rem;
    background-color: #f9f9f9;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

h2 {
    margin-bottom: 1.5rem;
    font-size: 1.5rem;
    color: #333;
}

.section {
    padding: 1rem;
    margin-bottom: 1.5rem;
    border: 1px solid #ddd;
    border-radius: 4px;
    background-color: #fff;
}

h3 {
    margin-bottom: 1rem;
    font-size: 1.25rem;
    color: #007bff;
}

.rule-item {
    padding: 0.5rem;
    border: 1px solid #ccc;
    border-radius: 4px;
    margin-bottom: 1rem;
}

h4 {
    margin-bottom: 0.5rem;
    font-size: 1.1rem;
    color: #555;
}

p,
strong,
div {
    font-size: 1rem;
    color: #333;
}

.validation-actions {
    display: flex;
    justify-content: flex-end;
    margin-top: 1rem;
}

.validation-actions button {
    margin-left: 0.5rem;
    padding: 0.5rem 1rem;
    border: none;
    border-radius: 4px;
    cursor: pointer;
}

.validation-actions button:first-child {
    background-color: #007bff;
    color: #fff;
}

.validation-actions button:last-child {
    background-color: #6c757d;
    color: #fff;
}

.loading-animation {
    position: fixed;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    background-color: rgba(0, 0, 0, 0.5);
    color: #fff;
    padding: 1rem;
    border-radius: 8px;
    z-index: 9999;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: column;
}

.loader {
    width: 48px;
    height: 48px;
    border: 5px solid #FFF;
    border-bottom-color: transparent;
    border-radius: 50%;
    margin-bottom: 0.5rem;
    display: inline-block;
    box-sizing: border-box;
    animation: rotation 1s linear infinite;
}

.loading-animation p {
    margin: 0;
    color: #fff;
}
</style>