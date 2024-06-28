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
            <p v-else-if="uploadingStep == 3">Uploading Plugin Scripts...</p>
        </div>
    </div>
</template>

<script>
import { validationPlugin, uploadPluginMetadata, uploadPluginScripts } from '@/api/index';

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
        };
    },
    mounted() {
        this.processDependencyFiles();
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
                    const rules = this.rules.map(rule => {
                        return {
                            name: rule.name,
                            input: rule.input,
                            output: rule.output,
                            script: rule.script.name,
                            parameters: rule.parameters,
                            nodeId: rule.nodeId,
                        };
                    });
                    this.uploadingStep += 1;
                    const validation_response = await validationPlugin(this.plugin, rules, this.drawflow);
                    console.log('Validation Response:', validation_response);

                    this.uploadingStep += 1;
                    const pluginCreate = validation_response.data.plugin;
                    const db_plugin = await uploadPluginMetadata(pluginCreate);
                    console.log("Plugin Data :", db_plugin.data);

                    this.uploadingStep += 1;
                    const scriptNameList = validation_response.data.scripts;
                    const scriptCreate = this.rules.map(rule => {
                        const scriptName = scriptNameList.find(s => s === rule.script.name);
                        if (scriptName) {
                            return {
                                name: scriptName,
                                content: rule.script  // File 객체 자체를 content로 저장
                            };
                        }
                    }).filter(s => s);

                    const formData = new FormData();

                    // 비동기 작업을 처리하기 위한 Promise 배열
                    const promises = scriptCreate.map((script, index) => {
                        return new Promise((resolve, reject) => {
                            const reader = new FileReader();
                            reader.onload = (event) => {
                                const fileContent = event.target.result;
                                console.log(`Script Name: ${script.name}`);
                                console.log(`Script Content: ${fileContent}`);

                                const file = new Blob([fileContent], { type: script.content.type });
                                formData.append('files', file, script.name);
                                resolve();
                            };
                            reader.onerror = (error) => {
                                reject(error);
                            };
                            reader.readAsText(script.content);  // File 객체를 텍스트로 읽기
                        });
                    });

                    // 모든 파일을 읽고 나서 FormData에 추가
                    Promise.all(promises).then(() => {
                        // plugin_name을 FormData에 추가
                        formData.append('plugin_name', db_plugin.data.plugin.name);

                        // FormData 내용 확인 (선택 사항)
                        for (let [key, value] of formData.entries()) {
                            console.log(`${key}: ${value}`);
                        }

                        // 스크립트 업로드
                        uploadPluginScripts(formData).then((scripts) => {
                            console.log(scripts);
                        }).catch((error) => {
                            console.error('Error while uploading plugin:', error);
                        });
                    }).catch((error) => {
                        console.error('Error reading files:', error);
                    });

                    this.uploadingStep = 0;
                    this.validationLoading = false;
                    alert('Plugin uploaded successfully!');

                    this.$emit('close');
                } catch (error) {
                    console.error('Error while uploading plugin:', error);
                    alert('Error while uploading plugin. Please try again.');
                    this.validationLoading = false;
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

            const scriptName = rule.script.name;
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