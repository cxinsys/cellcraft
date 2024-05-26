<template>
    <div class="validation-container">
        <h2>Validate and Upload Plugin</h2>
        <div class="section plugin-info">
            <h3>Plugin Information</h3>
            <p><strong>Name:</strong> {{ plugin.name }}</p>
            <p><strong>Description:</strong> {{ plugin.description }}</p>
            <p>
                <strong>Dependency Files:</strong>
            <ul>
                <li v-for="(dependency, index) in plugin.dependencyFiles" :key="index">
                    {{ dependency.type }} - {{ dependency.file }}
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
                    <div>{{ rule.script }}</div>
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
    </div>
</template>

<script>
import { uploadPlugin } from '@/api/index';

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
            reponse: null
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
                // 플러그인 업로드 로직을 구현
                console.log('Uploading plugin...', this.plugin, this.rules);
                const response = await uploadPlugin(this.plugin, this.rules, this.drawflow);
                console.log('Response:', response);
                this.response = response;
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

            const scriptName = rule.script ? rule.script : '';
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
</style>