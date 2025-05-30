<template>
  <div class="modal-overlay">
    <div class="modal-container">
      <div class="modal-header">
        <h2>Plugin Settings</h2>
        <button @click="$emit('close')" class="close-button">✕</button>
      </div>
      <div class="modal-info" v-if="currentStep === 1">
        <PluginInformation ref="pluginInfoComponent" :newPlugin="plugin" @update-plugin="updatePlugin" />
      </div>
      <div class="modal-flow" v-if="currentStep === 2">
        <PluginFlowchart ref="pluginFlowchartComponent" :newRules="rules" :newDrawflow="drawflow"
          @update-rules="updateRules" @update-drawflow="updateDrawflow" />
      </div>
      <div class="modal-val" v-if="currentStep === 3">
        <PluginValidation @close="close" :plugin="plugin" :rules="rules" :drawflow="drawflow" />
      </div>
      <div class="modal-actions">
        <button @click="prevStep" :disabled="currentStep === 1">Prev</button>
        <button @click="nextStep" :disabled="currentStep === 3">Next</button>
      </div>
    </div>
  </div>
</template>

<script>
import PluginInformation from "@/components/pluginComponents/PluginInformation.vue";
import PluginFlowchart from "@/components/pluginComponents/PluginFlowchart.vue";
import PluginValidation from "@/components/pluginComponents/PluginValidation.vue";

import { getPluginFile, getPluginReferenceFolders, getPluginPackageList } from "@/api/index";

export default {
  props: {
    editName: {
      type: String,
      required: true
    },
    editDescription: {
      type: String,
      required: true
    },
    editDependencies: {
      type: Object,
      required: true
    },
    editDrawflow: {
      type: Object,
      required: true
    },
    editRules: {
      type: Array,
      required: true
    }
  },
  data() {
    return {
      currentStep: 1,
      // plugin: {
      //   name: '',
      //   description: '',
      //   dependencyFiles: [],
      // },
      // rules: [],
      // drawflow: {},
      plugin: {
        name: this.editName,
        description: this.editDescription,
        dependencyFiles: this.editDependencies,
        referenceFolders: [],
        packageFiles: []
      },
      drawflow: this.editDrawflow,
      rules: this.editRules
    };
  },
  components: {
    PluginInformation,
    PluginFlowchart,
    PluginValidation,
  },
  async mounted() {
    const response = await getPluginReferenceFolders(this.editName);

    // Fetch file content for all files in the folder structure
    const reference_folders = await this.fetchReferenceFolderFiles(this.editName, response.data.reference_folders);
    console.log("Updated reference folders with file content:", reference_folders);

    this.plugin.referenceFolders = reference_folders;

    const dependencies = Object.keys(this.plugin.dependencyFiles);
    let dependencyFiles = [];
    if (dependencies.length > 0) {
      console.log(dependencies);
      for (let i = 0; i < dependencies.length; i++) {
        const dependency = dependencies[i];
        const file_info = {
          plugin_name: this.editName,
          file_name: dependency
        }
        const response = await getPluginFile(file_info);
        const fileBlob = response.data;
        const file = new File([fileBlob], dependency, { type: fileBlob.type });

        dependencyFiles.push({
          file: file,
          fileName: dependency,
          type: dependency
        });
      }

      this.plugin.dependencyFiles = dependencyFiles;
    }
    if (this.rules.length > 0) {
      console.log(this.rules);
      // this.rules를 순회하면서 script가 string type일 경우, getPluginFile을 활용해서 파일을 가져온다.
      for (let i = 0; i < this.rules.length; i++) {
        if (typeof this.rules[i].script === 'string') {
          const file_info = {
            plugin_name: this.editName,
            file_name: this.rules[i].script
          }
          const response = await getPluginFile(file_info);
          const fileBlob = response.data;
          const file = new File([fileBlob], this.rules[i].script, { type: fileBlob.type });

          this.rules[i].script = file;
        }
      }
    }

    const packageFiles = await getPluginPackageList(this.editName);

    if (packageFiles.data.package_files.length > 0) {
      console.log(packageFiles.data.package_files);
      let packageFilesList = [];
      for (let i = 0; i < packageFiles.data.package_files.length; i++) {
        const packageFile = packageFiles.data.package_files[i];
        const file_info = {
          plugin_name: this.editName,
          file_name: packageFile
        }
        const response = await getPluginFile(file_info);
        const fileBlob = response.data;
        const file = new File([fileBlob], packageFile, { type: fileBlob.type });

        packageFilesList.push({
          file: file,
          fileName: packageFile,
        });
      }

      this.plugin.packageFiles = packageFilesList;
    }
  },
  methods: {
    async fetchReferenceFolderFiles(pluginName, referenceFolders) {
      async function processFolders(pluginName, folders) {
        for (const folder of folders) {
          // Process files in the current folder
          for (const file of folder.files) {
            try {
              const file_info = {
                plugin_name: pluginName,
                file_name: file.name,
              };
              const response = await getPluginFile(file_info); // API call to get file content
              const fileBlob = response.data;

              // Create a File object and add it to the file object
              file.file = new File([fileBlob], file.name, { type: fileBlob.type || "text/plain" });
            } catch (error) {
              console.error(`Failed to fetch file: ${file.name}`, error);
            }
          }

          // Recursively process subfolders
          if (folder.subFolders && folder.subFolders.length > 0) {
            await processFolders(pluginName, folder.subFolders);
          }
        }
      }

      // Start processing the top-level reference folders
      await processFolders(pluginName, referenceFolders);
      return referenceFolders; // Return updated reference folders
    },

    prevStep() {
      if (this.currentStep > 1) {
        this.emitCurrentStepData();
        this.currentStep--;
      }
    },
    nextStep() {
      if (this.currentStep < 3) {
        this.emitCurrentStepData();
        this.currentStep++;
      }
    },
    emitCurrentStepData() {
      if (this.currentStep === 1) {
        this.$refs.pluginInfoComponent.emitPluginData();
      } else if (this.currentStep === 2) {
        this.$refs.pluginFlowchartComponent.emitFlowchartData();
      }
    },
    updatePlugin(pluginData) {
      this.plugin = pluginData;
    },
    updateRules(rulesData) {
      this.rules = rulesData;
    },
    updateDrawflow(drawflowData) {
      this.drawflow = drawflowData;
    },
    close() {
      this.$emit('close');
    }
  },
};
</script>

<style scoped>
.modal-overlay {
  position: fixed;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  background-color: rgba(0, 0, 0, 0.5);
  display: flex;
  justify-content: center;
  align-items: center;
  z-index: 1000;
}

.modal-container {
  background-color: white;
  padding: 1rem;
  border-radius: 8px;
  box-shadow: 0px 0px 10px rgba(0, 0, 0, 0.1);
  width: fit-content;
  height: 50rem;
  overflow-y: auto;
  display: flex;
  flex-direction: column;
  box-sizing: border-box;
}

.modal-header {
  display: flex;
  justify-content: center;
  align-items: center;
  border-bottom: 1px solid #eee;
  padding-bottom: 10px;
  margin-bottom: 20px;
  position: relative;
}

.modal-header h2 {
  margin: 0;
  flex-grow: 1;
  text-align: center;
}

.close-button {
  /* red round close button */
  background-color: #ff5c5c;
  color: white;
  border: none;
  border-radius: 50%;
  width: 1.2rem;
  height: 1.2rem;
  cursor: pointer;
  display: flex;
  justify-content: center;
  align-items: center;
  font-size: 0.8rem;
  font-weight: bold;
  position: absolute;
  right: 0;
}

.modal-info {
  width: 25rem;
}

.modal-flow {
  width: 70rem;
}

.modal-val {
  width: 50rem;
}

.modal-actions {
  display: flex;
  justify-content: flex-end;
  border-top: 1px solid #eee;
  padding-top: 10px;
  margin-top: 1rem;
}

.modal-actions button {
  background-color: #007BFF;
  color: white;
  border: none;
  padding: 0.6rem 1rem;
  border-radius: 5px;
  cursor: pointer;
  margin-left: 10px;
}

.modal-actions button:disabled {
  background-color: #ccc;
}

.modal-actions button:hover:not(:disabled) {
  background-color: #0056b3;
}
</style>
