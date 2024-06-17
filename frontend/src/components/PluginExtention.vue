<template>
  <div class="modal-overlay">
    <div class="modal-container">
      <div class="modal-header">
        <h2>Plugin Extension</h2>
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

export default {
  data() {
    return {
      currentStep: 1,
      plugin: {
        name: '',
        description: '',
        dependencyFiles: [],
      },
      rules: [],
      drawflow: {},
    };
  },
  components: {
    PluginInformation,
    PluginFlowchart,
    PluginValidation,
  },
  methods: {
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
