<template>
  <div class="plugin-container">
    <!-- 플러그인 이름 입력 필드 -->
      <div class="input-group">
        <label for="pluginName">Plugin Name:</label>
        <input type="text" id="pluginName" v-model="plugin.name" />
      </div>

      <!-- 플러그인 설명 입력 필드 -->
      <div class="input-group">
        <label for="pluginDescription">Plugin Description:</label>
        <textarea id="pluginDescription" v-model="plugin.description" rows="4"></textarea>
      </div>

      <!-- 의존성 파일 타입 선택 드롭다운 -->
      <div class="input-group">
        <label for="dependencyType">Select Dependency File Type:</label>
        <select id="dependencyType" v-model="selectedDependencyType" @change="addDependencyType">
          <option value="" disabled>Select a file type</option>
          <option value="requirements.txt" :disabled="isFileTypeAdded('requirements.txt')">requirements.txt</option>
          <option value="environment.yml" :disabled="isFileTypeAdded('environment.yml')">environment.yml</option>
          <option value="renv.lock" :disabled="isFileTypeAdded('renv.lock')">renv.lock</option>
        </select>
      </div>

      <!-- 의존성 파일 업로드 필드 -->
      <div v-for="(file, index) in plugin.dependencyFiles" :key="index" class="input-group">
        <label :for="file.type">Upload {{ file.type }}:</label>
        <div class="file-upload">
          <input type="file" :id="file.type" @change="handleFileUpload($event, file.type)" />
          <button type="button" @click="removeDependencyFile(file.type)">Remove</button>
        </div>
      </div>
  </div>
</template>

<script>
export default {
  data() {
    return {
      plugin: {
        name: '',
        description: '',
        dependencyFiles: [],
      },
      selectedDependencyType: '', // 선택한 의존성 파일 타입
    };
  },
  methods: {
    addDependencyType() {
      if (this.selectedDependencyType && !this.isFileTypeAdded(this.selectedDependencyType)) {
        this.plugin.dependencyFiles.push({ type: this.selectedDependencyType, file: null });
        this.selectedDependencyType = '';
      }
    },
    isFileTypeAdded(type) {
      return this.plugin.dependencyFiles.some(file => file.type === type);
    },
    handleFileUpload(event, type) {
      const file = event.target.files[0];
      const index = this.plugin.dependencyFiles.findIndex(file => file.type === type);
      if (index !== -1) {
        this.plugin.dependencyFiles[index].file = file;
      }
    },
    removeDependencyFile(type) {
      this.plugin.dependencyFiles = this.plugin.dependencyFiles.filter(file => file.type !== type);
    },
    emitPluginData() {
      this.$emit('update-plugin', this.plugin);
    },
  },
  watch: {
    plugin: {
      handler() {
        this.emitPluginData();
      },
      deep: true
    }
  },
};
</script>


<style scoped>
.plugin-container {
  width: 100%;
  background-color: white;
}

.input-group {
  width: 100%;
  margin-bottom: 1rem;
}

.input-group label {
  display: block;
  margin-bottom: 5px;
  font-weight: bold;
  color: #333;
}

.input-group input,
.input-group textarea,
.input-group select {
  width: 100%;
  padding: 10px;
  border: 1px solid #ccc;
  border-radius: 4px;
  font-size: 1rem;
  box-sizing: border-box;
}

.file-upload{
  display: flex;
  align-items: center;
}

input[type="file"] {
  border: none;
}

button {
  background-color: #ff0000;
  color: white;
  border: none;
  padding: 0.5rem;
  border-radius: 5px;
  cursor: pointer;
  font-size: 0.8rem;
  transition: background-color 0.3s ease;
}

button:hover {
  background-color: #b10101;
}
</style>
