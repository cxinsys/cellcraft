<template>
  <div class="plugin-container">
    <!-- 플러그인 이름 입력 필드 -->
    <div class="input-group">
      <label class="input-group__label" for="pluginName">Plugin Name</label>
      <input type="text" id="pluginName" v-model="plugin.name" :readonly="readOnly" />
    </div>

    <!-- 플러그인 설명 입력 필드 -->
    <div class="input-group">
      <label class="input-group__label" for="pluginDescription">Plugin Description</label>
      <textarea id="pluginDescription" v-model="plugin.description" rows="4" :readonly="readOnly"></textarea>
    </div>

    <!-- 플러그인 타입 선택 -->
    <div class="input-group" v-if="!readOnly">
      <label class="input-group__label">Plugin Type</label>
      <div class="plugin-type-selection">
        <label class="radio-option">
          <input type="radio" v-model="plugin.pluginType" value="analysis" />
          <span class="radio-label">
            <i class="fas fa-calculator"></i>
            Analysis Plugin
          </span>
          <p class="type-description">For data processing and computational analysis</p>
          <p class="type-constraints">• Multiple outputs allowed</p>
        </label>
        <label class="radio-option">
          <input type="radio" v-model="plugin.pluginType" value="visualization" />
          <span class="radio-label">
            <i class="fas fa-chart-bar"></i>
            Visualization Plugin
          </span>
          <p class="type-description">For data visualization using Plotly</p>
          <p class="type-constraints">• Must have exactly 1 JSON output per node</p>
        </label>
      </div>
    </div>

    <!-- 참조 스크립트 폴더 업로드 -->
    <div class="input-group" v-if="!readOnly">
      <label class="input-group__label">Custom Script Folder</label>
      <div class="file-upload">
        <input type="file" id="scriptFolder" webkitdirectory directory @change="handleScriptFolderUpload"
          class="file-input" />
        <label class="file-label" for="scriptFolder">
          Click to upload a folder
        </label>
      </div>
    </div>

    <div v-if="plugin.referenceFolders.length" class="input-group folder-tree">
      <nav class="tree-nav">
        <details v-for="(folder, idx) in plugin.referenceFolders" :key="folder.folderName"
          class="tree-nav__item is-expandable" :open="toggleFolder === idx" @click="toggleFolderState(idx)">
          <summary class="tree-nav__item-title">
            <img class="folder__item--icon" src="@/assets/open-folder.png" v-if="toggleFolder === idx"
              alt="Open Folder" />
            <img class="folder__item--icon" src="@/assets/folder.png" v-else alt="Closed Folder" />
            {{ folder.folderName }}
          </summary>
          <div class="tree-nav__item">
            <!-- 하위 폴더 재귀 렌더링 -->
            <details v-for="(subFolder, subIdx) in folder.subFolders" :key="subIdx" class="tree-nav__item">
              <summary class="tree-nav__item-title" :style="{ paddingLeft: (subIdx + 1) * 20 + 'px' }">
                <img class="folder__item--icon" src="@/assets/open-folder.png"
                  v-if="toggleFolder === idx + '-' + subIdx" alt="Open SubFolder" />
                <img class="folder__item--icon" src="@/assets/folder.png" v-else alt="Closed SubFolder" />
                {{ subFolder.folderName }}
              </summary>
              <div>
                <a v-for="file in subFolder.files" :key="file.name" class="tree-nav__item-title"
                  :style="{ paddingLeft: (subIdx + 2) * 20 + 'px' }">
                  {{ file.name }}
                </a>
              </div>
            </details>
            <a v-for="file in folder.files" :key="file.name" class="tree-nav__item-title"
              :style="{ paddingLeft: '20px' }">
              {{ file.name }}
            </a>
          </div>
          <button v-if="!readOnly" class="tree-nav__item-remove" type="button"
            @click="removeReferenceFolder(folder.folderName)">
            Remove
          </button>
        </details>
      </nav>
    </div>

    <!-- 로컬 의존성 패키지 파일 업로드 -->
    <div class="input-group" v-if="!readOnly">
      <h3 class="input-group__label">Local Dependency Files</h3>
      <div v-for="(file, index) in plugin.packageFiles" :key="index" class="input-group">
        <div class="file-upload">
          <label :for="'dependency-' + index" class="file-label">
            <input type="file" :id="'dependency-' + index" @change="handlePackageFileUpload($event, index)"
              class="file-input" accept=".whl,.gz" />
            {{ file.fileName || "Click to upload a .whl or .tar.gz file" }}
          </label>
          <button type="button" @click="removePackageFile(index)">Remove</button>
        </div>
      </div>
      <button class="add-button" type="button" @click="addPackageFile">Add Another File</button>
    </div>

    <!-- 의존성 파일 타입 선택 드롭다운 -->
    <div class="input-group" v-if="!readOnly">
      <label class="input-group__label" for="dependencyType">Standard Dependency Files</label>
      <select id="dependencyType" v-model="selectedDependencyType" @change="addDependencyType">
        <option value="" disabled>Select a file type</option>
        <option value="requirements.txt" :disabled="isFileTypeAdded('requirements.txt')">requirements.txt</option>
        <option value="environment.yml" :disabled="isFileTypeAdded('environment.yml')">environment.yml</option>
        <option value="renv.lock" :disabled="isFileTypeAdded('renv.lock')">renv.lock</option>
      </select>
    </div>

    <div v-for="(file, index) in plugin.dependencyFiles" :key="index" class="input-group">
      <label class="input-group__label" :for="file.type">{{ file.type }}</label>
      <div class="file-upload">
        <label :for="file.type" class="file-label">
          <input type="file" :id="file.type" @change="handleFileUpload($event, file.type)" class="file-input" />
          {{ file.fileName || 'Click to upload a file' }}
        </label>
        <button v-if="!readOnly" type="button" @click="removeDependencyFile(file.type)">Remove</button>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  props: {
    newPlugin: {
      type: Object,
      required: true
    },
    readOnly: {
      type: Boolean,
      default: false
    }
  },
  data() {
    return {
      plugin: {
        name: '',
        description: '',
        pluginType: 'analysis', // 기본값을 analysis로 설정
        referenceFolders: [],
        dependencyFiles: [],
        packageFiles: [],
      },
      referenceFolderName: '', // 업로드된 참조 폴더 이름 저장
      selectedDependencyType: '', // 선택한 의존성 파일 타입
      toggleFolder: null,
    };
  },
  watch: {
    'newPlugin.name': {
      handler(newValue) {
        this.plugin.name = newValue;
      },
      immediate: true
    },
    'newPlugin.description': {
      handler(newValue) {
        this.plugin.description = newValue;
      },
      immediate: true
    },
    'newPlugin.pluginType': {
      handler(newValue) {
        this.plugin.pluginType = newValue || 'analysis';
      },
      immediate: true
    },
    newPlugin: {
      handler(newValue) {
        // plugin.dependencyFiles을 직접 할당하기 전에 변화가 있을 경우에만 업데이트
        if (JSON.stringify(this.plugin.dependencyFiles) !== JSON.stringify(newValue.dependencyFiles)) {
          this.plugin.dependencyFiles = [...newValue.dependencyFiles];
        }
        if (JSON.stringify(this.plugin.referenceFolders) !== JSON.stringify(newValue.referenceFolders)) {
          this.plugin.referenceFolders = [...newValue.referenceFolders];
        }
        if (JSON.stringify(this.plugin.packageFiles) !== JSON.stringify(newValue.packageFiles)) {
          this.plugin.packageFiles = [...newValue.packageFiles];
        }
      },
      deep: true,
      immediate: true
    }
  },
  methods: {
    handlePackageFileUpload(event, index) {
      const file = event.target.files[0];
      if (file) {
        this.$set(this.plugin.packageFiles, index, {
          file,
          fileName: file.name,
        });
      }
    },
    removePackageFile(index) {
      this.plugin.packageFiles.splice(index, 1); // Remove the file from the list
      this.emitDependencyFiles();
    },
    addPackageFile() {
      this.plugin.packageFiles.push({ file: null, fileName: "" }); // Add a new empty file slot
    },
    toggleFolderState(idx) {
      // 열려 있으면 닫고, 닫혀 있으면 엽니다.
      this.toggleFolder = this.toggleFolder === idx ? null : idx;
    },
    handleScriptFolderUpload(event) {
      const files = Array.from(event.target.files); // FileList를 배열로 변환

      // 중첩된 디렉터리 구조로 변환
      const buildFolderStructure = (files) => {
        const root = {};
        files.forEach((file) => {
          const parts = file.webkitRelativePath.split("/");
          let current = root;
          parts.forEach((part, idx) => {
            if (idx === parts.length - 1) {
              if (!current.files) current.files = [];
              current.files.push({
                name: file.name,
                file: file,
                type: file.type,
              });
            } else {
              if (!current[part]) current[part] = { subFolders: [], files: [] };
              current = current[part];
            }
          });
        });
        return root;
      };

      const folderStructure = buildFolderStructure(files);

      // ReferenceFolder 스키마로 변환
      const convertToReferenceFolders = (folderName, folderData) => {
        return {
          folderName,
          files: folderData.files || [],
          subFolders: Object.entries(folderData)
            .filter(([key]) => key !== "files" && key !== "subFolders")
            .map(([subFolderName, subFolderData]) =>
              convertToReferenceFolders(subFolderName, subFolderData)
            ),
        };
      };

      const referenceFolders = Object.entries(folderStructure).map(([folderName, folderData]) =>
        convertToReferenceFolders(folderName, folderData)
      );

      this.plugin.referenceFolders = referenceFolders;

      this.emitPluginData();
    },
    removeReferenceFolder(folderName) {
      this.plugin.referenceFolders = this.plugin.referenceFolders.filter(folder => folder.folderName !== folderName);
      this.referenceFolderName = '';
      this.emitPluginData();
    },
    addDependencyType() {
      if (!this.plugin.dependencyFiles) {
        this.plugin.dependencyFiles = [{ type: this.selectedDependencyType, file: null }];
        this.selectedDependencyType = '';
        return true;
      }
      if (this.selectedDependencyType && !this.isFileTypeAdded(this.selectedDependencyType)) {
        this.plugin.dependencyFiles.push({ type: this.selectedDependencyType, file: null });
        this.selectedDependencyType = '';
      }
      this.emitPluginData();
    },
    isFileTypeAdded(type) {
      if (!this.plugin.dependencyFiles) return false;
      else return this.plugin.dependencyFiles.some(file => file.type === type);
    },
    handleFileUpload(event, type) {
      const file = event.target.files[0];
      const index = this.plugin.dependencyFiles.findIndex(file => file.type === type);
      if (index !== -1) {
        this.$set(this.plugin.dependencyFiles, index, { ...this.plugin.dependencyFiles[index], file, fileName: file.name });
        this.emitPluginData();
      }
    },
    removeDependencyFile(type) {
      this.plugin.dependencyFiles = this.plugin.dependencyFiles.filter(file => file.type !== type);
      this.emitPluginData();
    },
    emitPluginData() {
      this.$emit('update-plugin', this.plugin);
    }
  }
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

.input-group__label {
  display: block;
  margin-bottom: 1rem;
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

.file-label {
  display: inline-block;
  padding: 0.5rem 1rem;
  cursor: pointer;
  background-color: #f0f0f0;
  border: 1px solid #ccc;
  border-radius: 4px;
  margin-bottom: 0;
  margin-right: 1rem;
}

.file-upload {
  display: flex;
  align-items: center;
}

.file-input {
  display: none;
}

.folder-tree h3 {
  color: #2d2d2d;
  font-size: 1.1rem;
  margin-bottom: 1rem;
}

.tree-nav__item {
  display: block;
  white-space: nowrap;
  color: #ccc;
  position: relative;
}

.tree-nav details {
  position: relative;
}

.tree-nav__item-remove {
  position: absolute;
  right: 1rem;
  top: 0.25rem;
  cursor: pointer;
}

/* .tree-nav__item.is-expandable::before {
  border-left: 1px solid #333;
  content: "";
  height: 100%;
  left: 0.8rem;
  position: absolute;
  top: 2.4rem;
  height: calc(100% - 2.4rem);
}
 */

/* .tree-nav__item.is-expandable[open]>.tree-nav__item-title::before {
  font-family: "ionicons";
  transform: rotate(90deg);
} */

.tree-nav__item-title {
  cursor: pointer;
  display: block;
  outline: 0;
  color: #2d2d2d;
  font-size: 0.9rem;
  line-height: 2.5rem;
  padding: 0 1rem;
}

.tree-nav__item-title::-webkit-details-marker {
  display: none;
}

.tree-nav__item-title:hover {
  background-color: #f0f0f0;
}

.folder__item--icon {
  width: 16px;
  height: 16px;
  vertical-align: middle;
  margin-right: 8px;
}

.folder__item--icon.large {
  width: 24px;
  height: 24px;
  margin-right: 8px;
}

.add-button {
  background-color: #06bb00;
  color: white;
  border: none;
  padding: 0.5rem;
  border-radius: 5px;
  cursor: pointer;
  font-size: 0.8rem;
  transition: background-color 0.3s ease;
}

.add-button:hover {
  background-color: #009a00;
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

/* 읽기 전용 스타일 */
input[readonly],
textarea[readonly] {
  background-color: #f5f5f5;
  cursor: not-allowed;
  opacity: 0.8;
}

/* 플러그인 타입 선택 스타일 */
.plugin-type-selection {
  display: flex;
  gap: 1rem;
  margin-top: 0.5rem;
}

.radio-option {
  flex: 1;
  border: 2px solid #e0e0e0;
  border-radius: 12px;
  padding: 1.25rem;
  cursor: pointer;
  transition: all 0.3s ease;
  display: block;
  position: relative;
  background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
}

.radio-option:hover {
  border-color: #007BFF;
  background: linear-gradient(135deg, #f8f9fa 0%, #e3f2fd 100%);
  box-shadow: 0 4px 8px rgba(0, 123, 255, 0.1);
  transform: translateY(-1px);
}

.radio-option input[type="radio"] {
  display: none;
}

.radio-option input[type="radio"]:checked+.radio-label {
  color: #007BFF;
}

.radio-option input[type="radio"]:checked~.type-description,
.radio-option input[type="radio"]:checked~.type-constraints {
  color: #495057;
}

.radio-option:has(input[type="radio"]:checked) {
  border-color: #007BFF;
  background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%);
  box-shadow: 0 4px 12px rgba(0, 123, 255, 0.2);
  transform: translateY(-2px);
}

/* Analysis 플러그인 타입에 특별한 색상 */
.radio-option:has(input[value="analysis"]:checked) {
  border-color: #28a745;
  background: linear-gradient(135deg, #e8f5e8 0%, #d4edda 100%);
  box-shadow: 0 4px 12px rgba(40, 167, 69, 0.2);
}

.radio-option:has(input[value="analysis"]:checked) .radio-label {
  color: #28a745;
}

.radio-option:has(input[value="analysis"]) .radio-label i {
  color: #28a745;
}

/* Visualization 플러그인 타입에 특별한 색상 */
.radio-option:has(input[value="visualization"]:checked) {
  border-color: #6f42c1;
  background: linear-gradient(135deg, #f3e5f5 0%, #e8d4f8 100%);
  box-shadow: 0 4px 12px rgba(111, 66, 193, 0.2);
}

.radio-option:has(input[value="visualization"]:checked) .radio-label {
  color: #6f42c1;
}

.radio-option:has(input[value="visualization"]) .radio-label i {
  color: #6f42c1;
}

.radio-label {
  display: flex;
  align-items: center;
  gap: 0.75rem;
  font-weight: 600;
  font-size: 1.1rem;
  margin-bottom: 0.75rem;
  color: #333;
}

.radio-label i {
  font-size: 1.2rem;
}

.type-description {
  font-size: 0.95rem;
  color: #555;
  margin: 0.5rem 0;
  line-height: 1.5;
  font-weight: 400;
}

.type-constraints {
  font-size: 0.85rem;
  color: #777;
  margin: 0.5rem 0 0 0;
  font-style: normal;
  background-color: rgba(0, 0, 0, 0.03);
  padding: 0.5rem 0.75rem;
  border-radius: 6px;
  border-left: 3px solid #dee2e6;
}

/* Analysis 플러그인의 제약사항 스타일 */
.radio-option:has(input[value="analysis"]) .type-constraints {
  border-left-color: #28a745;
  background-color: rgba(40, 167, 69, 0.05);
}

/* Visualization 플러그인의 제약사항 스타일 */
.radio-option:has(input[value="visualization"]) .type-constraints {
  border-left-color: #6f42c1;
  background-color: rgba(111, 66, 193, 0.05);
}
</style>
