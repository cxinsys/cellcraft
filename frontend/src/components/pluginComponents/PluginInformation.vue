<template>
  <div class="plugin-container">
    <!-- 플러그인 이름 입력 필드 -->
    <div class="input-group">
      <label class="input-group__label" for="pluginName">Plugin Name</label>
      <input type="text" id="pluginName" v-model="plugin.name" />
    </div>

    <!-- 플러그인 설명 입력 필드 -->
    <div class="input-group">
      <label class="input-group__label" for="pluginDescription">Plugin Description</label>
      <textarea id="pluginDescription" v-model="plugin.description" rows="4"></textarea>
    </div>

    <!-- 참조 스크립트 폴더 업로드 -->
    <div class="input-group">
      <label class="input-group__label">Upload Reference Script Folder</label>
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
          <button class="tree-nav__item-remove" type="button" @click="removeReferenceFolder(folder.folderName)">
            Remove
          </button>
        </details>
      </nav>
    </div>

    <!-- 로컬 의존성 패키지 파일 업로드 -->
    <div class="input-group">
      <h3 class="input-group__label">Upload Local Dependency Package Files</h3>
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
    <div class="input-group">
      <label class="input-group__label" for="dependencyType">Select Dependency File Type</label>
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
        <button type="button" @click="removeDependencyFile(file.type)">Remove</button>
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
    }
  },
  data() {
    return {
      plugin: {
        name: '',
        description: '',
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
</style>
