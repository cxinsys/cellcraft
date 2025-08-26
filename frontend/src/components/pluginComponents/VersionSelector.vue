<template>
  <div class="version-selector">
    <!-- Version Selection Dropdown -->
    <div class="version-dropdown-container">
      <label class="version-label" :for="`version-selector-${pluginName}`">
        Plugin Version
      </label>
      <div class="dropdown-wrapper">
        <select 
          :id="`version-selector-${pluginName}`"
          v-model="selectedVersion"
          :disabled="disabled || loading"
          class="version-select"
          :class="{ 'loading': loading, 'disabled': disabled || loading }"
          @change="onVersionChange"
          :aria-label="`Select version for ${pluginName} plugin`"
        >
          <option v-if="loading" value="">Loading versions...</option>
          <option 
            v-else
            v-for="version in availableVersions" 
            :key="version" 
            :value="version"
          >
            {{ version }}{{ version === currentVersion ? ' (Current)' : '' }}
          </option>
        </select>
        
        <!-- Loading spinner -->
        <div v-if="loading" class="loading-spinner">
          <div class="spinner"></div>
        </div>
        
        <!-- Update indicator -->
        <div 
          v-if="!loading && hasNewerVersion" 
          class="update-indicator"
          title="New version available"
        >
          <span class="update-badge">Update Available</span>
        </div>
      </div>
    </div>

    <!-- Current Version Info -->
    <div v-if="!loading && currentVersion" class="version-info">
      <div class="current-version">
        <span class="info-label">Current:</span>
        <span class="version-tag current">{{ currentVersion }}</span>
      </div>
      <div v-if="selectedVersion && selectedVersion !== currentVersion" class="selected-version">
        <span class="info-label">Selected:</span>
        <span class="version-tag selected">{{ selectedVersion }}</span>
      </div>
    </div>

    <!-- Update Button -->
    <div 
      v-if="!loading && selectedVersion && selectedVersion !== currentVersion" 
      class="update-button-container"
    >
      <button 
        @click="requestUpdate"
        :disabled="disabled || updating"
        class="update-button"
        :class="{ 'updating': updating }"
        :aria-label="`Update ${pluginName} from ${currentVersion} to ${selectedVersion}`"
      >
        <div v-if="updating" class="updating-content">
          <div class="update-spinner"></div>
          <span>Updating...</span>
        </div>
        <span v-else>Update to {{ selectedVersion }}</span>
      </button>
    </div>

    <!-- Error State -->
    <div v-if="error" class="error-container">
      <div class="error-message">
        <img src="@/assets/error.png" alt="Error" class="error-icon" />
        <span>{{ error }}</span>
      </div>
      <button @click="retryFetch" class="retry-button">
        Retry
      </button>
    </div>
  </div>
</template>

<script>
import { getPluginVersions, updatePluginVersion } from "@/api/index";

export default {
  name: "VersionSelector",
  props: {
    pluginName: {
      type: String,
      required: true
    },
    currentVersion: {
      type: String,
      default: "latest"
    },
    disabled: {
      type: Boolean,
      default: false
    },
    value: {
      type: String,
      default: ""
    }
  },
  // Note: Vue 2.6 doesn't support emits option - this is handled through $emit calls
  data() {
    return {
      availableVersions: [],
      selectedVersion: this.value || this.currentVersion,
      loading: false,
      updating: false,
      error: null
    };
  },
  computed: {
    hasNewerVersion() {
      if (!this.availableVersions.length || !this.currentVersion) return false;
      
      // Simple check: if current version is not the first in the list
      // Note: This assumes versions are sorted with newest first from the API
      // For more robust version comparison, consider using semver library
      return this.currentVersion !== this.availableVersions[0] && 
             this.availableVersions[0] !== 'latest';
    }
  },
  watch: {
    pluginName: {
      handler: 'fetchVersions',
      immediate: true
    },
    currentVersion: {
      handler(newVersion) {
        if (newVersion && !this.value) {
          this.selectedVersion = newVersion;
        }
      },
      immediate: true
    },
    value: {
      handler(newValue) {
        if (newValue !== this.selectedVersion) {
          this.selectedVersion = newValue || this.currentVersion;
        }
      }
    }
  },
  methods: {
    async fetchVersions() {
      if (!this.pluginName) return;
      
      this.loading = true;
      this.error = null;
      
      try {
        const response = await getPluginVersions(this.pluginName);
        const data = response.data;
        
        this.availableVersions = data.available_versions || ['latest'];
        
        // Set current version if not already set
        if (data.current_version && !this.currentVersion) {
          this.currentVersion = data.current_version;
        }
        
        // Set selected version to current if no value provided
        if (!this.value && this.currentVersion) {
          this.selectedVersion = this.currentVersion;
          this.$emit('input', this.selectedVersion);
        }
        
      } catch (error) {
        console.error(`Failed to fetch versions for ${this.pluginName}:`, error);
        this.error = this.getErrorMessage(error);
        this.availableVersions = [this.currentVersion || 'latest'];
      } finally {
        this.loading = false;
      }
    },
    
    onVersionChange() {
      this.$emit('input', this.selectedVersion);
      this.$emit('version-changed', {
        pluginName: this.pluginName,
        selectedVersion: this.selectedVersion,
        currentVersion: this.currentVersion
      });
    },
    
    async requestUpdate() {
      if (!this.selectedVersion || this.selectedVersion === this.currentVersion) return;
      
      this.updating = true;
      
      try {
        await updatePluginVersion(this.pluginName, this.selectedVersion);
        
        this.$emit('update-requested', {
          pluginName: this.pluginName,
          oldVersion: this.currentVersion,
          newVersion: this.selectedVersion
        });
        
        // Update current version after successful update
        this.currentVersion = this.selectedVersion;
        
      } catch (error) {
        console.error(`Failed to update ${this.pluginName} to version ${this.selectedVersion}:`, error);
        this.error = this.getErrorMessage(error);
      } finally {
        this.updating = false;
      }
    },
    
    retryFetch() {
      this.fetchVersions();
    },
    
    getErrorMessage(error) {
      if (error.response?.data?.detail) {
        return error.response.data.detail;
      } else if (error.message) {
        return `Network error: ${error.message}`;
      } else {
        return 'Failed to load plugin versions';
      }
    }
  }
};
</script>

<style scoped>
.version-selector {
  display: flex;
  flex-direction: column;
  gap: 1rem;
  padding: 1rem;
  background: #f8f9fa;
  border-radius: 8px;
  border: 1px solid #e1e5e9;
  transition: all 0.3s ease;
}

.version-selector:hover {
  border-color: #c3c9d0;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
}

/* Version Dropdown */
.version-dropdown-container {
  display: flex;
  flex-direction: column;
  gap: 0.5rem;
}

.version-label {
  font-weight: 600;
  font-size: 0.9rem;
  color: #374151;
  margin: 0;
}

.dropdown-wrapper {
  position: relative;
  display: flex;
  align-items: center;
  gap: 0.75rem;
}

.version-select {
  flex: 1;
  padding: 0.75rem 1rem;
  border: 1px solid #d1d5db;
  border-radius: 6px;
  background: white;
  font-size: 0.9rem;
  color: #374151;
  cursor: pointer;
  transition: all 0.2s ease;
  min-width: 150px;
}

.version-select:hover:not(:disabled) {
  border-color: #2196f3;
  box-shadow: 0 0 0 2px rgba(33, 150, 243, 0.1);
}

.version-select:focus {
  outline: none;
  border-color: #2196f3;
  box-shadow: 0 0 0 2px rgba(33, 150, 243, 0.2);
}

.version-select.loading,
.version-select.disabled {
  background-color: #f3f4f6;
  color: #9ca3af;
  cursor: not-allowed;
  opacity: 0.7;
}

.version-select option {
  padding: 0.5rem;
  color: #374151;
}

/* Loading Spinner */
.loading-spinner {
  display: flex;
  align-items: center;
  justify-content: center;
}

.spinner {
  width: 20px;
  height: 20px;
  border: 2px solid #e5e7eb;
  border-top: 2px solid #2196f3;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}

/* Update Indicator */
.update-indicator {
  display: flex;
  align-items: center;
}

.update-badge {
  padding: 0.25rem 0.75rem;
  background: linear-gradient(135deg, #ff9800, #f57c00);
  color: white;
  font-size: 0.75rem;
  font-weight: 600;
  border-radius: 12px;
  box-shadow: 0 2px 4px rgba(255, 152, 0, 0.3);
  animation: pulse 2s infinite;
}

@keyframes pulse {
  0%, 100% { transform: scale(1); opacity: 1; }
  50% { transform: scale(1.05); opacity: 0.9; }
}

/* Version Info */
.version-info {
  display: flex;
  flex-direction: column;
  gap: 0.5rem;
  padding: 0.75rem;
  background: rgba(33, 150, 243, 0.05);
  border-radius: 6px;
  border: 1px solid rgba(33, 150, 243, 0.1);
}

.current-version,
.selected-version {
  display: flex;
  align-items: center;
  gap: 0.5rem;
}

.info-label {
  font-size: 0.8rem;
  color: #6b7280;
  font-weight: 500;
  min-width: 60px;
}

.version-tag {
  padding: 0.25rem 0.75rem;
  border-radius: 16px;
  font-size: 0.8rem;
  font-weight: 600;
  letter-spacing: 0.3px;
}

.version-tag.current {
  background: #e3f2fd;
  color: #1565c0;
  border: 1px solid #bbdefb;
}

.version-tag.selected {
  background: #e8f5e8;
  color: #2e7d32;
  border: 1px solid #c8e6c9;
}

/* Update Button */
.update-button-container {
  display: flex;
  justify-content: center;
  margin-top: 0.5rem;
}

.update-button {
  padding: 0.75rem 1.5rem;
  background: #2196f3;
  color: white;
  border: none;
  border-radius: 8px;
  font-size: 0.9rem;
  font-weight: 600;
  cursor: pointer;
  transition: all 0.3s ease;
  box-shadow: 0 2px 4px rgba(33, 150, 243, 0.2);
  min-width: 120px;
}

.update-button:hover:not(:disabled) {
  background: #1976d2;
  transform: translateY(-2px);
  box-shadow: 0 4px 8px rgba(33, 150, 243, 0.3);
}

.update-button:disabled {
  background: #9ca3af;
  cursor: not-allowed;
  transform: none;
  box-shadow: none;
}

.update-button.updating {
  background: #ff9800;
}

.updating-content {
  display: flex;
  align-items: center;
  gap: 0.5rem;
  justify-content: center;
}

.update-spinner {
  width: 16px;
  height: 16px;
  border: 2px solid rgba(255, 255, 255, 0.3);
  border-top: 2px solid white;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

/* Error State */
.error-container {
  padding: 0.75rem;
  background: #fef2f2;
  border: 1px solid #fecaca;
  border-radius: 6px;
  display: flex;
  flex-direction: column;
  gap: 0.75rem;
}

.error-message {
  display: flex;
  align-items: center;
  gap: 0.5rem;
  color: #dc2626;
  font-size: 0.9rem;
}

.error-icon {
  width: 16px;
  height: 16px;
  filter: invert(15%) sepia(92%) saturate(6151%) hue-rotate(357deg) brightness(89%) contrast(114%);
}

.retry-button {
  align-self: flex-start;
  padding: 0.5rem 1rem;
  background: #dc2626;
  color: white;
  border: none;
  border-radius: 6px;
  font-size: 0.8rem;
  font-weight: 600;
  cursor: pointer;
  transition: all 0.2s ease;
}

.retry-button:hover {
  background: #b91c1c;
  transform: translateY(-1px);
}

/* Responsive Design */
@media (max-width: 768px) {
  .version-selector {
    padding: 0.75rem;
    gap: 0.75rem;
  }
  
  .version-info {
    flex-direction: column;
    gap: 0.75rem;
  }
  
  .current-version,
  .selected-version {
    flex-direction: column;
    align-items: flex-start;
    gap: 0.25rem;
  }
  
  .info-label {
    min-width: auto;
  }
}

/* Accessibility */
.version-select:focus-visible {
  outline: 2px solid #2196f3;
  outline-offset: 2px;
}

.update-button:focus-visible {
  outline: 2px solid white;
  outline-offset: 2px;
}

/* High contrast mode support */
@media (prefers-contrast: high) {
  .version-selector {
    border-color: #000;
  }
  
  .version-select {
    border-color: #000;
  }
  
  .update-badge {
    background: #000;
    color: #fff;
  }
}

/* Reduced motion support */
@media (prefers-reduced-motion: reduce) {
  * {
    animation-duration: 0.01ms !important;
    animation-iteration-count: 1 !important;
    transition-duration: 0.01ms !important;
  }
}
</style>