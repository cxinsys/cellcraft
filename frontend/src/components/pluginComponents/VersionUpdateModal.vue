<template>
  <div v-if="show" class="version-update-modal-overlay" @click.self="handleClose" @keydown.esc="handleClose">
    <div class="version-update-modal" role="dialog" aria-labelledby="modal-title" aria-modal="true">
      <!-- Modal Header -->
      <div class="modal-header">
        <div class="modal-title-container">
          <h3 id="modal-title" class="modal-title">Update Plugin Version</h3>
          <div class="plugin-info">
            <span class="plugin-name">{{ pluginName }}</span>
            <span class="plugin-badge" :class="getBadgeClass()">
              {{ getPluginTypeLabel() }}
            </span>
          </div>
        </div>
        <button 
          @click="handleClose" 
          class="close-btn" 
          aria-label="Close modal"
          :disabled="loading"
        >
          <img src="@/assets/close.png" alt="Close" class="close-icon" />
        </button>
      </div>

      <!-- Modal Content -->
      <div class="modal-content">
        <!-- Version Comparison -->
        <div class="version-comparison">
          <div class="version-item current">
            <span class="version-label">Current Version</span>
            <span class="version-value">{{ currentVersion || 'Unknown' }}</span>
          </div>
          
          <div class="arrow-container">
            <div class="update-arrow">
              <img src="@/assets/arrow-right.png" alt="Update to" class="arrow-icon" />
            </div>
          </div>
          
          <div class="version-item target">
            <span class="version-label">Target Version</span>
            <span class="version-value">{{ targetVersion || 'Unknown' }}</span>
          </div>
        </div>

        <!-- Update Information -->
        <div v-if="updateInfo && hasUpdateDetails" class="update-info-section">
          <h4 class="update-info-title">Update Information</h4>
          <div class="update-details">
            <!-- Changelog -->
            <div v-if="updateInfo.changelog" class="changelog-container">
              <h5 class="changelog-title">What's New:</h5>
              <ul class="changelog-list">
                <li v-for="(change, index) in getChangelogItems" :key="index">
                  {{ change }}
                </li>
              </ul>
            </div>
            
            <!-- Breaking Changes Warning -->
            <div v-if="updateInfo.breaking_changes" class="breaking-changes-warning">
              <div class="warning-icon">⚠️</div>
              <div class="warning-content">
                <h5 class="warning-title">Breaking Changes</h5>
                <p class="warning-text">{{ updateInfo.breaking_changes }}</p>
              </div>
            </div>
            
            <!-- Additional Notes -->
            <div v-if="updateInfo.notes" class="update-notes">
              <p class="notes-text">{{ updateInfo.notes }}</p>
            </div>
          </div>
        </div>

        <!-- Status Messages -->
        <div v-if="statusMessage" class="status-container" :class="getStatusClass()">
          <div class="status-content">
            <img :src="getStatusIcon()" :alt="statusType" class="status-icon" />
            <span class="status-text">{{ statusMessage }}</span>
          </div>
        </div>

        <!-- Progress Indicator -->
        <div v-if="loading && progress !== null" class="progress-container">
          <div class="progress-bar-container">
            <div class="progress-bar" :style="{ width: `${progress}%` }"></div>
          </div>
          <span class="progress-text">{{ progressText || `${progress}% complete` }}</span>
        </div>
      </div>

      <!-- Modal Footer -->
      <div class="modal-footer">
        <button 
          @click="handleCancel" 
          class="cancel-btn"
          :disabled="loading"
          type="button"
        >
          {{ loading ? 'Updating...' : 'Cancel' }}
        </button>
        
        <button 
          @click="handleConfirm" 
          class="confirm-btn"
          :disabled="loading || !canUpdate"
          :class="{ 'loading': loading }"
          type="button"
          aria-describedby="update-description"
        >
          <div v-if="loading" class="loading-content">
            <div class="loading-spinner"></div>
            <span>Updating...</span>
          </div>
          <span v-else>{{ getConfirmButtonText() }}</span>
        </button>
      </div>

      <!-- Screen reader description -->
      <div id="update-description" class="visually-hidden">
        Update {{ pluginName }} from version {{ currentVersion }} to {{ targetVersion }}
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: "VersionUpdateModal",
  props: {
    show: {
      type: Boolean,
      default: false
    },
    pluginName: {
      type: String,
      required: true
    },
    currentVersion: {
      type: String,
      required: true
    },
    targetVersion: {
      type: String,
      required: true
    },
    updateInfo: {
      type: Object,
      default: () => ({})
    },
    loading: {
      type: Boolean,
      default: false
    },
    progress: {
      type: Number,
      default: null // null means no progress bar
    },
    progressText: {
      type: String,
      default: ""
    },
    statusMessage: {
      type: String,
      default: ""
    },
    statusType: {
      type: String,
      default: "info", // info, success, warning, error
      validator: (value) => ['info', 'success', 'warning', 'error'].includes(value)
    }
  },
  emits: ['close', 'confirm', 'cancel'],
  computed: {
    hasUpdateDetails() {
      return this.updateInfo && (
        this.updateInfo.changelog ||
        this.updateInfo.breaking_changes ||
        this.updateInfo.notes
      );
    },
    
    getChangelogItems() {
      if (!this.updateInfo.changelog) return [];
      
      if (Array.isArray(this.updateInfo.changelog)) {
        return this.updateInfo.changelog;
      }
      
      if (typeof this.updateInfo.changelog === 'string') {
        return this.updateInfo.changelog
          .split('\n')
          .filter(line => line.trim())
          .map(line => line.replace(/^[-*]\s*/, ''));
      }
      
      return [];
    },
    
    canUpdate() {
      return this.pluginName && 
             this.currentVersion && 
             this.targetVersion && 
             this.currentVersion !== this.targetVersion;
    },
    
    isUpgrade() {
      return this.compareVersions(this.targetVersion, this.currentVersion) > 0;
    },
    
    isDowngrade() {
      return this.compareVersions(this.targetVersion, this.currentVersion) < 0;
    }
  },
  watch: {
    show: {
      handler(newValue) {
        if (newValue) {
          this.$nextTick(() => {
            this.focusModal();
            document.addEventListener('keydown', this.handleKeydown);
          });
        } else {
          document.removeEventListener('keydown', this.handleKeydown);
        }
      },
      immediate: true
    }
  },
  beforeDestroy() {
    document.removeEventListener('keydown', this.handleKeydown);
  },
  methods: {
    handleClose() {
      if (!this.loading) {
        this.$emit('close');
      }
    },
    
    handleCancel() {
      if (!this.loading) {
        this.$emit('cancel');
      }
    },
    
    handleConfirm() {
      if (this.canUpdate && !this.loading) {
        this.$emit('confirm', {
          pluginName: this.pluginName,
          fromVersion: this.currentVersion,
          toVersion: this.targetVersion,
          updateInfo: this.updateInfo
        });
      }
    },
    
    handleKeydown(event) {
      if (event.key === 'Escape' && this.show) {
        this.handleClose();
      }
    },
    
    focusModal() {
      // Focus the modal for accessibility
      const modal = this.$el.querySelector('.version-update-modal');
      if (modal) {
        modal.focus();
      }
    },
    
    getBadgeClass() {
      // Determine plugin type for styling
      const pluginSource = this.updateInfo?.source || 'local';
      return pluginSource === 'official' ? 'badge-official' : 'badge-local';
    },
    
    getPluginTypeLabel() {
      const pluginSource = this.updateInfo?.source || 'local';
      return pluginSource === 'official' ? 'OFFICIAL' : 'LOCAL';
    },
    
    getConfirmButtonText() {
      if (this.isUpgrade) {
        return `Upgrade to ${this.targetVersion}`;
      } else if (this.isDowngrade) {
        return `Downgrade to ${this.targetVersion}`;
      } else {
        return `Update to ${this.targetVersion}`;
      }
    },
    
    getStatusClass() {
      return {
        'status-info': this.statusType === 'info',
        'status-success': this.statusType === 'success',
        'status-warning': this.statusType === 'warning',
        'status-error': this.statusType === 'error'
      };
    },
    
    getStatusIcon() {
      switch (this.statusType) {
        case 'success':
          return require('@/assets/check.png');
        case 'warning':
          return require('@/assets/error.png'); // Fallback to error icon
        case 'error':
          return require('@/assets/error.png');
        default:
          return require('@/assets/check.png'); // Fallback to check icon for info
      }
    },
    
    compareVersions(version1, version2) {
      // Simple version comparison - handles semantic versioning
      if (version1 === version2) return 0;
      
      const v1Parts = version1.split('.').map(n => parseInt(n, 10) || 0);
      const v2Parts = version2.split('.').map(n => parseInt(n, 10) || 0);
      
      const maxLength = Math.max(v1Parts.length, v2Parts.length);
      
      for (let i = 0; i < maxLength; i++) {
        const v1Part = v1Parts[i] || 0;
        const v2Part = v2Parts[i] || 0;
        
        if (v1Part > v2Part) return 1;
        if (v1Part < v2Part) return -1;
      }
      
      return 0;
    }
  }
};
</script>

<style scoped>
/* Modal Overlay */
.version-update-modal-overlay {
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: rgba(0, 0, 0, 0.8);
  display: flex;
  justify-content: center;
  align-items: center;
  z-index: 10000;
  backdrop-filter: blur(5px);
  animation: fadeIn 0.2s ease-out;
}

@keyframes fadeIn {
  from { opacity: 0; }
  to { opacity: 1; }
}

/* Modal Container */
.version-update-modal {
  background: #2c3e50;
  border-radius: 16px;
  max-width: 90vw;
  max-height: 90vh;
  width: 600px;
  min-height: 400px;
  display: flex;
  flex-direction: column;
  box-shadow: 0px 10px 30px rgba(0, 0, 0, 0.5);
  border: 1px solid rgba(255, 255, 255, 0.1);
  animation: slideIn 0.3s ease-out;
  outline: none;
  overflow: hidden;
}

@keyframes slideIn {
  from { 
    opacity: 0;
    transform: translateY(-20px) scale(0.95);
  }
  to { 
    opacity: 1;
    transform: translateY(0) scale(1);
  }
}

/* Modal Header */
.modal-header {
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
  padding: 1.5rem;
  border-bottom: 1px solid rgba(255, 255, 255, 0.1);
  background: #1f2a38;
  border-radius: 16px 16px 0 0;
}

.modal-title-container {
  flex: 1;
  display: flex;
  flex-direction: column;
  gap: 0.75rem;
}

.modal-title {
  margin: 0;
  color: #ecf0f1;
  font-weight: 600;
  font-size: 1.2rem;
}

.plugin-info {
  display: flex;
  align-items: center;
  gap: 0.75rem;
}

.plugin-name {
  color: #bdc3c7;
  font-size: 1rem;
  font-weight: 500;
}

.plugin-badge {
  padding: 0.25rem 0.75rem;
  border-radius: 12px;
  font-size: 0.7rem;
  font-weight: 600;
  letter-spacing: 0.5px;
  text-transform: uppercase;
}

.badge-official {
  background: linear-gradient(135deg, #3498db, #2980b9);
  color: white;
  box-shadow: 0 2px 4px rgba(52, 152, 219, 0.3);
}

.badge-local {
  background: linear-gradient(135deg, #95a5a6, #7f8c8d);
  color: white;
  box-shadow: 0 2px 4px rgba(149, 165, 166, 0.3);
}

.close-btn {
  background: #e74c3c;
  color: white;
  border: none;
  padding: 0.75rem;
  border-radius: 8px;
  cursor: pointer;
  transition: all 0.2s ease;
  display: flex;
  align-items: center;
  justify-content: center;
  width: 44px;
  height: 44px;
}

.close-btn:hover:not(:disabled) {
  background: #c0392b;
  transform: scale(1.05);
}

.close-btn:disabled {
  opacity: 0.5;
  cursor: not-allowed;
  transform: none;
}

.close-icon {
  width: 16px;
  height: 16px;
  filter: brightness(0) invert(1);
}

/* Modal Content */
.modal-content {
  flex: 1;
  padding: 1.5rem;
  overflow-y: auto;
  display: flex;
  flex-direction: column;
  gap: 1.5rem;
}

/* Version Comparison */
.version-comparison {
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 1.5rem;
  padding: 1.5rem;
  background: rgba(255, 255, 255, 0.05);
  border-radius: 12px;
  border: 1px solid rgba(255, 255, 255, 0.1);
}

.version-item {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 0.5rem;
  flex: 1;
  max-width: 150px;
}

.version-label {
  font-size: 0.9rem;
  color: #bdc3c7;
  font-weight: 500;
  text-align: center;
}

.version-value {
  padding: 0.75rem 1rem;
  border-radius: 8px;
  font-size: 1rem;
  font-weight: 600;
  text-align: center;
  min-width: 100px;
  word-break: break-all;
}

.version-item.current .version-value {
  background: #e3f2fd;
  color: #1565c0;
  border: 2px solid #bbdefb;
}

.version-item.target .version-value {
  background: #e8f5e8;
  color: #2e7d32;
  border: 2px solid #c8e6c9;
}

.arrow-container {
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 0 1rem;
}

.update-arrow {
  display: flex;
  align-items: center;
  justify-content: center;
  width: 40px;
  height: 40px;
  background: #3498db;
  border-radius: 50%;
  animation: pulse 2s infinite;
}

.arrow-icon {
  width: 20px;
  height: 20px;
  filter: brightness(0) invert(1);
}

@keyframes pulse {
  0%, 100% { transform: scale(1); }
  50% { transform: scale(1.1); }
}

/* Update Information */
.update-info-section {
  background: rgba(255, 255, 255, 0.03);
  border-radius: 12px;
  padding: 1.5rem;
  border: 1px solid rgba(255, 255, 255, 0.08);
}

.update-info-title {
  margin: 0 0 1rem 0;
  color: #ecf0f1;
  font-size: 1rem;
  font-weight: 600;
}

.update-details {
  display: flex;
  flex-direction: column;
  gap: 1.25rem;
}

/* Changelog */
.changelog-container {
  display: flex;
  flex-direction: column;
  gap: 0.75rem;
}

.changelog-title {
  margin: 0;
  color: #bdc3c7;
  font-size: 0.9rem;
  font-weight: 600;
}

.changelog-list {
  margin: 0;
  padding-left: 1.25rem;
  color: #95a5a6;
  font-size: 0.9rem;
  line-height: 1.5;
}

.changelog-list li {
  margin-bottom: 0.5rem;
}

/* Breaking Changes Warning */
.breaking-changes-warning {
  display: flex;
  gap: 1rem;
  padding: 1rem;
  background: rgba(231, 76, 60, 0.1);
  border: 1px solid rgba(231, 76, 60, 0.2);
  border-radius: 8px;
}

.warning-icon {
  font-size: 1.5rem;
  flex-shrink: 0;
}

.warning-content {
  flex: 1;
}

.warning-title {
  margin: 0 0 0.5rem 0;
  color: #e74c3c;
  font-size: 0.9rem;
  font-weight: 600;
}

.warning-text {
  margin: 0;
  color: #e74c3c;
  font-size: 0.85rem;
  line-height: 1.4;
}

/* Update Notes */
.update-notes {
  padding: 1rem;
  background: rgba(52, 152, 219, 0.1);
  border: 1px solid rgba(52, 152, 219, 0.2);
  border-radius: 8px;
}

.notes-text {
  margin: 0;
  color: #3498db;
  font-size: 0.9rem;
  line-height: 1.4;
}

/* Status Messages */
.status-container {
  padding: 1rem;
  border-radius: 8px;
  border: 1px solid;
}

.status-content {
  display: flex;
  align-items: center;
  gap: 0.75rem;
}

.status-icon {
  width: 20px;
  height: 20px;
}

.status-text {
  font-size: 0.9rem;
  font-weight: 500;
}

.status-info {
  background: rgba(52, 152, 219, 0.1);
  border-color: rgba(52, 152, 219, 0.2);
  color: #3498db;
}

.status-success {
  background: rgba(39, 174, 96, 0.1);
  border-color: rgba(39, 174, 96, 0.2);
  color: #27ae60;
}

.status-warning {
  background: rgba(241, 196, 15, 0.1);
  border-color: rgba(241, 196, 15, 0.2);
  color: #f1c40f;
}

.status-error {
  background: rgba(231, 76, 60, 0.1);
  border-color: rgba(231, 76, 60, 0.2);
  color: #e74c3c;
}

/* Progress Indicator */
.progress-container {
  display: flex;
  flex-direction: column;
  gap: 0.75rem;
}

.progress-bar-container {
  width: 100%;
  height: 8px;
  background: rgba(255, 255, 255, 0.1);
  border-radius: 4px;
  overflow: hidden;
}

.progress-bar {
  height: 100%;
  background: linear-gradient(90deg, #3498db, #2ecc71);
  border-radius: 4px;
  transition: width 0.3s ease;
  animation: progressPulse 2s infinite;
}

@keyframes progressPulse {
  0%, 100% { opacity: 1; }
  50% { opacity: 0.8; }
}

.progress-text {
  text-align: center;
  color: #bdc3c7;
  font-size: 0.9rem;
  font-weight: 500;
}

/* Modal Footer */
.modal-footer {
  display: flex;
  justify-content: flex-end;
  gap: 1rem;
  padding: 1.5rem;
  border-top: 1px solid rgba(255, 255, 255, 0.1);
  background: rgba(255, 255, 255, 0.02);
  border-radius: 0 0 16px 16px;
}

.cancel-btn,
.confirm-btn {
  padding: 0.75rem 1.5rem;
  border: none;
  border-radius: 8px;
  font-size: 0.9rem;
  font-weight: 600;
  cursor: pointer;
  transition: all 0.2s ease;
  min-width: 100px;
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 0.5rem;
}

.cancel-btn {
  background: #7f8c8d;
  color: white;
}

.cancel-btn:hover:not(:disabled) {
  background: #6c757d;
  transform: translateY(-1px);
}

.confirm-btn {
  background: #2ecc71;
  color: white;
}

.confirm-btn:hover:not(:disabled) {
  background: #27ae60;
  transform: translateY(-1px);
}

.confirm-btn:disabled,
.cancel-btn:disabled {
  background: #6c757d;
  cursor: not-allowed;
  transform: none;
  opacity: 0.6;
}

.confirm-btn.loading {
  background: #f39c12;
}

.loading-content {
  display: flex;
  align-items: center;
  gap: 0.5rem;
}

.loading-spinner {
  width: 16px;
  height: 16px;
  border: 2px solid rgba(255, 255, 255, 0.3);
  border-top: 2px solid white;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}

/* Accessibility */
.visually-hidden {
  position: absolute;
  width: 1px;
  height: 1px;
  padding: 0;
  margin: -1px;
  overflow: hidden;
  clip: rect(0, 0, 0, 0);
  white-space: nowrap;
  border: 0;
}

/* Focus Management */
.version-update-modal:focus {
  outline: none;
}

.close-btn:focus-visible,
.cancel-btn:focus-visible,
.confirm-btn:focus-visible {
  outline: 2px solid #3498db;
  outline-offset: 2px;
}

/* Responsive Design */
@media (max-width: 768px) {
  .version-update-modal {
    width: 95vw;
    max-height: 95vh;
    margin: 1rem;
  }
  
  .modal-header {
    padding: 1rem;
  }
  
  .modal-content {
    padding: 1rem;
    gap: 1rem;
  }
  
  .modal-footer {
    padding: 1rem;
    flex-direction: column;
  }
  
  .version-comparison {
    flex-direction: column;
    gap: 1rem;
    text-align: center;
  }
  
  .arrow-container {
    transform: rotate(90deg);
  }
  
  .cancel-btn,
  .confirm-btn {
    width: 100%;
    min-width: auto;
  }
}

/* High Contrast Mode */
@media (prefers-contrast: high) {
  .version-update-modal {
    border: 2px solid white;
  }
  
  .plugin-badge {
    border: 1px solid white;
  }
  
  .version-value {
    border-width: 2px;
  }
}

/* Reduced Motion */
@media (prefers-reduced-motion: reduce) {
  .version-update-modal-overlay,
  .version-update-modal,
  .update-arrow,
  .progress-bar,
  .loading-spinner {
    animation: none;
  }
  
  * {
    transition-duration: 0.01ms !important;
  }
}
</style>