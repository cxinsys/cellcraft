<template>
  <div class="recommendations-container">
    <div class="recommendations-header">
      <div class="header-content">
        <h3>
          <img src="@/assets/magic_wand.png" alt="Recommendations" class="header-icon" />
          Recommended Visualizations
        </h3>
        <p class="header-description">
          Based on your data types and analysis results
        </p>
      </div>
      <div class="header-actions">
        <button @click="refreshRecommendations" class="refresh-btn" :disabled="isLoading">
          <img src="@/assets/refresh.png" alt="Refresh" class="btn-icon" />
          Refresh
        </button>
        <button @click="showAdvancedFilters = !showAdvancedFilters" class="filter-btn">
          <img src="@/assets/filter.png" alt="Filter" class="btn-icon" />
          Filters
        </button>
      </div>
    </div>

    <!-- Advanced Filters -->
    <div v-if="showAdvancedFilters" class="advanced-filters">
      <div class="filter-section">
        <label class="filter-label">Data Type:</label>
        <div class="filter-chips">
          <button
            v-for="type in dataTypes"
            :key="type"
            @click="toggleDataTypeFilter(type)"
            :class="['filter-chip', { active: selectedDataTypes.includes(type) }]"
          >
            {{ type }}
          </button>
        </div>
      </div>

      <div class="filter-section">
        <label class="filter-label">Visualization Category:</label>
        <div class="filter-chips">
          <button
            v-for="category in visualizationCategories"
            :key="category"
            @click="toggleCategoryFilter(category)"
            :class="['filter-chip', { active: selectedCategories.includes(category) }]"
          >
            {{ category }}
          </button>
        </div>
      </div>

      <div class="filter-section">
        <label class="filter-label">Complexity:</label>
        <div class="complexity-slider">
          <input
            v-model="complexityRange"
            type="range"
            min="1"
            max="5"
            class="slider"
          />
          <div class="complexity-labels">
            <span>Simple</span>
            <span>Advanced</span>
          </div>
        </div>
      </div>
    </div>

    <!-- Loading State -->
    <div v-if="isLoading" class="loading-state">
      <div class="loading-animation">
        <div class="loading-spinner"></div>
        <p>Analyzing your data and generating recommendations...</p>
      </div>
    </div>

    <!-- Recommendations Grid -->
    <div v-else class="recommendations-grid">
      <div
        v-for="recommendation in filteredRecommendations"
        :key="recommendation.id"
        class="recommendation-card"
        :class="{ 'highly-recommended': recommendation.score >= 0.8 }"
      >
        <!-- Card Header -->
        <div class="card-header">
          <div class="plugin-info">
            <img :src="recommendation.icon" :alt="recommendation.name" class="plugin-icon" />
            <div class="plugin-details">
              <h4 class="plugin-name">{{ recommendation.name }}</h4>
              <p class="plugin-category">{{ recommendation.category }}</p>
            </div>
          </div>
          <div class="recommendation-score">
            <div class="score-circle" :class="getScoreClass(recommendation.score)">
              {{ Math.round(recommendation.score * 100) }}%
            </div>
          </div>
        </div>

        <!-- Card Content -->
        <div class="card-content">
          <p class="plugin-description">{{ recommendation.description }}</p>

          <!-- Why Recommended -->
          <div class="recommendation-reasons">
            <h5>Why recommended:</h5>
            <ul class="reasons-list">
              <li v-for="reason in recommendation.reasons" :key="reason" class="reason-item">
                <span class="reason-icon">{{ getReasonIcon(reason) }}</span>
                <span class="reason-text">{{ reason }}</span>
              </li>
            </ul>
          </div>

          <!-- Compatible Data Types -->
          <div class="data-compatibility">
            <h5>Works with:</h5>
            <div class="data-type-tags">
              <span
                v-for="dataType in recommendation.compatibleDataTypes"
                :key="dataType"
                class="data-type-tag"
                :class="getDataTypeClass(dataType)"
              >
                {{ dataType }}
              </span>
            </div>
          </div>

          <!-- Required Files -->
          <div v-if="recommendation.requiredFiles.length > 0" class="required-files">
            <h5>Required files:</h5>
            <div class="file-requirements">
              <div
                v-for="file in recommendation.requiredFiles"
                :key="file.name"
                class="file-requirement"
                :class="{ available: isFileAvailable(file) }"
              >
                <img :src="getFileTypeIcon(file.type)" :alt="file.type" class="file-icon" />
                <div class="file-info">
                  <span class="file-name">{{ file.name }}</span>
                  <span class="file-description">{{ file.description }}</span>
                </div>
                <div class="availability-status">
                  {{ isFileAvailable(file) ? '✓' : '!' }}
                </div>
              </div>
            </div>
          </div>
        </div>

        <!-- Card Actions -->
        <div class="card-actions">
          <button
            @click="previewVisualization(recommendation)"
            class="action-btn preview-btn"
            :disabled="!canPreview(recommendation)"
          >
            <img src="@/assets/eye.png" alt="Preview" class="btn-icon" />
            Quick Preview
          </button>
          <button
            @click="useVisualization(recommendation)"
            class="action-btn use-btn"
            :disabled="!canUse(recommendation)"
          >
            <img src="@/assets/plus.png" alt="Use" class="btn-icon" />
            Add to Workflow
          </button>
        </div>

        <!-- Highly Recommended Badge -->
        <div v-if="recommendation.score >= 0.8" class="recommended-badge">
          <img src="@/assets/star.png" alt="Highly Recommended" class="badge-icon" />
          <span>Highly Recommended</span>
        </div>
      </div>
    </div>

    <!-- No Recommendations -->
    <div v-if="!isLoading && filteredRecommendations.length === 0" class="no-recommendations">
      <img src="@/assets/no_results.png" alt="No recommendations" class="no-results-icon" />
      <h4>No matching visualizations found</h4>
      <p>Try adjusting your filters or upload more data files.</p>
      <button @click="resetFilters" class="reset-filters-btn">
        Reset Filters
      </button>
    </div>

    <!-- Preview Modal -->
    <div v-if="showPreview" class="preview-modal-overlay" @click.self="closePreview">
      <div class="preview-modal">
        <div class="preview-header">
          <h3>{{ selectedRecommendation.name }} Preview</h3>
          <button @click="closePreview" class="close-btn">
            <img src="@/assets/close.png" alt="Close" class="close-icon" />
          </button>
        </div>
        <div class="preview-content">
          <div class="preview-image">
            <img :src="selectedRecommendation.previewImage" :alt="selectedRecommendation.name" />
          </div>
          <div class="preview-details">
            <h4>Features:</h4>
            <ul class="feature-list">
              <li v-for="feature in selectedRecommendation.features" :key="feature">
                {{ feature }}
              </li>
            </ul>
            <h4>Configuration Options:</h4>
            <div class="config-options">
              <div v-for="option in selectedRecommendation.configOptions" :key="option.name" class="config-option">
                <label>{{ option.label }}:</label>
                <select v-if="option.type === 'select'" v-model="option.value">
                  <option v-for="choice in option.choices" :key="choice" :value="choice">
                    {{ choice }}
                  </option>
                </select>
                <input v-else :type="option.type" v-model="option.value" />
              </div>
            </div>
          </div>
        </div>
        <div class="preview-actions">
          <button @click="closePreview" class="cancel-btn">Cancel</button>
          <button @click="useFromPreview" class="use-btn">Add to Workflow</button>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: "VisualizationRecommendations",
  props: {
    availableFiles: {
      type: Array,
      default: () => []
    },
    workflowData: {
      type: Object,
      default: () => ({})
    },
    algorithmResults: {
      type: Object,
      default: () => ({})
    }
  },
  data() {
    return {
      isLoading: false,
      showAdvancedFilters: false,
      showPreview: false,
      selectedRecommendation: null,
      recommendations: [],
      dataTypes: ['Network', 'Expression', 'Metadata', 'Statistics', 'Trajectory'],
      visualizationCategories: ['Network', 'Heatmap', 'Scatter Plot', 'Bar Chart', 'Time Series'],
      selectedDataTypes: [],
      selectedCategories: [],
      complexityRange: 3,
      // Sample recommendations data
      sampleRecommendations: [
        {
          id: 'network_vis',
          name: 'Interactive Network Visualization',
          category: 'Network',
          description: 'Interactive force-directed network graph with node filtering and clustering.',
          score: 0.95,
          icon: require('@/assets/network_icon.png'),
          compatibleDataTypes: ['Network', 'Expression'],
          reasons: [
            'Perfect match for network data',
            'Supports large networks (>1000 nodes)',
            'Interactive filtering capabilities'
          ],
          requiredFiles: [
            {
              name: 'Network Edges',
              type: 'csv',
              description: 'Source, target, weight columns',
              pattern: /edges?\.csv$/i
            },
            {
              name: 'Node Attributes',
              type: 'csv',
              description: 'Node properties and metadata',
              pattern: /nodes?\.csv$/i
            }
          ],
          features: [
            'Force-directed layout',
            'Node size by degree/expression',
            'Edge filtering by weight',
            'Community detection',
            'Export to multiple formats'
          ],
          configOptions: [
            {
              name: 'layout',
              label: 'Layout Algorithm',
              type: 'select',
              value: 'force-directed',
              choices: ['force-directed', 'circular', 'hierarchical']
            },
            {
              name: 'nodeSize',
              label: 'Node Size Attribute',
              type: 'select',
              value: 'degree',
              choices: ['degree', 'expression', 'betweenness']
            }
          ],
          previewImage: require('@/assets/preview_network.png')
        },
        {
          id: 'expression_heatmap',
          name: 'Expression Heatmap',
          category: 'Heatmap',
          description: 'Clustered heatmap for gene expression data with dendrograms.',
          score: 0.88,
          icon: require('@/assets/heatmap_icon.png'),
          compatibleDataTypes: ['Expression', 'Metadata'],
          reasons: [
            'Excellent for expression patterns',
            'Built-in clustering algorithms',
            'Handles missing values'
          ],
          requiredFiles: [
            {
              name: 'Expression Matrix',
              type: 'csv',
              description: 'Genes × samples expression matrix',
              pattern: /expression|matrix\.csv$/i
            }
          ],
          features: [
            'Hierarchical clustering',
            'Custom color schemes',
            'Gene/sample filtering',
            'Annotation tracks',
            'Interactive zoom'
          ],
          configOptions: [
            {
              name: 'clustering',
              label: 'Clustering Method',
              type: 'select',
              value: 'ward',
              choices: ['ward', 'complete', 'average']
            },
            {
              name: 'colorScale',
              label: 'Color Scale',
              type: 'select',
              value: 'RdBu',
              choices: ['RdBu', 'viridis', 'plasma']
            }
          ],
          previewImage: require('@/assets/preview_heatmap.png')
        },
        {
          id: 'trajectory_plot',
          name: 'Trajectory Analysis',
          category: 'Time Series',
          description: 'Pseudotime trajectory visualization with branching detection.',
          score: 0.75,
          icon: require('@/assets/trajectory_icon.png'),
          compatibleDataTypes: ['Trajectory', 'Expression', 'Metadata'],
          reasons: [
            'Designed for trajectory data',
            'Automatic branching detection',
            'Gene expression overlay'
          ],
          requiredFiles: [
            {
              name: 'Trajectory Data',
              type: 'csv',
              description: 'Cell pseudotime and coordinates',
              pattern: /trajectory|pseudotime\.csv$/i
            }
          ],
          features: [
            'Branching trajectory paths',
            'Gene expression overlay',
            'Cell type annotation',
            'Interactive path selection',
            '3D trajectory option'
          ],
          configOptions: [
            {
              name: 'dimension',
              label: 'Dimensions',
              type: 'select',
              value: '2D',
              choices: ['2D', '3D']
            }
          ],
          previewImage: require('@/assets/preview_trajectory.png')
        }
      ]
    };
  },
  computed: {
    filteredRecommendations() {
      let filtered = [...this.recommendations];

      // Filter by data types
      if (this.selectedDataTypes.length > 0) {
        filtered = filtered.filter(rec =>
          rec.compatibleDataTypes.some(type => this.selectedDataTypes.includes(type))
        );
      }

      // Filter by categories
      if (this.selectedCategories.length > 0) {
        filtered = filtered.filter(rec =>
          this.selectedCategories.includes(rec.category)
        );
      }

      // Filter by complexity
      const maxComplexity = this.complexityRange / 5; // Convert to 0-1 scale
      filtered = filtered.filter(rec => {
        const recComplexity = rec.features.length / 10; // Rough complexity estimate
        return recComplexity <= maxComplexity;
      });

      // Sort by recommendation score
      return filtered.sort((a, b) => b.score - a.score);
    }
  },
  mounted() {
    this.generateRecommendations();
  },
  watch: {
    availableFiles: {
      handler() {
        this.generateRecommendations();
      },
      deep: true
    },
    algorithmResults: {
      handler() {
        this.generateRecommendations();
      },
      deep: true
    }
  },
  methods: {
    async generateRecommendations() {
      this.isLoading = true;
      
      try {
        // Simulate API call for generating recommendations
        await new Promise(resolve => setTimeout(resolve, 1500));
        
        // Analyze available files and algorithm results
        const analysisResults = this.analyzeDataContext();
        
        // Score and rank recommendations based on analysis
        this.recommendations = this.scoreRecommendations(analysisResults);
        
      } catch (error) {
        console.error('Failed to generate recommendations:', error);
      } finally {
        this.isLoading = false;
      }
    },
    analyzeDataContext() {
      const context = {
        hasNetworkData: this.availableFiles.some(f => 
          /network|edge/i.test(f.name) && f.type === 'csv'
        ),
        hasExpressionData: this.availableFiles.some(f => 
          /expression|matrix/i.test(f.name) && ['csv', 'h5ad'].includes(f.type)
        ),
        hasTrajectoryData: this.availableFiles.some(f => 
          /trajectory|pseudotime/i.test(f.name)
        ),
        hasMetadata: this.availableFiles.some(f => 
          /metadata|sample/i.test(f.name)
        ),
        algorithmType: this.algorithmResults.type || 'unknown',
        dataSize: this.availableFiles.reduce((sum, f) => sum + (f.size || 0), 0)
      };
      
      return context;
    },
    scoreRecommendations(context) {
      return this.sampleRecommendations.map(rec => {
        let score = rec.score; // Base score
        
        // Boost score based on available data
        if (rec.id === 'network_vis' && context.hasNetworkData) score += 0.1;
        if (rec.id === 'expression_heatmap' && context.hasExpressionData) score += 0.1;
        if (rec.id === 'trajectory_plot' && context.hasTrajectoryData) score += 0.15;
        
        // Reduce score if required files are missing
        const missingFiles = rec.requiredFiles.filter(req => !this.isFileAvailable(req));
        score -= missingFiles.length * 0.1;
        
        return {
          ...rec,
          score: Math.max(0, Math.min(1, score)) // Clamp between 0 and 1
        };
      });
    },
    isFileAvailable(requiredFile) {
      return this.availableFiles.some(file => 
        requiredFile.pattern.test(file.name) && 
        file.type === requiredFile.type
      );
    },
    canPreview(recommendation) {
      // Check if we have enough data to show a meaningful preview
      const hasRequiredFiles = recommendation.requiredFiles.every(req => 
        this.isFileAvailable(req)
      );
      return hasRequiredFiles || recommendation.score > 0.7;
    },
    canUse(recommendation) {
      // Check if we can actually use this visualization
      const hasRequiredFiles = recommendation.requiredFiles.every(req => 
        this.isFileAvailable(req)
      );
      return hasRequiredFiles;
    },
    previewVisualization(recommendation) {
      this.selectedRecommendation = { ...recommendation };
      this.showPreview = true;
    },
    closePreview() {
      this.showPreview = false;
      this.selectedRecommendation = null;
    },
    useVisualization(recommendation) {
      this.$emit('visualization-selected', {
        pluginId: recommendation.id,
        pluginName: recommendation.name,
        configuration: recommendation.configOptions?.reduce((acc, opt) => {
          acc[opt.name] = opt.value;
          return acc;
        }, {}) || {}
      });
    },
    useFromPreview() {
      this.useVisualization(this.selectedRecommendation);
      this.closePreview();
    },
    async refreshRecommendations() {
      await this.generateRecommendations();
    },
    toggleDataTypeFilter(type) {
      const index = this.selectedDataTypes.indexOf(type);
      if (index === -1) {
        this.selectedDataTypes.push(type);
      } else {
        this.selectedDataTypes.splice(index, 1);
      }
    },
    toggleCategoryFilter(category) {
      const index = this.selectedCategories.indexOf(category);
      if (index === -1) {
        this.selectedCategories.push(category);
      } else {
        this.selectedCategories.splice(index, 1);
      }
    },
    resetFilters() {
      this.selectedDataTypes = [];
      this.selectedCategories = [];
      this.complexityRange = 3;
    },
    getScoreClass(score) {
      if (score >= 0.8) return 'score-high';
      if (score >= 0.6) return 'score-medium';
      return 'score-low';
    },
    getReasonIcon(reason) {
      if (reason.includes('Perfect match') || reason.includes('Excellent')) return '🎯';
      if (reason.includes('Supports') || reason.includes('Built-in')) return '⚡';
      if (reason.includes('Interactive') || reason.includes('capabilities')) return '🔧';
      return '✓';
    },
    getDataTypeClass(dataType) {
      const classMap = {
        'Network': 'network',
        'Expression': 'expression',
        'Metadata': 'metadata',
        'Statistics': 'stats',
        'Trajectory': 'trajectory'
      };
      return classMap[dataType] || 'default';
    },
    getFileTypeIcon(fileType) {
      const iconMap = {
        'csv': require('@/assets/csv_icon.png'),
        'h5ad': require('@/assets/h5ad_icon.png'),
        'json': require('@/assets/json_icon.png')
      };
      return iconMap[fileType] || require('@/assets/file_icon.png');
    }
  }
};
</script>

<style scoped>
.recommendations-container {
  background: white;
  border-radius: 16px;
  padding: 2rem;
  box-shadow: 0 8px 32px rgba(0, 0, 0, 0.08);
  margin: 2rem 0;
}

.recommendations-header {
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
  margin-bottom: 1.5rem;
}

.header-content h3 {
  display: flex;
  align-items: center;
  gap: 0.75rem;
  margin: 0 0 0.5rem 0;
  font-size: 1.5rem;
  font-weight: 600;
  color: #333;
}

.header-icon {
  width: 1.75rem;
  height: 1.75rem;
}

.header-description {
  margin: 0;
  color: #666;
  font-size: 1rem;
}

.header-actions {
  display: flex;
  gap: 1rem;
}

.refresh-btn,
.filter-btn {
  display: flex;
  align-items: center;
  gap: 0.5rem;
  padding: 0.75rem 1rem;
  border: 1px solid #e0e0e0;
  border-radius: 8px;
  background: white;
  cursor: pointer;
  font-size: 0.875rem;
  font-weight: 500;
  transition: all 0.2s ease;
}

.refresh-btn:hover,
.filter-btn:hover {
  border-color: #2196f3;
  background: #f3f8ff;
}

.btn-icon {
  width: 1rem;
  height: 1rem;
}

/* Advanced Filters */
.advanced-filters {
  background: #f8f9fa;
  border-radius: 12px;
  padding: 1.5rem;
  margin-bottom: 2rem;
  border: 1px solid #e9ecef;
}

.filter-section {
  margin-bottom: 1.5rem;
}

.filter-section:last-child {
  margin-bottom: 0;
}

.filter-label {
  display: block;
  font-weight: 500;
  color: #333;
  margin-bottom: 0.75rem;
}

.filter-chips {
  display: flex;
  gap: 0.5rem;
  flex-wrap: wrap;
}

.filter-chip {
  padding: 0.5rem 1rem;
  border: 2px solid #e0e0e0;
  border-radius: 20px;
  background: white;
  cursor: pointer;
  font-size: 0.875rem;
  font-weight: 500;
  transition: all 0.2s ease;
}

.filter-chip:hover {
  border-color: #2196f3;
  background: #f3f8ff;
}

.filter-chip.active {
  background: #2196f3;
  color: white;
  border-color: #2196f3;
}

.complexity-slider {
  max-width: 300px;
}

.slider {
  width: 100%;
  height: 6px;
  border-radius: 5px;
  background: #d3d3d3;
  outline: none;
  -webkit-appearance: none;
}

.slider::-webkit-slider-thumb {
  -webkit-appearance: none;
  appearance: none;
  width: 20px;
  height: 20px;
  border-radius: 50%;
  background: #2196f3;
  cursor: pointer;
}

.complexity-labels {
  display: flex;
  justify-content: space-between;
  margin-top: 0.5rem;
  font-size: 0.8rem;
  color: #666;
}

/* Loading State */
.loading-state {
  text-align: center;
  padding: 4rem 2rem;
}

.loading-animation {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 1.5rem;
}

.loading-spinner {
  width: 3rem;
  height: 3rem;
  border: 4px solid #e3f2fd;
  border-top: 4px solid #2196f3;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}

/* Recommendations Grid */
.recommendations-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
  gap: 2rem;
}

.recommendation-card {
  background: white;
  border: 2px solid #e9ecef;
  border-radius: 16px;
  overflow: hidden;
  transition: all 0.3s ease;
  position: relative;
}

.recommendation-card:hover {
  transform: translateY(-4px);
  box-shadow: 0 12px 32px rgba(0, 0, 0, 0.15);
  border-color: #2196f3;
}

.recommendation-card.highly-recommended {
  border-color: #4caf50;
  background: linear-gradient(135deg, #ffffff 0%, #f8fff8 100%);
}

.card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 1.5rem;
  background: #f8f9fa;
  border-bottom: 1px solid #e9ecef;
}

.plugin-info {
  display: flex;
  align-items: center;
  gap: 1rem;
}

.plugin-icon {
  width: 3rem;
  height: 3rem;
  border-radius: 8px;
  object-fit: cover;
}

.plugin-name {
  margin: 0 0 0.25rem 0;
  font-size: 1.2rem;
  font-weight: 600;
  color: #333;
}

.plugin-category {
  margin: 0;
  font-size: 0.875rem;
  color: #666;
  text-transform: uppercase;
  letter-spacing: 0.5px;
}

.recommendation-score {
  display: flex;
  align-items: center;
}

.score-circle {
  width: 3rem;
  height: 3rem;
  border-radius: 50%;
  display: flex;
  align-items: center;
  justify-content: center;
  font-weight: 600;
  font-size: 0.875rem;
  color: white;
}

.score-circle.score-high {
  background: #4caf50;
}

.score-circle.score-medium {
  background: #ff9800;
}

.score-circle.score-low {
  background: #f44336;
}

.card-content {
  padding: 1.5rem;
}

.plugin-description {
  margin: 0 0 1.5rem 0;
  color: #666;
  line-height: 1.5;
}

.recommendation-reasons,
.data-compatibility,
.required-files {
  margin-bottom: 1.5rem;
}

.recommendation-reasons h5,
.data-compatibility h5,
.required-files h5 {
  margin: 0 0 0.75rem 0;
  font-size: 1rem;
  font-weight: 600;
  color: #333;
}

.reasons-list {
  list-style: none;
  padding: 0;
  margin: 0;
}

.reason-item {
  display: flex;
  align-items: center;
  gap: 0.75rem;
  margin-bottom: 0.5rem;
  padding: 0.5rem;
  background: #f8f9fa;
  border-radius: 6px;
}

.reason-icon {
  font-size: 1rem;
}

.reason-text {
  flex: 1;
  font-size: 0.875rem;
  color: #333;
}

.data-type-tags {
  display: flex;
  gap: 0.5rem;
  flex-wrap: wrap;
}

.data-type-tag {
  padding: 0.25rem 0.75rem;
  border-radius: 20px;
  font-size: 0.8rem;
  font-weight: 500;
}

.data-type-tag.network {
  background: #e3f2fd;
  color: #1565c0;
}

.data-type-tag.expression {
  background: #e8f5e9;
  color: #2e7d32;
}

.data-type-tag.metadata {
  background: #fff3e0;
  color: #e65100;
}

.data-type-tag.stats {
  background: #f3e5f5;
  color: #7b1fa2;
}

.data-type-tag.trajectory {
  background: #ffebee;
  color: #c62828;
}

.data-type-tag.default {
  background: #f5f5f5;
  color: #616161;
}

.file-requirements {
  display: flex;
  flex-direction: column;
  gap: 0.75rem;
}

.file-requirement {
  display: flex;
  align-items: center;
  gap: 0.75rem;
  padding: 0.75rem;
  background: #f8f9fa;
  border-radius: 8px;
  border-left: 4px solid #d1d5db;
}

.file-requirement.available {
  border-left-color: #4caf50;
  background: #f8fff8;
}

.file-icon {
  width: 1.5rem;
  height: 1.5rem;
  flex-shrink: 0;
}

.file-info {
  flex: 1;
  display: flex;
  flex-direction: column;
}

.file-name {
  font-weight: 500;
  color: #333;
}

.file-description {
  font-size: 0.8rem;
  color: #666;
  margin-top: 0.25rem;
}

.availability-status {
  font-size: 1.2rem;
  font-weight: 600;
}

.file-requirement.available .availability-status {
  color: #4caf50;
}

.file-requirement:not(.available) .availability-status {
  color: #f44336;
}

.card-actions {
  display: flex;
  gap: 1rem;
  padding: 1.5rem;
  background: #f8f9fa;
  border-top: 1px solid #e9ecef;
}

.action-btn {
  flex: 1;
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 0.5rem;
  padding: 0.75rem 1rem;
  border: none;
  border-radius: 8px;
  font-weight: 500;
  cursor: pointer;
  transition: all 0.2s ease;
}

.preview-btn {
  background: #f5f5f5;
  color: #666;
}

.preview-btn:hover:not(:disabled) {
  background: #e9ecef;
}

.use-btn {
  background: #2196f3;
  color: white;
}

.use-btn:hover:not(:disabled) {
  background: #1976d2;
  transform: translateY(-1px);
}

.action-btn:disabled {
  opacity: 0.5;
  cursor: not-allowed;
  transform: none;
}

.recommended-badge {
  position: absolute;
  top: 1rem;
  right: 1rem;
  background: linear-gradient(135deg, #4caf50 0%, #45a049 100%);
  color: white;
  padding: 0.5rem 1rem;
  border-radius: 20px;
  font-size: 0.8rem;
  font-weight: 500;
  display: flex;
  align-items: center;
  gap: 0.5rem;
  box-shadow: 0 4px 12px rgba(76, 175, 80, 0.3);
}

.badge-icon {
  width: 1rem;
  height: 1rem;
  filter: brightness(0) invert(1);
}

/* No Recommendations */
.no-recommendations {
  text-align: center;
  padding: 4rem 2rem;
  color: #666;
}

.no-results-icon {
  width: 4rem;
  height: 4rem;
  opacity: 0.5;
  margin-bottom: 1.5rem;
}

.no-recommendations h4 {
  margin: 0 0 1rem 0;
  font-size: 1.2rem;
  color: #333;
}

.reset-filters-btn {
  margin-top: 1.5rem;
  padding: 0.75rem 1.5rem;
  background: #2196f3;
  color: white;
  border: none;
  border-radius: 8px;
  cursor: pointer;
  font-weight: 500;
  transition: all 0.2s ease;
}

.reset-filters-btn:hover {
  background: #1976d2;
}

/* Preview Modal */
.preview-modal-overlay {
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
}

.preview-modal {
  background: white;
  border-radius: 16px;
  width: 90vw;
  max-width: 900px;
  height: 80vh;
  max-height: 700px;
  display: flex;
  flex-direction: column;
  box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
  overflow: hidden;
}

.preview-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 1.5rem;
  background: #f8f9fa;
  border-bottom: 1px solid #e9ecef;
}

.preview-header h3 {
  margin: 0;
  font-size: 1.25rem;
  font-weight: 600;
  color: #333;
}

.close-btn {
  background: none;
  border: none;
  cursor: pointer;
  padding: 0.5rem;
  border-radius: 50%;
  transition: background-color 0.2s;
}

.close-btn:hover {
  background: #e9ecef;
}

.close-icon {
  width: 1.2rem;
  height: 1.2rem;
}

.preview-content {
  flex: 1;
  display: grid;
  grid-template-columns: 1fr 1fr;
  overflow: hidden;
}

.preview-image {
  background: #f8f9fa;
  display: flex;
  align-items: center;
  justify-content: center;
  overflow: hidden;
}

.preview-image img {
  max-width: 100%;
  max-height: 100%;
  object-fit: contain;
}

.preview-details {
  padding: 2rem;
  overflow-y: auto;
}

.preview-details h4 {
  margin: 0 0 1rem 0;
  font-size: 1.1rem;
  font-weight: 600;
  color: #333;
}

.feature-list {
  list-style: none;
  padding: 0;
  margin: 0 0 2rem 0;
}

.feature-list li {
  padding: 0.5rem 0;
  border-bottom: 1px solid #f0f0f0;
}

.config-options {
  display: flex;
  flex-direction: column;
  gap: 1rem;
}

.config-option {
  display: flex;
  flex-direction: column;
  gap: 0.5rem;
}

.config-option label {
  font-weight: 500;
  color: #333;
}

.config-option select,
.config-option input {
  padding: 0.5rem;
  border: 1px solid #d1d5db;
  border-radius: 6px;
  font-size: 0.875rem;
}

.preview-actions {
  display: flex;
  justify-content: flex-end;
  gap: 1rem;
  padding: 1.5rem;
  background: #f8f9fa;
  border-top: 1px solid #e9ecef;
}

.cancel-btn {
  padding: 0.75rem 1.5rem;
  background: #f5f5f5;
  color: #666;
  border: none;
  border-radius: 8px;
  cursor: pointer;
  font-weight: 500;
  transition: all 0.2s ease;
}

.cancel-btn:hover {
  background: #e9ecef;
}

.preview-actions .use-btn {
  padding: 0.75rem 1.5rem;
  background: #4caf50;
  color: white;
  border: none;
  border-radius: 8px;
  cursor: pointer;
  font-weight: 500;
  transition: all 0.2s ease;
}

.preview-actions .use-btn:hover {
  background: #45a049;
  transform: translateY(-1px);
}

/* Responsive Design */
@media (max-width: 768px) {
  .recommendations-grid {
    grid-template-columns: 1fr;
  }
  
  .recommendations-header {
    flex-direction: column;
    gap: 1rem;
    align-items: flex-start;
  }
  
  .header-actions {
    width: 100%;
    justify-content: flex-start;
  }
  
  .preview-content {
    grid-template-columns: 1fr;
    grid-template-rows: 1fr 1fr;
  }
  
  .filter-chips {
    justify-content: center;
  }
}
</style>