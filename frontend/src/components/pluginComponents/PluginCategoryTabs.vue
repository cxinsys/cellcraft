<template>
  <div class="category-tabs-container">
    <!-- Main Category Tabs -->
    <div class="main-category-tabs">
      <button v-for="category in mainCategories" :key="category.key" @click="setActiveMainCategory(category.key)"
        :class="['main-tab', { active: activeMainCategory === category.key }]">
        <div class="tab-content">
          <img :src="category.icon" :alt="category.label" class="tab-icon" />
          <span class="tab-label">{{ category.label }}</span>
          <span class="tab-count">{{ getCategoryCount(category.key) }}</span>
        </div>
      </button>
    </div>


    <!-- Source Filter (Official/Local) -->
    <div class="source-filter">
      <div class="source-tabs">
        <button v-for="source in sourceFilters" :key="source.key" @click="setActiveSource(source.key)"
          :class="['source-tab', { active: activeSource === source.key }]">
          <span class="source-label">{{ source.label }}</span>
          <span class="source-count">{{ getSourceCount(source.key) }}</span>
        </button>
      </div>
    </div>

    <!-- Resource Filter (GPU/CPU) -->
    <div class="resource-filter">
      <div class="resource-tabs">
        <button v-for="resource in resourceFilters" :key="resource.key" @click="setActiveResource(resource.key)"
          :class="['resource-tab', { active: activeResource === resource.key }]">
          <span class="resource-label">{{ resource.label }}</span>
          <span class="resource-count">{{ getResourceCount(resource.key) }}</span>
        </button>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: "PluginCategoryTabs",
  props: {
    plugins: {
      type: Array,
      default: () => []
    }
  },
  data() {
    return {
      activeMainCategory: 'all',
      activeSource: 'all',
      activeResource: 'all',
      mainCategories: [
        {
          key: 'all',
          label: 'All Plugins',
          icon: require('@/assets/plugins.png')
        },
        {
          key: 'algorithm',
          label: 'Analysis Algorithms',
          icon: require('@/assets/Algorithm_logo.png')
        },
        {
          key: 'visualization',
          label: 'Visualization',
          icon: require('@/assets/Visualization_logo.png')
        }
      ],
      sourceFilters: [
        { key: 'all', label: 'All Sources' },
        { key: 'official', label: 'Official' },
        { key: 'local', label: 'Local' }
      ],
      resourceFilters: [
        { key: 'all', label: 'All Resources' },
        { key: 'gpu', label: 'GPU Required' },
        { key: 'cpu', label: 'CPU Only' }
      ]
    };
  },
  methods: {
    setActiveMainCategory(category) {
      this.activeMainCategory = category;
      this.emitFilters();
    },
    setActiveSource(source) {
      this.activeSource = source;
      this.emitFilters();
    },
    setActiveResource(resource) {
      this.activeResource = resource;
      this.emitFilters();
    },
    getCategoryCount(category) {
      if (category === 'all') return this.plugins.length;

      return this.plugins.filter(plugin => {
        const pluginCategory = this.getPluginCategory(plugin);
        return pluginCategory === category;
      }).length;
    },
    getSourceCount(source) {
      if (source === 'all') return this.plugins.length;

      return this.plugins.filter(plugin => {
        const pluginSource = plugin.source || 'local';
        return pluginSource === source;
      }).length;
    },
    getResourceCount(resource) {
      if (resource === 'all') return this.plugins.length;

      return this.plugins.filter(plugin => {
        const useGpu = plugin.use_gpu || false;
        if (resource === 'gpu') return useGpu;
        if (resource === 'cpu') return !useGpu;
        return true;
      }).length;
    },
    getPluginCategory(plugin) {
      // First check plugin_type field from backend (ANALYSIS/VISUALIZATION)
      if (plugin.plugin_type) {
        const type = plugin.plugin_type.toLowerCase();
        if (type === 'analysis' || type === 'algorithm') {
          return 'algorithm';
        } else if (type === 'visualization') {
          return 'visualization';
        }
      }

      // Legacy support: check category field
      if (plugin.category) {
        return plugin.category;
      }

      // Fallback: analyze plugin name/description for categorization
      const name = plugin.name.toLowerCase();
      const desc = (plugin.description || '').toLowerCase();

      const visualizationKeywords = ['plot', 'chart', 'graph', 'visual', 'heatmap', 'scatter'];
      const isVisualization = visualizationKeywords.some(keyword =>
        name.includes(keyword) || desc.includes(keyword)
      );

      return isVisualization ? 'visualization' : 'algorithm';
    },
    emitFilters() {
      this.$emit('filters-changed', {
        mainCategory: this.activeMainCategory,
        source: this.activeSource,
        resource: this.activeResource
      });
    }
  },
  mounted() {
    // Emit initial filters
    this.emitFilters();
  }
};
</script>

<style scoped>
.category-tabs-container {
  margin: 1rem 0 2rem 0;
  display: flex;
  flex-direction: column;
  gap: 1rem;
}

/* Main Category Tabs */
.main-category-tabs {
  display: flex;
  gap: 1rem;
  padding: 0 5px;
}

.main-tab {
  flex: 1;
  background: white;
  border: 2px solid #e0e0e0;
  border-radius: 12px;
  padding: 1.25rem;
  cursor: pointer;
  transition: all 0.3s ease;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
}

.main-tab:hover {
  border-color: #2196f3;
  box-shadow: 0 4px 16px rgba(33, 150, 243, 0.15);
  transform: translateY(-2px);
}

.main-tab.active {
  background: linear-gradient(135deg, #2196f3 0%, #1976d2 100%);
  border-color: #1976d2;
  color: white;
  box-shadow: 0 6px 20px rgba(33, 150, 243, 0.3);
}

.tab-content {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 0.75rem;
}

.tab-icon {
  width: 2.5rem;
  height: 2.5rem;
  object-fit: contain;
  transition: transform 0.3s ease;
}

.main-tab:hover .tab-icon {
  transform: scale(1.1);
}

.main-tab.active .tab-icon {
  filter: brightness(0) invert(1);
}

.tab-label {
  font-weight: 600;
  font-size: 1rem;
  text-align: center;
  line-height: 1.2;
}

.tab-count {
  background: rgba(0, 0, 0, 0.1);
  color: #666;
  font-size: 0.875rem;
  font-weight: 500;
  padding: 0.25rem 0.75rem;
  border-radius: 20px;
  min-width: 2rem;
  text-align: center;
}

.main-tab.active .tab-count {
  background: rgba(255, 255, 255, 0.2);
  color: white;
}


/* Source Filter */
.source-filter {
  display: flex;
  justify-content: center;
  padding: 0 5px;
}

.source-tabs {
  display: flex;
  gap: 0.5rem;
  background: #f8f9fa;
  padding: 0.25rem;
  border-radius: 8px;
  border: 1px solid #e9ecef;
}

.source-tab {
  background: transparent;
  border: none;
  border-radius: 6px;
  padding: 0.5rem 1rem;
  cursor: pointer;
  transition: all 0.2s ease;
  display: flex;
  align-items: center;
  gap: 0.5rem;
}

.source-tab:hover {
  background: #e3f2fd;
}

.source-tab.active {
  background: white;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

.source-label {
  font-weight: 500;
  font-size: 0.9rem;
}

.source-count {
  background: rgba(0, 0, 0, 0.1);
  color: #666;
  font-size: 0.8rem;
  padding: 0.15rem 0.4rem;
  border-radius: 10px;
  min-width: 1.2rem;
  text-align: center;
}

.source-tab.active .source-count {
  background: #2196f3;
  color: white;
}

/* Resource Filter */
.resource-filter {
  display: flex;
  justify-content: center;
  padding: 0 5px;
}

.resource-tabs {
  display: flex;
  gap: 0.5rem;
  background: #f0f8ff;
  padding: 0.25rem;
  border-radius: 8px;
  border: 1px solid #d0e7ff;
}

.resource-tab {
  background: transparent;
  border: none;
  border-radius: 6px;
  padding: 0.5rem 1rem;
  cursor: pointer;
  transition: all 0.2s ease;
  display: flex;
  align-items: center;
  gap: 0.5rem;
}

.resource-tab:hover {
  background: #e3f2fd;
}

.resource-tab.active {
  background: white;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

.resource-label {
  font-weight: 500;
  font-size: 0.9rem;
}

.resource-count {
  background: rgba(0, 0, 0, 0.1);
  color: #666;
  font-size: 0.8rem;
  padding: 0.15rem 0.4rem;
  border-radius: 10px;
  min-width: 1.2rem;
  text-align: center;
}

.resource-tab.active .resource-count {
  background: #ff9800;
  color: white;
}

/* Responsive Design */
@media (max-width: 768px) {
  .main-category-tabs {
    flex-direction: column;
    gap: 0.5rem;
  }

  .main-tab {
    padding: 1rem;
  }

  .tab-content {
    flex-direction: row;
    justify-content: space-between;
  }

  .tab-icon {
    width: 2rem;
    height: 2rem;
  }

  .sub-category-tabs {
    justify-content: center;
  }

  .sub-tab {
    padding: 0.5rem 1rem;
  }
}
</style>