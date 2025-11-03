// frontend/tests/e2e/fixtures/workflow-api.js

/**
 * Mock API responses for workflow configuration tests
 * These fixtures provide deterministic data for E2E tests
 */

/**
 * Mock response for /api/files/* endpoints
 * Returns folder structure with files
 */
export const mockFilesResponse = {
  folders: [
    ['PBMCLight1000', 1],
    ['TestData', 2],
    ['Results', 3],
  ],
  success: true,
};

/**
 * Mock response for specific folder files
 */
export const mockFolderFiles = {
  PBMCLight1000: {
    files: [
      {
        file_name: 'PBMCLight1000.h5ad',
        file_size: 1048576, // 1MB
        file_type: 'h5ad',
        created_at: '2024-11-01T00:00:00',
      },
      {
        file_name: 'gene_list.txt',
        file_size: 2048,
        file_type: 'txt',
        created_at: '2024-11-01T00:00:00',
      },
    ],
    success: true,
  },
  TestData: {
    files: [
      {
        file_name: 'test_matrix.csv',
        file_size: 524288, // 512KB
        file_type: 'csv',
        created_at: '2024-11-01T00:00:00',
      },
    ],
    success: true,
  },
};

/**
 * Mock response for /api/routes/plugin/list
 * Returns available plugins with their metadata
 */
export const mockPluginsResponse = {
  plugins: [
    {
      name: 'TENET',
      version: '1.0.0',
      description: 'Gene regulatory network reconstruction using TENET',
      inputs: [
        {
          name: 'expression_matrix',
          type: 'h5ad',
          required: true,
          description: 'Single-cell expression matrix',
        },
        {
          name: 'gene_list',
          type: 'txt',
          required: false,
          description: 'Optional gene list',
        },
      ],
      outputs: [
        {
          name: 'grn_results',
          type: 'csv',
          description: 'GRN reconstruction results',
        },
      ],
      parameters: [
        {
          name: 'FDR',
          type: 'number',
          default: 0.05,
          min: 0,
          max: 1,
          description: 'False discovery rate threshold',
        },
        {
          name: 'min_genes',
          type: 'integer',
          default: 10,
          min: 1,
          max: 1000,
          description: 'Minimum number of genes',
        },
        {
          name: 'use_gpu',
          type: 'boolean',
          default: false,
          description: 'Enable GPU acceleration',
        },
      ],
    },
    {
      name: 'GENIE3',
      version: '1.0.0',
      description: 'GRN inference using random forests',
      inputs: [
        {
          name: 'expression_data',
          type: 'h5ad',
          required: true,
          description: 'Expression data matrix',
        },
      ],
      outputs: [
        {
          name: 'network_results',
          type: 'csv',
          description: 'Network inference results',
        },
      ],
      parameters: [
        {
          name: 'n_trees',
          type: 'integer',
          default: 1000,
          min: 100,
          max: 10000,
          description: 'Number of trees in random forest',
        },
      ],
    },
  ],
  success: true,
};

/**
 * Mock response for /api/routes/datatable/load_data
 * Returns paginated table data
 */
export const mockDataTableResponse = {
  columns: [
    { label: 'Gene', field: 'gene', type: 'string' },
    { label: 'Cell_001', field: 'cell_001', type: 'number' },
    { label: 'Cell_002', field: 'cell_002', type: 'number' },
    { label: 'Cell_003', field: 'cell_003', type: 'number' },
  ],
  rows: [
    { gene: 'GENE1', cell_001: 2.5, cell_002: 3.1, cell_003: 1.8 },
    { gene: 'GENE2', cell_001: 0.5, cell_002: 1.2, cell_003: 2.3 },
    { gene: 'GENE3', cell_001: 4.2, cell_002: 3.8, cell_003: 4.5 },
    { gene: 'GENE4', cell_001: 1.1, cell_002: 0.9, cell_003: 1.5 },
    { gene: 'GENE5', cell_001: 3.3, cell_002: 2.7, cell_003: 3.9 },
  ],
  totalRecords: 5,
  success: true,
};

/**
 * Mock response for /api/routes/files/data/:name (ScatterPlot data)
 * Returns UMAP coordinates and metadata
 */
export const mockScatterPlotData = {
  data: {
    X: [1.5, 2.3, -1.2, 0.8, 3.1, -2.5, 1.0, -0.5, 2.8, -1.8],
    Y: [0.5, 1.8, -2.1, 1.5, 0.3, -1.5, 2.3, -0.8, 1.2, -2.8],
    cluster: [0, 0, 1, 0, 2, 1, 2, 1, 2, 1],
    cell_type: [
      'T Cell',
      'T Cell',
      'B Cell',
      'T Cell',
      'Monocyte',
      'B Cell',
      'Monocyte',
      'B Cell',
      'Monocyte',
      'B Cell',
    ],
  },
  columns: ['X', 'Y', 'cluster', 'cell_type'],
  success: true,
};

/**
 * Mock response for /api/routes/files/columns (H5AD column names)
 * Returns available columns for parameter configuration
 */
export const mockH5ADColumns = {
  obs_columns: ['cell_type', 'cluster', 'batch', 'sample_id'],
  var_columns: ['gene_name', 'highly_variable', 'mt', 'ribosomal'],
  obsm_keys: ['X_pca', 'X_umap', 'X_tsne'],
  uns_keys: ['neighbors', 'pca', 'umap'],
  success: true,
};

/**
 * Mock response for ResultFiles modal (workflow/results)
 */
export const mockResultFilesResponse = [
  {
    name: 'tenet_primary_network.csv',
    size: 24576,
    type: 'csv',
  },
  {
    name: 'tenet_primary_scores.json',
    size: 16384,
    type: 'json',
  },
  {
    name: 'tenet_intermediate_heatmap.png',
    size: 40960,
    type: 'png',
  },
];

/**
 * Mock execution manifest payload
 */
export const mockExecutionManifest = {
  manifest_info: {
    format_version: '1.0',
    generated_at: '2025-10-15T09:30:00Z',
    generated_by: 'mock-user',
    description: 'Mock execution manifest for GRNViz visualization test',
  },
  task_metadata: {
    task_id: 'mock-task-1234',
    workflow_id: 42,
    algorithm_id: '12',
    plugin_name: 'TENET',
    task_type: 'compile',
    status: 'SUCCESS',
    start_time: '2025-10-15 09:00:00',
    end_time: '2025-10-15 09:15:00',
    plugin_image_uri: 'ghcr.io/cellcraft/tenet:1.0.0',
  },
  plugin_metadata: {
    name: 'TENET',
    description: 'Mock TENET metadata',
    author: 'mock-user',
    version: '1.0.0',
    plugin_type: 'analysis',
    source: 'mock',
  },
  workflow_metadata: {
    title: 'SUCCESS',
    created_at: '2025-10-15T08:30:00Z',
    updated_at: '2025-10-15T09:15:00Z',
    node_count: 6,
  },
  execution_files: [
    {
      name: 'tenet_primary_network.csv',
      type: 'primary',
      size: 24576,
    },
    {
      name: 'tenet_primary_scores.json',
      type: 'primary',
      size: 16384,
    },
    {
      name: 'tenet_intermediate_heatmap.png',
      type: 'intermediate',
      size: 40960,
    },
  ],
};

/**
 * Mock response for GRNViz visualization execution/result
 */
export const mockVisualizationRunResponse = {
  success: true,
  task_id: 'mock-visualization-task',
  status: 'SUCCESS',
};

export const mockVisualizationResultResponse = {
  data: [
    {
      type: 'bar',
      x: ['Cluster A', 'Cluster B', 'Cluster C'],
      y: [32, 18, 27],
      marker: {
        color: ['#1f77b4', '#ff7f0e', '#2ca02c'],
      },
      name: 'Regulatory score',
    },
  ],
  layout: {
    title: 'Mock GRNViz Bar Plot',
    xaxis: { title: 'Cluster' },
    yaxis: { title: 'Score' },
    showlegend: true,
  },
  success: true,
};

/**
 * Mock response for /routes/task/monitoring (job list)
 */
export const mockTaskMonitoringResponse = {
  data: [
    {
      task_id: 'mock-task-1234',
      workflow_title: 'SUCCESS',
      plugin_name: 'TENET/MockPlugin : v1.0.0',
      plugin_type: 'Analysis',
      start_time: '2025-10-15 09:00:00',
      end_time: '2025-10-15 09:15:00',
      running_time: '00:15:00',
      status: 'SUCCESS',
    },
  ],
  success: true,
};

/**
 * Setup route handlers for workflow APIs
 * @param {Page} page - Playwright page object
 * @param {object} options - Configuration options
 * @param {boolean} options.useDefaultFixtures - Use default mock responses
 * @param {object} options.customResponses - Override specific responses
 */
export async function setupWorkflowRoutes(page, options = {}) {
  const {
    useDefaultFixtures = true,
    customResponses = {},
    enableResultMocks = false,
  } = options;

  if (!useDefaultFixtures && Object.keys(customResponses).length === 0) {
    return; // No mocking needed
  }

  // Mock files endpoint
  await page.route('**/api/routes/files/*', async (route) => {
    const url = route.request().url();

    if (customResponses.files) {
      await route.fulfill({ json: customResponses.files });
    } else if (useDefaultFixtures) {
      // Determine which folder is being requested
      const folderMatch = url.match(/folder=([^&]+)/);
      if (folderMatch) {
        const folderName = decodeURIComponent(folderMatch[1]);
        const response =
          mockFolderFiles[folderName] || { files: [], success: true };
        await route.fulfill({ json: response });
      } else {
        await route.fulfill({ json: mockFilesResponse });
      }
    } else {
      await route.continue();
    }
  });

  // Mock plugin list endpoint
  await page.route('**/api/routes/plugin/list', async (route) => {
    if (customResponses.plugins) {
      await route.fulfill({ json: customResponses.plugins });
    } else if (useDefaultFixtures) {
      await route.fulfill({ json: mockPluginsResponse });
    } else {
      await route.continue();
    }
  });

  // Mock datatable endpoint
  await page.route('**/api/routes/datatable/load_data', async (route) => {
    if (customResponses.datatable) {
      await route.fulfill({ json: customResponses.datatable });
    } else if (useDefaultFixtures) {
      await route.fulfill({ json: mockDataTableResponse });
    } else {
      await route.continue();
    }
  });

  // Mock scatter plot data endpoint
  await page.route('**/api/routes/files/data/*', async (route) => {
    if (customResponses.scatterData) {
      await route.fulfill({ json: customResponses.scatterData });
    } else if (useDefaultFixtures) {
      await route.fulfill({ json: mockScatterPlotData });
    } else {
      await route.continue();
    }
  });

  // Mock H5AD columns endpoint
  await page.route('**/api/routes/files/columns', async (route) => {
    if (customResponses.h5adColumns) {
      await route.fulfill({ json: customResponses.h5adColumns });
    } else if (useDefaultFixtures) {
      await route.fulfill({ json: mockH5ADColumns });
    } else {
      await route.continue();
    }
  });

  // Optional mocks for completed workflow artefacts
  if (enableResultMocks || customResponses.results) {
    await page.route('**/api/routes/workflow/results', async (route) => {
      const payload = customResponses.results || mockResultFilesResponse;
      await route.fulfill({ json: payload });
    });

    await page.route('**/api/routes/workflow/result', async (route) => {
      const bodyData = route.request().postDataJSON?.() || {};
      const fileName = bodyData?.filename || 'mock-result.txt';
      const fileContent = customResponses.resultFileContent || `mock content for ${fileName}`;

      await route.fulfill({
        body: fileContent,
        headers: {
          'Content-Type': 'application/octet-stream',
          'Content-Disposition': `attachment; filename="${fileName}"`,
        },
      });
    });
  }

  if (enableResultMocks || customResponses.manifest) {
    await page.route('**/api/routes/task/*/execution-manifest', async (route) => {
      const manifestPayload = customResponses.manifest || mockExecutionManifest;
      const manifestJson = JSON.stringify(manifestPayload, null, 2);

      await route.fulfill({
        body: manifestJson,
        headers: {
          'Content-Type': 'application/json',
          'Content-Disposition': 'attachment; filename="execution_manifest.json"',
        },
      });
    });
  }

  if (enableResultMocks || customResponses.taskMonitoring) {
    await page.route('**/api/routes/task/monitoring', async (route) => {
      const monitoringPayload = customResponses.taskMonitoring || mockTaskMonitoringResponse;
      await route.fulfill({ json: monitoringPayload });
    });
  }

  if (enableResultMocks || customResponses.visualizationRun || customResponses.visualizationResult) {
    await page.route('**/api/routes/workflow/visualization', async (route) => {
      const runPayload = customResponses.visualizationRun || mockVisualizationRunResponse;
      await route.fulfill({ json: runPayload });
    });

    await page.route('**/api/routes/workflow/visualization/result', async (route) => {
      const resultPayload = customResponses.visualizationResult || mockVisualizationResultResponse;
      await route.fulfill({ json: resultPayload });
    });
  }

  // Mock workflow save/update endpoints (respond with success)
  await page.route('**/api/routes/workflow/*', async (route) => {
    if (route.request().method() === 'POST' || route.request().method() === 'PUT') {
      await route.fulfill({
        json: { success: true, message: 'Workflow updated' },
      });
    } else {
      await route.continue();
    }
  });
}

/**
 * Clear all route handlers
 * @param {Page} page - Playwright page object
 */
export async function clearWorkflowRoutes(page) {
  await page.unroute('**/api/routes/files/*');
  await page.unroute('**/api/routes/plugin/list');
  await page.unroute('**/api/routes/datatable/load_data');
  await page.unroute('**/api/routes/files/data/*');
  await page.unroute('**/api/routes/files/columns');
  await page.unroute('**/api/routes/workflow/results');
  await page.unroute('**/api/routes/workflow/result');
  await page.unroute('**/api/routes/task/*/execution-manifest');
  await page.unroute('**/api/routes/task/monitoring');
  await page.unroute('**/api/routes/workflow/visualization');
  await page.unroute('**/api/routes/workflow/visualization/result');
  await page.unroute('**/api/routes/workflow/*');
}
