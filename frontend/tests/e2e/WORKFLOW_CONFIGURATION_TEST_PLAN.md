# Workflow Configuration Test Plan (Section 2.2)

## Overview

This document outlines the implementation plan for E2E tests covering workflow configuration scenarios (Section 2.2 of the test strategy). The tests verify node configuration, file assignment, parameter settings, and data visualization within the CellCraft workflow editor.

## Test Priority Matrix

Based on QA analysis (Codex review), tests are prioritized by:
- **Ease of validation** (DOM accessibility, stable selectors)
- **Test independence** (minimal dependencies)
- **Execution speed**
- **Critical user journey impact**

| Priority | Test Scenario | Status | Dependencies | Estimated Effort |
|----------|---------------|--------|--------------|------------------|
| 1 | InputFile Assignment | ✅ COMPLETE | None | 3h |
| 2 | Dropdown File Selection | 🔜 NEXT | Priority 1 | 2h |
| 3 | Algorithm Parameter Configuration | 📋 PLANNED | Priority 2 | 3h |
| 4 | DataTable Matrix Display | 📋 PLANNED | Priority 1 | 4h |
| 5 | ScatterPlot UMAP Visualization | 📋 PLANNED | Priority 1 | 5h |

**Legend:**
- ✅ COMPLETE: Implemented and passing
- 🔜 NEXT: Ready to implement
- 📋 PLANNED: Scheduled for implementation

---

## Implemented Infrastructure

### Page Objects

#### 1. WorkflowPage (`pages/WorkflowPage.js`)
Main page object for the Drawflow workflow editor.

**Key Methods:**
- `goto(workflowId)` - Navigate to workflow by ID
- `verifyPageLoaded()` - Verify editor is ready
- `getCanvasNodes()` - Get all nodes on canvas
- `findNodeByType(nodeType)` - Find specific node
- `openNodeModal(nodeType)` - Open node configuration modal
- `getCurrentTab()` - Get active tab information
- `switchTab(index)` - Switch between tabs
- `closeTab(name)` - Close a tab
- `dragNodeToCanvas(nodeType, position)` - Add node to canvas
- `connectNodes(source, target)` - Create node connections
- `verifyConnection(source, target)` - Verify connection exists
- `waitForMessage(text, timeout)` - Wait for notification
- `getTabCount()` - Get number of open tabs
- `getAvailableNodeTypes()` - Get node types from sidebar

**Usage Example:**
```javascript
await workflowPage.goto(workflowId);
await workflowPage.verifyPageLoaded();
await workflowPage.openNodeModal('InputFile');
```

#### 2. InputFileModal (`pages/modals/InputFileModal.js`)
Page object for InputFile node configuration modal.

**Key Methods:**
- `verifyModalOpen()` - Verify modal is visible
- `getFolders()` - Get available folders
- `selectFolder(name)` - Select a folder
- `verifyFolderSelected(name)` - Verify folder highlighted
- `getFiles()` - Get files in current folder
- `selectFile(name)` - Select a file
- `verifyFileSelected(name)` - Verify file highlighted
- `getCurrentFileInfo()` - Get selected file details
- `verifyCurrentFile(name)` - Verify current file display
- `clickApply()` - Apply file assignment
- `verifyApplyButtonState(applied)` - Verify button state
- `waitForLoading(timeout)` - Wait for loading indicators
- `assignFile(folder, file)` - Complete assignment workflow

**Usage Example:**
```javascript
await inputFileModal.verifyModalOpen();
await inputFileModal.assignFile('PBMCLight1000', 'PBMCLight1000.h5ad');
await inputFileModal.verifyApplyButtonState(true);
```

#### 3. AlgorithmModal (`pages/modals/AlgorithmModal.js`)
Page object for Algorithm node configuration modal.

**Key Methods:**
- `verifyModalOpen()` - Verify modal is visible
- `getAvailablePlugins()` - Get plugin list
- `selectPlugin(name)` - Select a plugin
- `getInputDropdown(identifier)` - Get input file dropdown
- `selectInputFile(dropdown, file)` - Assign file to input
- `getSelectedInputFile(dropdown)` - Get selected file
- `verifyFileInDropdown(dropdown, file)` - Verify file available
- `getAvailableFiles(dropdown)` - Get all available files
- `getParameterInput(identifier)` - Get parameter input
- `setParameter(param, value)` - Set parameter value
- `getParameterValue(param)` - Get parameter value
- `verifyDefaultValue(param, value)` - Verify default value
- `getParameterNames()` - Get all parameter names
- `clickApply()` - Apply configuration

**Usage Example:**
```javascript
await algorithmModal.verifyModalOpen();
await algorithmModal.selectInputFile(0, 'PBMCLight1000.h5ad');
await algorithmModal.setParameter('FDR', 0.01);
await algorithmModal.clickApply();
```

### Helper Utilities (`utils/workflow.js`)

**Vuex Store Utilities:**
- `waitForVuexState(page, condition, timeout)` - Wait for store update
- `getWorkflowMetadata(page)` - Get workflow state from store
- `getNodeFileAssignment(page, nodeId)` - Get node file assignment
- `waitForFilePropagation(page, fileName, timeout)` - Wait for file sharing

**Drawflow Utilities:**
- `getDrawflowNodeData(page, nodeId)` - Get node data
- `getDrawflowConnections(page)` - Get all connections
- `createDrawflowNode(page, type, position)` - Programmatically create node
- `connectDrawflowNodes(page, from, to)` - Programmatically connect nodes
- `nodeExistsOnCanvas(page, type)` - Check if node exists
- `countNodesOnCanvas(page, type)` - Count nodes by type

**Component Utilities:**
- `waitForPlotlyRender(page, containerId, timeout)` - Wait for Plotly chart
- `waitForTableLoad(page, timeout)` - Wait for vue-good-table load
- `getTableData(page)` - Extract table data
- `waitForModalAnimation(page, duration)` - Wait for animations
- `clearWorkflowStorage(page)` - Clear localStorage

### Network Fixtures (`fixtures/workflow-api.js`)

Mock API responses for deterministic testing:

**Available Fixtures:**
- `mockFilesResponse` - Folder structure
- `mockFolderFiles` - Files in folders
- `mockPluginsResponse` - Available plugins with metadata
- `mockDataTableResponse` - Paginated table data
- `mockScatterPlotData` - UMAP coordinates and clusters
- `mockH5ADColumns` - Available H5AD columns

**Setup Function:**
```javascript
await setupWorkflowRoutes(page, {
  useDefaultFixtures: true,
  customResponses: {
    files: customFileData,
    plugins: customPluginData,
  }
});
```

**Supported Endpoints:**
- `/api/routes/files/*` - File/folder listing
- `/api/routes/plugin/list` - Plugin metadata
- `/api/routes/datatable/load_data` - Table data
- `/api/routes/files/data/*` - Scatter plot data
- `/api/routes/files/columns` - H5AD columns
- `/api/routes/workflow/*` - Workflow save/update

---

## Test Implementation Status

### ✅ Priority 1: InputFile Assignment (COMPLETE)

**File:** `workflow-configuration.spec.js`

**Test Cases:**
1. ✅ Should assign file to InputFile node
   - Opens InputFile modal
   - Selects folder and file
   - Applies assignment
   - Verifies Apply button state change
   - Verifies persistence by reopening modal

2. ✅ Should change assigned file to different file
   - Assigns initial file
   - Changes to different file
   - Verifies new file persisted

3. ✅ Should handle empty folder gracefully
   - Selects empty folder
   - Verifies appropriate handling

**Key Selectors:**
- `.folder__item` - Folder list items
- `.file__item` - File list items
- `.form__button--apply` - Apply button
- `.form__info--name` - Current file name
- `.toggleFolder` - Selected folder class
- `.toggleFile` - Selected file class

**API Interactions:**
- `GET /api/routes/files/*` - Fetch folders/files
- `POST /api/routes/workflow/*` - Save file assignment

**Success Criteria:** ✅ All passing
- Folder selection highlights correctly
- File details display accurately
- Apply button state changes to "Applied"
- File assignment persists in Vuex store
- File can be changed after initial assignment

---

## Remaining Test Scenarios

### 🔜 Priority 2: Dropdown File Selection (NEXT)

**Objective:** Verify file propagation from InputFile to Algorithm node via dropdowns

**Prerequisites:**
- Priority 1 tests passing
- InputFile → Algorithm connection established

**Test Cases to Implement:**

1. **Should propagate file to Algorithm dropdown after connection**
   ```javascript
   test('Should propagate file to Algorithm dropdown', async ({ page }) => {
     // Setup: Create workflow with InputFile and Algorithm nodes
     await workflowPage.dragNodeToCanvas('Algorithm', { x: 500, y: 300 });

     // Connect InputFile to Algorithm
     await workflowPage.connectNodes('InputFile', 'Algorithm');

     // Assign file to InputFile
     await workflowPage.openNodeModal('InputFile');
     await inputFileModal.assignFile('PBMCLight1000', 'PBMCLight1000.h5ad');

     // Open Algorithm modal
     await workflowPage.openNodeModal('Algorithm');

     // Verify file appears in dropdown
     await algorithmModal.verifyFileInDropdown(0, 'PBMCLight1000.h5ad');

     // Select the file
     await algorithmModal.selectInputFile(0, 'PBMCLight1000.h5ad');

     // Verify selection
     const selected = await algorithmModal.getSelectedInputFile(0);
     expect(selected).toContain('PBMCLight1000.h5ad');
   });
   ```

2. **Should enable dropdown only when upstream node connected**
   - Verify dropdown disabled when no connection
   - Connect upstream node
   - Verify dropdown enabled
   - Verify files populated

3. **Should update dropdown when upstream file changes**
   - Assign initial file to InputFile
   - Verify file in Algorithm dropdown
   - Change InputFile assignment
   - Verify Algorithm dropdown updates

**Key Selectors:**
- `.plugin-input-dropdown` - Input file dropdowns
- `.plugin-input-dropdown option` - Dropdown options
- `.plugin-input-dropdown:disabled` - Disabled state

**API Interactions:**
- `GET /api/routes/plugin/list` - Plugin metadata (inputs/outputs)
- Vuex mutations for file propagation

**Estimated Time:** 2 hours

---

### 📋 Priority 3: Algorithm Parameter Configuration

**Objective:** Verify parameter default values and custom input handling

**Prerequisites:**
- Priority 2 tests passing

**Test Cases to Implement:**

1. **Should display default parameter values**
   - Open Algorithm modal
   - Select plugin (e.g., TENET)
   - Verify default values for all parameters
   - Check numeric, boolean, and text parameters

2. **Should accept custom parameter values**
   - Set custom numeric value
   - Toggle boolean checkbox
   - Input custom text value
   - Verify values persist after blur/change

3. **Should validate parameter constraints**
   - Test min/max bounds for numeric inputs
   - Test required field validation
   - Test data type validation

4. **Should persist parameter values**
   - Set custom parameters
   - Close and reopen modal
   - Verify parameters retained

**Key Selectors:**
- `.parameter__input` - Parameter inputs
- `.parameter__label` - Parameter labels
- `input[type="checkbox"]` - Boolean parameters
- `.parameter__error` - Validation errors

**API Interactions:**
- `GET /api/routes/plugin/list` - Parameter metadata

**Estimated Time:** 3 hours

---

### 📋 Priority 4: DataTable Matrix Display

**Objective:** Verify matrix data loads and displays correctly in DataTable node

**Prerequisites:**
- Priority 1 tests passing (file assigned to InputFile)

**Test Cases to Implement:**

1. **Should display matrix data from assigned file**
   - Connect InputFile → DataTable
   - Assign H5AD file to InputFile
   - Open DataTable modal
   - Verify table renders with columns and rows
   - Verify data matches expected values

2. **Should handle pagination**
   - Load large dataset
   - Verify pagination controls appear
   - Click next page
   - Verify new data loads

3. **Should handle sorting**
   - Click column header
   - Verify data sorts ascending
   - Click again
   - Verify data sorts descending

4. **Should show loading state**
   - Trigger data load
   - Verify loading overlay appears
   - Wait for load completion
   - Verify loading overlay disappears

**Key Selectors:**
- `.vgt-table` - Table container
- `.vgt-table thead th` - Table headers
- `.vgt-table tbody tr` - Table rows
- `.vgt-loading` - Loading overlay
- `.vgt-pagination` - Pagination controls

**API Interactions:**
- `POST /api/routes/datatable/load_data` - Table data with pagination

**Helper Functions:**
- `waitForTableLoad(page, timeout)`
- `getTableData(page)` - Extract table data

**Estimated Time:** 4 hours

---

### 📋 Priority 5: ScatterPlot UMAP Visualization

**Objective:** Verify UMAP visualization renders correctly with proper controls

**Prerequisites:**
- Priority 1 tests passing

**Test Cases to Implement:**

1. **Should render UMAP scatter plot**
   - Connect InputFile → ScatterPlot
   - Assign H5AD file
   - Open ScatterPlot modal
   - Verify Plotly chart renders
   - Verify scatter layer appears

2. **Should display axis selectors**
   - Verify X-axis dropdown populated
   - Verify Y-axis dropdown populated
   - Change axis selection
   - Verify plot updates

3. **Should display grouping selector**
   - Verify grouping dropdown populated (e.g., cluster, cell_type)
   - Change grouping
   - Verify colors update

4. **Should handle lasso selection**
   - Verify lasso tool available
   - Perform lasso selection (if possible in headless)
   - Verify selection persisted

5. **Should show loading state**
   - Trigger data load
   - Verify loading indicator
   - Wait for chart render
   - Verify chart interactive

**Key Selectors:**
- `#plotly__scatter` - Plotly container
- `#plotly__scatter .scatterlayer` - Scatter plot layer
- `#plotly__scatter svg.main-svg` - Main SVG element
- `.options__item select` - Axis/grouping selectors

**API Interactions:**
- `GET /api/routes/files/data/<file>` - Scatter plot data
- `POST /api/routes/workflow/node/data` - Save lasso selections

**Helper Functions:**
- `waitForPlotlyRender(page, containerId, timeout)`

**Estimated Time:** 5 hours

---

## Running the Tests

### Run All Workflow Configuration Tests
```bash
cd frontend
npm run test:e2e -- workflow-configuration.spec.js
```

### Run Specific Priority Tests
```bash
# Priority 1 only
npm run test:e2e -- workflow-configuration.spec.js -g "Priority 1"

# Priority 2 only (when implemented)
npm run test:e2e -- workflow-configuration.spec.js -g "Priority 2"
```

### Run with UI Mode (for debugging)
```bash
npm run test:e2e -- workflow-configuration.spec.js --ui
```

### Run with Headed Browser
```bash
npm run test:e2e -- workflow-configuration.spec.js --headed
```

---

## Test Maintenance

### When to Update Tests

1. **Frontend Component Changes**
   - Update selectors in page objects
   - Verify interaction methods still work
   - Update assertions if UI behavior changes

2. **API Changes**
   - Update fixtures with new response structure
   - Modify route handlers in `workflow-api.js`
   - Update assertions for new data fields

3. **New Features**
   - Add new test cases for new functionality
   - Extend page objects with new methods
   - Create new fixtures as needed

### Selector Stability Strategy

**Best Practices:**
- Use data attributes when available (`data-node`, `data-testid`)
- Prefer class-based selectors for stable elements
- Avoid position-based selectors (`:nth-child`)
- Use text-based selectors for user-facing labels
- Combine selectors for specificity (`.parent > .child`)

**Current Selector Patterns:**
```javascript
// Good: Stable, semantic
.locator('.file__item', { has: page.locator(`.folder__name:has-text("${name}")`) })

// Good: Data attribute
.locator(`[data-node="${nodeType}"]`)

// Avoid: Position-dependent
.locator('li:nth-child(3)')

// Avoid: Generic
.locator('button')
```

---

## Next Steps

### Immediate (This Sprint)
1. ✅ Complete Priority 1 implementation
2. 🔜 Implement Priority 2 (Dropdown file selection)
3. 🔜 Implement Priority 3 (Parameter configuration)

### Next Sprint
1. 📋 Implement Priority 4 (DataTable)
2. 📋 Implement Priority 5 (ScatterPlot)
3. 📋 Add cross-browser testing (Firefox, WebKit)

### Future Enhancements
1. Add visual regression testing for charts
2. Add performance benchmarks for data loading
3. Add accessibility testing for modals
4. Create reusable test fixtures for common scenarios
5. Add API contract testing with schema validation

---

## Risk Assessment & Mitigation

### Known Risks

1. **Drawflow Drag Interactions (Medium Risk)**
   - **Issue:** Drag-and-drop can be flaky in headless mode
   - **Mitigation:** Use programmatic helpers when possible, increase timeout for drag operations
   - **Fallback:** Use `createDrawflowNode()` and `connectDrawflowNodes()` helpers

2. **Plotly Async Rendering (Medium Risk)**
   - **Issue:** Chart rendering timing varies
   - **Mitigation:** Use `waitForPlotlyRender()` helper with adequate timeout
   - **Verification:** Check for SVG elements and scatter layer

3. **Vuex State Synchronization (Low Risk)**
   - **Issue:** State updates may lag behind UI
   - **Mitigation:** Use `waitForVuexState()` helper for critical assertions
   - **Verification:** Query store directly via `page.evaluate()`

4. **Backend API Dependencies (High Risk)**
   - **Issue:** Tests depend on backend availability
   - **Mitigation:** Use network fixtures for deterministic testing
   - **Toggle:** Enable/disable fixtures via `setupWorkflowRoutes()` options

5. **localStorage Persistence (Low Risk)**
   - **Issue:** State leaks between tests
   - **Mitigation:** Clear storage in `beforeEach` using `clearWorkflowStorage()`
   - **Verification:** Check for clean state before each test

---

## Performance Benchmarks

### Target Metrics
- **Test Execution Time:** <30s per test case
- **Page Load Time:** <3s for workflow page
- **Modal Open Time:** <1s for any modal
- **File List Load:** <2s for folder with <100 files
- **Table Load:** <3s for 1000 rows
- **Chart Render:** <2s for 10000 points

### Actual Metrics (to be measured)
- Priority 1 tests: TBD
- Priority 2 tests: TBD
- Priority 3 tests: TBD
- Priority 4 tests: TBD
- Priority 5 tests: TBD

---

## References

### Documentation
- [Playwright Documentation](https://playwright.dev/)
- [Vue.js Testing Guide](https://vuejs.org/guide/scaling-up/testing.html)
- [Page Object Pattern](https://playwright.dev/docs/pom)

### Related Files
- `frontend/src/views/WorkFlowPage.vue` - Main workflow editor
- `frontend/src/components/modals/InputFile.vue` - InputFile modal
- `frontend/src/components/modals/Algorithm.vue` - Algorithm modal
- `frontend/src/store/workflow/` - Vuex workflow module
- `frontend/tests/e2e/authenticated-workflow.spec.js` - Project initialization tests

### Codex Analysis Report
- See Codex QA analysis output for detailed component review
- Priority matrix and rationale documented
- Selector identification and validation strategies
- Risk assessment and mitigation recommendations

---

## Changelog

### 2024-11-01
- ✅ Initial test plan created
- ✅ Priority 1 tests implemented
- ✅ Page objects created (WorkflowPage, InputFileModal, AlgorithmModal)
- ✅ Helper utilities implemented
- ✅ Network fixtures established
- 📋 Remaining priorities planned (2-5)

### Future Updates
- Document Priority 2 implementation
- Document Priority 3 implementation
- Update performance benchmarks
- Add cross-browser test results
