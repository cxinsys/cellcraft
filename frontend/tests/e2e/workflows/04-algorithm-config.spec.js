// frontend/tests/e2e/workflows/04-algorithm-config.spec.js
import { test, expect } from '../fixtures/auth.js';
import { testWorkflow } from './support/workflow-constants.js';
import { setupPageObjects, setupTestFiles, cleanupTestFiles, createWorkflowWithUniqueTitle, cleanupWorkflows } from './support/workflow-setup.js';
import { getNodeFileAssignment, getWorkflowMetadata } from '../utils/workflow.js';

/**
 * Test Suite: Algorithm Node - Parameter Configuration
 *
 * This test suite verifies the Algorithm node parameter configuration:
 * - File propagation from InputFile to Algorithm node
 * - Plugin logo and name display
 * - Input file dropdown verification
 * - Parameter types (numeric, string, boolean, dropdown) modification
 * - Vuex store synchronization for parameter changes
 *
 * Success criteria:
 * - File assignment propagates to Algorithm node
 * - Plugin logo displays correct plugin name
 * - Input file dropdown shows assigned file
 * - All parameter types can be modified
 * - Parameter changes are reflected in Vuex selectedPluginRules
 */
test.describe('Algorithm Node - Parameter Configuration', () => {
  test.describe.configure({ mode: 'serial' });

  let pageObjects;
  const uploadedFiles = [];
  const createdWorkflows = [];
  let currentTestFileName = null;
  let currentWorkflowTitle = null;

  test.beforeEach(async ({ page }) => {
    pageObjects = setupPageObjects(page);
    const { uploadedFileName } = await setupTestFiles(pageObjects.filesPage, testWorkflow, uploadedFiles);
    currentTestFileName = uploadedFileName;
    await pageObjects.projectsPage.goto();
    await pageObjects.projectsPage.verifyPageLoaded();
  });

  test.afterEach(async ({ page }) => {
    await cleanupTestFiles(pageObjects.filesPage, uploadedFiles);
    await cleanupWorkflows(pageObjects.projectsPage, createdWorkflows);
  });

  /**
   * Test: Should configure TENET algorithm parameters
   *
   * Verifies Algorithm node parameter configuration:
   * - File assignment to InputFile node
   * - File propagation to Algorithm node
   * - Plugin logo verification
   * - Input file dropdown verification
   * - Parameter modification (numeric, string, boolean, dropdown)
   * - Vuex store update verification
   */
  test('Should configure TENET algorithm parameters', async ({ page }) => {
    test.setTimeout(60000);

    // Create workflow with unique title
    currentWorkflowTitle = await createWorkflowWithUniqueTitle(
      pageObjects.projectsPage,
      pageObjects.workflowPage,
      testWorkflow,
      createdWorkflows
    );

    await page.waitForSelector('.drawflow-node', { timeout: 10000 });

    await pageObjects.workflowPage.openNodeModal(testWorkflow.inputNodeName);
    await pageObjects.inputFileModal.assignFile(testWorkflow.folder, currentTestFileName);

    await pageObjects.workflowPage.closeTab(testWorkflow.inputNodeTabName);
    await page.waitForTimeout(300);

    const algorithmLocator = await pageObjects.workflowPage.findNodeByType('Algorithm');
    const algorithmNodeIdAttr = await algorithmLocator.getAttribute('id');
    const algorithmNodeId = algorithmNodeIdAttr?.replace('node-', '') ?? '12';

    await expect
      .poll(async () => {
        const assignment = await getNodeFileAssignment(page, algorithmNodeId);
        if (!assignment) return null;

        if (typeof assignment === 'string') {
          return assignment.includes(currentTestFileName) ? assignment : null;
        }

        if (Array.isArray(assignment)) {
          return assignment.find((value) =>
            typeof value === 'string' && value.includes(currentTestFileName)
          ) ?? null;
        }

        if (typeof assignment === 'object') {
          return (
            Object.values(assignment).find(
              (value) => typeof value === 'string' && value.includes(currentTestFileName)
            ) ?? null
          );
        }

        return null;
      }, {
        message: `Expected Algorithm node (${algorithmNodeId}) to receive file ${currentTestFileName}`,
        timeout: 10000,
      })
      .not.toBeNull();

    await pageObjects.workflowPage.openNodeModal('Algorithm');
    await pageObjects.algorithmModal.verifyModalOpen();

    await expect
      .poll(async () => {
        const logoText = await pageObjects.algorithmModal.getPluginLogoText();
        return logoText;
      }, {
        message: `Waiting for algorithm logo to display plugin name ${testWorkflow.name}`,
        timeout: 15000,
      })
      .toContain(testWorkflow.name);

    const pluginLogoText = await pageObjects.algorithmModal.getPluginLogoText();
    console.log('🔖 Algorithm logo text:', pluginLogoText);

    await pageObjects.algorithmModal.verifyFileInDropdown(0, currentTestFileName);
    const selectedInputFile = await pageObjects.algorithmModal.getSelectedInputFile(0);
    expect(selectedInputFile).toContain(currentTestFileName);

    const metadataBefore = await getWorkflowMetadata(page);
    console.log(
      '📦 Vuex drawflow data before Algorithm parameter changes:',
      JSON.stringify(metadataBefore?.workflowInfo?.drawflow?.Home?.data ?? null, null, 2)
    );

    const algorithmNodeDataBefore = metadataBefore?.workflowInfo?.drawflow?.Home?.data?.[algorithmNodeId];
    const pluginRulesBefore = algorithmNodeDataBefore?.data?.selectedPluginRules ?? [];
    console.log('🧮 Algorithm selectedPluginRules (before):', JSON.stringify(pluginRulesBefore, null, 2));

    const flattenedParamsBefore = pluginRulesBefore.flatMap((rule) => rule.parameters ?? []);

    const numericParam = flattenedParamsBefore.find((param) =>
      ['int', 'float', 'number'].includes(param?.type)
    );
    let numericParamName = null;
    let numericNewValue = null;
    if (numericParam) {
      numericParamName = numericParam.name;
      const numericInitialValue = Number(numericParam.defaultValue ?? 0);
      numericNewValue = String(numericInitialValue + 1);
      await pageObjects.algorithmModal.setParameterValueByName(numericParamName, numericNewValue);
      const numericUiValue = await pageObjects.algorithmModal.getParameterValueByName(numericParamName);
      expect(numericUiValue).toBe(numericNewValue);
    }

    const stringParam = flattenedParamsBefore.find((param) => param?.type === 'string');
    let stringParamName = null;
    let stringNewValue = null;
    if (!numericParam && stringParam) {
      stringParamName = stringParam.name;
      const stringInitialValue = stringParam.defaultValue ?? '';
      stringNewValue = stringInitialValue === '' ? 'test-value' : `${stringInitialValue}-updated`;
      await pageObjects.algorithmModal.setParameterValueByName(stringParamName, stringNewValue);
      const stringUiValue = await pageObjects.algorithmModal.getParameterValueByName(stringParamName);
      expect(stringUiValue).toBe(stringNewValue);
    }

    const booleanParam = flattenedParamsBefore.find((param) => param?.type === 'boolean');
    let booleanParamName = null;
    let booleanNewValue = null;
    if (booleanParam) {
      booleanParamName = booleanParam.name;
      const booleanInitialValue = booleanParam.defaultValue === true || booleanParam.defaultValue === 'true';
      booleanNewValue = !booleanInitialValue;
      await pageObjects.algorithmModal.setParameterValueByName(booleanParamName, booleanNewValue);
      const booleanUiValue = await pageObjects.algorithmModal.getParameterValueByName(booleanParamName);
      expect(booleanUiValue).toBe(booleanNewValue);
    }

    const dropdownParam = flattenedParamsBefore.find(
      (param) => param?.type === 'h5adParameter' && param?.name !== 'clusters'
    );
    let dropdownParamName = null;
    let dropdownNewValue = null;
    if (dropdownParam) {
      dropdownParamName = dropdownParam.name;
      const dropdownOptions = await pageObjects.algorithmModal.getParameterDropdownOptions(dropdownParamName);
      console.log(`📝 Dropdown options for ${dropdownParamName}:`, dropdownOptions);
      const preferredOption = dropdownOptions.find((opt) => opt && opt !== dropdownParam.defaultValue) ?? dropdownOptions[0];

      if (preferredOption && preferredOption !== dropdownParam.defaultValue) {
        const { next } = await pageObjects.algorithmModal.selectParameterDropdownOption(
          dropdownParamName,
          preferredOption
        );
        dropdownNewValue = next;
        const dropdownUiValue = await pageObjects.algorithmModal.getParameterValueByName(dropdownParamName);
        expect(dropdownUiValue).toBe(dropdownNewValue);
      }
    }

    await page.waitForTimeout(500);

    const metadataAfter = await getWorkflowMetadata(page);
    console.log(
      '📦 Vuex drawflow data after Algorithm polling:',
      JSON.stringify(metadataAfter?.workflowInfo?.drawflow?.Home?.data ?? null, null, 2)
    );

    const algorithmNodeDataAfter = metadataAfter?.workflowInfo?.drawflow?.Home?.data?.[algorithmNodeId];
    const pluginRulesAfter = algorithmNodeDataAfter?.data?.selectedPluginRules ?? [];
    const flattenedParamsAfter = pluginRulesAfter.flatMap((rule) => rule.parameters ?? []);

    if (numericParamName) {
      const updatedNumericParam = flattenedParamsAfter.find((param) => param?.name === numericParamName);
      expect(updatedNumericParam?.defaultValue).toBe(numericNewValue);
    }

    if (stringParamName) {
      const updatedStringParam = flattenedParamsAfter.find((param) => param?.name === stringParamName);
      expect(updatedStringParam?.defaultValue).toBe(stringNewValue);
    }

    if (booleanParamName) {
      const updatedBooleanParam = flattenedParamsAfter.find((param) => param?.name === booleanParamName);
      const booleanAfterValue = updatedBooleanParam?.defaultValue;
      const booleanAfterNormalized = booleanAfterValue === true || booleanAfterValue === 'true';
      expect(booleanAfterNormalized).toBe(booleanNewValue);
    }

    if (dropdownParamName && dropdownNewValue) {
      const updatedDropdownParam = flattenedParamsAfter.find((param) => param?.name === dropdownParamName);
      expect(updatedDropdownParam?.defaultValue).toBe(dropdownNewValue);
    }
  });
});
