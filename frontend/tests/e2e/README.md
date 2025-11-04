# CellCraft E2E Testing Guide

This directory contains End-to-End (E2E) tests for the CellCraft application using Playwright.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Running Tests](#running-tests)
- [Environment Variables](#environment-variables)
- [Test Scenarios](#test-scenarios)
- [Page Object Models](#page-object-models)
- [Utilities](#utilities)
- [Debugging](#debugging)
- [CI/CD Integration](#cicd-integration)

## Prerequisites

Before running E2E tests, ensure you have:

1. **Node.js** (v16 or later)
2. **Docker and Docker Compose** (for running the full application stack)
3. **Playwright browsers** installed

### Install Playwright Browsers

```bash
npx playwright install
```

## Quick Start

### Option 1: Test Against Running Docker Compose Stack (Recommended)

1. Start the CellCraft application:
```bash
docker compose -f docker-compose.dev.yml up
```

2. Wait for all services to be ready (frontend at http://localhost:8080)

3. Run the tests:
```bash
# Run all tests (skip webServer since app is already running)
PLAYWRIGHT_SKIP_WEBSERVER=true npm run test:e2e

# Run specific test file
PLAYWRIGHT_SKIP_WEBSERVER=true npx playwright test authenticated-workflow.spec.js

# Run with UI mode for debugging
PLAYWRIGHT_SKIP_WEBSERVER=true npx playwright test --ui
```

### Option 2: Let Playwright Start the Server

```bash
# Playwright will automatically start the dev server
npm run test:e2e
```

**Note**: This option only starts the frontend. You'll need the backend running separately for full functionality.

## Project Structure

```
tests/e2e/
├── .auth/                              # Authentication state (gitignored)
│   └── test-user.json                 # Stored login session
├── fixtures/
│   ├── auth.js                        # Authentication fixture with session management
│   ├── workflow-api.js                # Mock API responses for workflow tests
│   └── files/                         # Test data files
│       ├── test_data.h5ad             # H5AD format test file
│       ├── test_sample.csv            # CSV test file
│       ├── test_genes.txt             # TXT test file
│       └── README.md                  # Fixture files documentation
├── pages/                              # Page Object Models
│   ├── DatasetsPage.js                # Datasets page POM
│   ├── FilesPage.js                   # Files page POM
│   ├── ProjectsPage.js                # Projects page POM
│   ├── WorkflowPage.js                # Workflow editor page POM (NEW)
│   └── modals/                        # Modal POMs (NEW)
│       ├── InputFileModal.js          # InputFile node configuration
│       └── AlgorithmModal.js          # Algorithm node configuration
├── utils/                              # Utility functions
│   ├── files.js                       # File upload/download utilities
│   ├── navigation.js                  # Navigation helpers
│   └── workflow.js                    # Workflow/Drawflow helpers (NEW)
├── authenticated-workflow.spec.js      # Project initialization tests
├── workflow-configuration.spec.js      # Workflow configuration tests (NEW)
├── login.spec.js                      # Login tests
├── WORKFLOW_CONFIGURATION_TEST_PLAN.md # Detailed test plan (NEW)
├── .gitignore                         # Git ignore for test artifacts
└── README.md                          # This file
```

## Running Tests

### Run All Tests

```bash
npm run test:e2e
```

### Run Specific Browser

```bash
# Chromium only
npx playwright test --project=chromium

# Firefox
npx playwright test --project=firefox

# WebKit (Safari)
npx playwright test --project=webkit
```

### Run Specific Test File

```bash
# Project initialization tests
npx playwright test authenticated-workflow.spec.js

# Workflow configuration tests (Priority 1)
npx playwright test workflow-configuration.spec.js --project=chromium -g "Priority 1"
```

### Run in Headed Mode (See Browser)

```bash
PLAYWRIGHT_HEADLESS=false npm run test:e2e
```

### Run with UI Mode (Interactive Debugging)

```bash
npx playwright test --ui
```

### Run in Debug Mode

```bash
npx playwright test --debug
```

## Environment Variables

Configure tests using environment variables:

| Variable | Description | Default |
|----------|-------------|---------|
| `PLAYWRIGHT_BASE_URL` | Application URL | `http://localhost:8080` |
| `PLAYWRIGHT_HEADLESS` | Run in headless mode | `true` |
| `PLAYWRIGHT_USER` | Test user email | `test1234@test.com` |
| `PLAYWRIGHT_PASS` | Test user password | `test1234` |
| `PLAYWRIGHT_SKIP_WEBSERVER` | Skip auto-starting webServer | `false` |
| `CI` | CI environment flag | `false` |

### Example Usage

```bash
# Test against different URL
PLAYWRIGHT_BASE_URL=http://localhost:8081 npm run test:e2e

# Use different test credentials
PLAYWRIGHT_USER=myuser@test.com PLAYWRIGHT_PASS=mypass123 npm run test:e2e

# Run with browser visible
PLAYWRIGHT_HEADLESS=false npm run test:e2e

# Test with Docker Compose running
PLAYWRIGHT_SKIP_WEBSERVER=true npm run test:e2e
```

## Test Scenarios

### 1. Authenticated Workflow Test Suite

Located in: `authenticated-workflow.spec.js`

**Status**: ✅ Complete

#### Main Test: Full Project Initialization Flow

**Scenario Steps:**
1. **Download Tutorial Dataset**
   - Navigate to Datasets page
   - Search for `pbmc_light_1000.h5ad`
   - Download the dataset
   - Verify download completion

2. **Upload Files**
   - Navigate to Files page
   - Upload H5AD file (`test_data.h5ad`)
   - Upload CSV file (`test_sample.csv`)
   - Upload TXT file (`test_genes.txt`)
   - Verify all uploads completed

3. **Verify File List**
   - Check all uploaded files appear in the table
   - Verify file count is correct

4. **Create TENET Workflow**
   - Navigate to Projects page
   - Click "New Workflow"
   - Select TENET plugin template
   - Verify workflow creation and redirect

#### Additional Tests

- **File Operations**: Upload and delete file lifecycle
- **Dataset Search**: Search and filter functionality
- **Plugin Templates**: Display and selection of templates

---

### 2. Workflow Configuration Test Suite

Located in: `workflow-configuration.spec.js`

**Status**: 🚧 Priority 1 Complete, Priority 2-5 Planned

**Authentication**: Uses shared auth fixture (`.auth/test-user.json`) for optimal performance

#### Test Priority Matrix

| Priority | Test Scenario | Status | Dependencies |
|----------|---------------|--------|--------------|
| 1 | InputFile Assignment | ✅ Complete | None |
| 2 | Dropdown File Selection | 📋 Planned | Priority 1 |
| 3 | Algorithm Parameter Configuration | 📋 Planned | Priority 2 |
| 4 | DataTable Matrix Display | 📋 Planned | Priority 1 |
| 5 | ScatterPlot UMAP Visualization | 📋 Planned | Priority 1 |

#### Priority 1: InputFile Assignment (✅ Complete)

**Scenario Steps:**
1. Create workflow from TENET template
2. Verify InputFile node exists on canvas
3. Open InputFile configuration modal
4. Select folder (e.g., `PBMCLight1000`)
5. Select file (e.g., `PBMCLight1000.h5ad`)
6. Apply file assignment
7. Verify "Apply" button changes to "Applied"
8. Verify persistence by reopening modal

**Test Cases:**
- ✅ Should assign file to InputFile node
- ✅ Should change assigned file to different file
- ✅ Should handle empty folder gracefully

**Run Command:**
```bash
# Headless chromium (default, fastest)
npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 1"

# With visible browser (debugging)
PLAYWRIGHT_HEADLESS=false npm run test:e2e -- workflow-configuration.spec.js --project=chromium -g "Priority 1"

# With UI mode (recommended for development)
npm run test:e2e -- workflow-configuration.spec.js --project=chromium --ui
```

#### Priority 2-5: Upcoming Tests

See `WORKFLOW_CONFIGURATION_TEST_PLAN.md` for detailed implementation plan

## Page Object Models

Page Object Models (POMs) encapsulate page interactions and provide a clean API for tests.

### DatasetsPage

```javascript
import { DatasetsPage } from './pages/DatasetsPage.js';

const datasetsPage = new DatasetsPage(page);
await datasetsPage.goto();
await datasetsPage.searchDataset('pbmc');
const download = await datasetsPage.downloadDataset('pbmc_light_1000.h5ad');
```

**Methods:**
- `goto()`: Navigate to page
- `searchDataset(searchTerm)`: Filter datasets
- `downloadDataset(title)`: Download a dataset
- `verifyDownload(download, filename)`: Verify download
- `getVisibleDatasets()`: Get list of datasets
- `verifyPageLoaded()`: Verify page is displayed

### FilesPage

```javascript
import { FilesPage } from './pages/FilesPage.js';

const filesPage = new FilesPage(page);
await filesPage.goto();
await filesPage.uploadFile('test_data.h5ad');
await filesPage.verifyFileExists('test_data.h5ad');
await filesPage.deleteFile('test_data.h5ad');
```

**Methods:**
- `goto()`: Navigate to page
- `uploadFile(fileName)`: Upload a file
- `verifyFileExists(fileName)`: Check file presence
- `deleteFile(fileName)`: Delete a file
- `getFileList()`: Get all files with metadata
- `getFileCount()`: Get number of files
- `waitForUploadComplete()`: Wait for upload to finish

### ProjectsPage

```javascript
import { ProjectsPage } from './pages/ProjectsPage.js';

const projectsPage = new ProjectsPage(page);
await projectsPage.goto();
await projectsPage.clickNewWorkflow();
await projectsPage.selectPluginTemplate('TENET');
```

**Methods:**
- `goto()`: Navigate to page
- `clickNewWorkflow()`: Open plugin selection modal
- `selectPluginTemplate(name)`: Select a template
- `selectDefaultTemplate()`: Create blank workflow
- `getAvailablePlugins()`: Get plugin list
- `verifyPluginAvailable(name)`: Check plugin exists
- `openWorkflow(title)`: Open existing workflow
- `deleteWorkflow(title)`: Delete a workflow

## Utilities

### File Utilities (`utils/files.js`)

```javascript
import { uploadFixture, verifyFileInList, deleteFile } from './utils/files.js';

// Upload a file
await uploadFixture(page, 'test_data.h5ad');

// Verify file exists
await verifyFileInList(page, 'test_data.h5ad', { shouldExist: true });

// Delete a file
await deleteFile(page, 'test_data.h5ad');
```

### Navigation Utilities (`utils/navigation.js`)

```javascript
import { goToDatasets, goToFiles, goToProjects, verifyAuthenticated } from './utils/navigation.js';

// Navigate to pages
await goToDatasets(page);
await goToFiles(page);
await goToProjects(page);

// Verify authentication
await verifyAuthenticated(page);
```

## Authentication

Tests use a shared authentication state to avoid logging in for every test.

**How it works:**
1. The `auth.js` fixture logs in once per worker
2. Session state is saved to `.auth/test-user.json`
3. Subsequent tests reuse the stored session
4. Session file is gitignored for security

**Custom credentials:**
```bash
PLAYWRIGHT_USER=custom@test.com PLAYWRIGHT_PASS=custompass npm run test:e2e
```

## Debugging

### View Test Reports

After running tests, view the HTML report:

```bash
npx playwright show-report
```

### Debug Failed Tests

```bash
# Run in debug mode
npx playwright test --debug

# Run with traces
npx playwright test --trace on

# View traces after failure
npx playwright show-trace trace.zip
```

### Screenshots and Videos

- **Screenshots**: Captured on failure (saved to `test-results/`)
- **Videos**: Recorded on failure (saved to `test-results/`)
- **Traces**: Captured on first retry (saved to `test-results/`)

### Common Issues

#### 1. "Target closed" or connection errors

**Solution**: Ensure the application is running before tests:
```bash
docker compose -f docker-compose.dev.yml up
```

#### 2. Authentication failures

**Solution**: Verify test credentials are correct:
- Check `.env` file or environment variables
- Default: `test1234@test.com` / `test1234`

#### 3. File upload failures

**Solution**: Ensure fixture files exist:
```bash
ls -la tests/e2e/fixtures/files/
```

See `tests/e2e/fixtures/files/README.md` for setup instructions.

#### 4. Port conflicts

**Solution**: Check if another service is using port 8080:
```bash
lsof -i :8080
```

## CI/CD Integration

### GitHub Actions Example

```yaml
name: E2E Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Setup Node.js
        uses: actions/setup-node@v3
        with:
          node-version: '18'

      - name: Install dependencies
        run: |
          cd frontend
          npm ci

      - name: Install Playwright Browsers
        run: npx playwright install --with-deps

      - name: Start Docker Compose
        run: docker compose -f docker-compose.dev.yml up -d

      - name: Wait for services
        run: |
          npx wait-on http://localhost:8080
          npx wait-on http://localhost:8000

      - name: Run E2E tests
        run: |
          cd frontend
          PLAYWRIGHT_SKIP_WEBSERVER=true npm run test:e2e

      - name: Upload test results
        if: always()
        uses: actions/upload-artifact@v3
        with:
          name: playwright-report
          path: frontend/playwright-report/
          retention-days: 30
```

### Docker Container Testing

Run tests inside Docker:

```bash
docker run --rm \
  --network=host \
  -v $(pwd):/app \
  -w /app/frontend \
  mcr.microsoft.com/playwright:v1.40.0-focal \
  /bin/bash -c "npm ci && PLAYWRIGHT_SKIP_WEBSERVER=true npm run test:e2e"
```

## Best Practices

1. **Use Page Object Models**: Encapsulate page interactions
2. **Use test.step()**: Make test reports more readable
3. **Avoid hardcoded waits**: Use Playwright's auto-waiting
4. **Reuse authentication**: Use the auth fixture
5. **Clean up test data**: Delete created resources after tests
6. **Use descriptive test names**: Make failures easy to understand
7. **Isolate tests**: Each test should be independent

## Additional Resources

- [Playwright Documentation](https://playwright.dev)
- [Playwright Best Practices](https://playwright.dev/docs/best-practices)
- [Page Object Model Pattern](https://playwright.dev/docs/pom)
- [Test Fixtures](https://playwright.dev/docs/test-fixtures)

## Support

For issues or questions:
1. Check this README
2. Review test code and comments
3. Check Playwright documentation
4. Contact the development team
