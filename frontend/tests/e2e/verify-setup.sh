#!/bin/bash
# E2E Test Setup Verification Script

set -e

echo "🔍 Verifying E2E Test Setup..."
echo

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if we're in the correct directory
if [ ! -f "playwright.config.js" ]; then
    echo -e "${RED}❌ Error: Please run this script from the frontend directory${NC}"
    exit 1
fi

# 1. Check directory structure
echo "📁 Checking directory structure..."
REQUIRED_DIRS=(
    "tests/e2e/fixtures"
    "tests/e2e/fixtures/files"
    "tests/e2e/pages"
    "tests/e2e/utils"
)

for dir in "${REQUIRED_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        echo -e "${GREEN}✓${NC} $dir exists"
    else
        echo -e "${RED}✗${NC} $dir missing"
        exit 1
    fi
done
echo

# 2. Check required files
echo "📄 Checking required files..."
REQUIRED_FILES=(
    "tests/e2e/fixtures/auth.js"
    "tests/e2e/fixtures/files/test_data.h5ad"
    "tests/e2e/fixtures/files/test_sample.csv"
    "tests/e2e/fixtures/files/test_genes.txt"
    "tests/e2e/pages/DatasetsPage.js"
    "tests/e2e/pages/FilesPage.js"
    "tests/e2e/pages/ProjectsPage.js"
    "tests/e2e/utils/files.js"
    "tests/e2e/utils/navigation.js"
    "tests/e2e/authenticated-workflow.spec.js"
    "tests/e2e/README.md"
)

for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo -e "${GREEN}✓${NC} $file exists"
    else
        echo -e "${RED}✗${NC} $file missing"
        exit 1
    fi
done
echo

# 3. Check Playwright installation
echo "🎭 Checking Playwright installation..."
if npm list @playwright/test >/dev/null 2>&1; then
    echo -e "${GREEN}✓${NC} Playwright is installed"
else
    echo -e "${YELLOW}⚠${NC} Playwright not found in node_modules"
    echo "  Run: npm install"
    exit 1
fi
echo

# 4. Check if Playwright browsers are installed
echo "🌐 Checking Playwright browsers..."
if npx playwright --version >/dev/null 2>&1; then
    echo -e "${GREEN}✓${NC} Playwright CLI is available"
    PLAYWRIGHT_VERSION=$(npx playwright --version)
    echo "  Version: $PLAYWRIGHT_VERSION"
else
    echo -e "${YELLOW}⚠${NC} Playwright browsers may not be installed"
    echo "  Run: npx playwright install"
fi
echo

# 5. Check if Docker Compose is running
echo "🐳 Checking Docker Compose services..."
if docker compose -f ../docker-compose.dev.yml ps | grep -q "Up"; then
    echo -e "${GREEN}✓${NC} Docker Compose services are running"
else
    echo -e "${YELLOW}⚠${NC} Docker Compose services not detected"
    echo "  Start with: docker compose -f docker-compose.dev.yml up -d"
fi
echo

# 6. Check if application is accessible
echo "🌍 Checking application accessibility..."
if curl -s -o /dev/null -w "%{http_code}" http://localhost:8080 | grep -q "200\|301\|302"; then
    echo -e "${GREEN}✓${NC} Application is accessible at http://localhost:8080"
else
    echo -e "${YELLOW}⚠${NC} Application not accessible at http://localhost:8080"
    echo "  Make sure Docker Compose is running"
fi
echo

# 7. Test credentials information
echo "🔑 Test Credentials:"
echo "  Email: ${PLAYWRIGHT_USER:-test1234@test.com}"
echo "  Password: ${PLAYWRIGHT_PASS:-test1234}"
echo

# Summary
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}✅ E2E Test Setup Verification Complete!${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo
echo "📚 Next Steps:"
echo "  1. Ensure Docker Compose is running:"
echo "     docker compose -f ../docker-compose.dev.yml up -d"
echo
echo "  2. Install Playwright browsers (if not done):"
echo "     npx playwright install"
echo
echo "  3. Run the E2E tests:"
echo "     PLAYWRIGHT_SKIP_WEBSERVER=true npm run test:e2e"
echo
echo "  4. Run with UI mode for debugging:"
echo "     PLAYWRIGHT_SKIP_WEBSERVER=true npx playwright test --ui"
echo
echo "  5. Run specific test file:"
echo "     PLAYWRIGHT_SKIP_WEBSERVER=true npx playwright test authenticated-workflow.spec.js"
echo
echo "📖 For more information, see tests/e2e/README.md"
