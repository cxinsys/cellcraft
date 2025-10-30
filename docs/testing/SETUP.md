# CellCraft Backend Testing - Environment Setup Guide

> Step-by-step guide to configure and run the backend testing environment

**Last Updated**: 2025-01-27
**Difficulty**: Intermediate
**Estimated Setup Time**: 20-30 minutes

---

## Prerequisites

### System Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| **OS** | Linux/macOS | Ubuntu 20.04+ / macOS 12+ |
| **RAM** | 8GB | 16GB |
| **Disk Space** | 5GB | 10GB |
| **Python** | 3.9+ | 3.10+ |
| **Docker** | 20.10+ | Latest |
| **Conda/Miniconda** | Latest | Latest |

### Required Software

- ✅ **Docker** & **Docker Compose** (for PostgreSQL test database)
- ✅ **Conda** or **Miniconda** (for Python environment management)
- ✅ **Git** (for version control)

### Verify Installations

```bash
# Check Docker
docker --version
docker compose version

# Check Conda
conda --version

# Check Python (should be 3.9+)
python3 --version
```

**Expected Output**:
```
Docker version 24.0.7
Docker Compose version v2.23.0
conda 23.7.4
Python 3.10.13
```

---

## Step 1: Clone Repository and Navigate

```bash
# Clone repository (if not already cloned)
git clone https://github.com/cxinsys/cellcraft.git
cd cellcraft/backend
```

---

## Step 2: PostgreSQL Test Database Setup

### 2.1 Start PostgreSQL Test Container

```bash
# Navigate to project root (if in backend/)
cd ..

# Start test-db container
docker compose -f docker-compose.dev.yml up test-db -d
```

**What this does**:
- Creates isolated PostgreSQL 15 container
- Exposes port **5433** (not 5432, to avoid conflicts)
- Creates database: `cellcraft_test`
- Sets credentials: `cellcraft` / `cellcraft`

### 2.2 Verify Database is Running

```bash
# Check container status
docker ps | grep test-db

# Test connection (optional)
docker exec -it cellcraft-test-db psql -U cellcraft -d cellcraft_test -c "SELECT version();"
```

**Expected Output**:
```
cellcraft-test-db    postgres:15    Up    0.0.0.0:5433->5432/tcp
```

### 2.3 Troubleshooting Database Issues

**Issue**: Port 5433 already in use
```bash
# Find process using port 5433
lsof -i :5433

# Stop conflicting process or change port in docker-compose.dev.yml
```

**Issue**: Container fails to start
```bash
# Check logs
docker logs cellcraft-test-db

# Remove and recreate
docker compose -f docker-compose.dev.yml down test-db
docker compose -f docker-compose.dev.yml up test-db -d
```

---

## Step 3: Conda Environment Setup

### 3.1 Create Test Environment

```bash
# Navigate to backend directory
cd backend

# Create conda environment from environment-test.yml
conda env create -f environment-test.yml
```

**What this installs**:
- Python 3.10
- FastAPI 0.78.0
- SQLAlchemy 1.4.41
- pytest 7.4.3
- pytest-asyncio 0.21.1
- pytest-cov 4.1.0
- All backend dependencies

**Expected Output**:
```
Collecting package metadata: done
Solving environment: done
Preparing transaction: done
Verifying transaction: done
Executing transaction: done
#
# To activate this environment, use
#
#     $ conda activate cellcraft-test
```

### 3.2 Activate Environment

```bash
conda activate cellcraft-test
```

**Verify activation**:
```bash
# Check Python path (should include cellcraft-test)
which python

# Check installed packages
conda list | grep pytest
```

**Expected Output**:
```
/home/user/miniconda3/envs/cellcraft-test/bin/python
pytest                    7.4.3
pytest-asyncio            0.21.1
pytest-cov                4.1.0
```

### 3.3 Troubleshooting Conda Issues

**Issue**: Conda environment creation fails
```bash
# Update conda
conda update -n base -c defaults conda

# Retry with verbose output
conda env create -f environment-test.yml -v
```

**Issue**: Wrong Python version
```bash
# Remove and recreate environment
conda env remove -n cellcraft-test
conda env create -f environment-test.yml
```

---

## Step 4: Environment Variables Configuration

### 4.1 Set Required Environment Variables

The test environment requires `TESTING=1` to activate test-specific configurations.

**Option A: Export in Terminal (Temporary)**
```bash
export TESTING=1
```

**Option B: Add to Shell Profile (Permanent)**
```bash
# For bash
echo 'export TESTING=1' >> ~/.bashrc
source ~/.bashrc

# For zsh
echo 'export TESTING=1' >> ~/.zshrc
source ~/.zshrc
```

**Option C: Use .env File (Project-specific)**
```bash
# Create .env file in backend/ directory
echo "TESTING=1" > .env

# pytest-dotenv will auto-load this file
```

### 4.2 Verify Environment Variables

```bash
# Check TESTING variable
echo $TESTING
```

**Expected Output**: `1`

### 4.3 Database Connection Variables

These are configured automatically by `app/common/config.py` for test environment:

```python
# Test database configuration (when TESTING=1)
POSTGRES_HOST=localhost
POSTGRES_PORT=5433
POSTGRES_USER=cellcraft
POSTGRES_PASSWORD=cellcraft
POSTGRES_DB=cellcraft_test
```

---

## Step 5: Run First Test

### 5.1 Execute Test Suite

```bash
# Ensure you're in backend/ directory and cellcraft-test is activated
pytest tests/unit/ -v
```

**What this does**:
- Discovers all test files in `tests/unit/`
- Runs 93 tests across 7 test modules
- Shows verbose output with test names

**Expected Output** (abbreviated):
```
============================= test session starts ==============================
platform linux -- Python 3.10.13, pytest-7.4.3
collected 93 items

tests/unit/test_auth.py::test_register_new_user PASSED                    [ 1%]
tests/unit/test_auth.py::test_register_duplicate_email_fails PASSED       [ 2%]
...
tests/unit/test_schemas.py::test_user_create_valid_data PASSED            [98%]
tests/unit/test_schemas.py::test_user_create_missing_username PASSED      [99%]
tests/unit/test_security.py::test_password_hashing_verification PASSED    [100%]

======================== 93 passed, 2 xpassed in 45.23s ========================
```

### 5.2 Run Specific Test Module

```bash
# Run only authentication tests
pytest tests/unit/test_auth.py -v

# Run only schema validation tests
pytest tests/unit/test_schemas.py -v
```

### 5.3 Run Tests with Coverage Report

```bash
# Generate HTML coverage report
pytest tests/unit/ --cov=app --cov-report=html

# View report (opens in browser)
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
```

**Expected Coverage**:
```
Name                              Stmts   Miss  Cover
-----------------------------------------------------
app/common/config.py                 45      3    93%
app/common/security.py               12      0   100%
app/database/crud/crud_user.py       35      0   100%
app/database/models.py               67     15    78%
-----------------------------------------------------
TOTAL                              3456   2697    22%
```

---

## Step 6: Verify Setup Completeness

### 6.1 Run Full Test Suite with Markers

```bash
# Run only unit tests
pytest tests/unit/ -m unit -v

# Run only authentication tests
pytest tests/unit/ -m auth -v

# Run only CRUD tests
pytest tests/unit/ -m crud -v
```

### 6.2 Check Test Statistics

```bash
# Count total tests
pytest tests/unit/ --collect-only | grep "test session starts" -A 1

# List all test files
pytest tests/unit/ --collect-only -q
```

**Expected Statistics**:
```
collected 93 items

7 test files
93 test functions
```

### 6.3 Verify Database Connection

```bash
# Run database-dependent tests
pytest tests/unit/test_crud_user.py -v
```

**If tests pass**: ✅ Database connection is working
**If tests fail**: ⚠️ Check database configuration (see Troubleshooting)

---

## Common Issues and Solutions

### Issue 1: ModuleNotFoundError

**Symptom**:
```
ModuleNotFoundError: No module named 'app'
```

**Solution**:
```bash
# Ensure you're in backend/ directory
cd /path/to/cellcraft/backend

# Verify PYTHONPATH includes current directory
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Or run pytest with pythonpath
pytest tests/unit/ --pythonpath=.
```

---

### Issue 2: Database Connection Refused

**Symptom**:
```
sqlalchemy.exc.OperationalError: could not connect to server: Connection refused
```

**Solution**:
```bash
# Check if test-db is running
docker ps | grep test-db

# If not running, start it
docker compose -f docker-compose.dev.yml up test-db -d

# Wait 10 seconds for database to initialize
sleep 10

# Retry tests
pytest tests/unit/ -v
```

---

### Issue 3: Permission Denied on Database

**Symptom**:
```
psycopg2.OperationalError: FATAL: role "cellcraft" does not exist
```

**Solution**:
```bash
# Recreate database container
docker compose -f docker-compose.dev.yml down test-db
docker volume rm cellcraft_test_db_data  # Remove old data
docker compose -f docker-compose.dev.yml up test-db -d
```

---

### Issue 4: Tests Fail with "TESTING environment variable not set"

**Symptom**:
```
AssertionError: TESTING environment variable must be set
```

**Solution**:
```bash
# Set environment variable
export TESTING=1

# Verify
echo $TESTING

# Retry tests
pytest tests/unit/ -v
```

---

### Issue 5: Alembic Migration Errors

**Symptom**:
```
alembic.util.exc.CommandError: Can't locate revision identified by 'xyz'
```

**Solution**:
```bash
# Reset test database
docker compose -f docker-compose.dev.yml down test-db
docker volume rm cellcraft_test_db_data

# Restart and run migrations
docker compose -f docker-compose.dev.yml up test-db -d
sleep 10

# Apply migrations
alembic upgrade head
```

---

### Issue 6: pytest Command Not Found

**Symptom**:
```
bash: pytest: command not found
```

**Solution**:
```bash
# Ensure conda environment is activated
conda activate cellcraft-test

# Verify pytest installation
which pytest

# If not installed, reinstall dependencies
conda env update -f environment-test.yml
```

---

## Daily Development Workflow

### Starting a Test Session

```bash
# 1. Navigate to project
cd /path/to/cellcraft/backend

# 2. Activate conda environment
conda activate cellcraft-test

# 3. Ensure test database is running
docker ps | grep test-db
# If not running: docker compose -f ../docker-compose.dev.yml up test-db -d

# 4. Set environment variable (if not permanent)
export TESTING=1

# 5. Run tests
pytest tests/unit/ -v
```

### Running Specific Tests

```bash
# Test a single function
pytest tests/unit/test_auth.py::test_register_new_user -v

# Test a class
pytest tests/unit/test_crud_user.py::TestCreateUser -v

# Test with keyword matching
pytest tests/unit/ -k "auth" -v
```

### Watching for Changes (Optional)

```bash
# Install pytest-watch
conda install -c conda-forge pytest-watch

# Auto-run tests on file changes
ptw tests/unit/ -- -v
```

---

## Advanced Configuration

### Custom Test Database Port

If port 5433 conflicts with other services:

1. Edit `docker-compose.dev.yml`:
```yaml
test-db:
  ports:
    - "5434:5432"  # Change 5433 to 5434
```

2. Update `app/common/config.py` (DevelopmentConfig):
```python
POSTGRES_PORT: int = 5434  # Change from 5433
```

3. Restart database:
```bash
docker compose -f docker-compose.dev.yml restart test-db
```

---

### Running Tests in Parallel

```bash
# Install pytest-xdist
conda install -c conda-forge pytest-xdist

# Run tests with 4 workers
pytest tests/unit/ -n 4 -v
```

**Note**: Some tests may fail due to database transaction conflicts. Use with caution.

---

## Cleanup and Reset

### Stop Test Database

```bash
# Stop without removing data
docker compose -f docker-compose.dev.yml stop test-db

# Stop and remove container (keeps data)
docker compose -f docker-compose.dev.yml down test-db

# Remove everything including data
docker compose -f docker-compose.dev.yml down test-db -v
```

### Remove Conda Environment

```bash
# Deactivate environment
conda deactivate

# Remove environment
conda env remove -n cellcraft-test
```

### Fresh Start

```bash
# Complete reset
docker compose -f docker-compose.dev.yml down -v
conda env remove -n cellcraft-test
rm -rf htmlcov/ .coverage .pytest_cache/

# Recreate from scratch (follow Steps 2-5)
```

---

## Next Steps

Once your environment is set up and tests are passing:

1. **Read Test Writing Guide**: See [GUIDE.md](./GUIDE.md) for best practices
2. **Review Current Progress**: See [PROGRESS.md](./PROGRESS.md) for Phase 2 status
3. **Understand Codebase**: See [CODEBASE-ANALYSIS.md](./CODEBASE-ANALYSIS.md) for complexity analysis
4. **Start Contributing**: Write new tests following the AAA pattern

---

## Quick Reference Commands

```bash
# Start test environment
cd cellcraft/backend
conda activate cellcraft-test
export TESTING=1
docker compose -f ../docker-compose.dev.yml up test-db -d

# Run all tests
pytest tests/unit/ -v

# Run with coverage
pytest tests/unit/ --cov=app --cov-report=html

# Run specific test file
pytest tests/unit/test_auth.py -v

# Run tests matching keyword
pytest tests/unit/ -k "user" -v

# Run with markers
pytest tests/unit/ -m unit -v

# Stop test database
docker compose -f ../docker-compose.dev.yml stop test-db
```

---

## Support and Resources

- **Main Documentation**: [README.md](./README.md)
- **Progress Tracking**: [PROGRESS.md](./PROGRESS.md)
- **Test Writing Guide**: [GUIDE.md](./GUIDE.md)
- **Codebase Analysis**: [CODEBASE-ANALYSIS.md](./CODEBASE-ANALYSIS.md)
- **GitHub Issues**: Report problems at https://github.com/cxinsys/cellcraft/issues

---

**Happy Testing! 🧪**
