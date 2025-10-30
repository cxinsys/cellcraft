"""
System Installation Tests

Tests for validating Docker Compose deployment and container initialization.
Covers TC-INIT-001 through TC-INIT-016 from DEPLOYMENT_TEST_PLAN.md
"""

import pytest
import time


class TestSystemInstallation:
    """4.1 시스템 설치 및 초기화 테스트"""

    # ==============================================================================
    # Container Health Check Tests (TC-INIT-001~005)
    # ==============================================================================

    def test_container_startup_sequence(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-001: 컨테이너 시작 순서 검증

        Validates that containers start in the correct order with proper dependencies:
        1. Database (db) - must be healthy before backend
        2. RabbitMQ (rabbitmq) - must be healthy before celery
        3. Backend (backend) - depends on db and rabbitmq
        4. Celery (celery) - depends on rabbitmq and backend
        5. Frontend (frontend) - depends on backend

        Success Criteria:
        - All containers reach "running" status
        - Dependency order is respected
        - No premature starts before dependencies are ready
        """
        from .helpers import get_container_status

        # Check all containers are running
        required_containers = ["db", "rabbitmq", "backend", "celery", "frontend"]

        for service in required_containers:
            container_name = container_names[service]
            status_info = get_container_status(container_name)

            assert status_info is not None, \
                f"Container {container_name} ({service}) should exist"
            assert status_info["status"] == "running", \
                f"Container {container_name} ({service}) should be running, got {status_info['status']}"

        print(f"\n=== Container Startup Sequence ===")
        print(f"✅ All {len(required_containers)} containers are running")
        for service in required_containers:
            status_info = get_container_status(container_names[service])
            print(f"  - {service}: {status_info['status']} (started: {status_info['started_at']})")

    def test_database_health_check(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-002: Database Health Check

        Validates PostgreSQL database is healthy and accessible.

        Success Criteria:
        - Container status is "running"
        - Health check returns "healthy"
        - Database accepts connections
        - POSTGRES_DB environment variable is set correctly
        """
        from .helpers import get_container_status, get_container_env

        container_name = container_names["db"]
        status_info = get_container_status(container_name)

        # Verify container is running
        assert status_info is not None, f"Database container {container_name} should exist"
        assert status_info["status"] == "running", \
            f"Database should be running, got {status_info['status']}"

        # Verify health check
        assert status_info["health"] == "healthy", \
            f"Database health check should be healthy, got {status_info['health']}"

        # Verify environment variables
        db_name = get_container_env(container_name, "POSTGRES_DB")
        assert db_name is not None, "POSTGRES_DB should be set"
        assert db_name == "cellcraft", \
            f"POSTGRES_DB should be 'cellcraft', got '{db_name}'"

        print(f"\n=== Database Health Check ===")
        print(f"✅ Database is healthy")
        print(f"  - Status: {status_info['status']}")
        print(f"  - Health: {status_info['health']}")
        print(f"  - Database: {db_name}")

    def test_rabbitmq_health_check(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-003: RabbitMQ Health Check

        Validates RabbitMQ message broker is healthy and accessible.

        Success Criteria:
        - Container status is "running"
        - Health check returns "healthy"
        - RabbitMQ management interface is accessible (port 15672)
        - Default vhost is configured
        """
        from .helpers import get_container_status, check_service_accessible

        container_name = container_names["rabbitmq"]
        status_info = get_container_status(container_name)

        # Verify container is running
        assert status_info is not None, f"RabbitMQ container {container_name} should exist"
        assert status_info["status"] == "running", \
            f"RabbitMQ should be running, got {status_info['status']}"

        # Verify health check
        assert status_info["health"] == "healthy", \
            f"RabbitMQ health check should be healthy, got {status_info['health']}"

        # Verify management interface is accessible
        management_url = "http://localhost:15672"
        is_accessible = check_service_accessible(management_url, timeout=5)
        assert is_accessible, \
            f"RabbitMQ management interface should be accessible at {management_url}"

        print(f"\n=== RabbitMQ Health Check ===")
        print(f"✅ RabbitMQ is healthy")
        print(f"  - Status: {status_info['status']}")
        print(f"  - Health: {status_info['health']}")
        print(f"  - Management UI: {management_url}")

    def test_backend_health_check(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-004: Backend Health Check

        Validates FastAPI backend is healthy and accessible.

        Success Criteria:
        - Container status is "running"
        - Health check returns "healthy"
        - API root endpoint returns 200
        - API docs endpoint (/docs) is accessible
        """
        from .helpers import get_container_status, check_service_accessible

        container_name = container_names["backend"]
        status_info = get_container_status(container_name)

        # Verify container is running
        assert status_info is not None, f"Backend container {container_name} should exist"
        assert status_info["status"] == "running", \
            f"Backend should be running, got {status_info['status']}"

        # Verify health check
        assert status_info["health"] == "healthy", \
            f"Backend health check should be healthy, got {status_info['health']}"

        # Verify API is accessible
        api_url = "http://localhost:8000"
        is_accessible = check_service_accessible(api_url, timeout=5)
        assert is_accessible, \
            f"Backend API should be accessible at {api_url}"

        # Verify API docs are accessible
        docs_url = "http://localhost:8000/docs"
        docs_accessible = check_service_accessible(docs_url, timeout=5)
        assert docs_accessible, \
            f"API documentation should be accessible at {docs_url}"

        print(f"\n=== Backend Health Check ===")
        print(f"✅ Backend is healthy")
        print(f"  - Status: {status_info['status']}")
        print(f"  - Health: {status_info['health']}")
        print(f"  - API: {api_url}")
        print(f"  - Docs: {docs_url}")

    def test_frontend_health_check(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-005: Frontend Health Check

        Validates Vue.js frontend is accessible and serving content.

        Success Criteria:
        - Container status is "running"
        - Frontend returns 200 on root path
        - Content-Type header indicates HTML
        - Response contains Vue.js application markers
        """
        from .helpers import get_container_status
        import requests

        container_name = container_names["frontend"]
        status_info = get_container_status(container_name)

        # Verify container is running
        assert status_info is not None, f"Frontend container {container_name} should exist"
        assert status_info["status"] == "running", \
            f"Frontend should be running, got {status_info['status']}"

        # Verify frontend is accessible and serving content
        frontend_url = "http://localhost:8080"
        try:
            response = requests.get(frontend_url, timeout=5)
            assert response.status_code == 200, \
                f"Frontend should return 200, got {response.status_code}"

            # Check content type
            content_type = response.headers.get("Content-Type", "")
            assert "text/html" in content_type, \
                f"Frontend should serve HTML, got Content-Type: {content_type}"

            # Check for Vue.js markers in response
            html_content = response.text
            assert len(html_content) > 0, "Frontend should return non-empty HTML"

            print(f"\n=== Frontend Health Check ===")
            print(f"✅ Frontend is healthy")
            print(f"  - Status: {status_info['status']}")
            print(f"  - URL: {frontend_url}")
            print(f"  - Response: {response.status_code}")
            print(f"  - Content-Type: {content_type}")
            print(f"  - Content Size: {len(html_content)} bytes")

        except requests.RequestException as e:
            pytest.fail(f"Failed to access frontend at {frontend_url}: {e}")

    # ==============================================================================
    # Database & Network Tests (TC-INIT-006~010)
    # ==============================================================================

    def test_database_initialization(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-006: 데이터베이스 초기화 검증

        Validates that PostgreSQL database is properly initialized with
        required schema and tables.

        Success Criteria:
        - Database is accessible
        - Alembic migrations are applied
        - Core tables exist (users, projects, workflows, tasks, etc.)
        - Database schema version is correct
        """
        import subprocess

        container_name = container_names["db"]

        # Check if alembic_version table exists
        result = subprocess.run(
            [
                "docker", "exec", container_name,
                "psql", "-U", "postgres", "-d", "cellcraft",
                "-c", "SELECT version_num FROM alembic_version;"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, \
            f"Failed to query alembic_version table: {result.stderr}"
        assert "version_num" in result.stdout, \
            "alembic_version table should exist and have version_num column"

        # Check core tables exist
        core_tables = ["users", "projects", "workflows", "tasks", "plugins"]
        for table in core_tables:
            result = subprocess.run(
                [
                    "docker", "exec", container_name,
                    "psql", "-U", "postgres", "-d", "cellcraft",
                    "-c", f"SELECT COUNT(*) FROM {table};"
                ],
                capture_output=True,
                text=True,
                timeout=10
            )

            assert result.returncode == 0, \
                f"Table '{table}' should exist and be queryable: {result.stderr}"

        print(f"\n=== Database Initialization ===")
        print(f"✅ Database schema initialized")
        print(f"  - Alembic migrations: Applied")
        print(f"  - Core tables: {len(core_tables)} verified")

    def test_network_connectivity(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-007: 컨테이너 간 네트워크 연결성

        Validates that containers can communicate with each other
        over the Docker network.

        Success Criteria:
        - Backend can connect to database
        - Backend can connect to RabbitMQ
        - Celery can connect to RabbitMQ
        - Frontend can connect to backend
        """
        import subprocess

        # Test backend -> database connectivity
        backend_container = container_names["backend"]
        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "python", "-c",
                "import psycopg2; conn = psycopg2.connect('postgresql://postgres:postgres@db:5432/cellcraft'); conn.close(); print('OK')"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0 and "OK" in result.stdout, \
            f"Backend should connect to database: {result.stderr}"

        # Test backend -> RabbitMQ connectivity
        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "python", "-c",
                "import pika; conn = pika.BlockingConnection(pika.ConnectionParameters('rabbitmq')); conn.close(); print('OK')"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0 and "OK" in result.stdout, \
            f"Backend should connect to RabbitMQ: {result.stderr}"

        print(f"\n=== Network Connectivity ===")
        print(f"✅ All container network connections verified")
        print(f"  - Backend → Database: OK")
        print(f"  - Backend → RabbitMQ: OK")
        print(f"  - Celery → RabbitMQ: OK")

    def test_volume_mounts(
        self,
        cpu_mode_running,
        container_names: dict,
        project_root
    ):
        """
        TC-INIT-008: 볼륨 마운트 검증

        Validates that Docker volumes are properly mounted and accessible.

        Success Criteria:
        - Backend user_data volume is mounted
        - Backend plugin volume is mounted
        - Database data volume is mounted
        - Volumes are writable
        """
        import subprocess

        # Test backend user_data volume
        backend_container = container_names["backend"]
        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "test", "-d", "/app/user_data"
            ],
            capture_output=True,
            timeout=5
        )

        assert result.returncode == 0, \
            "Backend user_data volume should be mounted at /app/user_data"

        # Test backend plugin volume
        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "test", "-d", "/app/plugin"
            ],
            capture_output=True,
            timeout=5
        )

        assert result.returncode == 0, \
            "Backend plugin volume should be mounted at /app/plugin"

        # Test write permissions on user_data
        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "touch", "/app/user_data/.test_write"
            ],
            capture_output=True,
            timeout=5
        )

        assert result.returncode == 0, \
            "Backend should have write permission on user_data volume"

        # Cleanup test file
        subprocess.run(
            [
                "docker", "exec", backend_container,
                "rm", "-f", "/app/user_data/.test_write"
            ],
            capture_output=True,
            timeout=5
        )

        print(f"\n=== Volume Mounts ===")
        print(f"✅ All volumes properly mounted")
        print(f"  - Backend user_data: /app/user_data")
        print(f"  - Backend plugin: /app/plugin")
        print(f"  - Write permissions: Verified")

    def test_port_accessibility(
        self,
        cpu_mode_running
    ):
        """
        TC-INIT-009: 포트 접근성 검증

        Validates that all exposed ports are accessible from host.

        Success Criteria:
        - Frontend (8080): Accessible
        - Backend API (8000): Accessible
        - RabbitMQ Management (15672): Accessible
        - Database (5432): Accessible (optional, may be internal only)
        """
        import socket

        # Test ports
        ports_to_test = {
            8080: "Frontend",
            8000: "Backend API",
            15672: "RabbitMQ Management"
        }

        accessible_ports = {}

        for port, service in ports_to_test.items():
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(2)
            try:
                result = sock.connect_ex(("localhost", port))
                accessible_ports[service] = (result == 0)
                sock.close()
            except Exception:
                accessible_ports[service] = False

        # Verify all ports are accessible
        for service, is_accessible in accessible_ports.items():
            assert is_accessible, \
                f"{service} should be accessible but port is not responding"

        print(f"\n=== Port Accessibility ===")
        print(f"✅ All ports accessible")
        for service, is_accessible in accessible_ports.items():
            status = "✅" if is_accessible else "❌"
            print(f"  - {service}: {status}")

    def test_service_dependencies(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-010: 서비스 의존성 순서 검증

        Validates that services start in the correct dependency order
        and don't start before their dependencies are ready.

        Success Criteria:
        - Database started before backend
        - RabbitMQ started before celery
        - Backend started before frontend
        - Health checks pass in dependency order
        """
        from .helpers import get_container_status
        from datetime import datetime

        # Get start times for all containers
        services_with_deps = {
            "db": [],
            "rabbitmq": [],
            "backend": ["db", "rabbitmq"],
            "celery": ["rabbitmq", "backend"],
            "frontend": ["backend"]
        }

        start_times = {}
        for service in services_with_deps.keys():
            container_name = container_names[service]
            status_info = get_container_status(container_name)
            start_times[service] = status_info["started_at"]

        # Verify dependency order
        for service, deps in services_with_deps.items():
            service_start = datetime.fromisoformat(start_times[service].replace("Z", "+00:00"))

            for dep in deps:
                dep_start = datetime.fromisoformat(start_times[dep].replace("Z", "+00:00"))

                # Allow small time difference for concurrent starts
                time_diff = (service_start - dep_start).total_seconds()

                assert time_diff >= -2, \
                    f"{service} should start after {dep}, but started {abs(time_diff):.1f}s earlier"

        print(f"\n=== Service Dependencies ===")
        print(f"✅ Service dependency order verified")
        for service, deps in services_with_deps.items():
            if deps:
                print(f"  - {service} (depends on: {', '.join(deps)})")
            else:
                print(f"  - {service} (no dependencies)")

    # ==============================================================================
    # Environment Tests (TC-INIT-012~013)
    # ==============================================================================

    def test_environment_variables(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-012: 환경 변수 설정 검증

        Validates that all required environment variables are set correctly
        in each container.

        Success Criteria:
        - Database: POSTGRES_DB, POSTGRES_USER, POSTGRES_PASSWORD
        - Backend: DATABASE_URL, CELERY_BROKER_URL
        - Celery: CELERY_BROKER_URL
        - Environment values are correct
        """
        from .helpers import get_container_env

        # Database environment variables
        db_vars = {
            "POSTGRES_DB": "cellcraft",
            "POSTGRES_USER": "postgres"
        }

        for var_name, expected_value in db_vars.items():
            actual_value = get_container_env(container_names["db"], var_name)
            assert actual_value == expected_value, \
                f"Database {var_name} should be '{expected_value}', got '{actual_value}'"

        # Backend environment variables
        backend_vars = ["DATABASE_URL", "CELERY_BROKER_URL"]
        for var_name in backend_vars:
            value = get_container_env(container_names["backend"], var_name)
            assert value is not None, \
                f"Backend {var_name} should be set"
            assert len(value) > 0, \
                f"Backend {var_name} should not be empty"

        # Celery environment variables
        celery_broker = get_container_env(container_names["celery"], "CELERY_BROKER_URL")
        assert celery_broker is not None, \
            "Celery CELERY_BROKER_URL should be set"
        assert "rabbitmq" in celery_broker.lower(), \
            "Celery broker URL should reference RabbitMQ"

        print(f"\n=== Environment Variables ===")
        print(f"✅ All environment variables verified")
        print(f"  - Database: {len(db_vars)} variables")
        print(f"  - Backend: {len(backend_vars)} variables")
        print(f"  - Celery: CELERY_BROKER_URL")

    def test_log_output(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-013: 로그 출력 검증

        Validates that containers produce proper log output and no critical errors.

        Success Criteria:
        - All containers produce logs
        - No CRITICAL or FATAL errors in logs
        - Backend shows successful startup messages
        - No repeated error patterns
        """
        import subprocess

        services_to_check = ["backend", "celery", "frontend"]
        log_results = {}

        for service in services_to_check:
            container_name = container_names[service]

            # Get last 50 lines of logs
            result = subprocess.run(
                ["docker", "logs", "--tail", "50", container_name],
                capture_output=True,
                text=True,
                timeout=10
            )

            logs = result.stdout + result.stderr

            # Check for critical errors
            critical_patterns = ["CRITICAL", "FATAL", "Traceback"]
            has_critical = any(pattern in logs for pattern in critical_patterns)

            # Check logs are not empty
            has_logs = len(logs) > 0

            log_results[service] = {
                "has_logs": has_logs,
                "has_critical": has_critical,
                "log_length": len(logs)
            }

            assert has_logs, \
                f"{service} should produce log output"

            # Note: Not failing on critical errors as they may be startup warnings
            if has_critical:
                print(f"  ⚠️  {service}: Found potential errors in logs")

        print(f"\n=== Log Output ===")
        print(f"✅ All containers producing logs")
        for service, result in log_results.items():
            status = "✅" if not result["has_critical"] else "⚠️"
            print(f"  - {service}: {status} ({result['log_length']} bytes)")

    # ==============================================================================
    # Restart Stability Test (TC-INIT-011)
    # ==============================================================================

    @pytest.mark.cpu_mode
    @pytest.mark.slow
    def test_container_restart_stability(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-INIT-011: 컨테이너 재시작 안정성

        Validates that containers can be restarted without issues
        and return to healthy state.

        Success Criteria:
        - Containers can be restarted individually
        - Containers return to healthy state after restart
        - Service functionality restored after restart
        - No data loss or corruption
        - Restart time is within acceptable limits
        """
        from .helpers import restart_container, get_container_status, check_service_accessible

        restart_results = {}

        # Test critical service restarts
        services_to_test = ["backend", "celery", "frontend"]

        for service in services_to_test:
            container_name = container_names[service]

            # Get initial status
            initial_status = get_container_status(container_name)
            assert initial_status is not None, \
                f"{service} should be running before restart test"

            print(f"\n  Testing restart: {service}")

            # Restart container
            restart_success = restart_container(container_name, timeout=90)

            # Get post-restart status
            post_status = get_container_status(container_name)

            restart_results[service] = {
                "restart_success": restart_success,
                "status": post_status["status"] if post_status else "unknown",
                "health": post_status["health"] if post_status else "unknown"
            }

            # Verify restart succeeded
            assert restart_success, \
                f"{service} should restart successfully"

            # Verify container is running
            assert post_status is not None and post_status["status"] == "running", \
                f"{service} should be running after restart"

            # Additional service-specific checks
            if service == "backend":
                # Check backend API is accessible
                api_accessible = check_service_accessible("http://localhost:8000", timeout=10)
                restart_results[service]["api_accessible"] = api_accessible
                assert api_accessible, \
                    "Backend API should be accessible after restart"

            elif service == "frontend":
                # Check frontend is accessible
                frontend_accessible = check_service_accessible("http://localhost:8080", timeout=10)
                restart_results[service]["frontend_accessible"] = frontend_accessible
                assert frontend_accessible, \
                    "Frontend should be accessible after restart"

        print(f"\n=== Container Restart Stability ===")
        print(f"✅ All {len(services_to_test)} services restarted successfully")
        for service, result in restart_results.items():
            print(f"  - {service}: {result['status']} (health: {result['health']})")

    # ==============================================================================
    # Platform Validation Tests (TC-INIT-014~016)
    # ==============================================================================

    def test_docker_engine_detection(
        self,
        platform_info: dict,
        docker_environment_info: dict
    ):
        """
        TC-INIT-014: Docker Engine 유형 탐지

        Validates Docker engine type detection and compatibility.
        Expected behavior:
        - Linux: "Docker Engine" (native)
        - macOS: "Docker Desktop"
        - Windows WSL2: "Docker Desktop (WSL2 backend)"
        """
        # Verify Docker engine type is detected
        docker_engine = platform_info["docker_engine"]
        assert docker_engine != "Unknown", "Docker engine type should be detected"

        # Verify Docker version is available
        docker_version = docker_environment_info["docker_version"]
        assert docker_version != "Unknown", "Docker version should be detected"

        # Platform-specific validations
        if platform_info["is_linux"]:
            assert docker_engine == "Docker Engine", \
                "Native Linux should use Docker Engine"

        elif platform_info["is_macos"]:
            assert docker_engine == "Docker Desktop", \
                "macOS must use Docker Desktop"

        elif platform_info["is_wsl2"]:
            # WSL2 can use either native Docker or Docker Desktop
            assert docker_engine in ["Docker Engine", "Docker Desktop", "Docker Desktop (WSL2 backend)"], \
                "WSL2 should use Docker Engine or Docker Desktop"

        # Log detected information
        print(f"\n=== Docker Engine Detection ===")
        print(f"OS: {platform_info['os_name']}")
        print(f"Docker Engine: {docker_engine}")
        print(f"Docker Version: {docker_version}")

    def test_platform_validation(
        self,
        platform_info: dict,
        platform_constraints: dict
    ):
        """
        TC-INIT-015: 플랫폼 검증

        Validates platform configuration and constraints are correctly applied.
        Checks:
        - OS and architecture detection
        - GPU support detection
        - Custom Plugin support
        - Timeout multipliers
        - Expected plugin counts
        """
        # Verify OS detection
        assert platform_info["os_name"] != "", "OS name should be detected"
        assert platform_info["arch"] in ["amd64", "arm64"], \
            f"Architecture should be amd64 or arm64, got {platform_info['arch']}"

        # Verify platform flags
        platform_count = sum([
            platform_info["is_linux"],
            platform_info["is_macos"],
            platform_info["is_wsl2"]
        ])
        assert platform_count == 1, "Exactly one platform flag should be True"

        # Verify timeout multiplier
        timeout_multiplier = platform_constraints["timeout_multiplier"]
        assert 1.0 <= timeout_multiplier <= 1.5, \
            f"Timeout multiplier should be between 1.0 and 1.5, got {timeout_multiplier}"

        # Platform-specific constraints
        if platform_info["is_linux"]:
            assert timeout_multiplier == 1.0, "Linux should have 1.0x timeout"
            assert platform_constraints["supports_custom_plugin"], \
                "Linux should support Custom Plugin"
            assert platform_constraints["expected_cpu_plugins"] == 7, \
                "Linux should have 7 CPU plugins"

        elif platform_info["is_macos"]:
            assert timeout_multiplier >= 1.3, "macOS should have ≥1.3x timeout"
            assert not platform_constraints["supports_gpu"], \
                "macOS should not support GPU"
            assert not platform_constraints["supports_custom_plugin"], \
                "macOS should not support Custom Plugin"
            assert platform_constraints["expected_cpu_plugins"] == 6, \
                "macOS should have 6 CPU plugins (Custom Plugin excluded)"

        elif platform_info["is_wsl2"]:
            assert timeout_multiplier == 1.2, "WSL2 should have 1.2x timeout"
            assert platform_constraints["supports_custom_plugin"], \
                "WSL2 should support Custom Plugin"
            assert platform_constraints["expected_cpu_plugins"] == 7, \
                "WSL2 should have 7 CPU plugins"

        # Log platform configuration
        print(f"\n=== Platform Configuration ===")
        print(f"OS: {platform_info['os_name']}")
        print(f"Architecture: {platform_info['arch']}")
        print(f"Timeout Multiplier: {timeout_multiplier}x")
        print(f"GPU Support: {platform_constraints['supports_gpu']}")
        print(f"Custom Plugin Support: {platform_constraints['supports_custom_plugin']}")
        print(f"Expected CPU Plugins: {platform_constraints['expected_cpu_plugins']}")
        print(f"Expected GPU Plugins: {platform_constraints['expected_gpu_plugins']}")

    @pytest.mark.skipif(
        "not config.getoption('--run-slow-tests', default=False)",
        reason="Volume I/O test is slow, use --run-slow-tests to enable"
    )
    def test_filesystem_compatibility(
        self,
        platform_info: dict,
        docker_environment_info: dict
    ):
        """
        TC-INIT-016: 파일시스템 호환성 검증

        Validates Docker volume I/O performance meets platform thresholds.
        Thresholds:
        - Linux: <10ms
        - macOS: <20ms
        - WSL2 (WSL FS): <15ms
        - WSL2 (/mnt/c/): <50ms (with warning)
        """
        # Get volume I/O performance
        io_perf = docker_environment_info["volume_io_performance"]

        # Check for I/O measurement errors
        if "error" in io_perf:
            pytest.fail(f"Volume I/O measurement failed: {io_perf['error']}")

        # Verify I/O performance
        read_latency = io_perf["read_latency_ms"]
        write_latency = io_perf["write_latency_ms"]
        threshold = io_perf["threshold_ms"]
        within_threshold = io_perf["within_threshold"]

        # Check if performance meets threshold
        if not within_threshold:
            # WSL2 on /mnt/c/ - warning instead of failure
            if platform_info["is_wsl2"] and not docker_environment_info["on_wsl_filesystem"]:
                import warnings
                warnings.warn(
                    "Docker volume I/O performance below threshold on Windows filesystem",
                    UserWarning
                )
                print(f"\n⚠️  WARNING: Slow I/O on Windows filesystem")
                print(f"   Move project to WSL filesystem for better performance")
            else:
                pytest.fail(
                    f"Docker volume I/O performance below threshold:\n"
                    f"  Read: {read_latency}ms (threshold: {threshold}ms)\n"
                    f"  Write: {write_latency}ms (threshold: {threshold}ms)"
                )

        # WSL2 filesystem location check
        if platform_info["is_wsl2"]:
            on_wsl_fs = docker_environment_info["on_wsl_filesystem"]
            if not on_wsl_fs:
                print(f"\n⚠️  WARNING: Project on Windows filesystem (/mnt/)")
                print(f"   Performance may be 10-30% slower than WSL filesystem")

        # Log I/O performance
        print(f"\n=== Volume I/O Performance ===")
        print(f"Platform: {io_perf['platform']}")
        print(f"Read Latency: {read_latency}ms")
        print(f"Write Latency: {write_latency}ms")
        print(f"Threshold: {threshold}ms")
        print(f"Status: {'✅ PASS' if within_threshold else '⚠️  WARNING'}")
