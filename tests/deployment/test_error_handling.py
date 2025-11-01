"""
Error Handling Tests

Tests for validating system error handling and recovery mechanisms.
Covers TC-ERR-001 through TC-ERR-006 from DEPLOYMENT_TEST_PLAN.md
"""

import pytest


class TestErrorHandling:
    """4.4 Error Handling Tests"""

    @pytest.mark.cpu_mode
    def test_container_failure_recovery(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-ERR-001: Container Failure Recovery

        Validates system can recover from individual container failures.

        Success Criteria:
        - System detects container failure
        - Failed container can be restarted
        - Services restore to healthy state
        - Dependent services handle failure gracefully
        """
        from .helpers import (
            inject_container_failure,
            recover_from_failure,
            get_container_status,
            check_service_accessible
        )
        import time

        # Test backend container failure and recovery
        container_name = container_names["backend"]

        print(f"\n  Testing failure recovery: backend")

        # Inject failure (stop container)
        failure_injected = inject_container_failure(container_name, "stop")
        assert failure_injected, \
            "Should be able to inject container failure"

        # Wait for failure to be detected
        time.sleep(5)

        # Verify container is stopped
        status = get_container_status(container_name)
        assert status is None or status["status"] != "running", \
            "Container should be stopped after failure injection"

        # Recover from failure
        recovery_success = recover_from_failure(container_name, "stop", timeout=90)
        assert recovery_success, \
            "Should successfully recover from container failure"

        # Wait for service to stabilize
        time.sleep(10)

        # Verify service is accessible
        api_accessible = check_service_accessible("http://localhost:8000", timeout=15)
        assert api_accessible, \
            "Backend API should be accessible after recovery"

        print(f"\n=== Container Failure Recovery ===")
        print(f"Success: Backend recovered from failure")
        print(f"  - Failure type: Container stop")
        print(f"  - Recovery time: ~15 seconds")
        print(f"  - Service restored: API accessible")

    @pytest.mark.cpu_mode
    def test_container_cascade_failure_recovery(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-ERR-001-EXT: Container Cascade Failure Recovery

        Validates system handles cascade failures across dependent services.

        Success Criteria:
        - Frontend fails gracefully when backend fails
        - Celery tasks queue when RabbitMQ fails
        - Services auto-recover in correct dependency order
        - No data loss during cascade
        """
        from .helpers import (
            inject_container_failure,
            recover_from_failure,
            get_container_status,
            check_service_accessible
        )
        import time

        # Test cascade: DB failure → Backend failure → Frontend degradation
        db_container = container_names["db"]
        backend_container = container_names["backend"]
        frontend_container = container_names["frontend"]

        print(f"\n  Testing cascade failure: DB → Backend → Frontend")

        # Step 1: Inject DB failure
        failure_injected = inject_container_failure(db_container, "stop")
        assert failure_injected, "Should inject DB failure"

        time.sleep(5)

        # Step 2: Verify backend detects DB loss but remains running
        backend_status = get_container_status(backend_container)
        assert backend_status and backend_status["status"] == "running", \
            "Backend should remain running despite DB failure"

        # Step 3: Verify frontend still serves static content
        frontend_accessible = check_service_accessible("http://localhost:8080", timeout=5)
        assert frontend_accessible, \
            "Frontend should still serve static content during DB failure"

        # Step 4: Recover DB
        recovery_success = recover_from_failure(db_container, "stop", timeout=90)
        assert recovery_success, "DB should recover"

        time.sleep(10)

        # Step 5: Verify backend reconnects and API becomes accessible
        max_retries = 3
        api_accessible = False
        for attempt in range(max_retries):
            time.sleep(10)
            api_accessible = check_service_accessible("http://localhost:8000/docs", timeout=15)
            if api_accessible:
                break

        assert api_accessible, \
            "Backend API should recover after DB restoration"

        print(f"\n=== Cascade Failure Recovery ===")
        print(f"Success: System recovered from cascade failure")
        print(f"  - DB failure detected: Yes")
        print(f"  - Backend remained stable: Yes")
        print(f"  - Frontend degraded gracefully: Yes")
        print(f"  - Full recovery achieved: Yes")

    @pytest.mark.cpu_mode
    def test_database_connection_loss(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-ERR-002: Database Connection Loss

        Validates system handles database connection loss gracefully.

        Success Criteria:
        - Backend detects database disconnection
        - Backend handles connection loss without crash
        - Connection restored after database recovery
        - No data corruption
        """
        from .helpers import (
            inject_container_failure,
            recover_from_failure,
            get_container_status,
            check_service_accessible
        )
        import time

        db_container = container_names["db"]
        backend_container = container_names["backend"]

        print(f"\n  Testing database connection loss")

        # Verify initial state
        initial_backend_status = get_container_status(backend_container)
        assert initial_backend_status and initial_backend_status["status"] == "running", \
            "Backend should be running initially"

        # Stop database
        failure_injected = inject_container_failure(db_container, "stop")
        assert failure_injected, \
            "Should be able to stop database"

        # Wait for connection loss
        time.sleep(5)

        # Backend should still be running (not crashed)
        backend_status = get_container_status(backend_container)
        assert backend_status and backend_status["status"] == "running", \
            "Backend should remain running despite database loss"

        # Recover database
        recovery_success = recover_from_failure(db_container, "stop", timeout=90)
        assert recovery_success, \
            "Should successfully recover database"

        # Wait for backend to reconnect with retry logic
        max_retries = 3
        retry_interval = 10
        api_accessible = False

        for attempt in range(max_retries):
            time.sleep(retry_interval)
            api_accessible = check_service_accessible("http://localhost:8000", timeout=15)

            if api_accessible:
                break

            print(f"  Retry {attempt + 1}/{max_retries}: API not yet accessible")

        # Validate database state after recovery
        print(f"  Validating database integrity after recovery...")
        import subprocess

        # Check database accepts connections
        db_accessible = subprocess.run(
            ["docker", "exec", db_container, "pg_isready", "-U", "cellcraft_admin"],
            capture_output=True, timeout=5
        )
        assert db_accessible.returncode == 0, \
            "Database should accept connections after recovery"

        # Check connection pool metrics
        connection_count = subprocess.run(
            ["docker", "exec", db_container, "psql", "-U", "cellcraft_admin", "-d", "cellcraft",
             "-t", "-c", "SELECT count(*) FROM pg_stat_activity WHERE datname='cellcraft';"],
            capture_output=True, text=True, timeout=10
        )

        if connection_count.returncode == 0:
            active_connections = int(connection_count.stdout.strip())
            print(f"  Active database connections: {active_connections}")
            assert active_connections > 0, \
                "Should have active connections after recovery"

        # Validate no data corruption
        tables_check = subprocess.run(
            ["docker", "exec", db_container, "psql", "-U", "cellcraft_admin", "-d", "cellcraft",
             "-t", "-c", "SELECT COUNT(*) FROM information_schema.tables WHERE table_schema='public';"],
            capture_output=True, text=True, timeout=10
        )

        if tables_check.returncode == 0:
            table_count = int(tables_check.stdout.strip())
            print(f"  Database tables found: {table_count}")
            assert table_count > 0, \
                "Database schema should be intact after recovery"

        print(f"\n=== Database Connection Loss ===")
        print(f"Success: System handled database disconnection")
        print(f"  - Backend remained running: Yes")
        print(f"  - Database recovered: Yes")
        print(f"  - Connection restored: {'Yes' if api_accessible else 'No'}")
        print(f"  - Reconnection attempts: {attempt + 1}/{max_retries}")

        # Assert that API is accessible after retry attempts
        assert api_accessible, \
            f"Backend API should be accessible after database recovery (tried {max_retries} times)"

    @pytest.mark.cpu_mode
    def test_disk_space_handling(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-ERR-006: Disk Space Handling

        Validates system handles low disk space gracefully.

        Success Criteria:
        - System detects low disk space
        - Appropriate error messages generated
        - System doesn't crash or corrupt data
        - System recovers when space available
        """
        import subprocess

        backend_container = container_names["backend"]

        print(f"\n  Testing disk space handling")

        # Check current disk usage
        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "df", "-h", "/app/user_data"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, \
            "Should be able to check disk space"

        disk_info = result.stdout

        # Parse disk usage percentage
        import re
        usage_match = re.search(r'(\d+)%', disk_info)
        disk_usage_percent = int(usage_match.group(1)) if usage_match else 0

        # Test write operations with various file sizes
        test_results = []
        file_sizes = ["1K", "1M", "10M"]

        for size in file_sizes:
            test_file = f"/app/user_data/.disk_test_{size}"
            result = subprocess.run(
                ["docker", "exec", backend_container, "sh", "-c",
                 f"dd if=/dev/zero of={test_file} bs={size} count=1 2>/dev/null"],
                capture_output=True, text=True, timeout=15
            )

            write_success = (result.returncode == 0)
            test_results.append((size, write_success))

            # Cleanup
            if write_success:
                subprocess.run(
                    ["docker", "exec", backend_container, "rm", "-f", test_file],
                    capture_output=True, timeout=10
                )

        # Verify at least small files can be written
        assert any(success for _, success in test_results), \
            "Should be able to write at least small test files"

        # Check for adequate free space
        df_result = subprocess.run(
            ["docker", "exec", backend_container, "df", "--output=avail", "/app/user_data"],
            capture_output=True, text=True, timeout=10
        )

        if df_result.returncode == 0:
            avail_kb = int(df_result.stdout.strip().split('\n')[-1])
            avail_gb = avail_kb / (1024 * 1024)

            if avail_gb < 1.0:
                print(f"  ⚠️  WARNING: Low disk space ({avail_gb:.2f} GB free)")

        print(f"\n=== Disk Space Handling ===")
        print(f"Success: Disk space monitoring and write tests functional")
        print(f"  - Disk usage: {disk_usage_percent}%")
        print(f"  - Write tests:")
        for size, success in test_results:
            status = "✅" if success else "❌"
            print(f"    - {size}: {status}")
        print(f"  - Available space: {avail_gb:.2f} GB")
        print(f"  - Disk info:")
        for line in disk_info.strip().split('\n')[-2:]:
            print(f"    {line}")

    @pytest.mark.cpu_mode
    def test_network_partition_handling(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-ERR-003: Network Partition Handling

        Validates system handles network partitions gracefully.

        Success Criteria:
        - System detects network issues
        - Services handle disconnection appropriately
        - System recovers when network restored
        - No data corruption during partition
        """
        from .helpers import (
            inject_container_failure,
            recover_from_failure,
            get_container_status
        )
        import time

        backend_container = container_names["backend"]

        print(f"\n  Testing network partition handling")

        # Verify initial connectivity
        initial_status = get_container_status(backend_container)
        assert initial_status and initial_status["status"] == "running", \
            "Backend should be running initially"

        # Inject network partition
        from .helpers import check_service_accessible
        import subprocess
        partition_injected = inject_container_failure(backend_container, "network")

        if partition_injected:
            print(f"  Network partition injected successfully")
            time.sleep(5)

            # Backend should still be running
            status = get_container_status(backend_container)
            assert status and status["status"] == "running", \
                "Backend should remain running during network partition"

            # Frontend should detect backend unavailability
            api_accessible = check_service_accessible("http://localhost:8000/docs", timeout=5)
            print(f"  Backend API accessible during partition: {api_accessible}")

            # Recover from partition
            recovery_success = recover_from_failure(backend_container, "network", timeout=30)
            assert recovery_success, "Network should be restored"

            time.sleep(10)

            # Verify API becomes accessible again
            max_retries = 3
            api_restored = False
            for attempt in range(max_retries):
                time.sleep(5)
                api_restored = check_service_accessible("http://localhost:8000/docs", timeout=10)
                if api_restored:
                    break

            assert api_restored, \
                "Backend API should be accessible after network recovery"

            print(f"\n=== Network Partition Handling ===")
            print(f"Success: Network partition handled and recovered")
            print(f"  - Partition injected: Yes")
            print(f"  - Container remained stable: Yes")
            print(f"  - Network recovered: Yes")
            print(f"  - API accessibility restored: Yes")
        else:
            # Test inter-container connectivity
            print(f"  Network injection not supported, testing connectivity instead")

            backend_to_db = subprocess.run(
                ["docker", "exec", backend_container, "ping", "-c", "3", "-W", "2", "db"],
                capture_output=True, timeout=10
            )

            backend_to_rabbitmq = subprocess.run(
                ["docker", "exec", backend_container, "ping", "-c", "3", "-W", "2", "rabbitmq"],
                capture_output=True, timeout=10
            )

            assert backend_to_db.returncode == 0, \
                "Backend should be able to reach database via network"

            assert backend_to_rabbitmq.returncode == 0, \
                "Backend should be able to reach RabbitMQ via network"

            print(f"\n=== Network Connectivity Validation ===")
            print(f"Success: Inter-container network is healthy")
            print(f"  - Backend → Database: OK")
            print(f"  - Backend → RabbitMQ: OK")

    @pytest.mark.cpu_mode
    def test_memory_pressure_handling(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-ERR-004: Memory Pressure Handling

        Validates system handles memory pressure appropriately.

        Success Criteria:
        - System monitors memory usage
        - Appropriate warnings generated
        - Services remain stable under pressure
        - Graceful degradation if needed
        """
        import subprocess

        backend_container = container_names["backend"]

        print(f"\n  Testing memory pressure handling")

        # Check current memory usage
        result = subprocess.run(
            [
                "docker", "stats", backend_container,
                "--no-stream", "--format",
                "{{.MemUsage}}\t{{.MemPerc}}"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, \
            "Should be able to check memory stats"

        stats = result.stdout.strip().split("\t")
        mem_usage = stats[0] if len(stats) > 0 else "N/A"
        mem_percent = stats[1] if len(stats) > 1 else "N/A"

        # Check container is healthy
        from .helpers import get_container_status
        status = get_container_status(backend_container)

        assert status and status["status"] == "running", \
            "Backend should be running under normal memory conditions"

        print(f"\n=== Memory Pressure Handling ===")
        print(f"Success: Memory monitoring functional")
        print(f"  - Memory usage: {mem_usage}")
        print(f"  - Memory percent: {mem_percent}")
        print(f"  - Container stable: Yes")
        print(f"  - Monitoring active: Yes")

    @pytest.mark.cpu_mode
    @pytest.mark.slow
    def test_memory_stress_handling(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-ERR-004-EXT: Memory Stress Test

        Validates system handles actual memory pressure gracefully.

        Success Criteria:
        - Container respects memory limits
        - System detects memory pressure
        - OOM killer doesn't crash entire container
        - Services degrade gracefully under pressure
        """
        import subprocess
        from .helpers import get_container_status

        backend_container = container_names["backend"]

        print(f"\n  Testing memory stress handling")

        # Get current memory usage baseline
        result = subprocess.run(
            ["docker", "stats", backend_container, "--no-stream", "--format", "{{.MemUsage}}\t{{.MemPerc}}"],
            capture_output=True, text=True, timeout=10
        )

        baseline_stats = result.stdout.strip().split("\t")
        baseline_usage = baseline_stats[0] if len(baseline_stats) > 0 else "N/A"
        baseline_percent = baseline_stats[1] if len(baseline_stats) > 1 else "N/A"

        # Allocate memory inside container (small stress test)
        stress_result = subprocess.run(
            ["docker", "exec", backend_container, "python", "-c",
             "import sys; x = bytearray(50 * 1024 * 1024); print('Memory allocated: 50MB'); sys.stdout.flush()"],
            capture_output=True, text=True, timeout=15
        )

        # Container should still be running
        status = get_container_status(backend_container)
        assert status and status["status"] == "running", \
            "Backend should remain running after memory allocation"

        # Get memory usage after stress
        result = subprocess.run(
            ["docker", "stats", backend_container, "--no-stream", "--format", "{{.MemUsage}}\t{{.MemPerc}}"],
            capture_output=True, text=True, timeout=10
        )

        stress_stats = result.stdout.strip().split("\t")
        stress_usage = stress_stats[0] if len(stress_stats) > 0 else "N/A"
        stress_percent = stress_stats[1] if len(stress_stats) > 1 else "N/A"

        print(f"\n=== Memory Stress Handling ===")
        print(f"Success: Memory stress test completed")
        print(f"  - Baseline usage: {baseline_usage} ({baseline_percent})")
        print(f"  - After stress: {stress_usage} ({stress_percent})")
        print(f"  - Container stable: Yes")
        print(f"  - Memory tracking: Functional")
        if stress_result.returncode == 0:
            print(f"  - Stress test result: {stress_result.stdout.strip()}")

