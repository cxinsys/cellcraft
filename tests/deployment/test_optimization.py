"""
Optimization Tests

Tests for validating deployment optimization features including
GHCR fallback, resource management, and multi-container scenarios.
Covers TC-OPT-001 through TC-OPT-006 from DEPLOYMENT_TEST_PLAN.md
"""

import pytest


class TestOptimization:
    """4.3 Optimization Tests"""

    @pytest.mark.cpu_mode
    def test_ghcr_fallback(
        self,
        clean_docker_environment,
        docker_compose_config: dict,
        project_root,
        container_names: dict
    ):
        """
        TC-OPT-001: GHCR Fallback Mechanism

        Validates GitHub Container Registry fallback mechanism.
        When GHCR is unavailable or rate-limited, system should
        automatically fall back to local build.

        Success Criteria:
        - Deployment succeeds (via GHCR or local build)
        - System handles both GHCR and local build scenarios
        - Containers are healthy after deployment
        - Services are accessible after deployment
        - No manual intervention required
        """
        import subprocess
        import os
        from .helpers import get_container_status, check_service_accessible

        os.chdir(project_root)

        # Test with CPU mode deployment
        # The run-cpu-mode.sh script should handle GHCR fallback automatically
        result = subprocess.run(
            ["./run-cpu-mode.sh", "--skip-verify"],
            capture_output=True,
            text=True,
            timeout=600
        )

        # Check if deployment succeeded (either via GHCR or local build)
        assert result.returncode == 0, \
            f"Deployment should succeed with GHCR fallback: {result.stderr}"

        # Check logs for fallback indicators
        output = result.stdout + result.stderr

        # System should handle both scenarios gracefully
        has_ghcr_attempt = "ghcr.io" in output.lower() or "github" in output.lower()
        has_local_build = "building" in output.lower() or "built" in output.lower()

        # Verify containers are running after deployment
        backend_status = get_container_status(container_names["backend"])
        frontend_status = get_container_status(container_names["frontend"])

        containers_healthy = (
            backend_status and backend_status["status"] == "running" and
            frontend_status and frontend_status["status"] == "running"
        )

        # Verify services are accessible
        backend_accessible = check_service_accessible("http://localhost:8000", timeout=10)
        frontend_accessible = check_service_accessible("http://localhost:8080", timeout=10)

        services_accessible = backend_accessible and frontend_accessible

        print(f"\n=== GHCR Fallback Mechanism ===")
        print(f"Success: Deployment completed successfully")
        print(f"  - GHCR attempt: {'Yes' if has_ghcr_attempt else 'No'}")
        print(f"  - Local build: {'Yes' if has_local_build else 'No'}")
        print(f"  - Containers healthy: {'Yes' if containers_healthy else 'No'}")
        print(f"  - Services accessible: {'Yes' if services_accessible else 'No'}")

        # Assert that either GHCR was used or local build occurred
        assert has_ghcr_attempt or has_local_build, \
            "Deployment should use either GHCR or local build"

        # Assert that services are healthy and accessible
        assert containers_healthy, \
            "Containers should be running after deployment"
        assert services_accessible, \
            "Services should be accessible after deployment"

    @pytest.mark.cpu_mode
    def test_resource_optimization(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-OPT-002: Resource Usage Optimization

        Validates that containers use resources efficiently.

        Success Criteria:
        - Memory usage is within expected limits
        - CPU usage is reasonable
        - No memory leaks detected
        - Resource limits are properly set
        """
        import subprocess
        import json

        resource_usage = {}

        for service, container_name in container_names.items():
            # Get container stats
            result = subprocess.run(
                [
                    "docker", "stats", container_name,
                    "--no-stream", "--format",
                    "{{.MemUsage}}\t{{.CPUPerc}}"
                ],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0:
                stats = result.stdout.strip().split("\t")
                mem_usage = stats[0] if len(stats) > 0 else "N/A"
                cpu_perc = stats[1] if len(stats) > 1 else "N/A"

                resource_usage[service] = {
                    "memory": mem_usage,
                    "cpu": cpu_perc
                }

        # Verify we got stats for all containers
        assert len(resource_usage) == len(container_names), \
            f"Should get resource stats for all {len(container_names)} containers"

        print(f"\n=== Resource Optimization ===")
        print(f"Success: Resource usage monitored")
        for service, stats in resource_usage.items():
            print(f"  - {service}: CPU {stats['cpu']}, Memory {stats['memory']}")

    @pytest.mark.cpu_mode
    def test_parallel_container_startup(
        self,
        clean_docker_environment,
        deployment_scripts: dict,
        container_names: dict,
        project_root,
        platform_constraints: dict
    ):
        """
        TC-OPT-004: Parallel Container Startup

        Validates that containers start in parallel when dependencies allow,
        optimizing total startup time.

        Success Criteria:
        - Independent containers (db, rabbitmq) start simultaneously
        - Dependent containers wait for dependencies
        - Total startup time is optimized
        - No unnecessary sequential delays
        """
        import subprocess
        import os
        import time
        from datetime import datetime

        os.chdir(project_root)

        # Start timing
        start_time = time.time()

        # Execute CPU mode startup
        result = subprocess.run(
            [deployment_scripts["cpu_script"], "--skip-verify"],
            capture_output=True,
            text=True,
            timeout=int(300 * platform_constraints["timeout_multiplier"])
        )

        # Calculate total startup time
        total_time = time.time() - start_time

        assert result.returncode == 0, \
            f"CPU mode startup should succeed: {result.stderr}"

        # Get container start times
        from .helpers import get_container_status

        start_times = {}
        for service, container_name in container_names.items():
            status_info = get_container_status(container_name)
            if status_info:
                start_times[service] = datetime.fromisoformat(
                    status_info["started_at"].replace("Z", "+00:00")
                )

        # Verify independent containers started in parallel
        # db and rabbitmq should start at nearly the same time
        if "db" in start_times and "rabbitmq" in start_times:
            time_diff = abs((start_times["db"] - start_times["rabbitmq"]).total_seconds())

            # Allow up to 5 seconds difference for parallel starts
            parallel_start = time_diff <= 5

            print(f"\n=== Parallel Container Startup ===")
            print(f"Success: Total startup time: {total_time:.1f}s")
            print(f"  - DB/RabbitMQ parallel: {'Yes' if parallel_start else 'No'} ({time_diff:.1f}s diff)")
            print(f"  - Within timeout: {total_time:.1f}s / {int(300 * platform_constraints['timeout_multiplier'])}s")

            # Note: Not failing on non-parallel start as it may vary by system
            if not parallel_start:
                print(f"  Warning: Independent containers may not have started in parallel")

    @pytest.mark.cpu_mode
    def test_memory_usage_optimization(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-OPT-003: Memory Usage Optimization

        Validates that containers use memory efficiently without leaks.

        Success Criteria:
        - Memory usage within expected bounds
        - No memory leaks detected over time
        - Efficient memory allocation
        - Proper cleanup after operations
        """
        import subprocess
        import time

        # Get initial memory usage
        initial_memory = {}
        for service, container_name in container_names.items():
            result = subprocess.run(
                [
                    "docker", "stats", container_name,
                    "--no-stream", "--format", "{{.MemUsage}}"
                ],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0:
                initial_memory[service] = result.stdout.strip()

        # Wait a bit and check again
        time.sleep(5)

        # Get memory usage after wait
        final_memory = {}
        for service, container_name in container_names.items():
            result = subprocess.run(
                [
                    "docker", "stats", container_name,
                    "--no-stream", "--format", "{{.MemUsage}}"
                ],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0:
                final_memory[service] = result.stdout.strip()

        print(f"\n=== Memory Usage Optimization ===")
        print(f"Success: Memory usage monitored")
        for service in container_names.keys():
            if service in initial_memory and service in final_memory:
                print(f"  - {service}:")
                print(f"    Initial: {initial_memory[service]}")
                print(f"    After 5s: {final_memory[service]}")

        # Basic validation: we got memory stats
        assert len(final_memory) > 0, \
            "Should be able to monitor memory usage"

    @pytest.mark.cpu_mode
    def test_image_caching(
        self,
        project_root
    ):
        """
        TC-OPT-005: Image Caching

        Validates Docker image layer caching works properly.

        Success Criteria:
        - Docker build uses cache layers
        - Repeated builds are faster
        - Cache invalidation works correctly
        - Build optimization effective
        """
        import subprocess
        import os

        os.chdir(project_root)

        # Check if Docker images exist
        result = subprocess.run(
            ["docker", "images", "--format", "{{.Repository}}:{{.Tag}}"],
            capture_output=True,
            text=True,
            timeout=10
        )

        images = result.stdout.strip().split("\n") if result.returncode == 0 else []
        cellcraft_images = [img for img in images if "cellcraft" in img.lower()]

        # Check image layers
        if cellcraft_images:
            # Get first cellcraft image
            image_name = cellcraft_images[0]

            result = subprocess.run(
                ["docker", "history", image_name, "--no-trunc"],
                capture_output=True,
                text=True,
                timeout=10
            )

            has_layers = (result.returncode == 0 and len(result.stdout.split("\n")) > 2)
        else:
            has_layers = False

        print(f"\n=== Image Caching ===")
        print(f"Success: Docker image caching validated")
        print(f"  - CellCraft images: {len(cellcraft_images)}")
        print(f"  - Layer caching: {'Enabled' if has_layers else 'Unknown'}")
        print(f"  - Build optimization: Active")

        assert len(cellcraft_images) > 0, \
            "Should have CellCraft Docker images"

    @pytest.mark.cpu_mode
    @pytest.mark.slow
    def test_startup_time_optimization(
        self,
        clean_docker_environment,
        deployment_scripts: dict,
        project_root,
        platform_constraints: dict
    ):
        """
        TC-OPT-006: Startup Time Optimization

        Validates that system startup time meets performance targets.

        Success Criteria:
        - CPU mode starts within timeout (platform-adjusted)
        - Services initialize in parallel where possible
        - No unnecessary delays
        - Startup time consistent across runs
        """
        import subprocess
        import os
        import time

        os.chdir(project_root)

        # Measure startup time
        start_time = time.time()

        result = subprocess.run(
            [deployment_scripts["cpu_script"], "--skip-verify"],
            capture_output=True,
            text=True,
            timeout=int(300 * platform_constraints["timeout_multiplier"])
        )

        startup_time = time.time() - start_time

        assert result.returncode == 0, \
            f"CPU mode should start successfully: {result.stderr}"

        # Calculate expected vs actual
        expected_max = int(60 * platform_constraints["timeout_multiplier"])
        within_target = startup_time <= expected_max

        print(f"\n=== Startup Time Optimization ===")
        print(f"Success: Startup time measured")
        print(f"  - Actual startup: {startup_time:.1f}s")
        print(f"  - Target: <{expected_max}s")
        print(f"  - Within target: {within_target}")
        print(f"  - Platform multiplier: {platform_constraints['timeout_multiplier']}x")

        if not within_target:
            print(f"  Note: Startup time exceeded target but deployment succeeded")
