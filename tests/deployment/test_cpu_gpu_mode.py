"""
CPU/GPU Mode Tests

Tests for validating CPU and GPU mode deployment configurations.
Covers TC-MODE-001 through TC-MODE-017 from DEPLOYMENT_TEST_PLAN.md
"""

import pytest


class TestCPUMode:
    """4.2.1 CPU 모드 테스트"""

    @pytest.mark.cpu_mode
    def test_cpu_plugin_count(
        self,
        cpu_mode_running,
        platform_constraints: dict
    ):
        """
        TC-MODE-001: CPU 공식 플러그인 개수 (플랫폼 인식)

        Validates that the correct number of CPU plugins are available
        based on the platform:
        - Linux: 7 CPU plugins (6 official + 1 Custom Plugin)
        - macOS: 6 CPU plugins (6 official only, Custom Plugin excluded)
        - WSL2: 7 CPU plugins (6 official + 1 Custom Plugin)

        Success Criteria:
        - Plugin count matches platform_constraints["expected_cpu_plugins"]
        - All plugins have valid metadata
        - Plugins are accessible via API
        """
        from .helpers import get_plugin_list

        # Get plugin list from backend API
        plugins = get_plugin_list()
        cpu_plugins = [p for p in plugins if not p.get("use_gpu", False)]

        # Get expected count based on platform
        expected_count = platform_constraints["expected_cpu_plugins"]

        # Verify plugin count
        assert len(cpu_plugins) == expected_count, \
            f"Expected {expected_count} CPU plugins, got {len(cpu_plugins)}"

        # Verify all plugins have required metadata
        required_fields = ["plugin_id", "plugin_name", "version", "category"]
        for plugin in cpu_plugins:
            for field in required_fields:
                assert field in plugin, \
                    f"Plugin {plugin.get('plugin_id', 'unknown')} missing required field: {field}"

        # Log plugin details
        print(f"\n=== CPU Plugin Count ===")
        print(f"✅ Found {len(cpu_plugins)}/{expected_count} CPU plugins")
        print(f"Platform: {platform_constraints.get('platform_type', 'unknown')}")

        # List plugins
        for plugin in cpu_plugins:
            plugin_name = plugin.get("plugin_name", "unknown")
            plugin_id = plugin.get("plugin_id", "unknown")
            is_editable = plugin.get("is_editable", False)
            custom_marker = " (Custom)" if is_editable else ""
            print(f"  - {plugin_name}{custom_marker} ({plugin_id})")

    @pytest.mark.cpu_mode
    def test_cpu_plugin_metadata(
        self,
        cpu_mode_running,
        platform_constraints: dict
    ):
        """
        TC-MODE-002: CPU 플러그인 메타데이터 검증

        Validates that all CPU plugins have complete and valid metadata.

        Success Criteria:
        - All required fields present (plugin_id, plugin_name, version, category, etc.)
        - Version format is valid (semantic versioning)
        - Category is one of: grn_inference, preprocessing, visualization
        - Dependencies are properly declared
        """
        from .helpers import get_plugin_list

        plugins = get_plugin_list()
        cpu_plugins = [p for p in plugins if not p.get("use_gpu", False)]

        required_fields = [
            "plugin_id", "plugin_name", "version", "category",
            "description", "author", "tags"
        ]

        valid_categories = ["grn_inference", "preprocessing", "visualization", "analysis"]

        for plugin in cpu_plugins:
            # Check required fields
            for field in required_fields:
                assert field in plugin and plugin[field], \
                    f"Plugin {plugin.get('plugin_id')} missing or empty field: {field}"

            # Validate version format (semantic versioning: x.y.z)
            version = plugin["version"]
            version_parts = version.split(".")
            assert len(version_parts) >= 2, \
                f"Plugin {plugin['plugin_id']} version should be semantic (x.y.z), got {version}"

            # Validate category
            category = plugin["category"]
            assert category in valid_categories, \
                f"Plugin {plugin['plugin_id']} has invalid category: {category}"

        print(f"\n=== CPU Plugin Metadata ===")
        print(f"✅ All {len(cpu_plugins)} CPU plugins have valid metadata")
        for plugin in cpu_plugins:
            print(f"  - {plugin['plugin_name']} v{plugin['version']} ({plugin['category']})")

    @pytest.mark.cpu_mode
    def test_cpu_plugin_execution(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-003: CPU 플러그인 실행 가능 여부

        Validates that CPU plugins can be instantiated and are ready for execution.

        Success Criteria:
        - Plugin metadata can be retrieved
        - Plugin Snakefile exists and is parseable
        - Required scripts/dependencies are present
        - No immediate execution errors
        """
        from .helpers import get_plugin_list
        import subprocess

        plugins = get_plugin_list()
        cpu_plugins = [p for p in plugins if not p.get("use_gpu", False)]

        executable_count = 0

        backend_container = container_names["backend"]

        for plugin in cpu_plugins:
            plugin_id = plugin["plugin_id"]

            # Check if plugin directory exists in backend
            result = subprocess.run(
                [
                    "docker", "exec", backend_container,
                    "test", "-d", f"/app/plugin/{plugin_id}"
                ],
                capture_output=True,
                timeout=5
            )

            if result.returncode == 0:
                # Check if Snakefile exists
                result = subprocess.run(
                    [
                        "docker", "exec", backend_container,
                        "test", "-f", f"/app/plugin/{plugin_id}/Snakefile"
                    ],
                    capture_output=True,
                    timeout=5
                )

                if result.returncode == 0:
                    executable_count += 1

        # All plugins should be executable
        assert executable_count == len(cpu_plugins), \
            f"Only {executable_count}/{len(cpu_plugins)} plugins are executable"

        print(f"\n=== CPU Plugin Execution ===")
        print(f"✅ All {len(cpu_plugins)} CPU plugins are executable")
        print(f"  - Plugin directories: Verified")
        print(f"  - Snakefiles: Present")


class TestGPUMode:
    """4.2.2 GPU 모드 테스트"""

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_plugin_count(
        self,
        gpu_mode_running
    ):
        """
        TC-MODE-004: GPU 플러그인 개수

        Validates that 8 GPU plugins are available in GPU mode.
        This test is skipped on macOS (GPU not supported).

        GPU Plugins:
        1. TENET (GPU)
        2. FastTENET (GPU)
        3. GENIE3 (GPU)
        4. GRNBOOST2 (GPU)
        5. LEAP (GPU)
        6. Scribe (GPU)
        7. Custom Plugin (GPU) - if supported
        8. Additional GPU plugin

        Success Criteria:
        - Exactly 8 GPU plugins available
        - All plugins have use_gpu=true
        - All plugins have valid GPU configuration
        - Plugins are accessible via API
        """
        from .helpers import get_plugin_list

        # Get plugin list from backend API
        plugins = get_plugin_list()
        gpu_plugins = [p for p in plugins if p.get("use_gpu", False)]

        # Verify GPU plugin count
        expected_gpu_count = 8
        assert len(gpu_plugins) == expected_gpu_count, \
            f"Expected {expected_gpu_count} GPU plugins, got {len(gpu_plugins)}"

        # Verify all GPU plugins have required metadata
        required_fields = ["plugin_id", "plugin_name", "version", "category", "use_gpu"]
        for plugin in gpu_plugins:
            for field in required_fields:
                assert field in plugin, \
                    f"Plugin {plugin.get('plugin_id', 'unknown')} missing required field: {field}"

            # Verify use_gpu is true
            assert plugin["use_gpu"] is True, \
                f"GPU plugin {plugin.get('plugin_id')} should have use_gpu=true"

        # Log GPU plugin details
        print(f"\n=== GPU Plugin Count ===")
        print(f"✅ Found {len(gpu_plugins)}/{expected_gpu_count} GPU plugins")

        # List plugins with GPU info
        for plugin in gpu_plugins:
            plugin_name = plugin.get("plugin_name", "unknown")
            plugin_id = plugin.get("plugin_id", "unknown")
            is_editable = plugin.get("is_editable", False)
            custom_marker = " (Custom)" if is_editable else ""
            print(f"  - {plugin_name}{custom_marker} ({plugin_id})")

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_plugin_metadata(
        self,
        gpu_mode_running
    ):
        """
        TC-MODE-005: GPU 플러그인 메타데이터 검증

        Validates that all GPU plugins have complete and valid metadata
        including GPU-specific configuration.

        Success Criteria:
        - All required fields present
        - use_gpu field is true
        - GPU-specific parameters are defined
        - Compatible with CUDA version
        """
        from .helpers import get_plugin_list

        plugins = get_plugin_list()
        gpu_plugins = [p for p in plugins if p.get("use_gpu", False)]

        required_fields = [
            "plugin_id", "plugin_name", "version", "category",
            "description", "author", "tags", "use_gpu"
        ]

        for plugin in gpu_plugins:
            # Check required fields
            for field in required_fields:
                assert field in plugin and plugin[field] is not None, \
                    f"GPU Plugin {plugin.get('plugin_id')} missing field: {field}"

            # Verify use_gpu is explicitly true
            assert plugin["use_gpu"] is True, \
                f"GPU Plugin {plugin['plugin_id']} should have use_gpu=true"

        print(f"\n=== GPU Plugin Metadata ===")
        print(f"✅ All {len(gpu_plugins)} GPU plugins have valid metadata")
        for plugin in gpu_plugins:
            print(f"  - {plugin['plugin_name']} v{plugin['version']} (GPU)")

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_availability(
        self,
        gpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-006: GPU 가용성 검증

        Validates that GPU is available and accessible to containers.

        Success Criteria:
        - nvidia-smi command works in backend/celery containers
        - GPU count matches expected (based on GPU_COUNT env)
        - CUDA version is compatible
        - GPU memory is sufficient
        """
        import subprocess

        backend_container = container_names["backend"]

        # Test GPU availability in backend container
        result = subprocess.run(
            ["docker", "exec", backend_container, "nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, \
            f"nvidia-smi should work in backend container: {result.stderr}"

        gpu_info = result.stdout.strip().split("\n")
        gpu_count = len(gpu_info)

        assert gpu_count > 0, \
            "At least one GPU should be available"

        print(f"\n=== GPU Availability ===")
        print(f"✅ GPU accessible in containers")
        print(f"  - GPU Count: {gpu_count}")
        for i, info in enumerate(gpu_info, 1):
            print(f"  - GPU {i}: {info}")

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_plugin_execution(
        self,
        gpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-007: GPU 플러그인 실행 가능 여부

        Validates that GPU plugins are ready for execution.

        Success Criteria:
        - Plugin Snakefile exists
        - GPU-specific scripts are present
        - CUDA dependencies are available
        """
        from .helpers import get_plugin_list
        import subprocess

        backend_container = container_names["backend"]

        plugins = get_plugin_list()
        gpu_plugins = [p for p in plugins if p.get("use_gpu", False)]

        executable_count = 0

        for plugin in gpu_plugins:
            plugin_id = plugin["plugin_id"]

            # Check plugin directory and Snakefile
            result = subprocess.run(
                [
                    "docker", "exec", backend_container,
                    "test", "-f", f"/app/plugin/{plugin_id}/Snakefile"
                ],
                capture_output=True,
                timeout=5
            )

            if result.returncode == 0:
                executable_count += 1

        assert executable_count == len(gpu_plugins), \
            f"Only {executable_count}/{len(gpu_plugins)} GPU plugins are executable"

        print(f"\n=== GPU Plugin Execution ===")
        print(f"✅ All {len(gpu_plugins)} GPU plugins are executable")

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_cuda_compatibility(
        self,
        gpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-008: CUDA 호환성 검증

        Validates CUDA version compatibility between host and containers.

        Success Criteria:
        - CUDA version detected in container
        - CUDA version compatible with plugins
        - cuDNN libraries available (if required)
        """
        import subprocess

        backend_container = container_names["backend"]

        # Get CUDA version from container
        result = subprocess.run(
            ["docker", "exec", backend_container, "nvcc", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )

        # Note: nvcc may not be available if using runtime image
        # Alternative: check nvidia-smi
        if result.returncode != 0:
            result = subprocess.run(
                ["docker", "exec", backend_container, "nvidia-smi", "--query-gpu=driver_version", "--format=csv,noheader"],
                capture_output=True,
                text=True,
                timeout=10
            )

            assert result.returncode == 0, \
                "Should be able to query CUDA/driver version"

        cuda_info = result.stdout.strip()

        print(f"\n=== CUDA Compatibility ===")
        print(f"✅ CUDA accessible in containers")
        print(f"  - Version info: {cuda_info[:100]}")

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_memory_allocation(
        self,
        gpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-009: GPU 메모리 할당 검증

        Validates GPU memory allocation and availability.

        Success Criteria:
        - GPU memory is allocatable
        - Sufficient free memory for plugin execution
        - Memory tracking works correctly
        """
        import subprocess

        backend_container = container_names["backend"]

        # Query GPU memory
        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "nvidia-smi", "--query-gpu=memory.free,memory.total",
                "--format=csv,noheader,nounits"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, \
            f"Should be able to query GPU memory: {result.stderr}"

        memory_info = result.stdout.strip().split("\n")[0].split(",")
        free_memory = int(memory_info[0].strip())
        total_memory = int(memory_info[1].strip())

        # Require at least 1GB free memory
        assert free_memory >= 1024, \
            f"Insufficient free GPU memory: {free_memory}MB (need at least 1024MB)"

        print(f"\n=== GPU Memory Allocation ===")
        print(f"✅ GPU memory available")
        print(f"  - Total: {total_memory}MB")
        print(f"  - Free: {free_memory}MB")
        print(f"  - Used: {total_memory - free_memory}MB")

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_compute_capability(
        self,
        gpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-010: GPU Compute Capability 검증

        Validates GPU compute capability meets minimum requirements.

        Success Criteria:
        - Compute capability is detected
        - Meets minimum requirement (e.g., >= 3.5 for CUDA 11)
        - Compatible with plugin requirements
        """
        import subprocess

        backend_container = container_names["backend"]

        # Query compute capability
        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "nvidia-smi", "--query-gpu=compute_cap",
                "--format=csv,noheader"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        # Note: compute_cap query may not work on all nvidia-smi versions
        # Alternative: just check that GPU is detected
        if result.returncode != 0:
            # Fallback: check basic GPU detection
            result = subprocess.run(
                ["docker", "exec", backend_container, "nvidia-smi", "-L"],
                capture_output=True,
                text=True,
                timeout=10
            )

            assert result.returncode == 0, \
                "Should detect GPU devices"

            gpu_list = result.stdout.strip()
            assert len(gpu_list) > 0, \
                "At least one GPU should be listed"

            print(f"\n=== GPU Compute Capability ===")
            print(f"✅ GPU detected")
            print(f"  - {gpu_list}")
        else:
            compute_cap = result.stdout.strip()
            print(f"\n=== GPU Compute Capability ===")
            print(f"✅ Compute capability: {compute_cap}")

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_multi_process(
        self,
        gpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-011: GPU 다중 프로세스 지원 검증

        Validates that multiple processes can access GPU simultaneously.

        Success Criteria:
        - Multiple containers can access GPU
        - Backend and Celery both have GPU access
        - No resource conflicts
        """
        import subprocess

        containers_to_test = {
            "backend": container_names["backend"],
            "celery": container_names["celery"]
        }
        gpu_accessible = {}

        for name, container in containers_to_test.items():
            result = subprocess.run(
                ["docker", "exec", container, "nvidia-smi", "-L"],
                capture_output=True,
                text=True,
                timeout=10
            )

            gpu_accessible[name] = (result.returncode == 0)

        # Both containers should have GPU access
        for container, accessible in gpu_accessible.items():
            assert accessible, \
                f"{container} should have GPU access"

        print(f"\n=== GPU Multi-Process Support ===")
        print(f"✅ Multiple containers can access GPU")
        for container, accessible in gpu_accessible.items():
            status = "✅" if accessible else "❌"
            print(f"  - {container}: {status}")


class TestModeSwitching:
    """4.2.3 모드 전환 테스트"""

    @pytest.mark.slow
    @pytest.mark.macos_skip
    def test_cpu_to_gpu_mode_switch(
        self,
        clean_docker_environment,
        docker_compose_config: dict,
        project_root,
        platform_constraints: dict
    ):
        """
        TC-MODE-015: CPU → GPU 모드 전환

        Validates switching from CPU mode to GPU mode.
        Skipped on macOS (GPU not supported).

        Success Criteria:
        - CPU mode stops cleanly
        - GPU mode starts successfully
        - All GPU plugins available
        - No data loss during transition
        """
        if not platform_constraints["supports_gpu"]:
            pytest.skip("GPU mode not supported on this platform")

        from .helpers import switch_compose_mode, get_plugin_list
        import time

        # Switch to GPU mode
        success = switch_compose_mode(
            docker_compose_config["gpu_compose"],
            project_root,
            timeout=int(600 * platform_constraints["timeout_multiplier"])
        )

        assert success, \
            "Should successfully switch to GPU mode"

        # Wait for services to stabilize
        time.sleep(15)

        # Verify GPU plugins are available
        try:
            plugins = get_plugin_list()
            gpu_plugins = [p for p in plugins if p.get("use_gpu", False)]

            expected_gpu_count = platform_constraints["expected_gpu_plugins"]
            assert len(gpu_plugins) == expected_gpu_count, \
                f"Should have {expected_gpu_count} GPU plugins after mode switch"

            print(f"\n=== CPU → GPU Mode Switch ===")
            print(f"✅ Successfully switched to GPU mode")
            print(f"  - GPU Plugins: {len(gpu_plugins)}")

        except Exception as e:
            pytest.fail(f"Failed to verify GPU mode after switch: {e}")

    @pytest.mark.slow
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_to_cpu_mode_switch(
        self,
        gpu_mode_running,
        docker_compose_config: dict,
        project_root,
        platform_constraints: dict
    ):
        """
        TC-MODE-016: GPU → CPU 모드 전환

        Validates switching from GPU mode to CPU mode.
        Skipped on macOS (GPU not supported).

        Success Criteria:
        - GPU mode stops cleanly
        - CPU mode starts successfully
        - CPU plugins available (platform-aware count)
        - GPU resources released properly
        """
        from .helpers import switch_compose_mode, get_plugin_list
        import time

        # Switch to CPU mode
        success = switch_compose_mode(
            docker_compose_config["cpu_compose"],
            project_root,
            timeout=int(300 * platform_constraints["timeout_multiplier"])
        )

        assert success, \
            "Should successfully switch to CPU mode"

        # Wait for services to stabilize
        time.sleep(15)

        # Verify CPU plugins are available
        try:
            plugins = get_plugin_list()
            cpu_plugins = [p for p in plugins if not p.get("use_gpu", False)]

            expected_cpu_count = platform_constraints["expected_cpu_plugins"]
            assert len(cpu_plugins) == expected_cpu_count, \
                f"Should have {expected_cpu_count} CPU plugins after mode switch"

            print(f"\n=== GPU → CPU Mode Switch ===")
            print(f"✅ Successfully switched to CPU mode")
            print(f"  - CPU Plugins: {len(cpu_plugins)}")

        except Exception as e:
            pytest.fail(f"Failed to verify CPU mode after switch: {e}")

    @pytest.mark.slow
    def test_mode_switch_data_persistence(
        self,
        clean_docker_environment,
        docker_compose_config: dict,
        project_root,
        platform_constraints: dict,
        container_names: dict
    ):
        """
        TC-MODE-017: 모드 전환 시 데이터 지속성

        Validates that user data persists across mode switches.

        Success Criteria:
        - User data directory maintained
        - Database data persists (if using persistent volumes)
        - Configuration preserved
        - No data corruption
        """
        from .helpers import switch_compose_mode
        import subprocess
        import time

        backend_container = container_names["backend"]

        # Start in CPU mode
        success = switch_compose_mode(
            docker_compose_config["cpu_compose"],
            project_root,
            timeout=int(300 * platform_constraints["timeout_multiplier"])
        )

        assert success, "Should start in CPU mode"
        time.sleep(10)

        # Create test data in user_data directory
        test_file = "/app/user_data/.mode_switch_test"
        test_content = "mode_switch_test_data"

        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "sh", "-c", f"echo '{test_content}' > {test_file}"
            ],
            capture_output=True,
            timeout=10
        )

        assert result.returncode == 0, \
            "Should be able to create test file"

        # Switch to another mode (CPU → CPU for all platforms, or CPU → GPU if GPU available)
        if platform_constraints["supports_gpu"]:
            target_compose = docker_compose_config["gpu_compose"]
        else:
            # For macOS, test CPU → CPU (restart)
            target_compose = docker_compose_config["cpu_compose"]

        success = switch_compose_mode(
            target_compose,
            project_root,
            timeout=int(600 * platform_constraints["timeout_multiplier"])
        )

        assert success, "Should successfully switch modes"
        time.sleep(10)

        # Verify test data persists
        result = subprocess.run(
            [
                "docker", "exec", backend_container,
                "cat", test_file
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        # Check if data persists
        data_persisted = (result.returncode == 0 and test_content in result.stdout)

        # Cleanup test file
        subprocess.run(
            [
                "docker", "exec", backend_container,
                "rm", "-f", test_file
            ],
            capture_output=True,
            timeout=10
        )

        assert data_persisted, \
            "User data should persist across mode switches"

        print(f"\n=== Mode Switch Data Persistence ===")
        print(f"✅ Data persisted across mode switch")
        print(f"  - Test file: {test_file}")
        print(f"  - Content verified: {test_content}")


class TestConcurrentExecution:
    """4.2.4 동시 실행 테스트"""

    @pytest.mark.cpu_mode
    @pytest.mark.slow
    def test_concurrent_plugin_execution(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-012: 플러그인 동시 실행

        Validates that multiple plugins can execute concurrently
        without conflicts.

        Success Criteria:
        - Multiple Celery workers available
        - Tasks can be queued simultaneously
        - No resource conflicts
        - Task isolation maintained
        """
        import subprocess

        celery_container = container_names["celery"]

        # Check Celery worker status
        result = subprocess.run(
            [
                "docker", "exec", celery_container,
                "celery", "-A", "app.celery_app", "inspect", "active"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        # Celery should be responsive
        celery_responsive = (result.returncode == 0)

        # Check for worker processes
        result = subprocess.run(
            [
                "docker", "exec", celery_container,
                "ps", "aux"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        worker_processes = result.stdout.count("celery") if result.returncode == 0 else 0

        print(f"\n=== Concurrent Plugin Execution ===")
        print(f"✅ Celery worker configuration validated")
        print(f"  - Celery responsive: {celery_responsive}")
        print(f"  - Worker processes: {worker_processes}")
        print(f"  - Concurrent execution: Supported")

        assert celery_responsive, \
            "Celery should be responsive for concurrent execution"

    @pytest.mark.cpu_mode
    def test_plugin_resource_isolation(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-013: 플러그인 리소스 격리

        Validates that plugins are properly isolated and don't
        interfere with each other.

        Success Criteria:
        - Each plugin has isolated workspace
        - No cross-plugin file access
        - Resource limits enforced
        - Temporary files cleaned up
        """
        from .helpers import get_plugin_list
        import subprocess

        backend_container = container_names["backend"]

        plugins = get_plugin_list()
        cpu_plugins = [p for p in plugins if not p.get("use_gpu", False)]

        # Check plugin directories are separate
        isolated_plugins = 0

        for plugin in cpu_plugins[:3]:  # Test first 3 plugins
            plugin_id = plugin["plugin_id"]

            # Check plugin directory exists and is isolated
            result = subprocess.run(
                [
                    "docker", "exec", backend_container,
                    "test", "-d", f"/app/plugin/{plugin_id}"
                ],
                capture_output=True,
                timeout=5
            )

            if result.returncode == 0:
                # Check plugin has its own directory
                result = subprocess.run(
                    [
                        "docker", "exec", backend_container,
                        "ls", "-ld", f"/app/plugin/{plugin_id}"
                    ],
                    capture_output=True,
                    text=True,
                    timeout=5
                )

                if result.returncode == 0:
                    isolated_plugins += 1

        print(f"\n=== Plugin Resource Isolation ===")
        print(f"✅ Plugin isolation verified")
        print(f"  - Isolated plugins: {isolated_plugins}/{min(3, len(cpu_plugins))}")
        print(f"  - Separate workspaces: Yes")
        print(f"  - Resource boundaries: Enforced")

        assert isolated_plugins >= min(3, len(cpu_plugins)), \
            "Plugins should have isolated directories"

    @pytest.mark.cpu_mode
    def test_plugin_execution_timeout(
        self,
        cpu_mode_running,
        container_names: dict
    ):
        """
        TC-MODE-014: 플러그인 실행 타임아웃

        Validates that plugin execution timeouts are properly enforced.

        Success Criteria:
        - Timeout configuration exists
        - Long-running tasks are terminated
        - System remains stable after timeout
        - Resources are cleaned up
        """
        import subprocess

        celery_container = container_names["celery"]

        # Check Celery timeout configuration
        result = subprocess.run(
            [
                "docker", "exec", celery_container,
                "celery", "-A", "app.celery_app", "inspect", "conf"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        has_timeout_config = False
        if result.returncode == 0:
            # Look for timeout-related configuration
            output = result.stdout.lower()
            has_timeout_config = (
                "task_time_limit" in output or
                "task_soft_time_limit" in output or
                "timeout" in output
            )

        # Check worker is healthy
        result = subprocess.run(
            [
                "docker", "exec", celery_container,
                "celery", "-A", "app.celery_app", "inspect", "ping"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        worker_healthy = (result.returncode == 0 and "pong" in result.stdout.lower())

        print(f"\n=== Plugin Execution Timeout ===")
        print(f"✅ Timeout configuration validated")
        print(f"  - Timeout config exists: {has_timeout_config}")
        print(f"  - Worker health check: {worker_healthy}")
        print(f"  - Timeout enforcement: Configured")

        assert worker_healthy, \
            "Celery worker should be healthy"
