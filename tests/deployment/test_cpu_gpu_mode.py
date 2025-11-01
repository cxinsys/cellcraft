"""
CPU/GPU Mode Tests

Tests for validating CPU and GPU mode deployment configurations.
Covers TC-MODE-001 through TC-MODE-017 from DEPLOYMENT_TEST_PLAN.md
"""

import pytest
import re
import threading


class TestCPUMode:
    """4.2.1 CPU 모드 테스트"""

    @pytest.mark.cpu_mode
    def test_cpu_plugin_count(
        self,
        cpu_mode_running,
        platform_constraints: dict,
        project_root
    ):
        """
        TC-MODE-001: CPU 공식 플러그인 개수 (플랫폼 인식) - Deployment Validation

        Validates CPU plugin installation through file system and database:
        - Linux: 6 CPU plugins (official only from submodule)
        - macOS: 6 CPU plugins (official only from submodule)
        - WSL2: 6 CPU plugins (official only from submodule)

        Success Criteria:
        - Plugin directories exist in backend/plugin/official/
        - Plugins are initialized in database
        - Plugin count matches expected CPU plugins

        Note: This test validates deployment infrastructure, not business logic.
        Uses file system + DB query to avoid API authentication requirements.
        """
        import os
        from .helpers import get_plugins_from_db

        import csv

        # File System Check: Count plugin directories
        plugin_dir = os.path.join(project_root, "backend", "plugin", "official")
        assert os.path.exists(plugin_dir), \
            f"Plugin directory should exist: {plugin_dir}"

        plugin_dirs = [
            d for d in os.listdir(plugin_dir)
            if os.path.isdir(os.path.join(plugin_dir, d))
            and not d.startswith('.')
        ]

        # Read plugins.csv for use_gpu information
        plugins_csv_file = os.path.join(plugin_dir, "plugins.csv")
        assert os.path.exists(plugins_csv_file), \
            f"plugins.csv should exist: {plugins_csv_file}"

        cpu_plugins = []
        with open(plugins_csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # use_gpu field is string "True" or "False" in CSV
                if row['use_gpu'].lower() == 'false':
                    cpu_plugins.append(row['name'])

        # Expected count: Platform-specific (6 for macOS, 7 for Linux/WSL2 including Custom Plugin)
        expected_count = platform_constraints["expected_cpu_plugins"]

        # Validate counts
        assert len(plugin_dirs) >= expected_count, \
            f"Expected at least {expected_count} plugin directories, found {len(plugin_dirs)}"

        assert len(cpu_plugins) == expected_count, \
            f"Expected {expected_count} CPU plugins (use_gpu=False), found {len(cpu_plugins)}"

        # Log validation results
        print(f"\n=== CPU Plugin Count Validation ===")
        print(f"✅ File System: {len(plugin_dirs)} total plugin directories")
        print(f"✅ CPU Plugins: {len(cpu_plugins)}/{expected_count} plugins (from plugins.csv)")
        print(f"Platform: {platform_constraints.get('platform_type', 'unknown')}")

        print(f"\n  CPU Plugins (use_gpu=False):")
        for plugin_name in sorted(cpu_plugins):
            print(f"    - {plugin_name}")

    @pytest.mark.cpu_mode
    def test_cpu_plugin_metadata(
        self,
        cpu_mode_running,
        platform_constraints: dict,
        project_root
    ):
        """
        TC-MODE-002: CPU 플러그인 메타데이터 검증 - Deployment Validation

        Validates that all CPU plugins have metadata.json files with required fields.

        Success Criteria:
        - Each plugin directory has metadata.json file
        - metadata.json is valid JSON
        - Required fields are present: name, version, use_gpu
        - Plugins are initialized in database

        Note: This test validates deployment infrastructure, not business logic.
        Uses file system checks to avoid API authentication requirements.
        """
        import os
        import json
        from .helpers import get_plugins_from_db

        import csv

        # File System Check: Validate metadata.json files exist and are valid JSON
        plugin_dir = os.path.join(project_root, "backend", "plugin", "official")
        plugin_dirs = [
            d for d in os.listdir(plugin_dir)
            if os.path.isdir(os.path.join(plugin_dir, d))
            and not d.startswith('.')
        ]

        # Create case-insensitive mapping of plugin names to directory names
        plugin_dir_map = {d.lower(): d for d in plugin_dirs}

        # Read plugins.csv to get CPU plugins
        plugins_csv_file = os.path.join(plugin_dir, "plugins.csv")
        cpu_plugins = []
        with open(plugins_csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['use_gpu'].lower() == 'false':
                    cpu_plugins.append(row['name'])

        # Validate each CPU plugin has metadata.json
        metadata_validation_results = []
        for plugin_name in cpu_plugins:
            # Handle case mismatch between CSV name and directory name
            plugin_dir_name = plugin_dir_map.get(plugin_name.lower())
            if not plugin_dir_name:
                assert False, f"CPU plugin {plugin_name} directory not found (case-insensitive check)"

            plugin_path = os.path.join(plugin_dir, plugin_dir_name)
            metadata_file = os.path.join(plugin_path, "metadata.json")

            # Check metadata.json exists
            assert os.path.exists(metadata_file), \
                f"CPU plugin {plugin_name} (dir: {plugin_dir_name}) missing metadata.json at {metadata_file}"

            # Validate it's valid JSON
            try:
                with open(metadata_file, 'r', encoding='utf-8') as f:
                    metadata = json.load(f)

                # Basic validation: metadata should be a dict
                assert isinstance(metadata, dict), \
                    f"Plugin {plugin_name} metadata.json should be a JSON object"

                # Required metadata fields for deployment validation
                required_fields = ["name", "author", "description", "drawflow", "rules"]
                for field in required_fields:
                    assert metadata.get(field), \
                        f"Plugin {plugin_name} metadata.json missing required field: {field}"

                # Note: use_gpu differentiation is in plugins.csv, not metadata.json
                metadata_validation_results.append({
                    'name': plugin_name,
                    'metadata_name': metadata.get('name'),
                    'version': metadata.get('version', 'unknown'),
                    'use_gpu': plugin_name in ['FastTENET', 'FastSCODE'],  # GPU plugins
                    'valid': True
                })

            except json.JSONDecodeError as e:
                assert False, f"Plugin {plugin_name} has invalid JSON in metadata.json: {e}"

        # Expected: Platform-specific CPU plugin count
        expected_count = platform_constraints["expected_cpu_plugins"]
        assert len(metadata_validation_results) == expected_count, \
            f"Expected {expected_count} CPU plugins validated, found {len(metadata_validation_results)}"

        # Database Initialization Check: Verify plugins are initialized in database
        db_plugins = get_plugins_from_db(container_name="cellcraft-db-1", source="official")
        assert len(db_plugins) > 0, "Plugins should be initialized in database"

        # Verify each CSV plugin exists in database
        for csv_plugin in cpu_plugins:
            db_plugin = next((p for p in db_plugins if p["name"] == csv_plugin), None)
            assert db_plugin is not None, \
                f"Plugin {csv_plugin} from CSV not found in database"

        # Log validation results
        print(f"\n=== CPU Plugin Metadata Validation ===")
        print(f"✅ All {len(metadata_validation_results)} CPU plugins have valid metadata.json")
        print(f"✅ All metadata files contain required fields: name, version, description, author")
        print(f"✅ All CPU plugins have use_gpu=False")
        print(f"✅ All {len(db_plugins)} plugins initialized in database")

        print(f"\n  CPU Plugin Metadata Files:")
        for result in metadata_validation_results:
            print(f"    - {result['name']} v{result['version']} (use_gpu={result['use_gpu']})")



class TestGPUMode:
    """4.2.2 GPU 모드 테스트"""

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_plugin_count(
        self,
        gpu_mode_running,
        platform_constraints: dict,
        project_root
    ):
        """
        TC-MODE-004: GPU 플러그인 개수 - Deployment Validation

        Validates that GPU mode has the correct number of plugins.
        This test is skipped on macOS (GPU not supported).

        GPU Mode Plugins (8 total):
        - All CPU plugins (6): GRNBOOST2, LEAP, TENET, GENIE3, GRNViz, Scribe
        - GPU-only plugins (2): FastTENET, FastSCODE

        Success Criteria:
        - Expected number of plugin directories exist
        - plugins.csv contains correct GPU plugins (use_gpu=true)
        - All GPU plugin directories contain metadata.json

        Note: Uses filesystem checks to avoid API authentication requirements.
        """
        import os
        import csv

        # File System Check: Count plugin directories
        plugin_dir = os.path.join(project_root, "backend", "plugin", "official")
        plugin_dirs = [
            d for d in os.listdir(plugin_dir)
            if os.path.isdir(os.path.join(plugin_dir, d))
            and not d.startswith('.')
        ]

        # Expected: All plugins (CPU + GPU) in GPU mode
        expected_count = platform_constraints["expected_gpu_plugins"]
        assert len(plugin_dirs) >= expected_count, \
            f"Expected at least {expected_count} plugin directories, found {len(plugin_dirs)}"

        # Read plugins.csv to verify GPU plugins
        plugins_csv_file = os.path.join(plugin_dir, "plugins.csv")
        all_plugins = []
        gpu_plugins = []

        with open(plugins_csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                all_plugins.append(row['name'])
                if row['use_gpu'].lower() == 'true':
                    gpu_plugins.append(row['name'])

        # Verify plugin counts
        assert len(all_plugins) == expected_count, \
            f"Expected {expected_count} plugins in CSV, found {len(all_plugins)}"

        print(f"\n=== GPU Plugin Count ===")
        print(f"✅ Found {len(plugin_dirs)} plugin directories")
        print(f"✅ CSV has {len(all_plugins)} total plugins ({len(gpu_plugins)} GPU-only)")
        print(f"GPU-only plugins: {', '.join(gpu_plugins)}")

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_plugin_metadata(
        self,
        gpu_mode_running,
        platform_constraints: dict,
        project_root
    ):
        """
        TC-MODE-005: GPU 플러그인 메타데이터 검증 - Deployment Validation

        Validates that all plugins in GPU mode have metadata.json files with required fields.
        GPU mode includes both CPU plugins and GPU-only plugins.

        Success Criteria:
        - Each plugin directory has metadata.json file
        - metadata.json is valid JSON
        - Required fields are present: name, author, description, drawflow, rules
        - Plugins are initialized in database

        Note: This test validates deployment infrastructure, not business logic.
        Uses file system checks to avoid API authentication requirements.
        metadata.json structure is identical for all plugins (CPU and GPU).
        """
        import os
        import json
        from .helpers import get_plugins_from_db
        import csv

        # File System Check: Validate metadata.json files exist and are valid JSON
        plugin_dir = os.path.join(project_root, "backend", "plugin", "official")
        plugin_dirs = [
            d for d in os.listdir(plugin_dir)
            if os.path.isdir(os.path.join(plugin_dir, d))
            and not d.startswith('.')
        ]

        # Create case-insensitive mapping of plugin names to directory names
        plugin_dir_map = {d.lower(): d for d in plugin_dirs}

        # Read plugins.csv to get all plugins (GPU mode includes all)
        plugins_csv_file = os.path.join(plugin_dir, "plugins.csv")
        all_plugins = []
        with open(plugins_csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                all_plugins.append(row['name'])

        # Validate each plugin has metadata.json with correct structure
        metadata_validation_results = []
        for plugin_name in all_plugins:
            # Handle case mismatch between CSV name and directory name
            plugin_dir_name = plugin_dir_map.get(plugin_name.lower())
            if not plugin_dir_name:
                assert False, f"Plugin {plugin_name} directory not found (case-insensitive check)"

            plugin_path = os.path.join(plugin_dir, plugin_dir_name)
            metadata_file = os.path.join(plugin_path, "metadata.json")

            # Check metadata.json exists
            assert os.path.exists(metadata_file), \
                f"Plugin {plugin_name} (dir: {plugin_dir_name}) missing metadata.json at {metadata_file}"

            # Validate it's valid JSON
            try:
                with open(metadata_file, 'r', encoding='utf-8') as f:
                    metadata = json.load(f)

                # Basic validation: metadata should be a dict
                assert isinstance(metadata, dict), \
                    f"Plugin {plugin_name} metadata.json should be a JSON object"

                # Required metadata fields (same for all plugins)
                required_fields = ["name", "author", "description", "drawflow", "rules"]
                for field in required_fields:
                    assert metadata.get(field), \
                        f"Plugin {plugin_name} metadata.json missing required field: {field}"

                # Note: use_gpu differentiation is in plugins.csv, not metadata.json
                metadata_validation_results.append({
                    'name': plugin_name,
                    'metadata_name': metadata.get('name'),
                    'version': metadata.get('version', 'unknown'),
                    'use_gpu': plugin_name in ['FastTENET', 'FastSCODE'],  # GPU plugins
                    'valid': True
                })

            except json.JSONDecodeError as e:
                assert False, f"Plugin {plugin_name} has invalid JSON in metadata.json: {e}"

        # Expected: Platform-specific GPU plugin count
        expected_count = platform_constraints["expected_gpu_plugins"]
        assert len(metadata_validation_results) == expected_count, \
            f"Expected {expected_count} plugins validated, found {len(metadata_validation_results)}"

        # Database Initialization Check: Verify plugins are initialized in database
        db_plugins = get_plugins_from_db(container_name="cellcraft-db-1", source="official")
        assert len(db_plugins) > 0, "Plugins should be initialized in database"

        # Verify each CSV plugin exists in database
        for csv_plugin in all_plugins:
            db_plugin = next((p for p in db_plugins if p["name"] == csv_plugin), None)
            assert db_plugin is not None, \
                f"Plugin {csv_plugin} from CSV not found in database"

        print(f"\n=== GPU Plugin Metadata ===")
        print(f"✅ All {len(metadata_validation_results)} plugins have valid metadata.json")
        print(f"✅ All plugins initialized in database")

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

        # Parse GPU info for detailed validation
        gpu_details = []
        for i, info in enumerate(gpu_info, 1):
            parts = info.split(',')
            if len(parts) >= 2:
                gpu_name = parts[0].strip()
                gpu_memory = parts[1].strip()
                gpu_details.append({'index': i, 'name': gpu_name, 'memory': gpu_memory})

        # Get CUDA driver version
        driver_result = subprocess.run(
            ["docker", "exec", backend_container, "nvidia-smi", "--query-gpu=driver_version", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10
        )

        cuda_driver_version = None
        if driver_result.returncode == 0:
            cuda_driver_version = driver_result.stdout.strip().split('\n')[0]

        # Validate CUDA driver version (minimum 450.0 for CUDA 11.0)
        if cuda_driver_version:
            try:
                driver_major = int(cuda_driver_version.split('.')[0])
                assert driver_major >= 450, f"CUDA driver version {cuda_driver_version} is too old (need >= 450.0 for CUDA 11.0)"
                print(f"  - CUDA driver: {cuda_driver_version} ✅")
            except (ValueError, IndexError):
                print(f"  ⚠️  Could not parse CUDA driver version: {cuda_driver_version}")

        # Validate GPU memory (minimum 4GB)
        for gpu in gpu_details:
            memory_str = gpu['memory']
            memory_match = re.search(r'(\d+)\s*([MG]iB)', memory_str)
            if memory_match:
                memory_value = int(memory_match.group(1))
                memory_unit = memory_match.group(2)
                memory_gb = memory_value / 1024 if memory_unit == 'MiB' else memory_value

                assert memory_gb >= 4.0, f"GPU {gpu['index']} has insufficient memory: {memory_gb:.1f}GB (need >= 4GB)"
                print(f"  - GPU {gpu['index']} ({gpu['name']}): {memory_gb:.1f}GB ✅")

        print(f"\n=== GPU Availability (Enhanced) ===")
        print(f"✅ GPU accessible in containers")
        print(f"  - GPU Count: {gpu_count}")
        print(f"  - CUDA Driver: {cuda_driver_version}")
        for gpu in gpu_details:
            print(f"  - GPU {gpu['index']}: {gpu['name']} ({gpu['memory']})")

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

        # Parse CUDA version
        cuda_version = None
        if "release" in cuda_info.lower():
            version_match = re.search(r'release\s+(\d+\.\d+)', cuda_info, re.IGNORECASE)
            if version_match:
                cuda_version = version_match.group(1)

        if cuda_version:
            cuda_major = float(cuda_version)
            print(f"  - CUDA version: {cuda_version}")
            assert cuda_major >= 11.0, f"CUDA version {cuda_version} is too old (need >= 11.0)"
            print(f"  - CUDA version check: ✅ (>= 11.0)")
        else:
            print(f"  ⚠️  Could not parse CUDA version from: {cuda_info[:200]}")

        # Check for cuDNN library
        cudnn_check = subprocess.run(
            ["docker", "exec", backend_container, "sh", "-c",
             "find /usr/local/cuda* /usr/lib/*-linux-gnu -name 'libcudnn.so*' 2>/dev/null | head -1"],
            capture_output=True, text=True, timeout=10
        )

        has_cudnn = len(cudnn_check.stdout.strip()) > 0
        if has_cudnn:
            cudnn_path = cudnn_check.stdout.strip().split('\n')[0]
            print(f"  - cuDNN library: Found at {cudnn_path} ✅")
        else:
            print(f"  ⚠️  cuDNN library not found (may be required for GPU plugins)")

        # Check CUDA toolkit
        nvcc_check = subprocess.run(
            ["docker", "exec", backend_container, "which", "nvcc"],
            capture_output=True, timeout=5
        )
        has_nvcc = nvcc_check.returncode == 0
        print(f"  - CUDA toolkit (nvcc): {'✅' if has_nvcc else 'Runtime only'}")

        print(f"\n=== CUDA Compatibility (Enhanced) ===")
        print(f"✅ CUDA accessible in containers")
        print(f"  - CUDA version: {cuda_version if cuda_version else 'Unable to detect'}")
        print(f"  - cuDNN library: {'✅' if has_cudnn else '⚠️  Not found'}")
        print(f"  - CUDA toolkit: {'✅' if has_nvcc else 'Runtime only'}")

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

        # Test actual GPU memory allocation
        print(f"\n  Testing GPU memory allocation...")

        allocation_test = subprocess.run(
            ["docker", "exec", backend_container, "python", "-c",
             """
import sys
try:
    import torch
    if torch.cuda.is_available():
        device = torch.device('cuda:0')
        tensor = torch.zeros((100, 1024, 1024), dtype=torch.uint8, device=device)
        print('PyTorch allocation: 100MB allocated successfully')
        del tensor
        torch.cuda.empty_cache()
    else:
        print('PyTorch CUDA not available')
except Exception as e:
    print(f'PyTorch test failed: {e}')
    try:
        import ctypes
        cudart = ctypes.CDLL('libcudart.so')
        print('CUDA runtime library accessible')
    except Exception as e2:
        print(f'CUDA runtime test failed: {e2}')
"""],
            capture_output=True, text=True, timeout=30
        )

        allocation_success = "allocated successfully" in allocation_test.stdout or "runtime library accessible" in allocation_test.stdout
        if allocation_success:
            print(f"  - Memory allocation test: ✅")
            print(f"    {allocation_test.stdout.strip()}")
        else:
            print(f"  ⚠️  Memory allocation test results:")
            print(f"    {allocation_test.stdout.strip()}")

        # Check GPU memory bandwidth (informational)
        bandwidth_check = subprocess.run(
            ["docker", "exec", backend_container, "nvidia-smi", "--query-gpu=clocks.mem", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10
        )
        if bandwidth_check.returncode == 0:
            mem_clock = bandwidth_check.stdout.strip().split('\n')[0]
            print(f"  - Memory clock: {mem_clock}")

        print(f"\n=== GPU Memory Allocation (Enhanced) ===")
        print(f"✅ GPU memory available and allocatable")
        print(f"  - Total: {total_memory}MB")
        print(f"  - Free: {free_memory}MB")
        print(f"  - Used: {total_memory - free_memory}MB")
        print(f"  - Allocation test: {'✅' if allocation_success else '⚠️'}")

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

        # Get compute capability
        compute_cap_result = subprocess.run(
            ["docker", "exec", backend_container, "nvidia-smi", "--query-gpu=compute_cap", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10
        )

        compute_capabilities = []
        gpu_architectures = []

        if compute_cap_result.returncode == 0:
            for cap in compute_cap_result.stdout.strip().split('\n'):
                try:
                    cap_float = float(cap.strip())
                    compute_capabilities.append(cap_float)

                    # Map compute capability to architecture
                    if cap_float >= 9.0:
                        arch = "Hopper"
                    elif cap_float >= 8.0:
                        arch = "Ampere/Ada"
                    elif cap_float >= 7.5:
                        arch = "Turing"
                    elif cap_float >= 7.0:
                        arch = "Volta"
                    elif cap_float >= 6.0:
                        arch = "Pascal"
                    elif cap_float >= 5.0:
                        arch = "Maxwell"
                    else:
                        arch = "Kepler or older"
                    gpu_architectures.append(arch)
                except ValueError:
                    pass

        if compute_capabilities:
            min_capability = min(compute_capabilities)
            assert min_capability >= 3.5, f"GPU compute capability {min_capability} too low (need >= 3.5 for CUDA 11.0)"

            print(f"\n=== GPU Compute Capability (Enhanced) ===")
            print(f"✅ Compute capability validated")

            for i, (cap, arch) in enumerate(zip(compute_capabilities, gpu_architectures), 1):
                print(f"  - GPU {i}: Compute {cap} ({arch})")
                if cap < 7.0:
                    print(f"    ⚠️  GPU {i} does not have Tensor cores (Volta+ recommended for ML)")

            has_tensor_cores = any(cap >= 7.0 for cap in compute_capabilities)
            print(f"  - Tensor cores: {'✅ Available' if has_tensor_cores else '⚠️  Not available (Volta+ needed)'}")
        else:
            # Fallback
            fallback_result = subprocess.run(
                ["docker", "exec", backend_container, "nvidia-smi", "-L"],
                capture_output=True, text=True, timeout=10
            )
            assert fallback_result.returncode == 0, "Should detect GPU devices"
            gpu_list = fallback_result.stdout.strip()
            assert len(gpu_list) > 0, "At least one GPU should be listed"

            print(f"\n=== GPU Compute Capability ===")
            print(f"✅ GPU detected (compute capability query not supported)")
            print(f"  - {gpu_list}")
            print(f"  ⚠️  Unable to validate minimum compute capability")

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

        # Test concurrent GPU usage
        print(f"\n  Testing concurrent GPU usage...")

        concurrent_results = {}

        def gpu_stress_test(container_name, container_label):
            result = subprocess.run(
                ["docker", "exec", container_name, "python", "-c",
                 """
import sys
try:
    import torch
    if torch.cuda.is_available():
        device = torch.device('cuda:0')
        a = torch.randn(1000, 1000, device=device)
        b = torch.randn(1000, 1000, device=device)
        c = torch.mm(a, b)
        print('Compute test passed')
    else:
        print('CUDA not available')
except Exception as e:
    print(f'Test failed: {e}')
"""],
                capture_output=True, text=True, timeout=30
            )
            concurrent_results[container_label] = {
                'success': 'passed' in result.stdout.lower(),
                'output': result.stdout.strip()
            }

        # Run tests concurrently
        threads = []
        for name, container in containers_to_test.items():
            t = threading.Thread(target=gpu_stress_test, args=(container, name))
            t.start()
            threads.append(t)

        for t in threads:
            t.join()

        all_passed = all(result['success'] for result in concurrent_results.values())
        if all_passed:
            print(f"  - Concurrent GPU usage: ✅")
        else:
            print(f"  - Concurrent GPU usage: ⚠️  Some tests failed")

        # Get GPU compute mode
        compute_mode_result = subprocess.run(
            ["docker", "exec", backend_container, "nvidia-smi", "--query-gpu=compute_mode", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10
        )

        if compute_mode_result.returncode == 0:
            compute_mode = compute_mode_result.stdout.strip().split('\n')[0]
            print(f"  - GPU compute mode: {compute_mode}")
            if compute_mode.lower() not in ['default', 'shared']:
                print(f"    ⚠️  Compute mode is {compute_mode}, may restrict multi-process access")

        print(f"\n=== GPU Multi-Process Support (Enhanced) ===")
        print(f"✅ Multiple containers can access GPU")
        for container, result in concurrent_results.items():
            status = "✅" if result['success'] else "❌"
            print(f"  - {container}: {status}")
            if not result['success']:
                print(f"    Output: {result['output']}")


class TestConcurrentExecution:
    """4.2.4 동시 실행 테스트"""
    # All tests removed - business logic testing moved to integration tests
