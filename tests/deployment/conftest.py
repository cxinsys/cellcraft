"""
Shared pytest fixtures for CellCraft deployment tests.

This module provides common fixtures for:
- Docker Compose configuration
- Deployment script paths
- Container management
- Environment setup and cleanup
- CPU/GPU mode system initialization
- Platform detection and constraints
"""

import pytest
import subprocess
import time
import os
import platform
from pathlib import Path
from typing import Generator, Dict

# Import platform detection helpers
from .helpers import (
    detect_docker_engine,
    get_os_name,
    check_nvidia_gpu,
    check_wsl2_gpu,
    detect_wsl_filesystem,
    measure_volume_io_latency
)


# ==============================================================================
# Configuration Fixtures
# ==============================================================================

@pytest.fixture(scope="session")
def project_root() -> Path:
    """
    Get project root directory.

    Returns:
        Path: Absolute path to project root
    """
    return Path(__file__).parent.parent.parent


@pytest.fixture(scope="session")
def docker_compose_config() -> Dict[str, str]:
    """
    Docker Compose configuration information.

    Returns:
        dict: Compose file paths and project name
    """
    return {
        "cpu_compose": "docker-compose.cpu.yml",
        "gpu_compose": "docker-compose.gpu.yml",
        "cpu_local_compose": "docker-compose.cpu.local.yml",
        "gpu_local_compose": "docker-compose.gpu.local.yml",
        "project_name": "cellcraft"
    }


@pytest.fixture(scope="session")
def deployment_scripts(project_root: Path) -> Dict[str, str]:
    """
    Deployment script paths.

    Args:
        project_root: Project root directory

    Returns:
        dict: Script paths (absolute)
    """
    return {
        "cpu_script": str(project_root / "run-cpu-mode.sh"),
        "gpu_script": str(project_root / "run-gpu-mode.sh"),
        "check_script": str(project_root / "check-installation.sh")
    }


@pytest.fixture(scope="session")
def container_names() -> Dict[str, str]:
    """
    Container name mappings.

    Returns:
        dict: Service name to container name mapping
    """
    return {
        "frontend": "cellcraft-frontend-1",
        "backend": "cellcraft-backend-1",
        "db": "cellcraft-db-1",
        "rabbitmq": "cellcraft-rabbitmq-1",
        "celery": "cellcraft-celery-1"
    }


# ==============================================================================
# Platform Detection Fixtures
# ==============================================================================

@pytest.fixture(scope="session")
def platform_info() -> Dict[str, any]:
    """
    Detect platform information (OS, architecture, Docker engine).

    Returns:
        dict: Platform information with keys:
            - os_type: "Linux", "Darwin", or "Windows"
            - os_name: Human-readable OS name (e.g., "Ubuntu 22.04")
            - arch: "amd64" or "arm64"
            - docker_engine: Docker engine type
            - is_wsl2: Boolean indicating WSL2 environment
            - is_macos: Boolean indicating macOS
            - is_linux: Boolean indicating native Linux
            - is_windows: Boolean indicating Windows (WSL2)
    """
    os_type = platform.system()
    os_name = get_os_name()
    arch_raw = platform.machine()

    # Normalize architecture
    if arch_raw in ["x86_64", "AMD64"]:
        arch = "amd64"
    elif arch_raw in ["arm64", "aarch64"]:
        arch = "arm64"
    else:
        arch = arch_raw

    # Detect Docker engine
    docker_engine = detect_docker_engine()

    # Detect WSL2
    is_wsl2 = False
    if os_type == "Linux" and os.path.exists("/proc/version"):
        with open("/proc/version", "r") as f:
            if "microsoft" in f.read().lower():
                is_wsl2 = True

    return {
        "os_type": os_type,
        "os_name": os_name,
        "arch": arch,
        "docker_engine": docker_engine,
        "is_wsl2": is_wsl2,
        "is_macos": os_type == "Darwin",
        "is_linux": os_type == "Linux" and not is_wsl2,
        "is_windows": is_wsl2
    }


@pytest.fixture(scope="session")
def platform_constraints(platform_info: Dict[str, any]) -> Dict[str, any]:
    """
    Calculate platform-specific constraints and adjustments.

    Args:
        platform_info: Platform information from platform_info fixture

    Returns:
        dict: Platform constraints with keys:
            - timeout_multiplier: Timeout adjustment factor
            - supports_gpu: Boolean indicating GPU support
            - supports_custom_plugin: Boolean indicating Custom Plugin support
            - expected_cpu_plugins: Expected number of CPU plugins
            - expected_gpu_plugins: Expected number of GPU plugins
            - max_startup_time_cpu: Maximum CPU mode startup time (seconds)
            - max_startup_time_gpu: Maximum GPU mode startup time (seconds)
    """
    # Determine timeout multiplier
    if platform_info["is_macos"]:
        # macOS: Docker Desktop with virtualization overhead
        if platform_info["arch"] == "arm64":
            timeout_multiplier = 1.5  # Apple Silicon +50%
        else:
            timeout_multiplier = 1.3  # Intel +30%
    elif platform_info["is_wsl2"]:
        # WSL2: Docker Desktop with WSL2 overhead
        timeout_multiplier = 1.2  # +20%
    else:
        # Native Linux: baseline performance
        timeout_multiplier = 1.0

    # GPU support
    supports_gpu = False
    if platform_info["is_macos"]:
        # macOS: No NVIDIA GPU support
        supports_gpu = False
    elif platform_info["is_wsl2"]:
        # WSL2: Check for WSL2 GPU support
        supports_gpu = check_wsl2_gpu()
    else:
        # Native Linux: Check for NVIDIA GPU
        supports_gpu = check_nvidia_gpu()

    # Custom Plugin support (is_editable=true plugins)
    # macOS Docker Desktop cannot mount host source code properly
    supports_custom_plugin = not platform_info["is_macos"]

    # Expected plugin counts
    if platform_info["is_macos"]:
        # macOS: CPU plugins only, exclude Custom Plugin
        expected_cpu_plugins = 6  # 7 - 1 (Custom Plugin)
        expected_gpu_plugins = 0  # GPU not supported
    else:
        # Linux/WSL2: All plugins supported
        expected_cpu_plugins = 6  # CPU mode: GRNBOOST2, LEAP, TENET, GENIE3, GRNViz, Scribe
        expected_gpu_plugins = 8 if supports_gpu else 0

    # Maximum startup times (baseline * multiplier)
    max_startup_time_cpu = int(60 * timeout_multiplier)
    max_startup_time_gpu = int(90 * timeout_multiplier)

    return {
        "platform": platform_info["os_type"],
        "timeout_multiplier": timeout_multiplier,
        "supports_gpu": supports_gpu,
        "supports_custom_plugin": supports_custom_plugin,
        "expected_cpu_plugins": expected_cpu_plugins,
        "expected_gpu_plugins": expected_gpu_plugins,
        "max_startup_time_cpu": max_startup_time_cpu,
        "max_startup_time_gpu": max_startup_time_gpu
    }


@pytest.fixture(scope="session")
def docker_environment_info(platform_info: Dict[str, any]) -> Dict[str, any]:
    """
    Get Docker environment information.

    Args:
        platform_info: Platform information from platform_info fixture

    Returns:
        dict: Docker environment information with keys:
            - docker_version: Docker version string
            - compose_version: Docker Compose version string
            - volume_io_performance: I/O performance metrics
            - on_wsl_filesystem: Boolean (WSL2 only)
    """
    # Get Docker version
    try:
        result = subprocess.run(
            ["docker", "version", "--format", "{{.Server.Version}}"],
            capture_output=True,
            text=True,
            timeout=5
        )
        docker_version = result.stdout.strip() if result.returncode == 0 else "Unknown"
    except Exception:
        docker_version = "Unknown"

    # Get Docker Compose version
    try:
        result = subprocess.run(
            ["docker", "compose", "version", "--short"],
            capture_output=True,
            text=True,
            timeout=5
        )
        compose_version = result.stdout.strip() if result.returncode == 0 else "Unknown"
    except Exception:
        compose_version = "Unknown"

    # Measure volume I/O performance
    volume_io_performance = measure_volume_io_latency()

    # WSL2 filesystem detection
    on_wsl_filesystem = None
    if platform_info["is_wsl2"]:
        on_wsl_filesystem = detect_wsl_filesystem()

    return {
        "docker_version": docker_version,
        "compose_version": compose_version,
        "volume_io_performance": volume_io_performance,
        "on_wsl_filesystem": on_wsl_filesystem
    }


# ==============================================================================
# Environment Setup Fixtures
# ==============================================================================

@pytest.fixture(scope="session")
def clean_docker_environment(
    docker_compose_config: Dict[str, str],
    project_root: Path
) -> Generator[None, None, None]:
    """
    Clean Docker environment at session start and end (session-scoped).

    NOTE: Changed to session scope for performance. Individual tests should
    not rely on container restart for state isolation.

    Args:
        docker_compose_config: Compose file configuration
        project_root: Project root directory

    Yields:
        None
    """
    # Change to project root directory
    original_dir = os.getcwd()
    os.chdir(project_root)

    # Pre-session cleanup: Remove all CellCraft containers and volumes
    for compose_file in [
        docker_compose_config["cpu_compose"],
        docker_compose_config["gpu_compose"]
    ]:
        subprocess.run(
            ["docker", "compose", "-f", compose_file, "down", "-v"],
            capture_output=True,
            timeout=60
        )

    # Wait for complete cleanup
    time.sleep(2)

    yield

    # Post-session cleanup
    for compose_file in [
        docker_compose_config["cpu_compose"],
        docker_compose_config["gpu_compose"]
    ]:
        subprocess.run(
            ["docker", "compose", "-f", compose_file, "down", "-v"],
            capture_output=True,
            timeout=60
        )

    # Restore original directory
    os.chdir(original_dir)


@pytest.fixture(scope="session")
def cpu_mode_running(
    deployment_scripts: Dict[str, str],
    clean_docker_environment: None,
    project_root: Path,
    platform_constraints: Dict[str, any]
) -> Generator[None, None, None]:
    """
    Start CellCraft in CPU mode (session-scoped for performance).

    Executes run-cpu-mode.sh once per test session and waits for system initialization.
    Timeout is adjusted based on platform (Linux: 300s, macOS Intel: 390s,
    macOS M1: 450s, WSL2: 360s).

    NOTE: Session-scoped to avoid 5-10 minute startup per test.
    Tests should ensure proper cleanup of their own state.

    Args:
        deployment_scripts: Script paths
        clean_docker_environment: Environment cleanup fixture
        project_root: Project root directory
        platform_constraints: Platform-specific constraints

    Yields:
        None

    Raises:
        pytest.Failed: If CPU mode startup fails
    """
    # Change to project root directory
    original_dir = os.getcwd()
    os.chdir(project_root)

    # Calculate dynamic timeout based on platform
    base_timeout = 300  # 5 minutes baseline
    dynamic_timeout = int(base_timeout * platform_constraints["timeout_multiplier"])

    # Execute CPU mode script
    result = subprocess.run(
        [deployment_scripts["cpu_script"], "--skip-verify"],
        capture_output=True,
        text=True,
        timeout=dynamic_timeout
    )

    # Check execution result
    if result.returncode != 0:
        print(f"=== STDOUT ===")
        print(result.stdout)
        print(f"=== STDERR ===")
        print(result.stderr)
        pytest.fail(f"CPU mode startup failed with exit code {result.returncode}")

    # Additional stabilization wait
    time.sleep(10)

    yield

    # Cleanup handled by clean_docker_environment
    os.chdir(original_dir)


@pytest.fixture(scope="session")
def gpu_mode_running(
    deployment_scripts: Dict[str, str],
    clean_docker_environment: None,
    project_root: Path,
    platform_constraints: Dict[str, any]
) -> Generator[None, None, None]:
    """
    Start CellCraft in GPU mode (session-scoped for performance).

    Executes run-gpu-mode.sh once per test session and waits for system initialization.
    Requires NVIDIA GPU and drivers. Timeout is adjusted based on platform
    (Linux: 600s, WSL2: 720s).

    NOTE: Session-scoped to avoid 10+ minute startup per test.
    Tests should ensure proper cleanup of their own state.

    Args:
        deployment_scripts: Script paths
        clean_docker_environment: Environment cleanup fixture
        project_root: Project root directory
        platform_constraints: Platform-specific constraints

    Yields:
        None

    Raises:
        pytest.Failed: If GPU mode startup fails
        pytest.Skip: If GPU not supported on platform
    """
    # Skip if GPU not supported
    if not platform_constraints["supports_gpu"]:
        pytest.skip("GPU mode not supported on this platform")

    # Change to project root directory
    original_dir = os.getcwd()
    os.chdir(project_root)

    # Calculate dynamic timeout based on platform
    base_timeout = 600  # 10 minutes baseline
    dynamic_timeout = int(base_timeout * platform_constraints["timeout_multiplier"])

    # Execute GPU mode script
    result = subprocess.run(
        [deployment_scripts["gpu_script"], "--skip-verify"],
        capture_output=True,
        text=True,
        timeout=dynamic_timeout
    )

    # Check execution result
    if result.returncode != 0:
        print(f"=== STDOUT ===")
        print(result.stdout)
        print(f"=== STDERR ===")
        print(result.stderr)
        pytest.fail(f"GPU mode startup failed with exit code {result.returncode}")

    # Additional stabilization wait
    time.sleep(10)

    yield

    # Cleanup handled by clean_docker_environment
    os.chdir(original_dir)


# ==============================================================================
# Test State Isolation Fixtures
# ==============================================================================

@pytest.fixture(autouse=True, scope="function")
def isolate_test_state():
    """
    Ensure test state isolation without restarting containers.

    This fixture runs before and after each test to provide minimal state
    cleanup while keeping containers running (session-scoped).

    Currently a placeholder for future state cleanup mechanisms such as:
    - Database cleanup (if needed)
    - File system cleanup in /app/user_data
    - Celery queue purging

    Yields:
        None
    """
    # Pre-test: Could add cleanup logic here if needed
    yield
    # Post-test: Could add cleanup logic here if needed
    pass


# ==============================================================================
# Pytest Configuration
# ==============================================================================

def pytest_configure(config):
    """
    Register custom pytest markers.

    Args:
        config: Pytest configuration object
    """
    # Register markers for test categorization
    config.addinivalue_line(
        "markers",
        "cpu_mode: CPU mode deployment tests"
    )
    config.addinivalue_line(
        "markers",
        "gpu_mode: GPU mode deployment tests"
    )
    config.addinivalue_line(
        "markers",
        "slow: Slow running tests (> 30 seconds)"
    )
    config.addinivalue_line(
        "markers",
        "requires_network: Tests requiring internet connection"
    )
    config.addinivalue_line(
        "markers",
        "requires_gpu: Tests requiring NVIDIA GPU"
    )
