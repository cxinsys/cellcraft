"""
Helper utilities for CellCraft deployment tests.

This module provides reusable functions for:
- Script execution
- Container status monitoring
- Service connectivity checks
- API interactions
- Git operations
- Environment variable access
- Performance measurement
- Platform detection (Linux/macOS/Windows WSL2)
- Docker engine detection
- Volume I/O performance measurement
"""

import subprocess
import time
import requests
import json
import os
import platform
from typing import Optional, List, Dict, Any, Callable, Tuple


# ==============================================================================
# Script Execution
# ==============================================================================

def run_script(
    script_path: str,
    *args: str,
    timeout: int = 300
) -> subprocess.CompletedProcess:
    """
    Execute deployment script and return result.

    Args:
        script_path: Path to script file
        *args: Script arguments
        timeout: Maximum execution time in seconds (default: 300)

    Returns:
        subprocess.CompletedProcess: Execution result with stdout/stderr

    Raises:
        subprocess.TimeoutExpired: If script execution exceeds timeout
    """
    result = subprocess.run(
        [script_path, *args],
        capture_output=True,
        text=True,
        timeout=timeout
    )
    return result


# ==============================================================================
# Container Status
# ==============================================================================

def wait_for_container_healthy(
    container_name: str,
    timeout: int = 60
) -> bool:
    """
    Wait for container to reach healthy status.

    Polls container health status every 2 seconds until healthy or timeout.

    Args:
        container_name: Name of Docker container
        timeout: Maximum wait time in seconds (default: 60)

    Returns:
        bool: True if container reached healthy status, False if timeout
    """
    start_time = time.time()

    while time.time() - start_time < timeout:
        result = subprocess.run(
            ["docker", "inspect", "--format", "{{.State.Health.Status}}", container_name],
            capture_output=True,
            text=True
        )

        if result.returncode == 0 and "healthy" in result.stdout:
            return True

        time.sleep(2)

    return False


def wait_for_container_running(
    container_name: str,
    timeout: int = 60
) -> bool:
    """
    Wait for container to reach running status.

    Polls container status every 2 seconds until running or timeout.

    Args:
        container_name: Name of Docker container
        timeout: Maximum wait time in seconds (default: 60)

    Returns:
        bool: True if container reached running status, False if timeout
    """
    start_time = time.time()

    while time.time() - start_time < timeout:
        result = subprocess.run(
            ["docker", "inspect", "--format", "{{.State.Status}}", container_name],
            capture_output=True,
            text=True
        )

        if result.returncode == 0 and "running" in result.stdout:
            return True

        time.sleep(2)

    return False


def get_container_status(container_name: str) -> Optional[Dict[str, Any]]:
    """
    Get detailed container status information.

    Args:
        container_name: Name of Docker container

    Returns:
        dict: Container status information with keys:
            - status: running/exited/etc
            - health: healthy/unhealthy/none
            - image: Docker image name
            - started_at: Container start timestamp
        None: If container does not exist or inspection fails
    """
    result = subprocess.run(
        ["docker", "inspect", container_name],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        return None

    try:
        inspect_data = json.loads(result.stdout)
        if not inspect_data:
            return None

        container = inspect_data[0]
        return {
            "status": container["State"]["Status"],
            "health": container["State"].get("Health", {}).get("Status", "none"),
            "image": container["Config"]["Image"],
            "started_at": container["State"]["StartedAt"]
        }
    except (json.JSONDecodeError, KeyError, IndexError):
        return None


# ==============================================================================
# Service Connectivity
# ==============================================================================

def check_service_accessible(url: str, timeout: int = 10) -> bool:
    """
    Check if service is accessible via HTTP.

    Args:
        url: Service URL to check
        timeout: Request timeout in seconds (default: 10)

    Returns:
        bool: True if service returns HTTP 200, False otherwise
    """
    try:
        response = requests.get(url, timeout=timeout)
        return response.status_code == 200
    except Exception:
        return False


def wait_for_service_ready(url: str, timeout: int = 60) -> bool:
    """
    Wait for service to become accessible.

    Polls service every 2 seconds until accessible or timeout.

    Args:
        url: Service URL to check
        timeout: Maximum wait time in seconds (default: 60)

    Returns:
        bool: True if service became accessible, False if timeout
    """
    start_time = time.time()

    while time.time() - start_time < timeout:
        if check_service_accessible(url, timeout=5):
            return True
        time.sleep(2)

    return False


# ==============================================================================
# API Interactions
# ==============================================================================

def get_plugin_list() -> List[Dict[str, Any]]:
    """
    Get plugin list from backend API.

    Returns:
        list: List of plugin dictionaries with metadata

    Raises:
        RuntimeError: If API call fails or returns non-200 status
    """
    try:
        response = requests.get(
            "http://localhost:8000/api/plugin/list",
            timeout=10
        )
        response.raise_for_status()
        return response.json()
    except Exception as e:
        raise RuntimeError(f"Failed to get plugin list: {e}")


# ==============================================================================
# Git Operations
# ==============================================================================

def get_current_git_branch(repo_path: str) -> str:
    """
    Get current Git branch name.

    Args:
        repo_path: Path to Git repository

    Returns:
        str: Current branch name (e.g., "main", "release/v1.0")

    Raises:
        subprocess.CalledProcessError: If git command fails
    """
    result = subprocess.run(
        ["git", "-C", repo_path, "symbolic-ref", "--short", "HEAD"],
        capture_output=True,
        text=True,
        check=True
    )
    return result.stdout.strip()


# ==============================================================================
# Environment Variables
# ==============================================================================

def get_container_env(container_name: str, var_name: str) -> Optional[str]:
    """
    Get environment variable from running container.

    Args:
        container_name: Name of Docker container
        var_name: Environment variable name (e.g., "POSTGRES_DB")

    Returns:
        str: Environment variable value
        None: If container not found or variable not set
    """
    result = subprocess.run(
        ["docker", "exec", container_name, "env"],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        return None

    # Parse env output and find matching variable
    for line in result.stdout.splitlines():
        if line.startswith(f"{var_name}="):
            return line.split("=", 1)[1]

    return None


# ==============================================================================
# Performance Measurement
# ==============================================================================

def measure_execution_time(
    func: Callable,
    *args,
    **kwargs
) -> Tuple[Any, float]:
    """
    Measure function execution time.

    Args:
        func: Function to execute
        *args: Positional arguments for function
        **kwargs: Keyword arguments for function

    Returns:
        tuple: (function_result, execution_time_in_seconds)

    Example:
        >>> result, exec_time = measure_execution_time(my_func, arg1, arg2)
        >>> print(f"Executed in {exec_time:.2f} seconds")
    """
    start_time = time.time()
    result = func(*args, **kwargs)
    execution_time = time.time() - start_time
    return result, execution_time


# ==============================================================================
# Platform Detection
# ==============================================================================

def detect_docker_engine() -> str:
    """
    Detect Docker engine type.

    Returns:
        str: Docker engine type
            - "Docker Engine": Native Docker on Linux
            - "Docker Desktop": Docker Desktop on macOS/Windows
            - "Docker Desktop (WSL2 backend)": Docker Desktop with WSL2 on Windows
            - "Unknown": Unable to detect

    Example:
        >>> engine = detect_docker_engine()
        >>> if engine == "Docker Engine":
        ...     print("Running on native Docker")
    """
    try:
        result = subprocess.run(
            ["docker", "info", "--format", "{{.OperatingSystem}}"],
            capture_output=True,
            text=True,
            timeout=5
        )

        if result.returncode != 0:
            return "Unknown"

        os_info = result.stdout.strip().lower()

        # Check for Docker Desktop patterns
        if "docker desktop" in os_info:
            # Check if running in WSL2
            if os.path.exists("/proc/version"):
                with open("/proc/version", "r") as f:
                    version_info = f.read().lower()
                    if "microsoft" in version_info or "wsl" in version_info:
                        return "Docker Desktop (WSL2 backend)"
            return "Docker Desktop"

        # Native Docker Engine (Linux)
        return "Docker Engine"

    except Exception:
        return "Unknown"


def get_os_name() -> str:
    """
    Get human-readable OS name and version.

    Returns:
        str: OS name and version (e.g., "Ubuntu 22.04", "macOS 14.2", "Windows 11")

    Example:
        >>> os_name = get_os_name()
        >>> print(f"Running on {os_name}")
    """
    system = platform.system()

    if system == "Linux":
        # Check if WSL2
        if os.path.exists("/proc/version"):
            with open("/proc/version", "r") as f:
                version_info = f.read().lower()
                if "microsoft" in version_info or "wsl" in version_info:
                    # WSL2 - try to get Windows version
                    try:
                        result = subprocess.run(
                            ["powershell.exe", "-Command", "(Get-WmiObject Win32_OperatingSystem).Caption"],
                            capture_output=True,
                            text=True,
                            timeout=3
                        )
                        if result.returncode == 0:
                            return f"WSL2 ({result.stdout.strip()})"
                    except Exception:
                        pass
                    return "WSL2 (Windows)"

        # Regular Linux - read /etc/os-release
        if os.path.exists("/etc/os-release"):
            os_info = {}
            with open("/etc/os-release", "r") as f:
                for line in f:
                    if "=" in line:
                        key, value = line.strip().split("=", 1)
                        os_info[key] = value.strip('"')

            name = os_info.get("NAME", "Linux")
            version = os_info.get("VERSION_ID", "")
            return f"{name} {version}".strip()

        return "Linux"

    elif system == "Darwin":
        # macOS
        mac_version = platform.mac_ver()[0]
        # Determine processor architecture
        arch = platform.machine()
        if arch == "arm64":
            return f"macOS {mac_version} (Apple Silicon)"
        else:
            return f"macOS {mac_version} (Intel)"

    elif system == "Windows":
        # Windows
        win_version = platform.win32_ver()[0]
        return f"Windows {win_version}"

    else:
        return f"{system} {platform.release()}"


def check_nvidia_gpu() -> bool:
    """
    Check if NVIDIA GPU is available.

    Returns:
        bool: True if nvidia-smi command succeeds, False otherwise

    Example:
        >>> if check_nvidia_gpu():
        ...     print("NVIDIA GPU detected")
    """
    try:
        result = subprocess.run(
            ["nvidia-smi"],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except Exception:
        return False


def check_wsl2_gpu() -> bool:
    """
    Check if WSL2 GPU support is available.

    WSL2 GPU requires:
    - /dev/dxg device file (DirectX Graphics)
    - nvidia-smi working

    Returns:
        bool: True if WSL2 GPU support is available, False otherwise

    Example:
        >>> if check_wsl2_gpu():
        ...     print("WSL2 GPU support available")
    """
    # Check if /dev/dxg exists (WSL2 GPU device)
    if not os.path.exists("/dev/dxg"):
        return False

    # Check if nvidia-smi works
    return check_nvidia_gpu()


def detect_wsl_filesystem() -> bool:
    """
    Detect if current directory is on WSL filesystem or Windows filesystem.

    WSL filesystem: /home/user/... (optimal performance)
    Windows filesystem: /mnt/c/Users/... (10-30% slower I/O)

    Returns:
        bool: True if on WSL filesystem, False if on Windows filesystem (/mnt/*)

    Raises:
        Warning: Prints warning if on Windows filesystem

    Example:
        >>> on_wsl_fs = detect_wsl_filesystem()
        >>> if not on_wsl_fs:
        ...     print("WARNING: Using Windows filesystem - expect slower performance")
    """
    current_dir = os.getcwd()

    # Check if on /mnt/ (Windows filesystem)
    on_windows_fs = current_dir.startswith("/mnt/")

    if on_windows_fs:
        print("⚠️  WARNING: Current directory is on Windows filesystem (/mnt/)")
        print("   Docker volume performance may be 10-30% slower.")
        print("   For optimal performance, move project to WSL filesystem (e.g., ~/cellcraft)")

    return not on_windows_fs


def measure_volume_io_latency() -> Dict[str, Any]:
    """
    Measure Docker volume I/O performance.

    Creates a test Docker volume, measures write and read latency,
    and compares against platform-specific thresholds.

    Returns:
        dict: I/O performance metrics with keys:
            - read_latency_ms: Read latency in milliseconds
            - write_latency_ms: Write latency in milliseconds
            - within_threshold: Boolean indicating if performance is acceptable
            - platform: Platform type (Linux/macOS/WSL2)
            - threshold_ms: Expected threshold for platform

    Platform Thresholds:
        - Linux: <10ms (native Docker)
        - macOS: <20ms (Docker Desktop)
        - WSL2 (WSL FS): <15ms
        - WSL2 (/mnt/c/): <50ms

    Example:
        >>> perf = measure_volume_io_latency()
        >>> if not perf['within_threshold']:
        ...     print(f"WARNING: Slow I/O detected ({perf['write_latency_ms']}ms)")
    """
    # Determine platform and threshold
    system = platform.system()
    is_wsl2 = False
    on_wsl_fs = True

    if system == "Linux" and os.path.exists("/proc/version"):
        with open("/proc/version", "r") as f:
            if "microsoft" in f.read().lower():
                is_wsl2 = True
                on_wsl_fs = not os.getcwd().startswith("/mnt/")

    # Set threshold
    if system == "Darwin":
        threshold_ms = 20.0
        platform_name = "macOS"
    elif is_wsl2:
        if on_wsl_fs:
            threshold_ms = 15.0
            platform_name = "WSL2 (WSL FS)"
        else:
            threshold_ms = 50.0
            platform_name = "WSL2 (/mnt/c/)"
    else:
        threshold_ms = 10.0
        platform_name = "Linux"

    # Create test volume
    volume_name = "cellcraft_io_test"

    try:
        # Clean up any existing test volume
        subprocess.run(
            ["docker", "volume", "rm", volume_name],
            capture_output=True,
            timeout=5
        )

        # Create test volume
        result = subprocess.run(
            ["docker", "volume", "create", volume_name],
            capture_output=True,
            text=True,
            timeout=5,
            check=True
        )

        # Measure write latency
        write_start = time.time()
        subprocess.run(
            [
                "docker", "run", "--rm",
                "-v", f"{volume_name}:/data",
                "alpine:latest",
                "sh", "-c", "dd if=/dev/zero of=/data/test bs=1M count=10 2>/dev/null"
            ],
            capture_output=True,
            timeout=10,
            check=True
        )
        write_latency_ms = (time.time() - write_start) * 1000

        # Measure read latency
        read_start = time.time()
        subprocess.run(
            [
                "docker", "run", "--rm",
                "-v", f"{volume_name}:/data",
                "alpine:latest",
                "sh", "-c", "cat /data/test > /dev/null"
            ],
            capture_output=True,
            timeout=10,
            check=True
        )
        read_latency_ms = (time.time() - read_start) * 1000

        # Clean up test volume
        subprocess.run(
            ["docker", "volume", "rm", volume_name],
            capture_output=True,
            timeout=5
        )

        # Check if within threshold
        within_threshold = (write_latency_ms <= threshold_ms and
                           read_latency_ms <= threshold_ms)

        return {
            "read_latency_ms": round(read_latency_ms, 2),
            "write_latency_ms": round(write_latency_ms, 2),
            "within_threshold": within_threshold,
            "platform": platform_name,
            "threshold_ms": threshold_ms
        }

    except Exception as e:
        # Clean up on error
        subprocess.run(
            ["docker", "volume", "rm", volume_name],
            capture_output=True
        )

        return {
            "read_latency_ms": -1,
            "write_latency_ms": -1,
            "within_threshold": False,
            "platform": platform_name,
            "threshold_ms": threshold_ms,
            "error": str(e)
        }


# ==============================================================================
# Container Restart and Mode Switching (Phase 2C)
# ==============================================================================

def restart_container(container_name: str, timeout: int = 60) -> bool:
    """
    Restart a Docker container and wait for it to be healthy.

    Args:
        container_name: Name of Docker container
        timeout: Maximum wait time in seconds (default: 60)

    Returns:
        bool: True if container restarted successfully and is healthy, False otherwise

    Example:
        >>> if restart_container("cellcraft-backend-1"):
        ...     print("Backend restarted successfully")
    """
    try:
        # Restart the container
        result = subprocess.run(
            ["docker", "restart", container_name],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode != 0:
            return False

        # Wait for container to be healthy
        return wait_for_container_healthy(container_name, timeout)

    except Exception:
        return False


def stop_all_containers(project_name: str = "cellcraft") -> bool:
    """
    Stop all containers for a Docker Compose project.

    Args:
        project_name: Docker Compose project name (default: "cellcraft")

    Returns:
        bool: True if all containers stopped successfully, False otherwise

    Example:
        >>> if stop_all_containers():
        ...     print("All containers stopped")
    """
    try:
        result = subprocess.run(
            ["docker", "compose", "down"],
            capture_output=True,
            text=True,
            timeout=120
        )

        return result.returncode == 0

    except Exception:
        return False


def switch_compose_mode(
    compose_file: str,
    project_root: str,
    timeout: int = 300
) -> bool:
    """
    Switch Docker Compose mode by bringing down current stack and starting new one.

    Args:
        compose_file: Docker Compose file to use (e.g., "docker-compose.cpu.yml")
        project_root: Path to project root directory
        timeout: Maximum wait time in seconds (default: 300)

    Returns:
        bool: True if mode switch succeeded, False otherwise

    Example:
        >>> if switch_compose_mode("docker-compose.gpu.yml", "/path/to/project"):
        ...     print("Switched to GPU mode")
    """
    import os

    try:
        original_dir = os.getcwd()
        os.chdir(project_root)

        # Stop current stack
        subprocess.run(
            ["docker", "compose", "down", "-v"],
            capture_output=True,
            timeout=120
        )

        # Wait for cleanup
        time.sleep(5)

        # Start new stack
        result = subprocess.run(
            ["docker", "compose", "-f", compose_file, "up", "-d"],
            capture_output=True,
            text=True,
            timeout=timeout
        )

        os.chdir(original_dir)

        return result.returncode == 0

    except Exception:
        try:
            os.chdir(original_dir)
        except Exception:
            pass
        return False


def inject_container_failure(
    container_name: str,
    failure_type: str = "stop"
) -> bool:
    """
    Inject a failure into a container for testing error handling.

    Args:
        container_name: Name of Docker container
        failure_type: Type of failure to inject:
            - "stop": Stop the container
            - "kill": Kill the container forcefully
            - "pause": Pause the container
            - "network": Disconnect container from network

    Returns:
        bool: True if failure injected successfully, False otherwise

    Example:
        >>> inject_container_failure("cellcraft-backend-1", "stop")
        >>> # Test error handling
        >>> restart_container("cellcraft-backend-1")
    """
    try:
        if failure_type == "stop":
            result = subprocess.run(
                ["docker", "stop", container_name],
                capture_output=True,
                timeout=30
            )
        elif failure_type == "kill":
            result = subprocess.run(
                ["docker", "kill", container_name],
                capture_output=True,
                timeout=10
            )
        elif failure_type == "pause":
            result = subprocess.run(
                ["docker", "pause", container_name],
                capture_output=True,
                timeout=10
            )
        elif failure_type == "network":
            result = subprocess.run(
                ["docker", "network", "disconnect", "cellcraft_default", container_name],
                capture_output=True,
                timeout=10
            )
        else:
            return False

        return result.returncode == 0

    except Exception:
        return False


def recover_from_failure(
    container_name: str,
    failure_type: str = "stop",
    timeout: int = 60
) -> bool:
    """
    Recover a container from an injected failure.

    Args:
        container_name: Name of Docker container
        failure_type: Type of failure to recover from (must match inject_container_failure)
        timeout: Maximum wait time in seconds (default: 60)

    Returns:
        bool: True if container recovered successfully, False otherwise

    Example:
        >>> inject_container_failure("cellcraft-backend-1", "stop")
        >>> recover_from_failure("cellcraft-backend-1", "stop")
    """
    try:
        if failure_type in ["stop", "kill"]:
            result = subprocess.run(
                ["docker", "start", container_name],
                capture_output=True,
                timeout=30
            )
        elif failure_type == "pause":
            result = subprocess.run(
                ["docker", "unpause", container_name],
                capture_output=True,
                timeout=10
            )
        elif failure_type == "network":
            result = subprocess.run(
                ["docker", "network", "connect", "cellcraft_default", container_name],
                capture_output=True,
                timeout=10
            )
        else:
            return False

        if result.returncode != 0:
            return False

        # Wait for container to be healthy (if applicable)
        if failure_type in ["stop", "kill"]:
            return wait_for_container_healthy(container_name, timeout)

        return True

    except Exception:
        return False
