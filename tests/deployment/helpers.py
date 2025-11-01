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


def wait_for_service_ready(url: str, timeout: int = 60, stability_checks: int = 3) -> bool:
    """
    Wait for service to become accessible with stability verification.

    Polls service every 2 seconds until accessible or timeout.
    Requires consecutive successful checks to ensure service stability.

    Args:
        url: Service URL to check
        timeout: Maximum wait time in seconds (default: 60)
        stability_checks: Number of consecutive successes required (default: 3)

    Returns:
        bool: True if service became accessible and stable, False if timeout
    """
    start_time = time.time()
    consecutive_successes = 0

    while time.time() - start_time < timeout:
        if check_service_accessible(url, timeout=5):
            consecutive_successes += 1
            if consecutive_successes >= stability_checks:
                return True  # Service is stable
        else:
            consecutive_successes = 0  # Reset on failure
        time.sleep(2)

    return False


def check_service_with_retry(url: str, retries: int = 3, retry_delay: int = 2) -> bool:
    """
    Check service accessibility with retry logic.

    Attempts multiple times to verify service is accessible,
    with delay between attempts.

    Args:
        url: Service URL to check
        retries: Number of retry attempts (default: 3)
        retry_delay: Delay between retries in seconds (default: 2)

    Returns:
        bool: True if service is accessible, False after all retries failed
    """
    for attempt in range(retries):
        if check_service_accessible(url, timeout=10):
            return True
        if attempt < retries - 1:  # Don't sleep after last attempt
            time.sleep(retry_delay)
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


def get_plugins_from_db(
    container_name: str = "cellcraft-db-1",
    source: str = "official"
) -> List[Dict[str, Any]]:
    """
    Query plugins table directly via Docker exec for deployment validation.

    This bypasses API authentication requirements and directly validates
    database initialization state.

    Args:
        container_name: Name of database container (default: cellcraft-db-1)
        source: Plugin source filter (default: "official")

    Returns:
        list: List of plugin dictionaries with fields:
            - id: Plugin ID
            - name: Plugin name
            - use_gpu: Boolean GPU requirement
            - source: Plugin source (official/custom)
            - created_at: Timestamp

    Raises:
        RuntimeError: If Docker exec fails or query returns error

    Example:
        >>> plugins = get_plugins_from_db()
        >>> cpu_plugins = [p for p in plugins if not p['use_gpu']]
        >>> len(cpu_plugins)  # Should be 6 or 7 depending on platform
    """
    import json

    # SQL query to get plugins with specified source
    query = f"""
    SELECT
        id,
        name,
        use_gpu,
        source,
        created_at,
        updated_at
    FROM plugins
    WHERE source = '{source}'
    ORDER BY name;
    """

    try:
        # Execute query via docker exec
        result = subprocess.run(
            [
                "docker", "exec", container_name,
                "psql", "-U", "cellcraft_admin", "-d", "cellcraft",
                "-t", "-A", "-F", ",",  # -t: tuples only, -A: unaligned, -F: field separator
                "-c", query
            ],
            capture_output=True,
            text=True,
            timeout=10,
            check=True
        )

        # Parse CSV output
        plugins = []
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            parts = line.split(',')
            if len(parts) >= 6:
                plugins.append({
                    'id': int(parts[0]),
                    'name': parts[1],
                    'use_gpu': parts[2].lower() == 't',  # PostgreSQL boolean: 't' or 'f'
                    'source': parts[3],
                    'created_at': parts[4],
                    'updated_at': parts[5]
                })

        return plugins

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Database query timed out on container '{container_name}'")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Failed to query plugins from database: {e.stderr}\n"
            f"Container: {container_name}, Query: {query}"
        )
    except Exception as e:
        raise RuntimeError(f"Unexpected error querying database: {e}")


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
    Restart a Docker container and wait for it to be healthy (or running if no health check).

    Args:
        container_name: Name of Docker container
        timeout: Maximum wait time in seconds (default: 60)

    Returns:
        bool: True if container restarted successfully and is healthy/running, False otherwise

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

        # Check if container has a health check
        status_info = get_container_status(container_name)
        if status_info and status_info["health"] != "none":
            # Container has health check, wait for it to be healthy
            return wait_for_container_healthy(container_name, timeout)
        else:
            # No health check, just wait for container to be running
            return wait_for_container_running(container_name, timeout)

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


# ==============================================================================
# Submodule Validation (Phase 2 - Submodule Tests)
# ==============================================================================

def get_submodule_status(submodule_path: str) -> Dict[str, Any]:
    """
    Get git submodule status including initialization state, branch, commit hash, and sync state.

    Args:
        submodule_path: Relative path to submodule from project root (e.g., "backend/plugin/official")

    Returns:
        dict: Submodule status information with keys:
            - initialized: bool - Whether submodule is initialized
            - commit_hash: str - Current commit hash (short)
            - commit_hash_full: str - Current commit hash (full)
            - branch: str - Current branch name (or None if detached HEAD)
            - remote_url: str - Remote repository URL
            - sync_status: str - "synced", "ahead", "behind", or "diverged"
            - has_uncommitted_changes: bool - Whether there are uncommitted changes
            - error: str - Error message if status check failed

    Example:
        >>> status = get_submodule_status("backend/plugin/official")
        >>> if status["initialized"]:
        ...     print(f"Branch: {status['branch']}, Commit: {status['commit_hash']}")
    """
    import subprocess
    import os

    try:
        status = {
            "initialized": False,
            "commit_hash": None,
            "commit_hash_full": None,
            "branch": None,
            "remote_url": None,
            "sync_status": "unknown",
            "has_uncommitted_changes": False,
            "error": None
        }

        # Check if submodule directory exists
        if not os.path.exists(submodule_path):
            status["error"] = f"Submodule path does not exist: {submodule_path}"
            return status

        # Check if directory is empty (not initialized)
        if not os.listdir(submodule_path):
            status["error"] = "Submodule directory is empty (not initialized)"
            return status

        status["initialized"] = True

        # Get current commit hash (short and full)
        result = subprocess.run(
            ["git", "-C", submodule_path, "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            status["commit_hash"] = result.stdout.strip()

        result = subprocess.run(
            ["git", "-C", submodule_path, "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            status["commit_hash_full"] = result.stdout.strip()

        # Get current branch (or detect detached HEAD)
        result = subprocess.run(
            ["git", "-C", submodule_path, "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            branch = result.stdout.strip()
            status["branch"] = None if branch == "HEAD" else branch

        # Get remote URL
        result = subprocess.run(
            ["git", "-C", submodule_path, "config", "--get", "remote.origin.url"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            status["remote_url"] = result.stdout.strip()

        # Check sync status with remote (if branch is known)
        if status["branch"]:
            # Fetch remote to get latest refs (quiet mode)
            subprocess.run(
                ["git", "-C", submodule_path, "fetch", "origin", "--quiet"],
                capture_output=True,
                timeout=30
            )

            # Check ahead/behind status
            result = subprocess.run(
                ["git", "-C", submodule_path, "rev-list", "--left-right", "--count",
                 f"HEAD...origin/{status['branch']}"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                counts = result.stdout.strip().split()
                if len(counts) == 2:
                    ahead, behind = int(counts[0]), int(counts[1])
                    if ahead == 0 and behind == 0:
                        status["sync_status"] = "synced"
                    elif ahead > 0 and behind == 0:
                        status["sync_status"] = "ahead"
                    elif ahead == 0 and behind > 0:
                        status["sync_status"] = "behind"
                    else:
                        status["sync_status"] = "diverged"

        # Check for uncommitted changes
        result = subprocess.run(
            ["git", "-C", submodule_path, "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            status["has_uncommitted_changes"] = len(result.stdout.strip()) > 0

        return status

    except subprocess.TimeoutExpired:
        return {
            **status,
            "error": "Git command timed out"
        }
    except Exception as e:
        return {
            **status,
            "error": f"Failed to get submodule status: {str(e)}"
        }


def get_submodule_branch(submodule_path: str) -> Optional[str]:
    """
    Get current branch name of submodule.

    Args:
        submodule_path: Relative path to submodule from project root

    Returns:
        str: Branch name, or None if detached HEAD or error

    Example:
        >>> branch = get_submodule_branch("backend/plugin/official")
        >>> print(f"Current branch: {branch}")
        Current branch: release/plugins-v1.0-cpu
    """
    import subprocess

    try:
        result = subprocess.run(
            ["git", "-C", submodule_path, "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5
        )

        if result.returncode == 0:
            branch = result.stdout.strip()
            return None if branch == "HEAD" else branch

        return None

    except Exception:
        return None


def validate_version_file(version_file_path: str) -> Dict[str, Any]:
    """
    Validate version.json file format and contents.

    Expected format:
    {
        "version": "1.0.0",
        "branch": "release/plugins-v1.0-cpu",
        "commit": "5e49d2b",
        "build_date": "2024-01-15",
        "plugins": ["TENET", "GENIE3", ...]
    }

    Args:
        version_file_path: Absolute path to version.json file

    Returns:
        dict: Validation result with keys:
            - valid: bool - Whether file is valid
            - version: str - Version string (if valid)
            - branch: str - Branch name (if valid)
            - commit: str - Commit hash (if valid)
            - plugins: list - Plugin list (if valid)
            - error: str - Error message if invalid

    Example:
        >>> result = validate_version_file("backend/plugin/official/version.json")
        >>> if result["valid"]:
        ...     print(f"Version: {result['version']}, Branch: {result['branch']}")
    """
    import json
    import os

    result = {
        "valid": False,
        "version": None,
        "branch": None,
        "commit": None,
        "plugins": [],
        "error": None
    }

    try:
        # Check if file exists
        if not os.path.exists(version_file_path):
            result["error"] = f"version.json file not found: {version_file_path}"
            return result

        # Read and parse JSON
        with open(version_file_path, 'r') as f:
            data = json.load(f)

        # Validate required fields
        required_fields = ["version", "branch", "commit"]
        missing_fields = [field for field in required_fields if field not in data]

        if missing_fields:
            result["error"] = f"Missing required fields: {', '.join(missing_fields)}"
            return result

        # Extract values
        result["version"] = data.get("version")
        result["branch"] = data.get("branch")
        result["commit"] = data.get("commit")
        result["plugins"] = data.get("plugins", [])
        result["valid"] = True

        return result

    except json.JSONDecodeError as e:
        result["error"] = f"Invalid JSON format: {str(e)}"
        return result
    except Exception as e:
        result["error"] = f"Failed to validate version file: {str(e)}"
        return result


def check_submodule_clean(submodule_path: str) -> bool:
    """
    Check if submodule has uncommitted changes (clean working directory).

    Args:
        submodule_path: Relative path to submodule from project root

    Returns:
        bool: True if working directory is clean, False if there are uncommitted changes

    Example:
        >>> if check_submodule_clean("backend/plugin/official"):
        ...     print("Submodule is clean")
        ... else:
        ...     print("Submodule has uncommitted changes")
    """
    import subprocess

    try:
        result = subprocess.run(
            ["git", "-C", submodule_path, "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=5
        )

        if result.returncode == 0:
            # Empty output means clean working directory
            return len(result.stdout.strip()) == 0

        return False

    except Exception:
        return False
