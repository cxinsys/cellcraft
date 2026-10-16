"""Plugin Docker build/tag operations & build orchestration.

Originally the router-facing orchestration lived here (``prepare_build_context``
/ ``dispatch_build_task``, PR-5). PR-9 folds in the Docker image build/tag/check
logic that used to sit in ``plugin/utils.py``:

- platform detection (``get_target_platform`` / ``get_optimal_platform`` / ...)
- ``build_plugin_docker_image`` / ``check_plugin_docker_image``
- ``generate_plugin_image_uri`` and ``setup_plugin_environments_via_docker``

Dockerfile *content* generation stays in ``plugin/runtime.py``
(``generate_plugin_dockerfile``); this module owns the imperative Docker calls.
Function signatures and behavior are unchanged from the original ``utils.py``.
"""
import os
import json
import platform
from datetime import datetime
from typing import Optional

import docker
from fastapi import HTTPException

from app.plugin.files import get_plugin_path
from app.plugin.runtime import generate_plugin_dockerfile
from app.worker.tasks import build_plugin_task

# Build log directory (Docker build logs)
BUILD_LOGS_DIR = "./plugin/build_logs"


def detect_target_platform():
    """
    시스템 아키텍처를 자동 감지하여 Docker 플랫폼 문자열 반환
    
    Returns:
        str: Docker 플랫폼 문자열 (예: "linux/amd64", "linux/arm64")
    """
    machine = platform.machine().lower()
    system = platform.system().lower()
    
    # 아키텍처 매핑 테이블
    arch_mapping = {
        'x86_64': 'amd64',
        'amd64': 'amd64',
        'aarch64': 'arm64',
        'arm64': 'arm64',
        'armv7l': 'arm/v7',
        'armv6l': 'arm/v6'
    }
    
    # 시스템이 Linux가 아닌 경우 기본값 반환
    if system != 'linux':
        system = 'linux'  # Docker 컨테이너는 Linux 기반
    
    arch = arch_mapping.get(machine, 'amd64')  # 알 수 없는 아키텍처는 amd64 기본값
    return f"{system}/{arch}"

def detect_platform_via_docker(client):
    """
    Docker 데몬에서 플랫폼 정보 추출
    
    Args:
        client: Docker 클라이언트 객체
        
    Returns:
        str: Docker 플랫폼 문자열, 실패 시 "linux/amd64"
    """
    try:
        info = client.info()
        architecture = info.get('Architecture', 'x86_64')
        os_type = info.get('OSType', 'linux').lower()
        
        arch_map = {
            'x86_64': 'amd64',
            'aarch64': 'arm64',
            'arm64': 'arm64'
        }
        
        docker_arch = arch_map.get(architecture, 'amd64')
        return f"{os_type}/{docker_arch}"
    except Exception:
        return "linux/amd64"  # fallback

def get_target_platform(client=None):
    """
    우선순위에 따른 플랫폼 감지
    
    우선순위:
    1. 환경변수 DOCKER_DEFAULT_PLATFORM
    2. Docker 데몬 정보 
    3. Python 시스템 정보
    
    Args:
        client: Docker 클라이언트 객체 (선택사항)
        
    Returns:
        str: 감지된 Docker 플랫폼 문자열
    """
    # 1순위: 환경변수 직접 지정
    if env_platform := os.environ.get('DOCKER_DEFAULT_PLATFORM'):
        return env_platform
    
    # 2순위: Docker 데몬 정보
    if client:
        docker_platform = detect_platform_via_docker(client)
        if docker_platform != "linux/amd64":  # 기본값이 아닌 경우 우선 사용
            return docker_platform
    
    # 3순위: Python 시스템 정보
    return detect_target_platform()

def check_gpu_requirement(plugin_path: str) -> bool:
    """
    플러그인의 GPU 요구사항 확인
    
    Args:
        plugin_path (str): 플러그인 폴더 경로
        
    Returns:
        bool: GPU 필요 여부
    """
    try:
        metadata_path = os.path.join(plugin_path, "metadata.json")
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
                return metadata.get('requiresGPU', False)
    except Exception:
        pass
    
    # metadata에서 확인 불가능한 경우, 기본값 False (CPU 전용)
    return False

def get_optimal_platform(plugin_path: str, client=None):
    """
    GPU 사용 여부에 따른 최적 플랫폼 선택
    
    Args:
        plugin_path (str): 플러그인 폴더 경로
        client: Docker 클라이언트 객체 (선택사항)
        
    Returns:
        str: 최적화된 Docker 플랫폼 문자열
    """
    use_gpu = check_gpu_requirement(plugin_path)
    
    if use_gpu:
        # GPU는 현재 linux/amd64만 지원 (NVIDIA CUDA)
        return "linux/amd64"
    
    # CPU 전용은 자동 감지
    return get_target_platform(client)


def ensure_build_logs_dir():
    """
    Ensure the build logs directory exists.
    """
    if not os.path.exists(BUILD_LOGS_DIR):
        os.makedirs(BUILD_LOGS_DIR, exist_ok=True)
        print(f"Created build logs directory: {BUILD_LOGS_DIR}")


def setup_plugin_environments_via_docker(plugin_name: str):
    """
    플러그인의 Python과 R 가상환경을 설치합니다.
    """
    try:
        client = docker.DockerClient(base_url="unix://var/run/docker.sock")
        celery_container = next(
            (c for c in client.containers.list() if 'celery' in c.name),
            None
        )
        if not celery_container:
            raise RuntimeError("Celery container not found")

        plugin_path = f"/app/plugin/{plugin_name}"
        dep_path = f"{plugin_path}/dependency"
        env_path = f"{plugin_path}/env"
        py_env = f"{env_path}/py_env"
        r_env = f"{env_path}/r_env"

        def run(cmd, error_msg, workdir=None):
            result = celery_container.exec_run(cmd, workdir=workdir, user="root")
            if result.exit_code != 0:
                raise RuntimeError(f"{error_msg}: {result.output.decode()}")
            return result

        run(f"mkdir -p {env_path} {r_env}", "Failed to create environment directories")
        run(f"chmod 777 {env_path} {r_env}", "Failed to set permissions")

        plugin_path, _ = get_plugin_path(plugin_name)
        plugin_dep_path = os.path.join(plugin_path, "dependency")
        
        if os.path.exists(os.path.join(plugin_dep_path, "environment.yml")):
            run(f"micromamba create -y -p {py_env} -f {dep_path}/environment.yml", "Failed to create Python env")
        elif os.path.exists(os.path.join(plugin_dep_path, "requirements.txt")):
            run(f"micromamba create -y -p {py_env} python=3.10", "Failed to create base Python env")
            run(f"micromamba run -p {py_env} pip install -r {dep_path}/requirements.txt", "Failed to install pip packages")

        if os.path.exists(os.path.join(plugin_dep_path, "renv.lock")):
            run(f"cp {plugin_dep_path}/renv.lock {r_env}/", "Failed to copy renv.lock")
            run(f"cp {plugin_dep_path}/*.tar.gz {r_env}/ || true", "Failed to copy tar.gz files")

            r_commands = [
                "R -e 'if (!requireNamespace(\"renv\", quietly=TRUE)) install.packages(\"renv\", repos=\"https://cloud.r-project.org\")'",
                "R -e 'renv::init(bare=TRUE)'",
                "R -e 'tryCatch({ renv::restore(lockfile=\"renv.lock\", prompt=FALSE) }, error=function(e) { message(\"Warning: partial restore\") })'",
                "R -e 'writeLines(c(\"Sys.setenv(RENV_PATHS_LIBRARY=\\\"renv_library\\\")\", \"source(\\\"renv/activate.R\\\")\"), \".Rprofile\")'"
            ]
            for cmd in r_commands:
                run(cmd, "Failed to setup R environment", workdir=r_env)

        run(f"chmod -R 777 {env_path}", "Failed to finalize permissions")

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to setup plugin environment: {str(e)}")


def build_plugin_docker_image(plugin_path: str, plugin_name: str) -> dict:
    """
    플러그인의 Docker 이미지를 빌드하고 로그를 저장합니다.

    Parameters:
        plugin_path (str): 플러그인 폴더 경로
        plugin_name (str): 플러그인 이름

    Returns:
        dict: {
            'success': bool,
            'message': str,
            'log_file': str,
            'image_tag': str
        }
    """
    client = None  # Docker 클라이언트 누수 방지를 위해 초기화
    try:
        client = docker.from_env()
        client.ping()  # Docker 데몬 연결 확인
    except Exception as e:
        error_msg = f"Failed to connect to Docker daemon: {str(e)}"
        return {
            'success': False,
            'message': error_msg,
            'log_file': None,
            'image_tag': None
        }

    # 로그 디렉토리 및 파일 설정
    ensure_build_logs_dir()
    log_file = os.path.join(BUILD_LOGS_DIR, f"{plugin_name.lower()}.log")
    image_tag = f"plugin-{plugin_name.lower()}"

    # 플랫폼 자동 감지
    target_platform = get_optimal_platform(plugin_path, client)
    build_platform = get_target_platform(client)
    use_gpu = check_gpu_requirement(plugin_path)

    try:
        with open(log_file, "w", encoding="utf-8") as f:
            f.write(f"Building Docker image for plugin: {plugin_name}\n")
            f.write(f"Start time: {datetime.now().isoformat()}\n")
            f.write(f"Docker version: {client.version().get('Version', 'unknown')}\n")
            f.write(f"Host platform: {platform.machine()} ({platform.system()})\n")
            f.write(f"Target platform: {target_platform}\n")
            f.write(f"Build platform: {build_platform}\n")
            f.write(f"GPU mode: {use_gpu}\n\n")
            f.flush()

        # 기존 이미지 제거 (있다면)
        try:
            client.images.remove(image=image_tag, force=True)
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"Removed existing image: {image_tag}\n\n")
        except docker.errors.ImageNotFound:
            pass

        # Dockerfile 존재 확인
        dockerfile_path = os.path.join(plugin_path, "Dockerfile")
        if not os.path.exists(dockerfile_path):
            error_msg = f"Dockerfile not found in {plugin_path}"
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"Error: {error_msg}\n")
            return {
                'success': False,
                'message': error_msg,
                'log_file': log_file,
                'image_tag': image_tag
            }

        # Dockerfile 내용 기록
        with open(log_file, "a", encoding="utf-8") as f:
            f.write("Dockerfile content:\n")
            with open(dockerfile_path, "r", encoding="utf-8") as df:
                f.write(df.read())
            f.write("\n\nStarting build...\n")
            f.flush()

        image = None
        build_output = []

        # 이미지 빌드
        try:
            # client.images.build()는 BuildKit을 올바르게 지원
            image, build_logs = client.images.build(
                path=plugin_path,
                tag=image_tag,
                rm=True,
                forcerm=True,
                buildargs={
                    'BUILDKIT_PROGRESS': 'plain',
                    'DOCKER_BUILDKIT': '1',  # BuildKit 명시적 활성화
                    'TARGETPLATFORM': target_platform,  # 자동 감지된 타겟 플랫폼
                    'BUILDPLATFORM': build_platform     # 자동 감지된 빌드 플랫폼
                },
                platform=None,  # 자동 플랫폼 감지
                pull=False,
                nocache=False
            )
            
            # 빌드 로그 처리
            for log_line in build_logs:
                with open(log_file, "a", encoding="utf-8") as f:
                    if 'stream' in log_line:
                        f.write(log_line['stream'])
                        build_output.append(log_line['stream'])
                    elif 'error' in log_line:
                        f.write(f"Build error: {log_line['error']}\n")
                        raise Exception(log_line['error'])
                    f.flush()
            
        except Exception as e:
            # Fallback: legacy API 사용 (BuildKit 없이)
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"Modern API failed, trying legacy API: {str(e)}\n")
            
            for chunk in client.api.build(
                path=plugin_path,
                tag=image_tag,
                rm=True,
                forcerm=True,
                decode=True,
                buildargs={
                    'BUILDKIT_PROGRESS': 'plain',
                    'DOCKER_BUILDKIT': '0',  # Legacy API에서는 BuildKit 비활성화
                    'TARGETPLATFORM': target_platform,  # 자동 감지된 타겟 플랫폼
                    'BUILDPLATFORM': build_platform     # 자동 감지된 빌드 플랫폼
                }
            ):
                with open(log_file, "a", encoding="utf-8") as f:
                    if 'stream' in chunk:
                        f.write(chunk['stream'])
                        build_output.append(chunk['stream'])
                    elif 'error' in chunk:
                        f.write(f"Build error: {chunk['error']}\n")
                        raise Exception(chunk['error'])
                    elif 'aux' in chunk and 'ID' in chunk['aux']:
                        image = client.images.get(chunk['aux']['ID'])
                    f.flush()

        if not image:
            image = client.images.get(image_tag)

        # 빌드 성공 기록
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(f"\nBuild completed successfully at {datetime.now().isoformat()}\n")
            if image:
                f.write(f"Image ID: {image.id}\n")
                f.write(f"Size: {image.attrs['Size'] / 1024 / 1024:.2f} MB\n")
                f.write(f"Created: {image.attrs['Created']}\n")
            f.flush()

        return {
            'success': True,
            'message': f"Successfully built image: {image_tag}",
            'log_file': log_file,
            'image_tag': image_tag
        }

    except Exception as e:
        error_msg = f"Unexpected error during Docker build: {str(e)}"
        try:
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"\nCritical Error: {error_msg}\n")
        except:
            pass
        return {
            'success': False,
            'message': error_msg,
            'log_file': log_file,
            'image_tag': image_tag
        }
    finally:
        # Docker 클라이언트 연결 해제 - TCP 소켓 누수 방지
        if client:
            try:
                client.close()
            except Exception:
                pass


def check_plugin_docker_image(plugin_name: str) -> bool:
    """
    플러그인의 Docker 이미지가 존재하는지 확인합니다.

    Parameters:
        plugin_name (str): 플러그인 이름

    Returns:
        bool: Docker 이미지 존재 여부
    """
    client = None  # Docker 클라이언트 누수 방지를 위해 초기화
    try:
        client = docker.from_env()
        image_tag = f"plugin-{plugin_name.lower()}"

        # 이미지 존재 여부 확인
        try:
            client.images.get(image_tag)
            return True
        except docker.errors.ImageNotFound:
            return False
        except Exception as e:
            print(f"Error checking Docker image: {str(e)}")
            return False

    except Exception as e:
        print(f"Failed to connect to Docker daemon: {str(e)}")
        return False
    finally:
        # Docker 클라이언트 연결 해제 - TCP 소켓 누수 방지
        if client:
            try:
                client.close()
            except Exception:
                pass

def generate_plugin_image_uri(plugin_name: str, source: str, version: Optional[str] = None) -> str:
    """
    Generate plugin image URI based on plugin type.
    
    Args:
        plugin_name: Name of the plugin
        source: Plugin source ("official" or "local")  
        version: Version for official plugins (optional)
        
    Returns:
        Plugin image URI string
    """
    if source == "official":
        # For official plugins, use GitHub Container Registry
        version = version or "latest"
        return f"ghcr.io/cxinsys/cellcraft-{plugin_name.lower()}:{version}"
    else:
        # For local plugins, use timestamp-based identifier
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        return f"local-build:{plugin_name.lower()}-{timestamp}"


def prepare_build_context(*, plugin_folder: str, script_folder: str, use_gpu: bool) -> str:
    """Ensure the scripts folder is build-ready and generate the Dockerfile.

    Mirrors the inline logic from ``build_plugin_docker``:
    - create the scripts folder if missing
    - drop a ``.gitkeep`` dummy file when the scripts folder is empty (Docker
      build needs a non-empty directory)
    - generate the Dockerfile at ``<plugin_folder>/Dockerfile``

    Returns the generated Dockerfile path.
    """
    # scripts 폴더 존재 확인 및 더미 파일 생성 (Docker 빌드를 위해 필요)
    print(f"Checking scripts folder before Docker build: {script_folder}")
    print(f"Scripts folder exists: {os.path.exists(script_folder)}")

    if not os.path.exists(script_folder):
        os.makedirs(script_folder)
        print(f"Created empty scripts folder at {script_folder}")

    # scripts 폴더가 비어있다면 더미 파일 생성
    scripts_contents = os.listdir(script_folder) if os.path.exists(script_folder) else []
    print(f"Scripts folder contents: {scripts_contents}")

    if not scripts_contents:
        dummy_file_path = os.path.join(script_folder, ".gitkeep")
        with open(dummy_file_path, 'w') as f:
            f.write("# This file ensures the scripts directory is not empty\n")
        print(f"Created dummy file at {dummy_file_path}")
        print(f"Scripts folder contents after dummy file: {os.listdir(script_folder)}")

    # Dockerfile 생성
    dockerfile_path = os.path.join(plugin_folder, "Dockerfile")
    generate_plugin_dockerfile(plugin_folder, dockerfile_path, use_gpu=use_gpu)
    print(f"Generated Dockerfile at: {dockerfile_path} (GPU: {use_gpu})")

    return dockerfile_path


def dispatch_build_task(*, plugin_name: str, user_id: int):
    """Dispatch the asynchronous plugin Docker build Celery task.

    Wire format (kwargs keys) is preserved exactly as the router used it.
    Returns the Celery ``AsyncResult``.
    """
    return build_plugin_task.apply_async(
        args=[],
        kwargs={
            'plugin_name': plugin_name,
            'user_id': user_id,
            'workflow_id': None,
            'algorithm_id': None,
            'task_type': "plugin_build"
        },
        ignore_result=False
    )
