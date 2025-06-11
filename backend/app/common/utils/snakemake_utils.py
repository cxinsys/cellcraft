import os
import docker
import yaml
from pathlib import Path
from typing import Dict, Any
import time

def get_log_path(snakefile_path: str) -> Path:
    """일관된 로그 파일 경로 반환"""
    snakefile_dir = Path(snakefile_path).parent
    log_dir = snakefile_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
    return log_dir / "run.log"

def wait_for_container_ready(container, max_retries=10):
    """컨테이너가 준비될 때까지 대기"""
    print(f"컨테이너 준비 대기 중 (최대 {max_retries}회 시도)...")
    
    for attempt in range(max_retries):
        try:
            container.reload()
            if container.status != 'running':
                print(f"컨테이너 상태: {container.status}")
                raise Exception(f"컨테이너가 종료됨: {container.status}")
            
            # 기본적인 명령 실행 가능 여부 확인
            result = container.exec_run("echo 'test'")
            if result.exit_code == 0:
                print(f"컨테이너 응답 확인 (시도 {attempt + 1})")
                
                # Python 환경 상세 확인
                python_path = "/opt/micromamba/envs/plugin_env/bin/python"
                python_checks = [
                    f"test -f {python_path} && echo 'Python executable exists' || echo 'Python executable not found'",
                    f"ls -l {python_path}",
                    f"{python_path} --version",
                    "ls -la /opt/micromamba/envs/plugin_env/bin/",
                    "echo '=== Python Environment Info ==='",
                    "which python",
                    "python --version",
                    "echo '=== Micromamba Environment Info ==='",
                    "micromamba env list",
                    "echo '=== Environment Variables ==='",
                    "env | grep -E '(PATH|MAMBA|PYTHON)'",
                    "echo '=== Plugin Environment Dependencies ==='",
                    "micromamba list -n plugin_env",
                    "echo '=== Pip Installed Packages ==='",
                    f"{python_path} -m pip list"
                ]
                
                for check in python_checks:
                    check_result = container.exec_run(check)
                    print(f"Check '{check}':\n{check_result.output.decode().strip()}")
                
                # Python 실행 파일 존재 여부 확인
                python_exists = container.exec_run(f"test -f {python_path}").exit_code == 0
                if not python_exists:
                    print(f"Python 실행 파일이 존재하지 않음: {python_path}")
                    continue
                
                # Python 버전 확인
                python_version = container.exec_run(f"{python_path} --version")
                if python_version.exit_code == 0:
                    print(f"Python 환경 확인됨: {python_version.output.decode().strip()}")
                    return True
                else:
                    print(f"Python 환경 확인 실패 (시도 {attempt + 1}): {python_version.output.decode().strip()}")
                
        except Exception as e:
            print(f"컨테이너 확인 실패 (시도 {attempt + 1}): {str(e)}")
        
        time.sleep(2)
    
    # 실패 시 디버그 정보 수집
    try:
        print("\n=== 컨테이너 디버그 정보 ===")
        # 환경 변수 확인
        env_result = container.exec_run("env | grep -E '(PATH|MAMBA|PYTHON)'")
        print(f"환경 변수:\n{env_result.output.decode()}")
        
        # 파일 시스템 확인
        ls_result = container.exec_run("ls -la /opt/micromamba/envs/")
        print(f"Micromamba 환경 디렉토리:\n{ls_result.output.decode()}")
        
        # Entrypoint 확인
        entrypoint_result = container.exec_run("cat /entrypoint.sh")
        print(f"Entrypoint 스크립트:\n{entrypoint_result.output.decode()}")
    except Exception as debug_error:
        print(f"디버그 정보 수집 실패: {debug_error}")
    
    raise Exception("모든 시도 후에도 컨테이너가 준비되지 않음")

def exec_in_plugin(plugin_name: str, snakefile_path: str, targets: list, workspace_path: str) -> Dict[str, Any]:
    """
    플러그인 컨테이너에서 Snakemake 작업을 실행하고 로그를 파일로 저장합니다.
    
    Parameters:
        plugin_name (str): 플러그인 이름
        snakefile_path (str): Snakefile 경로
        targets (list): Snakemake 타겟 목록
        workspace_path (str): 작업 공간 경로
    
    Returns:
        Dict[str, Any]: 실행 결과
    """
    container = None
    try:
        client = docker.from_env()
        
        # 도커 이미지 존재 여부 확인
        image_name = f'plugin-{plugin_name.lower()}'
        try:
            client.images.get(image_name)
        except docker.errors.ImageNotFound:
            return {
                "returncode": 1,
                "stdout": "",
                "stderr": f"Plugin image '{image_name}' not found. Please build the plugin first."
            }
        
        # 작업 디렉토리 설정
        workspace_path = Path(os.path.abspath(workspace_path))
        snakefile_dir = Path(os.path.dirname(snakefile_path))
        log_file = get_log_path(snakefile_path)
        
        # Snakefile 경로를 workspace 기준으로 수정
        container_snakefile_path = f"/workspace{snakefile_path[1:]}"
        
        # Celery 컨테이너 ID 가져오기
        celery_container_id = os.environ.get("HOSTNAME")
        if not celery_container_id:
            print("WARNING: HOSTNAME env var not found. Using default celery container name.")
            celery_container_id = "cellcraft_celery_1"  # docker-compose 서비스 이름
        
        print(f"Using celery container ID: {celery_container_id}")
        
        # Celery 컨테이너의 마운트 정보 조회
        try:
            celery_container = client.containers.get(celery_container_id)
            celery_mounts = celery_container.attrs.get("Mounts", [])
            
            # /app 마운트의 호스트 경로 찾기
            host_backend_path = None
            for mount in celery_mounts:
                if mount["Destination"] == "/app":
                    host_backend_path = mount["Source"]
                    break
            
            if not host_backend_path:
                print("WARNING: Could not find host path for /app mount in Celery container")
                host_backend_path = workspace_path  # 기본값으로 workspace_path 사용
        except Exception as e:
            print(f"WARNING: Error getting Celery container mounts: {e}")
            host_backend_path = workspace_path  # 기본값으로 workspace_path 사용
        
        print(f"Using host backend path: {host_backend_path}")
        
        # 컨테이너 실행 설정
        container_config = {
            'image': image_name,
            'volumes': {
                host_backend_path: {
                    "bind": "/workspace",
                    "mode": "rw"
                }
            },
            'volumes_from': [celery_container_id],  # 기존 volumes_from 유지
            'working_dir': '/workspace',
            'environment': {
                'MAMBA_ROOT_PREFIX': '/opt/micromamba',
                'PATH': '/opt/micromamba/envs/plugin_env/bin:/opt/micromamba/bin:/usr/local/bin:/usr/bin:/bin',
                'PYTHONUNBUFFERED': '1'
            },
            'user': f"{os.getuid()}:{os.getgid()}",  # 권한 문제 해결
            'detach': True,
            'network': 'cellcraft_app-network',
            'tty': True,
            'stdin_open': True
        }
        
        # Task마다 새로운 컨테이너 생성
        print(f"Starting new container with image: {image_name}")
        print(f"Using volumes from celery container: {celery_container_id}")
        container = client.containers.run(**container_config)
        print(f"Container started with ID: {container.id}")
        
        # 컨테이너가 준비될 때까지 대기
        wait_for_container_ready(container)
        print("Container is ready for execution")
        
        # Snakemake 실행 명령어
        bash_cmd = (
            f"cd /workspace && "
            f"snakemake {' '.join(targets)} "
            f"--snakefile {container_snakefile_path} "
            "-j 1 --printshellcmds "
            f"2>&1 | tee {log_file}"
        )
        
        print("Executing command in container...")
        # 컨테이너에서 명령어 실행
        exec_id = client.api.exec_create(
            container.id,
            ["bash", "-c", bash_cmd],
            workdir="/workspace"
        )["Id"]
        
        # 실행 결과 수집
        for chunk in client.api.exec_start(exec_id, stream=True):
            print(f"Container output: {chunk.decode().strip()}")
            
        exit_code = client.api.exec_inspect(exec_id)["ExitCode"]
        print(f"Container execution completed with exit code: {exit_code}")
        
        # 메타데이터 저장
        meta_file = snakefile_dir / "meta.yml"
        meta_file.write_text(
            yaml.safe_dump({
                "exit_code": exit_code,
                "plugin_name": plugin_name,
                "targets": targets,
                "snakefile_path": snakefile_path
            })
        )
        
        # 로그 파일 읽기
        stdout = ""
        stderr = ""
        if log_file.exists():
            with open(log_file, 'r') as f:
                log_content = f.read()
                stdout = log_content
                if exit_code != 0:
                    stderr = log_content
                print(f"Log file content length: {len(log_content)}")
        else:
            print("Warning: Log file not found")
            stderr = "Log file not created - Snakemake may not have run properly"
        
        return {
            "returncode": exit_code,
            "stdout": stdout,
            "stderr": stderr,
            "log_path": str(log_file)
        }
        
    except Exception as e:
        # 더 상세한 에러 정보 수집
        if container:
            try:
                logs = container.logs(tail=100).decode()
                return {
                    "returncode": 1,
                    "stdout": "",
                    "stderr": f"Error: {str(e)}\nContainer logs:\n{logs}"
                }
            except:
                pass
        raise
    finally:
        # 확실한 컨테이너 정리
        if container:
            try:
                container.stop(timeout=10)
                container.remove(force=True)
            except Exception as cleanup_error:
                print(f"Container cleanup error: {cleanup_error}")

def snakemakeProcess(targets, snakefile_path, plugin_name):
    """기존 snakemakeProcess 함수는 유지하되, 내부에서 run_plugin_container를 호출하도록 수정"""
    print("Targets:", targets)
    print("Snakefile path:", snakefile_path)
    
    # 호스트 시스템의 backend 폴더 절대 경로 지정
    workspace_path = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))))
    print(f"Workspace path: {workspace_path}")
    
    # Docker 컨테이너로 실행
    return exec_in_plugin(plugin_name, snakefile_path, targets, workspace_path)

def read_stream(stream, output):
    for line in iter(stream.readline, b''):
        output.append(line.decode())

def is_snakemake_installed() -> bool:
    import subprocess
    try:
        # Try running `snakemake --version`
        subprocess.run(["snakemake", "--version"], check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        return False
    
def list_conda_libraries():
    import subprocess
    try:
        result = subprocess.run(["conda", "list"], check=True, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error fetching installed libraries: {e}")
        return None
    
def list_pip_libraries():
    import subprocess
    try:
        result = subprocess.run(["pip", "list"], check=True, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error fetching installed libraries: {e}")
        return None

def get_conda_environment_name():
    import sys
    print(sys.executable)
    stat_info = os.stat('/app')
    owner_uid = stat_info.st_uid
    owner_gid = stat_info.st_gid

    print(f"Owner UID: {owner_uid}")
    print(f"Owner GID: {owner_gid}")

    return os.environ.get('CONDA_DEFAULT_ENV', None)

def filter_and_add_suffix(input_string):
    # input_string에 밑줄("_")이 포함되어 있는지 확인
    if "_" in input_string:
        # "_"로 구분된 문자열을 배열로 변환
        segments = input_string.split("_")
        # 마지막 두 요소를 제외한 나머지를 합침
        file_name = "_".join(segments[:-2]) + ".h5ad"
        return file_name
    # 밑줄이 없는 경우, 원래 문자열을 반환
    return input_string

def change_snakefile_parameter(snakefile_path, output_path, user_input, plugin_params):
    # 기존 Snakefile 읽기
    with open(snakefile_path, 'r') as file:
        content = file.read()

    # user_input의 값을 문자열로 변환해서 처리
    for key, value in user_input.items():
        if isinstance(value, str):
            content = content.replace(f"{{{key}}}", value)

    # user_input의 값을 문자열로 변환해서 처리
    for key, value in plugin_params.items():
        # 리스트인 경우, []로 감싸고 쉼표로 구분된 문자열로 변환
        if isinstance(value, list):
            value = "[" + ", ".join(f"\"{v}\"" if isinstance(v, str) else str(v) for v in value) + "]"
        # 딕셔너리인 경우, {}로 감싸고 키:값 쌍을 쉼표로 구분된 문자열로 변환
        elif isinstance(value, dict):
            value = '{' + ', '.join(f"{k}: \"{v}\"" if isinstance(v, str) else f"{k}: {v}" for k, v in value.items()) + '}'
        # 숫자형(int, float)은 그대로 사용
        elif isinstance(value, (int, float)):
            value = str(value)
        # 문자열인 경우, " "로 감싸줌
        elif isinstance(value, str):
            value = f"\"{value}\""

        # 변환된 값을 content에 반영
        content = content.replace(f"{{{key}}}", value)

    # 새로운 Snakefile 생성
    with open(output_path, 'w') as file:
        file.write(content)

    # 생성된 파일 경로를 반환
    return output_path
