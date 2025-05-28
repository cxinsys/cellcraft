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
    for attempt in range(max_retries):
        container.reload()
        if container.status != 'running':
            raise Exception(f"Container exited with status: {container.status}")
        
        # 더 신뢰할 수 있는 준비 상태 확인
        try:
            exit_code, output = container.exec_run(
                "test -f /opt/micromamba/envs/plugin_env/bin/python && echo 'ready'",
                demux=True
            )
            if exit_code == 0 and b'ready' in (output[0] or b''):
                return True
        except Exception as e:
            print(f"Container not ready yet: {e}")
        
        time.sleep(2)
    
    raise Exception("Container failed to become ready")

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
        
        # 더 안정적인 볼륨 마운트 설정
        host_workspace_path = os.environ.get("HOST_WORKSPACE_PATH", workspace_path)
        
        # 컨테이너 실행 설정
        container_config = {
            'image': image_name,
            'volumes': {
                host_workspace_path: {
                    "bind": "/workspace",
                    "mode": "rw"
                }
            },
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
        container = client.containers.run(**container_config)
        print(f"Container started with ID: {container.id}")
        
        # 컨테이너가 준비될 때까지 대기
        wait_for_container_ready(container)
        print("Container is ready for execution")
        
        # Snakemake 실행 명령어 (entrypoint가 환경을 설정하므로 단순화)
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
            "stderr": stderr
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
