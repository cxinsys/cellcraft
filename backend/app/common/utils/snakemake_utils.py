import os
from subprocess import Popen, PIPE

def snakemakeProcess(targets, snakefile_path):
    print("Targets:", targets)
    print("Snakefile path:", snakefile_path)
    print("Environment PATH:", os.environ["PATH"])

    # Snakemake 실행 명령어
    command = [
        '/opt/conda/envs/snakemake/bin/snakemake',
        *targets,  
        '--snakefile', snakefile_path,
        '-j',
    ]

    # Snakemake 프로세스 실행
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    print("STDOUT:", stdout.decode())
    print("STDERR:", stderr.decode())

    return {
        "returncode": process.returncode,
        "stdout": stdout.decode(),
        "stderr": stderr.decode()
    }

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
