import os

def snakemakeProcess(target, snakefile_path):
    from subprocess import Popen, PIPE
    print("Target:", target)
    print("Snakefile path:", snakefile_path)
    print("Environment PATH:", os.environ["PATH"])

    # Snakemake 명령어 구성
    command = [
        '/opt/conda/envs/snakemake/bin/snakemake',
        target,
        '--snakefile', snakefile_path,
        '-j',
    ]

    # Snakemake 프로세스 실행
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    print("STDOUT:", stdout.decode())
    print("STDERR:", stderr.decode())


    # 종료 요청을 감지하고 프로세스에 신호를 보냄
    # while True:
    #     if terminate_event.is_set():
    #         process.send_signal(signal.SIGTERM)
    #         break
    #     if process.poll() is not None:
    #         break
        
    return {"returncode": process.returncode, "stdout": stdout.decode(), "stderr": stderr.decode()}

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

def create_snakefile(file_path, user_input):
    # 기존 Snakefile 읽기
    with open("./workflow/Snakefile", 'r') as file:
        content = file.read()

    # {input...}과 {output...}을 제외한 {}로 감싸진 부분을 사용자 입력 값으로 대체
    for key, value in user_input.items():
        content = content.replace(f"{{{key}}}", value)

    # 새로운 Snakefile 생성
    with open(file_path, 'w') as file:
        file.write(content)

    # 생성된 파일 경로를 반환
    return file_path