"""
CellCraft Host GPU Environment Tests

올바른 GPU 테스트: 호스트 시스템의 GPU 환경 검증
- backend/celery 컨테이너가 아닌 호스트에서 GPU 환경 확인
- 실제 플러그인 컨테이너가 GPU에 접근할 수 있는지 확인
"""

import pytest
import subprocess
import re
import time


class TestHostGPUEnvironment:
    """호스트 시스템의 GPU 환경 검증"""

    @pytest.mark.gpu_mode
    def test_host_nvidia_driver(self):
        """
        호스트의 NVIDIA 드라이버 설치 및 동작 확인

        Success Criteria:
        - nvidia-smi 명령어 실행 가능
        - GPU 장치 감지됨
        - 드라이버 버전이 최소 요구사항 충족 (>= 450.0)
        """
        # nvidia-smi 실행 가능 확인
        result = subprocess.run(
            ["nvidia-smi"],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, \
            f"nvidia-smi not available on host: {result.stderr}"

        print(f"\n=== Host NVIDIA Driver ===")
        print(result.stdout)

        # 드라이버 버전 확인
        driver_result = subprocess.run(
            ["nvidia-smi", "--query-gpu=driver_version", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert driver_result.returncode == 0, "Cannot query driver version"

        driver_version = driver_result.stdout.strip().split('\n')[0]
        driver_major = int(driver_version.split('.')[0])

        assert driver_major >= 450, \
            f"Driver version {driver_version} too old (need >= 450.0 for CUDA 11.0)"

        print(f"✅ Host NVIDIA Driver: {driver_version}")

    @pytest.mark.gpu_mode
    def test_host_gpu_detection(self):
        """
        호스트의 GPU 장치 감지 및 정보 확인

        Success Criteria:
        - 최소 1개 이상의 GPU 감지
        - 각 GPU의 이름과 메모리 정보 획득
        - GPU 메모리가 최소 요구사항 충족 (>= 4GB)
        """
        # GPU 목록 및 정보 조회
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,name,memory.total", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, "Cannot query GPU information"

        gpu_lines = result.stdout.strip().split('\n')
        gpu_count = len(gpu_lines)

        assert gpu_count > 0, "No GPUs detected on host"

        print(f"\n=== Host GPU Detection ===")
        print(f"GPU Count: {gpu_count}")

        # 각 GPU 정보 파싱 및 검증
        for line in gpu_lines:
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 3:
                gpu_index = parts[0]
                gpu_name = parts[1]
                gpu_memory = parts[2]

                # 메모리 크기 파싱
                memory_match = re.search(r'(\d+)\s*([MG]iB)', gpu_memory)
                if memory_match:
                    memory_value = int(memory_match.group(1))
                    memory_unit = memory_match.group(2)
                    memory_gb = memory_value / 1024 if memory_unit == 'MiB' else memory_value

                    assert memory_gb >= 4.0, \
                        f"GPU {gpu_index} has insufficient memory: {memory_gb:.1f}GB (need >= 4GB)"

                    print(f"  GPU {gpu_index}: {gpu_name} ({memory_gb:.1f}GB)")

        print(f"✅ {gpu_count} GPU(s) detected on host")

    @pytest.mark.gpu_mode
    def test_host_nvidia_docker_runtime(self):
        """
        NVIDIA Docker Runtime 설치 및 설정 확인

        Success Criteria:
        - docker info에서 nvidia runtime 감지
        - nvidia runtime이 사용 가능 상태
        """
        # Docker info 조회
        result = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, "Cannot run docker info"

        docker_info = result.stdout.lower()

        # NVIDIA runtime 확인
        assert "nvidia" in docker_info, \
            "NVIDIA Docker runtime not found in docker info"

        print(f"\n=== NVIDIA Docker Runtime ===")
        print(f"✅ NVIDIA runtime detected in docker info")

        # Runtimes 섹션 출력
        if "runtimes:" in docker_info:
            lines = result.stdout.split('\n')
            in_runtime_section = False
            for line in lines:
                if 'Runtimes:' in line:
                    in_runtime_section = True
                    print(line)
                elif in_runtime_section:
                    if line.strip() and not line.startswith(' '):
                        break
                    print(line)

    @pytest.mark.gpu_mode
    def test_plugin_container_gpu_access(self):
        """
        플러그인 컨테이너가 GPU에 접근 가능한지 확인

        실제 CellCraft 동작:
        1. celery가 Docker-in-Docker로 플러그인 컨테이너 생성
        2. 플러그인 컨테이너 내부에서 GPU 사용
        3. 이 테스트는 그 과정을 시뮬레이션

        Success Criteria:
        - --gpus all 옵션으로 컨테이너 생성 가능
        - 컨테이너 내부에서 nvidia-smi 실행 가능
        - GPU 장치가 컨테이너 내부에서 감지됨
        """
        print(f"\n=== Plugin Container GPU Access Test ===")
        print(f"Testing: docker run --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi")

        # NVIDIA CUDA 베이스 이미지로 GPU 접근 테스트
        result = subprocess.run(
            [
                "docker", "run", "--rm", "--gpus", "all",
                "nvidia/cuda:11.8.0-base-ubuntu22.04",
                "nvidia-smi", "-L"
            ],
            capture_output=True,
            text=True,
            timeout=60
        )

        assert result.returncode == 0, \
            f"Plugin container cannot access GPU: {result.stderr}"

        gpu_output = result.stdout.strip()
        assert "GPU" in gpu_output, \
            f"No GPU detected in plugin container: {gpu_output}"

        gpu_count = len([line for line in gpu_output.split('\n') if 'GPU' in line])

        print(f"✅ Plugin container GPU access OK")
        print(f"  - {gpu_count} GPU(s) accessible from container")
        print(f"\nContainer output:")
        print(gpu_output)

    @pytest.mark.gpu_mode
    def test_host_gpu_memory_availability(self):
        """
        호스트의 GPU 메모리 상태 확인

        Success Criteria:
        - 각 GPU의 메모리 상태 조회 가능
        - 충분한 여유 메모리 확보 (최소 1GB)
        """
        # GPU 메모리 정보 조회
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=index,name,memory.total,memory.free,memory.used",
                "--format=csv,noheader,nounits"
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, "Cannot query GPU memory"

        print(f"\n=== Host GPU Memory ===")

        for line in result.stdout.strip().split('\n'):
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 5:
                gpu_index = parts[0]
                gpu_name = parts[1]
                total_mb = int(parts[2])
                free_mb = int(parts[3])
                used_mb = int(parts[4])

                # 최소 여유 메모리 확인
                assert free_mb >= 1024, \
                    f"GPU {gpu_index} has insufficient free memory: {free_mb}MB (need >= 1024MB)"

                print(f"GPU {gpu_index} ({gpu_name}):")
                print(f"  Total: {total_mb}MB ({total_mb/1024:.1f}GB)")
                print(f"  Free:  {free_mb}MB ({free_mb/1024:.1f}GB)")
                print(f"  Used:  {used_mb}MB ({used_mb/1024:.1f}GB)")

        print(f"✅ All GPUs have sufficient free memory")

    @pytest.mark.gpu_mode
    def test_host_cuda_driver_version(self):
        """
        호스트의 CUDA 드라이버 버전 및 호환성 확인

        Success Criteria:
        - CUDA 드라이버 버전 감지
        - CUDA 11.0 이상 지원 (드라이버 >= 450.0)
        """
        result = subprocess.run(
            ["nvidia-smi"],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, "nvidia-smi failed"

        output = result.stdout

        # CUDA 버전 파싱
        cuda_match = re.search(r'CUDA Version:\s+(\d+\.\d+)', output)

        if cuda_match:
            cuda_version = cuda_match.group(1)
            cuda_major = float(cuda_version)

            assert cuda_major >= 11.0, \
                f"CUDA version {cuda_version} too old (need >= 11.0)"

            print(f"\n=== Host CUDA Driver ===")
            print(f"✅ CUDA Version: {cuda_version}")
            print(f"  Compatible with CUDA 11.0+ requirements")
        else:
            print(f"\n⚠️  Could not parse CUDA version from nvidia-smi")
