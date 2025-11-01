# Deployment Tests

Docker Compose 배포 자동화 테스트

## 개요

이 디렉토리는 CellCraft의 Docker Compose 배포를 검증하는 통합 테스트를 포함합니다.

- **대상**: run-cpu-mode.sh, run-gpu-mode.sh, docker-compose.*.yml
- **범위**: 컨테이너 초기화, 플러그인 로드, CPU/GPU 모드 전환, 에러 복구
- **문서**: `docs/testing/DEPLOYMENT_TEST_PLAN.md`, `docs/testing/DEPLOYMENT_TEST_TASKS.md`

## 필수 요구사항

### 소프트웨어
- Docker 26.0+
- Docker Compose v2+
- Python 3.9+
- pytest 7.0+
- Git 2.0+

### 시스템 리소스
- **CPU 테스트**: 4 cores, 8GB RAM, 20GB 디스크
- **GPU 테스트**: NVIDIA GPU (CUDA 지원), 16GB RAM, 30GB 디스크

## 환경 설정

### 1. Git 서브모듈 초기화
```bash
git submodule update --init --recursive
```

### 2. 환경 변수 파일 준비
```bash
# .env.test 파일이 이미 프로젝트 루트에 있습니다
cat .env.test
```

### 3. 테스트 의존성 설치
```bash
pip install -r backend/requirements-test.txt
```

## 플랫폼별 실행 가이드

### 지원 플랫폼

CellCraft 배포 테스트는 3개 주요 플랫폼을 지원합니다:

| 플랫폼 | CPU 모드 | GPU 모드 | Custom Plugin | Docker 타입 |
|--------|----------|----------|---------------|-------------|
| **Linux (Ubuntu/Debian)** | ✅ 전체 지원 | ✅ 전체 지원 | ✅ 지원 | Native Docker Engine |
| **macOS (Intel/Apple Silicon)** | ✅ 지원 | ❌ **미지원** | ❌ **미지원** | Docker Desktop (필수) |
| **Windows 11 (WSL2)** | ✅ 전체 지원 | ✅ 지원 (조건부) | ✅ 지원 | Docker Desktop (WSL2 backend) |

### 플랫폼 감지

시스템 환경을 자동으로 감지하여 적절한 제약사항을 적용합니다:

```bash
# Docker Engine 타입 확인
docker info --format "{{.OperatingSystem}}"

# OS 버전 확인
# Linux
cat /etc/os-release

# macOS
sw_vers

# WSL2 확인
cat /proc/version | grep -i microsoft
```

### 플랫폼별 제약사항

#### ⚠️ macOS 제한사항

**GPU 모드 미지원**:
- macOS에는 NVIDIA GPU가 없어 GPU 모드 테스트 불가능
- GPU 관련 테스트 자동 스킵 (~15개 테스트)

**Custom Plugin 미지원**:
- Docker Desktop의 볼륨 마운트 제한으로 `is_editable=true` 플러그인 미지원
- CPU 모드 플러그인 수: 6개 (7개 중 Custom Plugin 제외)

**성능 오버헤드**:
- Docker Desktop 가상화로 인한 성능 저하
- Intel Mac: +30% 타임아웃 조정
- Apple Silicon (M1/M2): +50% 타임아웃 조정

#### ⚠️ Windows WSL2 고려사항

**프로젝트 위치 최적화**:
```bash
# ✅ 권장: WSL 파일시스템 사용
cd ~/cellcraft
pwd
# /home/username/cellcraft

# ❌ 비권장: Windows 파일시스템 (/mnt/c/)
cd /mnt/c/Users/username/cellcraft
# 볼륨 I/O 성능 30-50% 저하
```

**GPU 지원 요구사항**:
- Windows NVIDIA Driver 470.76+ 설치 (Windows에 설치, WSL 아님)
- WSL 커널 5.10.43+ (GPU 지원)
- DirectX Graphics 장치 확인:
  ```bash
  ls -l /dev/dxg  # 존재해야 함
  nvidia-smi      # WSL에서 작동해야 함
  ```

**성능 기대치**:
- WSL 파일시스템: +20% 타임아웃 조정
- /mnt/c/ 위치: 30-50% 추가 지연

### Docker Desktop 설정 (macOS/Windows)

Docker Desktop을 사용하는 플랫폼의 권장 리소스 설정:

**Settings → Resources**:
```
CPUs: 4+ cores (6+ 권장)
Memory: 6+ GB (8GB 권장)
Swap: 2+ GB
Disk: 20+ GB (30GB 권장)
```

**macOS 추가 설정**:
- Settings → General → "Use Virtualization framework" 활성화
- Settings → General → "Use Rosetta for x86/amd64 emulation" 활성화 (Apple Silicon)

**Windows WSL2 통합 확인**:
```powershell
# PowerShell에서 실행
wsl --version        # WSL 2 확인
wsl -l -v           # VERSION 열이 2여야 함
```

Docker Desktop → Settings → Resources → WSL Integration:
- "Enable integration with my default WSL distro" 활성화
- 사용하는 WSL 배포판 선택 (예: Ubuntu-22.04)

## 테스트 실행

### Quick Start (모든 플랫폼)

**전체 테스트 실행**:
```bash
pytest tests/deployment/ -v -s
```

**빠른 실행 (slow 테스트 제외)**:
```bash
pytest tests/deployment/ -m "not slow" -v
```

### 플랫폼별 실행 가이드

#### Linux (Ubuntu/Debian) - 기준 플랫폼

**전체 테스트 실행**:
```bash
# 모든 테스트 실행 (CPU + GPU 모드)
pytest tests/deployment/ -v -s
```

**성능 기대치**:
- CPU 모드 시작: <60초
- GPU 모드 시작: <90초
- 전체 테스트 소요시간: ~15-20분

**단계별 실행 전략**:
```bash
# Phase 1: 플랫폼 검증 (1-2분)
pytest tests/deployment/test_system_installation.py::TestSystemInstallation::test_docker_engine_detection -v
pytest tests/deployment/test_system_installation.py::TestSystemInstallation::test_platform_validation -v

# Phase 2: 시스템 설치 검증 (5-10분)
pytest tests/deployment/ -m "cpu_mode and not slow" -v

# Phase 3: CPU/GPU 모드 전환 (10-15분)
pytest tests/deployment/test_cpu_gpu_mode.py -v

# Phase 4: 스크립트 옵션 검증 (5-10분)
pytest tests/deployment/test_script_options.py -v

# Phase 5: 에러 복구 (10-15분)
pytest tests/deployment/test_error_recovery.py -v

# Phase 6: 성능 벤치마크 (5-10분)
pytest tests/deployment/test_performance.py -v
```

#### macOS (Intel/Apple Silicon)

**⚠️ 중요**: GPU 및 Custom Plugin 테스트는 자동으로 스킵됩니다.

**권장 실행 명령어**:
```bash
# GPU 및 macOS 비호환 테스트 제외
pytest tests/deployment/ -m "not gpu_mode and not macos_skip" -v
```

**Docker Desktop 설정 확인**:
```bash
# Docker Desktop 리소스 확인
docker info | grep -E "CPUs|Total Memory"

# 최소 요구사항 확인 (4 CPUs, 6GB RAM)
```

**Apple Silicon 고려사항**:
- Rosetta 2 활성화 필요
- 타임아웃 1.5배 증가 (90초 CPU 시작)
- 일부 테스트에서 아키텍처 경고 발생 가능

**예상 테스트 결과**:
- 전체 테스트: ~50개
- 스킵 (GPU): ~15개
- 스킵 (macOS 비호환): ~3개
- 실행: ~32개
- 목표 통과율: 90%+ (실행된 테스트 기준)

**플랫폼 검증 테스트**:
```bash
# macOS 플랫폼 및 Docker Desktop 감지 확인
pytest tests/deployment/test_system_installation.py::TestSystemInstallation::test_docker_engine_detection -v -s
pytest tests/deployment/test_system_installation.py::TestSystemInstallation::test_platform_validation -v -s
```

#### Windows 11 (WSL2)

**⚠️ 필수**: WSL2 터미널에서 실행 (PowerShell/Command Prompt 아님).

**사전 요구사항 확인**:
```bash
# WSL 버전 확인 (PowerShell에서)
wsl --version
wsl -l -v  # VERSION이 2여야 함

# WSL 터미널에서 Docker 확인
docker version
docker compose version

# 프로젝트 위치 확인 (WSL 터미널에서)
pwd
# ✅ /home/username/cellcraft (권장)
# ❌ /mnt/c/Users/... (성능 저하)
```

**테스트 실행** (WSL 터미널에서):
```bash
# 전체 테스트 실행
pytest tests/deployment/ -v -s
```

**GPU 지원 검증**:
```bash
# DirectX 장치 확인
ls -l /dev/dxg

# NVIDIA GPU 확인
nvidia-smi

# GPU 테스트 실행
pytest tests/deployment/ -m "gpu_mode" -v
```

**성능 기대치**:
- CPU 모드 시작: <72초 (WSL FS), <90초+ (/mnt/c/)
- GPU 모드 시작: <108초
- 전체 테스트: ~18-24분

**파일시스템 위치별 성능**:
```bash
# 볼륨 I/O 성능 테스트
pytest tests/deployment/test_system_installation.py::TestSystemInstallation::test_filesystem_compatibility -v -s --run-slow-tests

# 예상 결과:
# WSL FS: Read/Write <15ms ✅
# /mnt/c/: Read/Write <50ms (경고) ⚠️
```

**Docker Desktop WSL2 통합 확인**:
```bash
# Docker가 WSL2 백엔드를 사용하는지 확인
docker info --format "{{.OperatingSystem}}"
# 예상 출력: "Docker Desktop" 또는 "Docker Desktop (WSL2 backend)"
```

### 플랫폼 비교

| 플랫폼 | 권장 명령어 | 예상 소요시간 | 목표 통과율 |
|--------|------------|--------------|------------|
| **Linux** | `pytest tests/deployment/ -v` | ~15-20분 | 100% |
| **macOS** | `pytest tests/deployment/ -m "not gpu_mode and not macos_skip" -v` | ~10-15분 | 90%+ (실행 기준) |
| **WSL2** | `pytest tests/deployment/ -v` | ~18-24분 | 95%+ |

### 권장 실행 순서 (모든 플랫폼)

**Step 1: 플랫폼 검증** (~2분):
```bash
pytest tests/deployment/test_system_installation.py -k "test_docker_engine_detection or test_platform_validation" -v
```

**Step 2: CPU 모드 빠른 테스트** (~5분):
```bash
pytest tests/deployment/ -m "cpu_mode and not slow" -v
```

**Step 3: CPU 모드 전체** (~10분):
```bash
pytest tests/deployment/ -m "cpu_mode" -v
```

**Step 4: GPU 모드 (Linux/WSL2만)** (~15분):
```bash
pytest tests/deployment/ -m "gpu_mode" -v
```

**Step 5: 에러 복구 및 스크립트 옵션** (~10분):
```bash
pytest tests/deployment/test_error_recovery.py test_script_options.py -v
```

**Step 6: 성능 벤치마크** (~5분):
```bash
pytest tests/deployment/test_performance.py -v
```

### 카테고리별 실행

#### CPU 모드만
```bash
pytest tests/deployment/ -m cpu_mode -v
```

#### GPU 모드만 (GPU 환경 필요)
```bash
pytest tests/deployment/ -m gpu_mode -v
```

#### 빠른 테스트만 (slow 제외)
```bash
pytest tests/deployment/ -m "not slow" -v
```

#### 특정 테스트 파일
```bash
pytest tests/deployment/test_system_installation.py -v
pytest tests/deployment/test_cpu_gpu_mode.py -v
```

### 상세 출력 옵션

#### 실시간 출력 보기
```bash
pytest tests/deployment/ -v -s
```

#### 실패한 테스트만 재실행
```bash
pytest tests/deployment/ --lf
```

#### 커버리지 리포트
```bash
pytest tests/deployment/ --cov=tests/deployment --cov-report=html
```

## 테스트 마커

### 기본 마커 (5개)

| 마커 | 설명 | 사용 예 |
|------|------|---------|
| `cpu_mode` | CPU 모드 배포 테스트 | `@pytest.mark.cpu_mode` |
| `gpu_mode` | GPU 모드 배포 테스트 | `@pytest.mark.gpu_mode` |
| `slow` | 30초 이상 소요되는 테스트 | `@pytest.mark.slow` |
| `requires_network` | 인터넷 연결 필요 | `@pytest.mark.requires_network` |
| `requires_gpu` | NVIDIA GPU 필요 | `@pytest.mark.requires_gpu` |

### 플랫폼 마커 (8개)

| 마커 | 설명 | 사용 예 |
|------|------|---------|
| `windows` | Windows WSL2 전용 테스트 | `@pytest.mark.windows` |
| `macos` | macOS 전용 테스트 | `@pytest.mark.macos` |
| `linux` | Linux 전용 테스트 | `@pytest.mark.linux` |
| `macos_skip` | macOS 비호환 테스트 (GPU, Custom Plugin) | `@pytest.mark.macos_skip` |
| `docker_desktop` | Docker Desktop 전용 테스트 | `@pytest.mark.docker_desktop` |
| `native_docker` | Native Docker Engine 전용 테스트 | `@pytest.mark.native_docker` |
| `arm64` | ARM64/Apple Silicon 아키텍처 테스트 | `@pytest.mark.arm64` |
| `amd64` | AMD64/x86_64 아키텍처 테스트 | `@pytest.mark.amd64` |

### 마커 사용 예제

#### 기본 사용법

```python
import pytest

@pytest.mark.cpu_mode
def test_cpu_plugin_count(cpu_mode_running):
    """CPU 모드 플러그인 개수 테스트"""
    plugins = get_plugin_list()
    assert len(plugins) >= 6

@pytest.mark.gpu_mode
@pytest.mark.requires_gpu
def test_gpu_plugin_count(gpu_mode_running):
    """GPU 모드 플러그인 개수 테스트"""
    plugins = get_plugin_list()
    assert len(plugins) >= 8

@pytest.mark.macos_skip
def test_custom_plugin_loading(cpu_mode_running):
    """Custom Plugin 로드 테스트 (macOS 제외)"""
    plugins = get_plugin_list()
    custom_plugins = [p for p in plugins if p.get("is_editable")]
    assert len(custom_plugins) > 0
```

#### 복합 마커 조합

```bash
# CPU 모드 + macOS 호환 테스트만
pytest tests/deployment/ -m "cpu_mode and not macos_skip" -v

# GPU 모드 + Linux 전용
pytest tests/deployment/ -m "gpu_mode and linux" -v

# 빠른 테스트 + macOS 비호환 제외
pytest tests/deployment/ -m "not slow and not macos_skip" -v

# Docker Desktop 환경 테스트
pytest tests/deployment/ -m "docker_desktop" -v

# ARM64 아키텍처 테스트
pytest tests/deployment/ -m "arm64" -v
```

### 플랫폼별 권장 마커 조합

#### Linux (모든 테스트)
```bash
# 전체 실행
pytest tests/deployment/ -v

# 또는 Linux 전용만
pytest tests/deployment/ -m "linux" -v
```

#### macOS (GPU 및 비호환 제외)
```bash
# GPU 및 macOS 비호환 테스트 제외
pytest tests/deployment/ -m "not gpu_mode and not macos_skip" -v

# macOS 전용 테스트만
pytest tests/deployment/ -m "macos" -v

# Docker Desktop 테스트
pytest tests/deployment/ -m "docker_desktop" -v
```

#### Windows WSL2 (전체 또는 선택)
```bash
# 전체 실행
pytest tests/deployment/ -v

# WSL2 전용 테스트만
pytest tests/deployment/ -m "windows" -v

# GPU 테스트만 (GPU 드라이버 설치 시)
pytest tests/deployment/ -m "gpu_mode and windows" -v
```

### 마커 확인 및 미리보기

**등록된 모든 마커 확인**:
```bash
pytest --markers
```

**특정 마커 조합으로 테스트 수집 미리보기**:
```bash
# 실행하지 않고 테스트 목록만 확인
pytest tests/deployment/ -m "not gpu_mode" --collect-only

# 예상 출력:
# <Module test_system_installation.py>
#   <Class TestSystemInstallation>
#     <Function test_docker_engine_detection>
#     <Function test_platform_validation>
#     ...
```

**마커별 테스트 수 확인**:
```bash
# CPU 모드 테스트 개수
pytest tests/deployment/ -m "cpu_mode" --collect-only -q

# GPU 모드 테스트 개수
pytest tests/deployment/ -m "gpu_mode" --collect-only -q

# macOS 호환 테스트 개수
pytest tests/deployment/ -m "not gpu_mode and not macos_skip" --collect-only -q
```

## 디렉토리 구조

```
tests/deployment/
├── __init__.py                  # 패키지 초기화
├── conftest.py                  # pytest fixtures (11개)
├── helpers.py                   # 유틸리티 함수 (16개)
├── README.md                    # 이 파일
├── test_system_installation.py # 시스템 설치 테스트
├── test_cpu_gpu_mode.py         # CPU/GPU 모드 테스트
├── test_script_options.py       # 스크립트 옵션 테스트
├── test_error_recovery.py       # 에러 복구 테스트
└── test_performance.py          # 성능 벤치마크
```

## Fixtures

### Configuration Fixtures (4개)
- `project_root`: 프로젝트 루트 디렉토리
- `docker_compose_config`: Docker Compose 파일 설정
- `deployment_scripts`: 배포 스크립트 경로
- `container_names`: 컨테이너 이름 매핑

### Platform Detection Fixtures (3개 - NEW)
- `platform_info`: OS 타입, 아키텍처, Docker Engine 감지
- `platform_constraints`: 타임아웃 배수, GPU/Custom Plugin 지원, 예상 플러그인 수
- `docker_environment_info`: Docker 버전, I/O 성능, WSL 파일시스템 위치

### Environment Fixtures (4개)
- `clean_docker_environment`: 테스트 전/후 Docker 정리
- `cpu_mode_running`: CPU 모드로 시스템 자동 시작 (플랫폼 인식)
- `gpu_mode_running`: GPU 모드로 시스템 자동 시작 (플랫폼 인식, macOS 자동 스킵)
- `script_execution_env`: 스크립트 실행 환경 설정

## Helper Functions

### Platform Detection (6개 - NEW)
- `detect_docker_engine()`: Docker Engine 타입 감지 (Native/Desktop/WSL2)
- `get_os_name()`: 사람이 읽을 수 있는 OS 이름 (Ubuntu 22.04, macOS 14.2)
- `check_nvidia_gpu()`: NVIDIA GPU 가용성 확인
- `check_wsl2_gpu()`: WSL2 GPU 지원 확인 (/dev/dxg + nvidia-smi)
- `detect_wsl_filesystem()`: WSL vs Windows 파일시스템 감지
- `measure_volume_io_latency()`: Docker 볼륨 I/O 성능 측정

### Container Management (3개)
- `wait_for_container_healthy()`: 컨테이너 healthy 상태 대기
- `wait_for_container_running()`: 컨테이너 running 상태 대기
- `get_container_status()`: 컨테이너 상태 조회

### Service Connectivity (2개)
- `check_service_accessible()`: HTTP 서비스 접근 확인
- `wait_for_service_ready()`: 서비스 준비 완료 대기

### API & Data (4개)
- `get_plugin_list()`: 플러그인 API 조회
- `get_current_git_branch()`: Git 브랜치 확인
- `get_container_env()`: 컨테이너 환경 변수 조회
- `measure_execution_time()`: 함수 실행 시간 측정

### Utilities (1개)
- `run_script()`: 배포 스크립트 실행

## 트러블슈팅

### 컨테이너가 시작되지 않음
```bash
# 로그 확인
docker compose -f docker-compose.cpu.yml logs

# 컨테이너 상태 확인
docker compose -f docker-compose.cpu.yml ps

# 수동으로 정리
docker compose -f docker-compose.cpu.yml down -v
```

### 플러그인 로드 실패
```bash
# Git 서브모듈 상태 확인
git submodule status

# 서브모듈 재초기화
git submodule update --init --recursive --force

# 올바른 브랜치인지 확인
cd backend/plugin/official
git branch --show-current
```

### 테스트 타임아웃
```bash
# 타임아웃 시간 증가
pytest tests/deployment/ --timeout=600

# 또는 특정 테스트만 실행
pytest tests/deployment/test_system_installation.py::TestSystemInstallation::test_container_startup_sequence -v
```

## CI/CD 통합

GitHub Actions 워크플로우는 `.github/workflows/deployment-tests.yml`에 정의되어 있습니다.

### 로컬에서 CI 시뮬레이션
```bash
# 환경 변수 설정
cp .env.test .env

# 전체 테스트 실행 (CI와 동일)
pytest tests/deployment/ -v --tb=short
```

## 참고 자료

- **테스트 계획서**: `docs/testing/DEPLOYMENT_TEST_PLAN.md`
- **작업 리스트**: `docs/testing/DEPLOYMENT_TEST_TASKS.md`
- **메인 문서**: `README.md` (프로젝트 루트)
- **기존 테스트**: `backend/tests/` (단위 및 통합 테스트)

## 기여 가이드

새로운 테스트 추가 시:
1. 적절한 테스트 파일 선택 (또는 새로 생성)
2. 테스트 클래스 및 함수 작성
3. 적절한 마커 추가 (`@pytest.mark.cpu_mode` 등)
4. Docstring 작성 (TC-XXX-XXX 참조)
5. `DEPLOYMENT_TEST_PLAN.md` 업데이트

## 라이선스

프로젝트 라이선스 참조
