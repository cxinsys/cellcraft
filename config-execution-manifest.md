## **Execution Manifest 필수 구성 요소**

### **1. Core Metadata Section**

- **manifest_version**: Manifest 스키마 버전
- **execution_id**: 고유 실행 식별자 (UUID)
- **execution_timestamp**:
    - start: 시작 시간
    - end: 종료 시간
    - duration: 총 실행 시간
- **project_name**: 프로젝트명
- **description**: 실행 작업 설명

### **2. Environment Specification**

- **cellcraft_version**: CellCraft 플랫폼 버전
- **platform**: 운영체제 플랫폼 (linux/amd64, darwin/arm64 등)
- **os_version**: 운영체제 상세 버전
- **cpu_info**:
    - model: CPU 모델명
    - cores: CPU 코어 수
- **memory_total**: 전체 메모리 용량
- **container_info** (if applicable):
    - engine: 컨테이너 엔진 종류 (docker, singularity 등)
    - image: 사용된 이미지 이름 및 태그
    - image_digest: 이미지 체크섬

### **3. Data Provenance**

- **input_files**: 각 입력 파일별
    - filename: 파일명
    - filepath: 파일 경로
    - checksum: 체크섬 값 (SHA256)
    - size: 파일 크기
    - format: 파일 형식 (h5ad, csv, fastq 등)

### **4. Workflow Execution Details**

- **workflow_name**: 워크플로우 이름
- **stages**: 각 스테이지별
    - stage_name: 스테이지 이름
    - plugin_name: 사용된 플러그인 이름
    - plugin_version: 플러그인 버전
    - execution_status: 실행 상태 (completed, failed, skipped)
    - parameters: 사용된 모든 파라미터와 값
    - input_files: 입력 파일 리스트 및 체크섬
    - output_files: 출력 파일 리스트 및 체크섬
    - start_time: 스테이지 시작 시간
    - end_time: 스테이지 종료 시간

### **5. Quality Control & Validation**

- **computational_metrics**:
    - total_runtime: 전체 실행 시간
- **output_validation**:
    - completeness_check: 출력 파일 완전성 확인
    - integrity_check: 출력 파일 무결성 확인

### **6. Dependencies & Software Stack**

- **plugin_dependencies**: 각 플러그인별
    - plugin_name: 플러그인명
    - version: 플러그인 버전
    - repository: 소스 저장소 URL
    - commit_hash: Git 커밋 해시 (if applicable)
- **system_libraries** (if required):
    - library_name: 시스템 라이브러리명
    - version: 라이브러리 버전

### **선택적 구성 요소 (권장)**

#### **Environment**
- gpu_info: GPU 정보 (모델, 메모리, CUDA 버전)
