# CellCraft FastAPI 백엔드 테스트 코드 작성 & 리팩토링 계획서

## 📋 Executive Summary

**목표**:
1. FastAPI 백엔드의 API 단위 테스트 인프라 구축 및 핵심 기능 테스트 커버리지 확보
2. 테스트 작성 과정에서 발견되는 코드 품질 이슈 개선 및 리팩토링

**타겟 커버리지**: 초기 60% → 중기 75% → 장기 85%
**예상 기간**: 8-12주 (테스트 작성 + 리팩토링 포함)

**현재 상태** (2025-01-XX):
- ✅ 테스트 환경 구축 완료 (PostgreSQL test-db, Conda environment)
- ✅ 기본 픽스처 구현 완료 (conftest.py)
- ✅ 인증 API 샘플 테스트 작성 (11개 테스트)
- 📊 초기 테스트 결과: **8 passed, 3 failed** (격리 문제 해결 완료)
- 📈 현재 커버리지: **20%**

---

## 🎯 테스트 환경 구축 진행 상황

### ✅ 완료된 작업

#### 1. Conda 테스트 환경 생성 (environment-test.yml)
```yaml
✓ 환경 이름: cellcraft-test
✓ Python 버전: 3.10.6
✓ 포함 내용:
  - 기존 environment.yml의 모든 의존성
  - 테스트 도구: pytest, pytest-asyncio, pytest-cov, pytest-mock
  - HTTP 클라이언트: httpx==0.24.1 (starlette 0.25.0 호환)
  - 코드 품질 도구: black, isort, pylint, flake8
```

**생성 및 활성화**:
```bash
conda env create -f backend/environment-test.yml
conda activate cellcraft-test
```

#### 2. Docker Compose 테스트 데이터베이스 추가
```yaml
✓ 서비스 이름: test-db
✓ 이미지: postgres:15
✓ 포트: 5433 (호스트) → 5432 (컨테이너)
✓ 데이터베이스: cellcraft_test
✓ 사용자: test_user / test_pass
✓ 특징:
  - 메인 db 서비스와 포트 충돌 방지 (5433 사용)
  - 볼륨 미사용 (ephemeral, 테스트 간 격리)
  - 재시작 정책 없음 (수동 제어)
```

**실행 방법**:
```bash
docker compose -f docker-compose.dev.yml up test-db -d
```

#### 3. 데이터베이스 연결 설정 (app/database/conn.py)
```python
✓ TESTING 환경 변수 감지
✓ TESTING=1 시 PostgreSQL 테스트 DB 사용
  - URI: postgresql://test_user:test_pass@localhost:5433/cellcraft_test
  - Echo: False (SQL 로깅 비활성화)
  - Pool pre-ping: True (연결 유지)
✓ 프로덕션 스키마와 100% 동일 (JSONB 타입 보존)
```

#### 4. pytest 설정 파일 (pytest.ini)
```ini
✓ 테스트 경로: tests/
✓ 커버리지: app/ 디렉토리
✓ 리포트: HTML + 터미널
✓ 마커 정의:
  - unit: 단위 테스트
  - integration: 통합 테스트
  - auth: 인증 테스트
```

#### 5. 공통 픽스처 구현 (tests/conftest.py)
```python
✓ test_engine: PostgreSQL 세션 스코프 엔진
✓ db_session: 함수 스코프, 테스트 격리 보장
  - 각 테스트 전 테이블 drop/create
  - 테스트 간 데이터 격리 완벽 보장
✓ client: FastAPI TestClient with DB override
✓ sample_user: 테스트 사용자 픽스처
✓ sample_admin_user: 관리자 사용자 픽스처
✓ sample_plugin: 로컬 플러그인 픽스처
✓ sample_workflow: 샘플 워크플로우 픽스처
✓ sample_file: 파일 메타데이터 픽스처
✓ auth_headers: JWT 인증 헤더 픽스처
✓ admin_auth_headers: 관리자 인증 헤더 픽스처
```

#### 6. 샘플 테스트 작성 (tests/unit/test_auth.py)
```python
✓ 11개 인증 API 테스트:
  - 회원가입 (정상, 중복 이메일)
  - 로그인 (정상, 잘못된 자격증명)
  - 토큰 검증 (정상, 만료, 없음)
  - 사용자 정보 조회
✓ AAA 패턴 (Arrange-Act-Assert) 적용
```

### 📊 현재 테스트 실행 결과

**명령어**:
```bash
cd /home/dmshin/cellcraft/backend
conda run -n cellcraft-test pytest -v
```

**결과** (2025-01-XX):
```
================================== test session starts ==================================
collected 11 items

tests/unit/test_auth.py::test_register_success PASSED                            [  9%]
tests/unit/test_auth.py::test_register_duplicate_email PASSED                    [ 18%]
tests/unit/test_auth.py::test_login_success PASSED                               [ 27%]
tests/unit/test_auth.py::test_login_invalid_credentials PASSED                   [ 36%]
tests/unit/test_auth.py::test_login_nonexistent_user PASSED                      [ 45%]
tests/unit/test_auth.py::test_get_current_user_success PASSED                    [ 54%]
tests/unit/test_auth.py::test_get_current_user_without_token PASSED              [ 63%]
tests/unit/test_auth.py::test_get_current_user_with_invalid_token PASSED         [ 72%]
tests/unit/test_auth.py::test_refresh_token FAILED                               [ 81%]
tests/unit/test_auth.py::test_associate_plugin FAILED                            [ 90%]
tests/unit/test_auth.py::test_dissociate_plugin FAILED                           [100%]

======================== 8 passed, 3 failed in 2.34s ==========================

Coverage: 20%
```

**실패 원인 분석**:
- 3개 실패: API 응답 스키마 검증 오류 (is_active, hashed_password 필드 누락)
- 격리 문제 (UniqueViolation) 해결 완료
- 다음 단계: Pydantic 스키마 수정 필요

### ⚠️ 주요 의사결정 및 변경 사항

#### 1. SQLite → PostgreSQL 테스트 DB 전환
**원래 계획** (Section 4.1):
```python
# In-Memory SQLite 사용
engine = create_engine(
    "sqlite:///:memory:",
    connect_args={"check_same_thread": False},
    poolclass=StaticPool,
)
```

**변경 이유**:
- JSONB 타입이 SQLite에서 지원되지 않음
- 프로덕션 스키마 수정은 부적절함
- PostgreSQL 테스트 DB로 프로덕션 환경과 100% 동일하게 유지

**최종 구현**:
- Docker test-db 서비스 (PostgreSQL 15)
- 포트 5433으로 격리
- TESTING=1 환경 변수로 자동 전환

#### 2. 테스트 격리 전략 변경
**원래 계획** (Section 4.1):
```python
# 트랜잭션 롤백 방식
connection = test_db_engine.connect()
transaction = connection.begin()
# ...
transaction.rollback()  # 테스트 후 롤백
```

**변경 이유**:
- 픽스처에서 commit() 사용 시 롤백 불가
- UniqueViolation 에러 발생 (이메일 중복)

**최종 구현**:
```python
# 각 테스트 전 테이블 재생성
Base.metadata.drop_all(bind=test_engine)
Base.metadata.create_all(bind=test_engine)
# → 완벽한 격리 보장
```

#### 3. 의존성 관리 방식
**원래 계획** (Section 1.1):
```yaml
# environment.yml 업데이트
dependencies:
  - pytest>=8.3.0
  - pytest-asyncio>=0.23.0
  # ...
```

**최종 구현**:
- `environment-test.yml` 별도 생성
- environment.yml 전체 복사 + 테스트 도구 추가
- Conda 환경 격리로 충돌 방지

#### 4. httpx 버전 호환성 이슈
**문제**:
```
TypeError: Client.__init__() got an unexpected keyword argument 'app'
```

**해결**:
- httpx 0.28.1 → 0.24.1 다운그레이드
- starlette 0.25.0과 호환성 확보

### 🚧 다음 단계

#### 우선순위 P0 (즉시 처리)
- [ ] 실패한 3개 테스트 수정 (Pydantic 스키마)
- [ ] test_workflows.py 작성 (8개 테스트)
- [ ] test_models.py 작성 (모델 검증)

#### 우선순위 P1 (Week 3-4)
- [ ] test_plugins.py 작성 (25개 테스트)
- [ ] 플러그인 API 리팩토링 시작
- [ ] Repository 패턴 적용 준비

#### 인프라 개선 (선택)
- [ ] pytest-xdist로 병렬 실행 설정
- [ ] GitHub Actions CI 워크플로우 추가
- [ ] 커버리지 리포트 자동 생성

### 📝 테스트 실행 명령어

**전체 테스트 실행**:
```bash
conda run -n cellcraft-test pytest -v
```

**커버리지 포함**:
```bash
conda run -n cellcraft-test pytest -v --cov=app --cov-report=html
```

**특정 테스트만 실행**:
```bash
conda run -n cellcraft-test pytest tests/unit/test_auth.py -v
```

**마커 기반 실행**:
```bash
conda run -n cellcraft-test pytest -m unit -v
```

---

## 1. 테스트 환경 설정

### 1.1 의존성 패키지

**⚠️ 실제 구현**: `environment-test.yml` 별도 생성 (원래 계획과 다름)

**파일**: `backend/environment-test.yml`

```yaml
name: cellcraft-test
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  # [기존 environment.yml의 모든 의존성 포함]
  # ... (생략) ...

  - pip:
      # [기존 pip 의존성]

      # Testing dependencies (추가됨)
      - pytest-asyncio>=0.23.0
      - pytest-cov>=4.1.0
      - pytest-mock>=3.12.0
      - httpx==0.24.1               # starlette 0.25.0 호환
      - faker>=20.0.0
      - pytest-xdist>=3.5.0
      - pytest-timeout>=2.2.0

      # Code quality tools (추가됨)
      - black>=23.0.0
      - isort>=5.12.0
      - pylint>=3.0.0
      - flake8>=6.0.0
```

**생성 이유**:
- 기존 environment.yml 보존 (프로덕션 의존성)
- 테스트 전용 의존성 격리
- 전체 앱 의존성 + 테스트 도구 통합

**사용법**:
```bash
# 환경 생성
conda env create -f backend/environment-test.yml

# 환경 활성화
conda activate cellcraft-test

# 테스트 실행
pytest -v
```

### 1.2 pytest 설정
```ini
# backend/pytest.ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts = 
    -v
    --strict-markers
    --cov=app
    --cov-report=html
    --cov-report=term-missing
    -n auto

markers =
    unit: 단위 테스트
    integration: 통합 테스트
    slow: 느린 테스트 (5초 이상)
    celery: Celery 태스크 테스트
    refactor: 리팩토링 필요 항목
```

---

## 2. 현재 프로젝트 구조 분석

### 2.1 실제 백엔드 구조 (코드베이스 검토 결과)

```
backend/
├── app/
│   ├── routes/
│   │   ├── api.py                    # API 라우터 통합
│   │   ├── dep.py                    # 의존성 주입 (get_db, get_current_user)
│   │   ├── celery_tasks.py           # Celery 태스크 정의
│   │   └── endpoints/                # 도메인별 라우터
│   │       ├── auth.py               # 인증 (92줄) ✓ 단순함
│   │       ├── workflow.py           # 워크플로우 관리
│   │       ├── plugin.py             # 플러그인 관리 (1576줄) ⚠️ 매우 복잡
│   │       ├── files.py              # 파일 업로드/다운로드
│   │       ├── task.py               # 태스크 실행 및 모니터링
│   │       ├── admin.py              # 관리자 기능
│   │       └── datatable.py          # 데이터테이블 조회
│   │
│   ├── database/
│   │   ├── models.py                 # SQLAlchemy 모델 (User, Workflow, Plugin, Task, File)
│   │   ├── conn.py                   # 데이터베이스 연결 및 세션
│   │   ├── crud/                     # CRUD 레이어 (단순 쿼리 래퍼)
│   │   │   ├── base.py
│   │   │   ├── crud_user.py          # 74줄
│   │   │   ├── crud_workflow.py      # 42줄
│   │   │   ├── crud_plugin.py        # 376줄
│   │   │   ├── crud_task.py
│   │   │   ├── crud_file.py
│   │   │   └── crud_admin.py
│   │   └── schemas/                  # Pydantic 스키마
│   │       ├── user.py
│   │       ├── workflow.py
│   │       ├── plugin.py
│   │       ├── task.py
│   │       ├── file.py
│   │       └── admin.py
│   │
│   ├── common/
│   │   ├── config.py                 # 설정 관리 (Celery, DB, CORS)
│   │   ├── security.py               # JWT 토큰 생성/검증, 비밀번호 해싱
│   │   ├── enums.py                  # Enum 타입 정의
│   │   └── utils/                    # 유틸리티 함수들
│   │       ├── docker_utils.py       # Docker 컨테이너 관리
│   │       ├── snakemake_utils.py    # Snakemake 워크플로우 실행
│   │       ├── plugin_utils.py       # 플러그인 파일 처리
│   │       ├── celery_utils.py       # Celery 앱 생성
│   │       ├── workflow_utils.py     # 워크플로우 검증
│   │       ├── h5ad_utils.py         # H5AD 파일 처리
│   │       ├── cache_utils.py        # 결과 캐싱
│   │       ├── error_utils.py        # 에러 처리
│   │       └── plugin_*              # 플러그인 관련 유틸리티들
│   │
│   └── main.py                       # FastAPI 애플리케이션 엔트리포인트
│
├── plugin/                           # 플러그인 저장소
│   ├── local/                        # 사용자 생성 플러그인
│   └── official/                     # 공식 플러그인
│
├── user/                             # 사용자 데이터
├── alembic/                          # 데이터베이스 마이그레이션
├── requirements.txt
└── environment.yml
```

### 2.2 주요 발견사항

#### ✅ 잘 구성된 부분
1. **모델 계층**: SQLAlchemy ORM을 사용한 명확한 모델 정의
2. **스키마 검증**: Pydantic을 활용한 입출력 검증
3. **비동기 처리**: Celery를 통한 장기 실행 태스크 처리
4. **의존성 주입**: 기본적인 DI 패턴 (get_db, get_current_user)

#### ⚠️ 개선 필요 부분
1. **비즈니스 로직 위치**
   - **문제**: 라우터 함수에 비즈니스 로직이 직접 구현됨
   - **예시**: `plugin.py` 1576줄 중 대부분이 비즈니스 로직
   - **영향**: 테스트 어려움, 재사용성 낮음, 유지보수 복잡

2. **CRUD 레이어의 역할**
   - **현재**: 단순 쿼리 래퍼 함수 모음
   - **부족**: 도메인 로직, 복잡한 쿼리, 트랜잭션 관리
   - **예시**: `crud_workflow.py`는 단순 CRUD만 제공

3. **중복 코드**
   - 파일 업로드/백업/롤백 로직이 여러 엔드포인트에 반복
   - 에러 핸들링 패턴이 일관되지 않음
   - Docker 이미지 처리 로직 중복

4. **복잡도 문제**
   - `upload_scripts()`: 454줄의 단일 함수
   - `upload_plugin()`: 287줄의 단일 함수
   - 순환 복잡도(Cyclomatic Complexity) 높음

5. **누락된 계층**
   - **서비스 레이어**: 비즈니스 로직을 캡슐화할 계층 부재
   - **Repository 패턴**: 데이터 접근을 추상화할 패턴 미적용
   - **도메인 이벤트**: 플러그인 빌드 완료 등의 이벤트 처리 부재

### 2.3 제안하는 테스트 디렉토리 구조

```
backend/
├── app/                        # 기존 애플리케이션 코드
│   ├── routes/endpoints/       # 기존: 비즈니스 로직 포함 (리팩토링 대상)
│   ├── database/crud/          # 기존: 단순 쿼리 래퍼
│   ├── services/               # 🆕 신규: 비즈니스 로직 분리
│   │   ├── auth_service.py
│   │   ├── workflow_service.py
│   │   ├── plugin_service.py
│   │   ├── file_service.py
│   │   └── task_service.py
│   ├── repositories/           # 🆕 신규: Repository 패턴 적용
│   │   ├── base_repository.py
│   │   ├── user_repository.py
│   │   ├── workflow_repository.py
│   │   ├── plugin_repository.py
│   │   └── task_repository.py
│   └── exceptions.py           # 🆕 신규: 커스텀 예외 클래스
│
└── tests/
    ├── conftest.py             # 공통 픽스처 (DB, Client, Auth)
    ├── pytest.ini              # pytest 설정
    │
    ├── unit/                   # 단위 테스트 (격리된 테스트)
    │   ├── test_models.py      # SQLAlchemy 모델 테스트
    │   ├── test_schemas.py     # Pydantic 스키마 검증
    │   ├── services/           # 서비스 레이어 테스트
    │   │   ├── test_auth_service.py
    │   │   ├── test_workflow_service.py
    │   │   ├── test_plugin_service.py
    │   │   └── test_file_service.py
    │   ├── repositories/       # Repository 테스트
    │   │   ├── test_user_repository.py
    │   │   ├── test_workflow_repository.py
    │   │   └── test_plugin_repository.py
    │   └── utils/              # 유틸리티 함수 테스트
    │       ├── test_docker_utils.py
    │       ├── test_plugin_utils.py
    │       └── test_workflow_utils.py
    │
    ├── api/                    # API 엔드포인트 테스트 (통합)
    │   ├── test_auth.py        # 인증 API (회원가입, 로그인, JWT)
    │   ├── test_workflows.py   # 워크플로우 CRUD
    │   ├── test_plugins.py     # 플러그인 업로드, 빌드, 버전 관리
    │   ├── test_tasks.py       # 태스크 실행 및 모니터링
    │   ├── test_files.py       # 파일 업로드/다운로드
    │   ├── test_admin.py       # 관리자 기능
    │   └── test_datatable.py   # 데이터테이블 조회
    │
    ├── integration/            # 통합 테스트 (E2E 시나리오)
    │   ├── test_workflow_execution.py    # 워크플로우 전체 실행 흐름
    │   ├── test_plugin_lifecycle.py      # 플러그인 생명주기
    │   └── test_celery_tasks.py          # Celery 태스크 통합
    │
    └── fixtures/               # 테스트 데이터
        ├── sample_workflow.json
        ├── sample_plugin.json
        ├── sample_drawflow.json
        └── sample.h5ad
```

---

## 3. 코드 품질 이슈 상세 분석

### 3.1 복잡도 메트릭스 (실제 측정치)

| 파일 | 줄 수 | 함수 수 | 평균 함수 길이 | 최대 함수 길이 | 복잡도 점수 |
|------|------|---------|----------------|----------------|------------|
| `plugin.py` | 1576 | 30 | 52줄 | 454줄 (`upload_scripts`) | 🔴 매우 높음 |
| `auth.py` | 92 | 3 | 30줄 | 48줄 | 🟢 낮음 |
| `crud_plugin.py` | 376 | 12 | 31줄 | - | 🟡 중간 |
| `crud_workflow.py` | 42 | 6 | 7줄 | 10줄 | 🟢 낮음 |
| `crud_user.py` | 74 | 8 | 9줄 | 19줄 | 🟢 낮음 |
| `celery_tasks.py` | 200+ | 3 | 66줄 | 200줄+ | 🔴 높음 |

### 3.2 문제점 상세 분석

#### 🔴 **Critical: plugin.py 복잡도**

**문제 요약**:
- 1576줄의 단일 파일에 30개 함수
- 비즈니스 로직, 파일 I/O, Docker 관리, DB 작업이 혼재
- 테스트 불가능한 구조

**주요 함수 분석**:

```python
# 1. upload_scripts() - 454줄
# 문제: 파일 업로드 + 백업 + 검증 + 롤백 + 복원 로직이 하나의 함수에
# 책임: 7개 이상 (SRP 위반)
# - 파일 저장
# - 스테이징 디렉토리 관리
# - 백업 생성/복원
# - 에러 핸들링
# - 로깅
# - Reference 폴더 관리
# - 트랜잭션 관리

# 2. upload_plugin() - 287줄
# 문제: 플러그인 검증 + DB 저장 + 파일 시스템 관리
# 책임: 5개 이상
# - 플러그인 중복 확인
# - 폴더 생성
# - 메타데이터 생성
# - Snakefile 생성
# - DB 트랜잭션

# 3. build_plugin_docker() - 80줄
# 문제: Dockerfile 생성 + Celery 태스크 시작
# 중복: check_and_pull_official_plugin_images()와 유사한 로직
```

**리팩토링 우선순위**:
1. **High**: `upload_scripts`, `upload_plugin` 분리
2. **Medium**: Docker 관련 로직 추상화
3. **Low**: 버전 관리 엔드포인트 서비스화

#### 🟡 **Medium: 파일 처리 로직 중복**

**중복 패턴 발견**:
```python
# 패턴 1: 백업-실행-복원 (3곳에서 반복)
backup_folder = f"{target}_backup_{int(time.time())}"
shutil.copytree(target, backup_folder)
try:
    # 작업 수행
    shutil.rmtree(target)
    # 새로운 작업
    shutil.rmtree(backup_folder)  # 성공 시
except:
    shutil.copytree(backup_folder, target)  # 실패 시 복원
    shutil.rmtree(backup_folder)

# 위치:
# - upload_scripts()
# - upload_plugin()
# - upload_package()
```

**필요한 추상화**:
```python
# app/common/utils/file_transaction.py (신규 생성 필요)
class FileTransactionManager:
    def __init__(self, target_path: str):
        self.target = Path(target_path)
        self.backup = None

    def __enter__(self):
        # 백업 생성
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # 롤백 or 커밋
        pass
```

#### 🟡 **Medium: Celery 태스크 구조**

**현재 문제**:
```python
# celery_tasks.py의 MyTask 클래스
# - before_start: 55줄 (DB 조회 + 이미지 URI 생성)
# - on_success: 30줄 (캐싱 로직 포함)
# - on_failure: 16줄 (컨테이너 정리)
# - on_revoke: 33줄 (복잡한 컨테이너 패턴 매칭)
```

**개선 방향**:
- 각 lifecycle hook을 독립 함수로 분리
- 컨테이너 관리 로직을 `docker_utils.py`로 이동
- 캐싱 로직을 `cache_utils.py`로 이동

#### 🟢 **Low: CRUD 레이어 개선**

**현재 한계**:
```python
# crud_workflow.py
def get_user_workflows(db: Session, user_id: int):
    return db.query(models.Workflow).filter(
        models.Workflow.user_id == user_id
    ).all()

# 문제:
# - 단순 쿼리 래퍼
# - 페이지네이션 없음
# - 정렬 옵션 없음
# - 관계 로딩 최적화 없음 (N+1 쿼리 가능성)
```

**개선 목표**:
```python
# repositories/workflow_repository.py (신규 생성)
class WorkflowRepository(BaseRepository):
    def get_by_user(
        self,
        user_id: int,
        page: int = 1,
        limit: int = 10,
        sort_by: str = "created_at"
    ) -> Page[Workflow]:
        # 페이지네이션 + 정렬 + eager loading
        pass
```

### 3.3 에러 핸들링 패턴 분석

**발견된 패턴 (일관성 없음)**:

```python
# 패턴 A: 단순 HTTPException
if not user:
    raise HTTPException(status_code=404, detail="User not found")

# 패턴 B: try-except + HTTPException
try:
    result = operation()
except Exception as e:
    raise HTTPException(status_code=500, detail=f"Error: {str(e)}")

# 패턴 C: try-except + DB 롤백
try:
    db.commit()
except Exception as e:
    db.rollback()
    raise HTTPException(...)

# 패턴 D: 다단계 try-except (plugin.py)
try:
    try:
        # 실행
    except Exception as e_promote:
        try:
            # 롤백
        except Exception as e_rollback:
            # 롤백 실패 처리
```

**표준화 필요**:
```python
# exceptions.py (신규 생성)
class CellCraftException(Exception):
    """Base exception"""
    status_code = 500
    detail = "Internal server error"

class ResourceNotFound(CellCraftException):
    status_code = 404

class PluginBuildFailed(CellCraftException):
    status_code = 500

# 사용:
raise PluginBuildFailed(detail="Docker build failed", context={"plugin": name})
```

### 3.4 테스트 가능성 평가

| 컴포넌트 | 현재 테스트 가능성 | 주요 장애물 | 리팩토링 후 개선도 |
|----------|-------------------|-------------|------------------|
| `auth.py` | 🟢 높음 | 없음 | +10% |
| `crud_user.py` | 🟢 높음 | DB 의존성 | +20% (모킹 개선) |
| `plugin.py` | 🔴 매우 낮음 | 비즈니스 로직 혼재, 파일 I/O, Docker | +300% (서비스 분리) |
| `celery_tasks.py` | 🟡 중간 | Celery 의존성, 컨테이너 관리 | +150% (로직 분리) |
| `workflow.py` | 🟢 높음 | 없음 | +30% (검증 로직 분리) |

**테스트 불가능 코드 예시**:
```python
# plugin.py - upload_scripts() 중 일부
async def upload_scripts(...):
    # 454줄 중...

    # 파일 I/O (모킹 어려움)
    with open(file_path, "wb") as f:
        f.write(content)

    # 로깅 (부작용)
    logger.info(f"Saved {file.filename}")

    # 파일 시스템 작업 (통합 테스트 필요)
    shutil.move(str(scripts_staging_dir), str(scripts_dir))

    # DB 작업
    db_plugin.dependencies = dependencies_dict

    # 모든 것이 하나의 함수에!
```

**리팩토링 후 테스트 가능 코드**:
```python
# services/plugin_service.py
class PluginService:
    def __init__(
        self,
        plugin_repo: PluginRepository,
        file_manager: FileTransactionManager,
        logger: Logger
    ):
        # 의존성 주입

    def upload_scripts(self, plugin_name: str, files: List[UploadFile]):
        # 비즈니스 로직만 (테스트 가능)
        with self.file_manager.transaction(plugin_name):
            self.plugin_repo.update_scripts(...)

# 테스트 시:
# - plugin_repo: Mock
# - file_manager: Mock
# - logger: Mock
# -> 완전 격리된 단위 테스트 가능
```

---

## 4. 핵심 픽스처 정의 (conftest.py)

### 4.1 데이터베이스 픽스처

**⚠️ 실제 구현**: PostgreSQL test-db 사용 (원래 계획과 다름)

**변경 이유**:
- JSONB 타입이 SQLite에서 지원되지 않음
- 프로덕션 스키마 수정은 부적절
- PostgreSQL로 프로덕션 환경과 100% 동일하게 유지

**실제 구현** (`tests/conftest.py`):

```python
# tests/conftest.py
import os
import pytest
from typing import Generator
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session
from fastapi.testclient import TestClient

# IMPORTANT: Set TESTING environment variable BEFORE importing app
# This makes app/database/conn.py use PostgreSQL test DB (localhost:5433)
os.environ["TESTING"] = "1"

from app.main import app
from app.database.conn import Base
from app.database import models
from app.routes.dep import get_db
from app.common.security import get_password_hash
from app.common.enums import PluginType


@pytest.fixture(scope="session")
def test_engine():
    """
    Create PostgreSQL test database engine.

    Scope: session (shared across all tests in the session)

    Connects to: Docker test-db service (localhost:5433)
    Database: cellcraft_test
    Credentials: test_user / test_pass

    Note: Requires test-db service running:
        docker compose -f docker-compose.dev.yml up test-db -d
    """
    TEST_DATABASE_URI = "postgresql://test_user:test_pass@localhost:5433/cellcraft_test"

    engine = create_engine(
        TEST_DATABASE_URI,
        echo=False,  # Set to True for SQL debugging
        pool_pre_ping=True
    )

    # Create all tables
    Base.metadata.create_all(bind=engine)

    yield engine

    # Cleanup
    Base.metadata.drop_all(bind=engine)
    engine.dispose()


@pytest.fixture(scope="function")
def db_session(test_engine) -> Generator[Session, None, None]:
    """
    Create a new database session for each test with complete isolation.

    Scope: function (new session for each test)

    Isolation Strategy: Drop and recreate all tables before each test
    to ensure a clean state. This prevents data leakage between tests.

    Yields:
        Session: SQLAlchemy session for database operations
    """
    # Complete isolation: Drop and recreate all tables before each test
    Base.metadata.drop_all(bind=test_engine)
    Base.metadata.create_all(bind=test_engine)

    TestingSessionLocal = sessionmaker(
        autocommit=False,
        autoflush=False,
        bind=test_engine
    )

    session = TestingSessionLocal()

    yield session

    # Cleanup
    session.close()
    # Note: No rollback needed - tables will be dropped in next test


@pytest.fixture(scope="function")
def client(db_session: Session) -> Generator[TestClient, None, None]:
    """
    Create a FastAPI TestClient with test database dependency override.

    Scope: function (new client for each test)

    Args:
        db_session: Test database session

    Yields:
        TestClient: FastAPI test client for API testing
    """
    def override_get_db():
        try:
            yield db_session
        finally:
            pass  # Session cleanup handled by db_session fixture

    app.dependency_overrides[get_db] = override_get_db

    with TestClient(app) as test_client:
        yield test_client

    app.dependency_overrides.clear()
```

**Docker Compose 설정** (`docker-compose.dev.yml`):

```yaml
test-db:
  image: postgres:15
  environment:
    - TZ=Asia/Seoul
    - POSTGRES_DB=cellcraft_test
    - POSTGRES_USER=test_user
    - POSTGRES_PASSWORD=test_pass
  ports:
    - "5433:5432"  # Host 5433 -> Container 5432 (avoid conflict with db)
  networks:
    - app-network
  # No volumes - ephemeral test database
  # No restart - manual control for testing
```

**테스트 DB 실행**:
```bash
# 테스트 DB 시작
docker compose -f docker-compose.dev.yml up test-db -d

# 테스트 실행
conda activate cellcraft-test
pytest -v

# 테스트 DB 중지
docker compose -f docker-compose.dev.yml down test-db
```

**핵심 차이점**:

| 항목 | 원래 계획 (SQLite) | 실제 구현 (PostgreSQL) |
|------|-------------------|----------------------|
| 데이터베이스 | SQLite in-memory | PostgreSQL test-db |
| 포트 | N/A | localhost:5433 |
| 격리 전략 | Transaction rollback | Drop/create tables |
| JSONB 지원 | ❌ 불가 | ✅ 완전 지원 |
| 프로덕션 일치도 | 낮음 | 100% 동일 |
| 설정 복잡도 | 낮음 | 중간 (Docker 필요) |
```

### 4.2 인증 픽스처

```python
@pytest.fixture
def sample_user(db_session):
    """
    테스트 사용자 생성

    실제 models.User를 기반으로 생성
    """
    user = models.User(
        email="test@example.com",
        username="testuser",
        hashed_password=security.get_password_hash("testpassword123"),
        is_active=True,
        is_superuser=False
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    return user


@pytest.fixture
def superuser(db_session):
    """관리자 권한 사용자"""
    admin = models.User(
        email="admin@example.com",
        username="admin",
        hashed_password=security.get_password_hash("adminpassword123"),
        is_active=True,
        is_superuser=True
    )
    db_session.add(admin)
    db_session.commit()
    db_session.refresh(admin)
    return admin


@pytest.fixture
def user_token(sample_user):
    """
    JWT 액세스 토큰 생성

    실제 security.create_access_token 사용
    """
    from datetime import timedelta
    from app.common.config import settings

    access_token = security.create_access_token(
        subject=sample_user.id,
        expires_delta=timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    )
    return access_token


@pytest.fixture
def authenticated_client(client, user_token):
    """
    인증된 TestClient

    모든 요청에 Authorization 헤더 자동 추가
    """
    client.headers = {
        **client.headers,
        "Authorization": f"Bearer {user_token}",
    }
    return client


@pytest.fixture
def admin_client(client, superuser):
    """관리자 권한 TestClient"""
    from datetime import timedelta
    from app.common.config import settings

    admin_token = security.create_access_token(
        subject=superuser.id,
        expires_delta=timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    )
    client.headers = {
        **client.headers,
        "Authorization": f"Bearer {admin_token}",
    }
    return client
```

### 4.3 도메인 모델 픽스처

```python
@pytest.fixture
def sample_workflow(db_session, sample_user):
    """
    샘플 워크플로우 생성

    실제 Workflow 모델 기반
    """
    workflow_info = {
        "drawflow": {
            "Home": {
                "data": {
                    "1": {
                        "id": 1,
                        "name": "data",
                        "data": {"file_id": 1},
                        "class": "data",
                        "html": "",
                        "inputs": {},
                        "outputs": {"output_1": {"connections": [{"node": "2", "output": "input_1"}]}},
                        "pos_x": 100,
                        "pos_y": 100
                    },
                    "2": {
                        "id": 2,
                        "name": "TENET",
                        "data": {"parameters": {"FDR": 0.05}},
                        "class": "algorithm",
                        "html": "",
                        "inputs": {"input_1": {"connections": [{"node": "1", "input": "output_1"}]}},
                        "outputs": {},
                        "pos_x": 300,
                        "pos_y": 100
                    }
                }
            }
        }
    }

    workflow = models.Workflow(
        title="Test Workflow",
        thumbnail="base64_encoded_image",
        workflow_info=workflow_info,
        user_id=sample_user.id
    )
    db_session.add(workflow)
    db_session.commit()
    db_session.refresh(workflow)
    return workflow


@pytest.fixture
def sample_plugin(db_session):
    """
    샘플 플러그인 (Local)
    """
    plugin = models.Plugin(
        name="TestPlugin",
        description="Test plugin for unit testing",
        author="test_author",
        plugin_path="./plugin/local/TestPlugin/",
        plugin_type=None,
        dependencies={"requirements.txt": "numpy==1.24.0\npandas==2.0.0"},
        drawflow={
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "id": 1,
                            "name": "TestPlugin",
                            "data": {"parameters": {"param1": "value1"}},
                            "class": "algorithm",
                            "html": "",
                            "inputs": {"input_1": {"connections": []}},
                            "outputs": {"output_1": {"connections": []}},
                            "pos_x": 100,
                            "pos_y": 100
                        }
                    }
                }
            }
        },
        rules={
            "0": {
                "name": "rule_all",
                "input": ["input.h5ad"],
                "output": ["output.csv"],
                "script": "scripts/process.py",
                "parameters": {"param1": "default_value"},
                "nodeId": "1",
                "isVisualization": False
            }
        },
        use_gpu=False,
        source="local",
        is_editable=True,
        version="1.0.0"
    )
    db_session.add(plugin)
    db_session.commit()
    db_session.refresh(plugin)
    return plugin


@pytest.fixture
def official_plugin(db_session):
    """공식 플러그인 (편집 불가)"""
    plugin = models.Plugin(
        name="TENET",
        description="Official TENET plugin",
        author="CellCraft Team",
        plugin_path="./plugin/official/TENET/",
        plugin_type=None,
        dependencies=None,
        drawflow={"drawflow": {"Home": {"data": {}}}},
        rules={"0": {"name": "tenet_rule"}},
        use_gpu=True,
        source="official",
        is_editable=False,
        version="1.5.0",
        submodule_path="plugin/official/TENET"
    )
    db_session.add(plugin)
    db_session.commit()
    db_session.refresh(plugin)
    return plugin


@pytest.fixture
def sample_file(db_session, sample_user):
    """업로드된 파일 레코드"""
    file = models.File(
        file_name="test_data.h5ad",
        file_size="1024 KB",
        file_path="./user/testuser/data/test_data.h5ad",
        folder="data",
        user_id=sample_user.id
    )
    db_session.add(file)
    db_session.commit()
    db_session.refresh(file)
    return file


@pytest.fixture
def sample_task(db_session, sample_user, sample_workflow, sample_plugin):
    """실행 태스크 레코드"""
    from datetime import datetime

    task = models.Task(
        task_id="test-task-id-12345",
        start_time=datetime.now(),
        end_time=None,
        status="PENDING",
        user_id=sample_user.id,
        workflow_id=sample_workflow.id,
        algorithm_id="algorithm_1",
        plugin_name=sample_plugin.name,
        plugin_id=sample_plugin.id,
        task_type="compile",
        plugin_image_uri=f"ghcr.io/cxinsys/cellcraft-{sample_plugin.name.lower()}:1.0.0"
    )
    db_session.add(task)
    db_session.commit()
    db_session.refresh(task)
    return task
```

### 4.4 Celery 및 외부 의존성 모킹

```python
@pytest.fixture(autouse=True)
def mock_celery(monkeypatch):
    """
    Celery EAGER 모드 설정

    테스트 시 Celery 태스크를 동기적으로 실행
    RabbitMQ 없이 테스트 가능
    """
    monkeypatch.setenv("CELERY_ALWAYS_EAGER", "True")
    monkeypatch.setenv("CELERY_EAGER_PROPAGATES", "True")


@pytest.fixture
def mock_docker(mocker):
    """
    Docker 클라이언트 모킹

    실제 Docker 없이 테스트 가능
    """
    mock_client = mocker.patch("docker.from_env")
    mock_images = mocker.MagicMock()
    mock_containers = mocker.MagicMock()

    mock_client.return_value.images = mock_images
    mock_client.return_value.containers = mock_containers

    return {
        "client": mock_client,
        "images": mock_images,
        "containers": mock_containers
    }


@pytest.fixture
def mock_snakemake(mocker):
    """Snakemake 실행 모킹"""
    mock_snakemake = mocker.patch("app.common.utils.snakemake_utils.snakemakeProcess")
    mock_snakemake.return_value = {
        "returncode": 0,
        "stdout": "Snakemake execution completed",
        "stderr": ""
    }
    return mock_snakemake


@pytest.fixture
def tmp_plugin_dir(tmp_path):
    """
    임시 플러그인 디렉토리

    파일 시스템 테스트용
    """
    plugin_dir = tmp_path / "plugin" / "local" / "TestPlugin"
    plugin_dir.mkdir(parents=True)

    # 기본 폴더 구조 생성
    (plugin_dir / "scripts").mkdir()
    (plugin_dir / "dependency").mkdir()

    # 메타데이터 파일 생성
    import json
    metadata = {
        "name": "TestPlugin",
        "description": "Test",
        "drawflow": {},
        "rules": {}
    }
    (plugin_dir / "metadata.json").write_text(json.dumps(metadata))

    return plugin_dir
```

### 4.5 테스트 데이터 픽스처

```python
@pytest.fixture
def sample_h5ad_file(tmp_path):
    """
    임시 H5AD 파일 생성

    실제 AnnData 객체는 무겁기 때문에 모킹
    """
    import h5py
    import numpy as np

    h5ad_path = tmp_path / "test_data.h5ad"

    with h5py.File(h5ad_path, "w") as f:
        # 기본 구조만 생성
        f.create_dataset("obs", data=np.array([]))
        f.create_dataset("var", data=np.array([]))
        f.create_dataset("X", data=np.array([[]]))

    return h5ad_path


@pytest.fixture
def sample_workflow_json():
    """워크플로우 JSON 데이터"""
    return {
        "title": "Test Workflow",
        "thumbnail": "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg==",
        "workflow_info": {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "id": 1,
                            "name": "data",
                            "data": {"file_id": 1},
                            "class": "data",
                            "inputs": {},
                            "outputs": {"output_1": {"connections": [{"node": "2", "output": "input_1"}]}}
                        }
                    }
                }
            }
        }
    }


@pytest.fixture
def sample_plugin_metadata():
    """플러그인 메타데이터"""
    return {
        "name": "TestPlugin",
        "description": "A test plugin",
        "author": "test_author",
        "plugin_type": "analysis",
        "dependencies": {
            "requirements.txt": "numpy>=1.24.0\npandas>=2.0.0\nscipy>=1.10.0"
        },
        "drawflow": {
            "drawflow": {
                "Home": {
                    "data": {
                        "1": {
                            "id": 1,
                            "name": "TestPlugin",
                            "class": "algorithm",
                            "data": {"parameters": {"threshold": 0.05}},
                            "inputs": {"input_1": {"connections": []}},
                            "outputs": {"output_1": {"connections": []}}
                        }
                    }
                }
            }
        },
        "rules": {
            "0": {
                "name": "process_data",
                "input": ["{input_file}"],
                "output": ["{output_file}"],
                "script": "scripts/process.py",
                "parameters": {"threshold": 0.05},
                "nodeId": "1",
                "isVisualization": False
            }
        },
        "use_gpu": False
    }
```

---

## 4. 테스트 & 리팩토링 통합 워크플로우

### Phase 1: 기본 인프라 + 코드 구조 개선 (Week 1-2) 🔴

**테스트 작업**
```
✓ conftest.py 픽스처 구현
✓ 테스트 DB 설정
✓ TestClient 설정
✓ CI/CD 기본 설정
```

**리팩토링 작업**
```
✓ 의존성 주입(DI) 패턴 도입
  - app/dependencies.py 생성
  - get_db(), get_current_user() 등 의존성 함수 정의
  
✓ Repository 패턴 적용 준비
  - app/repositories/ 디렉토리 생성
  - 기본 Base Repository 클래스 정의
```

### Phase 2: 인증 & 모델 + 계층 분리 (Week 3-4) 🟡

**테스트 작업**
```
✓ test_auth.py - 로그인, 회원가입, JWT
✓ test_models.py - SQLAlchemy 모델 검증
✓ test_schemas.py - Pydantic 스키마 검증
```

**리팩토링 작업**
```
✓ 인증 로직 리팩토링
  - API 라우터와 비즈니스 로직 분리
  - services/auth_service.py 생성
  - JWT 토큰 생성/검증 유틸리티 분리
  
✓ 중복 코드 제거
  - 공통 예외 핸들러 통합
  - 공통 응답 스키마 정의
```

### Phase 3: 핵심 API + 서비스 레이어 강화 (Week 5-8) 🟢

**테스트 작업**
```
✓ test_workflows.py - 워크플로우 CRUD
✓ test_plugins.py - 플러그인 등록/검증
✓ test_files.py - 파일 업로드/다운로드
✓ test_tasks.py - 태스크 실행
```

**리팩토링 작업**
```
✓ Repository 패턴 완전 적용
  - WorkflowRepository, PluginRepository 구현
  - 데이터 접근 로직을 API 라우터에서 분리
  
✓ 서비스 레이어 강화
  - WorkflowService: 복잡한 비즈니스 로직 처리
  - ValidationService: 워크플로우/플러그인 검증
  - FileService: 파일 처리 로직 캡슐화
  
✓ 파일 처리 개선
  - 파일 저장 로직 추상화
  - 스토리지 전략 패턴 적용 (로컬/S3 대응)
  
✓ 에러 핸들링 일관성
  - 커스텀 Exception 클래스 정의
  - 전역 예외 핸들러 개선
```

### Phase 4: 통합 테스트 + 성능 최적화 (Week 9-12) 🔵

**테스트 작업**
```
✓ test_workflow_execution.py - E2E 워크플로우
✓ test_e2e_scenarios.py - 사용자 시나리오
✓ 성능 테스트 추가
```

**리팩토링 작업**
```
✓ Celery 태스크 구조 개선
  - 태스크 재시도 로직 추가
  - 태스크 상태 추적 개선
  - 에러 복구 메커니즘 강화
  
✓ 데이터베이스 쿼리 최적화
  - N+1 쿼리 문제 해결
  - Eager Loading 적용
  - 인덱스 추가
  
✓ API 응답 최적화
  - 페이지네이션 개선
  - 응답 캐싱 전략 적용 (선택)
```

---

## 5. 리팩토링 체크리스트

### 5.1 코드 구조 개선

#### ✅ 의존성 주입 (Dependency Injection)
```python
# Before (리팩토링 전)
@router.post("/workflows")
def create_workflow(workflow: WorkflowCreate, db: Session = Depends(get_db)):
    # DB 직접 조작
    db_workflow = Workflow(**workflow.dict())
    db.add(db_workflow)
    db.commit()

# After (리팩토링 후)
@router.post("/workflows")
def create_workflow(
    workflow: WorkflowCreate,
    workflow_service: WorkflowService = Depends(get_workflow_service)
):
    # 서비스 레이어로 위임
    return workflow_service.create(workflow)
```

#### ✅ Repository 패턴
```python
# app/repositories/base.py
class BaseRepository:
    def __init__(self, model, db: Session):
        self.model = model
        self.db = db
    
    def get_by_id(self, id: int):
        return self.db.query(self.model).filter(self.model.id == id).first()
    
    def create(self, data: dict):
        instance = self.model(**data)
        self.db.add(instance)
        self.db.commit()
        return instance

# app/repositories/workflow_repository.py
class WorkflowRepository(BaseRepository):
    def __init__(self, db: Session):
        super().__init__(Workflow, db)
    
    def get_by_user_id(self, user_id: int):
        return self.db.query(self.model).filter(
            self.model.user_id == user_id
        ).all()
```

#### ✅ 서비스 레이어
```python
# app/services/workflow_service.py
class WorkflowService:
    def __init__(
        self,
        workflow_repo: WorkflowRepository,
        validation_service: ValidationService
    ):
        self.workflow_repo = workflow_repo
        self.validation_service = validation_service
    
    def create(self, workflow_data: WorkflowCreate, user_id: int):
        # 비즈니스 로직
        self.validation_service.validate_workflow(workflow_data)
        
        # 데이터 저장
        return self.workflow_repo.create({
            **workflow_data.dict(),
            "user_id": user_id
        })
```

### 5.2 테스트 가능성 개선

| 개선 항목 | Before | After | 효과 |
|---------|--------|-------|------|
| 하드코딩된 의존성 | DB 직접 접근 | DI로 주입 | 모킹 용이 |
| 복잡한 함수 | 200줄 이상 | 20-30줄 단위 | 단위 테스트 가능 |
| 전역 상태 | 전역 변수 사용 | 함수 파라미터 | 테스트 격리 |
| 파일 I/O | 직접 파일 쓰기 | 추상화 계층 | 모킹 가능 |

### 5.3 코드 품질 지표

```python
# 리팩토링 전후 비교 목표
┌─────────────────────┬─────────┬─────────┐
│ 지표                │ Before  │ After   │
├─────────────────────┼─────────┼─────────┤
│ 평균 함수 길이      │ 150줄   │ 30줄    │
│ 순환 복잡도         │ 15-20   │ 5-10    │
│ 코드 중복률         │ 20%     │ <5%     │
│ 테스트 커버리지     │ 0%      │ 60%+    │
│ Import 의존성       │ 순환참조│ 단방향  │
└─────────────────────┴─────────┴─────────┘
```

---

## 6. 테스트 작성 가이드라인

### 6.1 명명 규칙
```python
# 함수명: test_{기능}_{시나리오}_{예상결과}
test_create_workflow_with_valid_data_returns_201()
test_login_with_invalid_credentials_returns_401()
test_upload_file_exceeds_size_limit_returns_413()
```

### 6.2 테스트 패턴 (AAA Pattern)
```python
def test_example(authenticated_client, sample_workflow):
    # Arrange - 테스트 데이터 준비
    workflow_data = {...}
    
    # Act - 실행
    response = authenticated_client.post("/api/workflows", json=workflow_data)
    
    # Assert - 검증
    assert response.status_code == 201
    assert response.json()["name"] == workflow_data["name"]
```

### 6.3 핵심 검증 항목
```python
# HTTP 응답 검증
✓ status_code
✓ response JSON schema
✓ error messages

# 데이터베이스 검증
✓ 레코드 생성/수정/삭제
✓ 관계(Relationship) 무결성

# 비즈니스 로직 검증
✓ 워크플로우 검증 로직
✓ 파일 포맷 검증
✓ 권한 검사

# 리팩토링 검증
✓ 서비스 레이어 호출 확인
✓ Repository 메서드 호출 확인
✓ 의존성 주입 정상 작동
```

---

## 7. 리팩토링 발견 항목 추적

### 7.1 리팩토링 이슈 마킹
```python
# 테스트 작성 중 발견된 리팩토링 필요 항목 마킹
@pytest.mark.refactor
def test_complex_workflow_validation():
    """
    TODO: workflow validation 로직이 너무 복잡함
    - ValidationService로 분리 필요
    - 단일 책임 원칙 위반
    """
    pass
```

### 7.2 리팩토링 백로그 관리
```markdown
# docs/refactoring-backlog.md

## High Priority 🔴
- [ ] API 라우터에서 비즈니스 로직 분리 (Week 3-4)
- [ ] Repository 패턴 적용 (Week 5-6)
- [ ] 에러 핸들링 일관성 개선 (Week 5)

## Medium Priority 🟡
- [ ] 중복 코드 제거 (Week 6-7)
- [ ] 파일 처리 로직 추상화 (Week 7)
- [ ] Celery 태스크 재시도 로직 (Week 8)

## Low Priority 🟢
- [ ] 캐싱 전략 적용 (Week 10+)
- [ ] 성능 최적화 (Week 11+)
```

### 7.3 Phase 2.1 테스트 중 발견된 버그 및 개선사항

#### 🔴 Critical (P0) - 즉시 수정 필요

**BUG-001: Pydantic Validation Gap - None Values Bypass Validation**
- **발견 위치**: `test_auth.py::TestUserRegistration::test_register_missing_required_fields`
- **증상**: None 값이 required fields에 대한 Pydantic validation을 우회하여 데이터베이스 레벨에서 IntegrityError 발생
- **기대 동작**: 422 Unprocessable Entity (Pydantic validation error)
- **실제 동작**: 500 Internal Server Error (Database IntegrityError)
- **영향 범위**: User registration endpoint (`/routes/auth/register`)
- **보안 영향**: API validation layer 우회 가능, 잘못된 에러 응답으로 내부 구조 노출
- **수정 방법**:
  - `app/database/schemas/user.py`의 `UserCreate` 스키마에 필드 validator 추가
  - 또는 Pydantic v2의 `Field(...)` 설정으로 None 불허 명시
- **수정 대상 파일**: `app/database/schemas/user.py`
- **테스트 파일**: `tests/unit/test_auth.py:105-138` (현재 @pytest.mark.xfail)
- **예상 수정 기간**: 0.5일
- **Phase 2.3에서 수정 예정**

```python
# 현재 (문제)
class UserCreate(BaseModel):
    username: str
    email: EmailStr
    password: str

# 수정 후 (권장)
from pydantic import Field, validator

class UserCreate(BaseModel):
    username: str = Field(..., min_length=1)
    email: EmailStr = Field(...)
    password: str = Field(..., min_length=8)

    @validator('username', 'email', 'password', pre=True)
    def prevent_none(cls, v):
        if v is None:
            raise ValueError('Field cannot be None')
        return v
```

#### 🟡 High Priority (P1) - Phase 2.3에서 수정

**BUG-002: Empty String Validation Gap**
- **발견 위치**: `test_auth.py::TestUserRegistration::test_register_empty_fields`
- **증상**: 빈 문자열("")이 username 필드에 허용됨
- **기대 동작**: 422 validation error로 거부
- **실제 동작**: 200 OK로 빈 username을 가진 사용자 생성
- **영향 범위**: User registration, data quality
- **데이터 무결성**: 빈 username으로 사용자 생성 가능
- **수정 방법**: `UserCreate` 스키마에 `min_length=1` constraint 추가
- **수정 대상 파일**: `app/database/schemas/user.py`
- **테스트 파일**: `tests/unit/test_auth.py:140-166` (현재 @pytest.mark.xfail)
- **예상 수정 기간**: 0.25일
- **Phase 2.3에서 수정 예정**

```python
# 수정 방법
class UserCreate(BaseModel):
    username: str = Field(..., min_length=1, max_length=50)
    email: EmailStr = Field(...)
    password: str = Field(..., min_length=8)
```

#### 🟢 Medium Priority (P2) - 보안 강화 권장

**LIMITATION-001: JWT Token Uniqueness Issue**
- **발견 위치**: `test_auth.py::TestUserLogin::test_concurrent_logins_generate_different_tokens`
- **증상**: 동일 초(second) 내 생성된 JWT 토큰이 완전히 동일함
- **원인**: JWT payload에 'iat' (issued at) claim이 없고, 'exp'만 초 단위 정밀도 사용
- **기대 동작**: 각 로그인마다 고유한 토큰 생성
- **실제 동작**: 1초 이내 연속 로그인 시 동일한 토큰 생성
- **보안 영향**:
  - 토큰 재사용 가능
  - 토큰 revocation 어려움
  - 세션 관리 복잡성 증가
- **영향 범위**: All authenticated endpoints
- **수정 방법 (3가지 옵션)**:
  1. **'iat' claim 추가** (권장): 마이크로초 정밀도로 발급 시간 기록
  2. **'jti' claim 추가**: UUID로 각 토큰 고유 식별
  3. **Session Management**: 데이터베이스 기반 세션 추적
- **수정 대상 파일**: `app/common/security.py::create_access_token()`
- **테스트 파일**: `tests/unit/test_auth.py:300-354` (현재 @pytest.mark.xfail)
- **예상 수정 기간**: 1일 (토큰 구조 변경 + 기존 토큰 호환성 고려)
- **Phase 3에서 검토 예정**

```python
# 현재 create_access_token()
def create_access_token(subject: Union[str, Any], expires_delta: timedelta = None) -> str:
    expire = datetime.utcnow() + (expires_delta or timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES))
    to_encode = {"exp": expire, "sub": str(subject)}
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt

# 수정 후 (Option 1: iat claim 추가)
def create_access_token(subject: Union[str, Any], expires_delta: timedelta = None) -> str:
    now = datetime.utcnow()
    expire = now + (expires_delta or timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES))
    to_encode = {
        "exp": expire,
        "iat": now,  # 발급 시간 추가
        "sub": str(subject)
    }
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt

# 수정 후 (Option 2: jti claim 추가 - 더 강력)
import uuid

def create_access_token(subject: Union[str, Any], expires_delta: timedelta = None) -> str:
    now = datetime.utcnow()
    expire = now + (expires_delta or timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES))
    to_encode = {
        "exp": expire,
        "iat": now,
        "jti": str(uuid.uuid4()),  # 고유 토큰 ID
        "sub": str(subject)
    }
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt
```

#### 📊 발견 통계 (Phase 2.1)

- **총 발견 이슈**: 3건
- **Critical (P0)**: 1건 - Validation bypass
- **High (P1)**: 1건 - Empty string acceptance
- **Medium (P2)**: 1건 - Token uniqueness
- **테스트로 문서화된 이슈**: 3건 (@pytest.mark.xfail)
- **즉시 수정 필요**: 2건 (P0, P1)

#### 🔄 이슈 추적

모든 발견된 버그는 다음과 같이 추적됩니다:
1. **테스트 코드**: `@pytest.mark.xfail`로 마킹 (실패가 예상됨을 명시)
2. **문서화**: 본 섹션에 상세 기록
3. **수정 계획**: 해당 Phase에 수정 작업 포함
4. **검증**: 수정 후 xfail 마크 제거 및 테스트 통과 확인

---

## 8. 테스트 실행 전략

### 로컬 개발
```bash
# 전체 테스트
pytest

# 리팩토링 필요 항목만 확인
pytest -m refactor

# 특정 파일/마커
pytest tests/api/test_auth.py
pytest -m "unit and not slow"

# 커버리지 확인
pytest --cov=app --cov-report=html
open htmlcov/index.html

# 병렬 실행
pytest -n auto
```

### 코드 품질 체크
```bash
# 리팩토링 후 코드 품질 검증
black app/ tests/                 # 포맷팅
isort app/ tests/                 # import 정렬
pylint app/                       # 정적 분석
pytest --cov=app --cov-fail-under=60  # 커버리지 최소 60%
```

### CI/CD (GitHub Actions)
```yaml
# .github/workflows/test.yml
name: Backend Tests & Code Quality
on: [push, pull_request]
jobs:
  test:
    steps:
      - Checkout
      - Setup Python 3.10
      - Install dependencies
      - Run black check
      - Run pylint
      - Run pytest with coverage
      - Upload coverage
      - Comment coverage on PR
```

---

## 9. 실제 엔드포인트 기반 테스트 케이스 상세 계획

> 이 섹션은 실제 백엔드 코드베이스를 분석하여 구체적인 API 엔드포인트와 테스트 시나리오를 문서화합니다.

### 9.1 인증 API (test_auth.py) - 우선순위: P0

**파일**: `app/routes/endpoints/auth.py` (92줄 - 단순)
**복잡도**: 🟢 낮음
**리팩토링 필요도**: 낮음

**엔드포인트 목록**:
- ✓ `POST /auth/register` - 회원가입
- ✓ `POST /auth/login/access-token` - 로그인
- ✓ `GET /auth/me` - 사용자 정보
- ✓ `POST /auth/plugins` - 플러그인 연결

**테스트 시나리오** (13개):

| No | 테스트 케이스 | HTTP Method | 예상 결과 | 우선순위 | 검증 항목 | 테스트 유형 |
|----|------------|-------------|----------|---------|----------|-----------|
| 1 | 유효한 정보로 회원가입 | POST /auth/register | 201 | P0 | User 생성, 폴더 생성 | 통합 |
| 2 | 중복 이메일로 회원가입 | POST /auth/register | 400 | P0 | 에러 메시지 | 단위 |
| 3 | 유효한 자격증명으로 로그인 | POST /auth/login/access-token | 200 + JWT | P0 | 토큰 구조 | 통합 |
| 4 | 잘못된 비밀번호로 로그인 | POST /auth/login/access-token | 400 | P0 | 에러 메시지 | 단위 |
| 5 | 존재하지 않는 이메일 로그인 | POST /auth/login/access-token | 400 | P0 | 에러 메시지 | 단위 |
| 6 | 비활성 사용자 로그인 | POST /auth/login/access-token | 400 | P1 | "please login to activate" | 단위 |
| 7 | 유효한 토큰으로 사용자 정보 조회 | GET /auth/me | 200 | P0 | email, username | 통합 |
| 8 | 만료된 토큰으로 접근 | GET /auth/me | 403 | P1 | "Could not validate credentials" | 단위 |
| 9 | 토큰 없이 접근 | GET /auth/me | 403 | P1 | "Not authenticated" | 단위 |

**예상 커버리지**: 95% (단순한 CRUD)

---

### 9.2 플러그인 API (test_plugins.py) - 우선순위: P0 (Critical)

**파일**: `app/routes/endpoints/plugin.py` (1576줄 - 매우 복잡!)
**복잡도**: 🔴 매우 높음
**리팩토링 필요도**: 🔴 Critical

**주요 엔드포인트** (15개):
- ⚠️ `POST /plugin/upload_scripts` - **454줄 함수!**
- ⚠️ `POST /plugin/upload` - 287줄
- ✓ `POST /plugin/validation`
- ✓ `POST /plugin/upload_package`
- ✓ `POST /plugin/build_docker/{name}`
- ✓ `GET /plugin/list`
- ✓ `GET /plugin/info/{name}`
- ✓ `POST /plugin/associate`
- ✓ `POST /plugin/dissociate`
- ✓ `GET /plugin/template/{id}`
- ✓ `GET /plugin/versions/{name}` (official 전용)
- ✓ `POST /plugin/versions/{name}/update`

**테스트 시나리오** (25개 이상):

#### A. 플러그인 생명주기 테스트

| No | 테스트 케이스 | HTTP Method | 예상 결과 | 우선순위 | 리팩토링 | 비고 |
|----|------------|-------------|----------|---------|----------|------|
| 1 | 플러그인 메타데이터 검증 (유효) | POST /plugin/validation | 200 | P0 | 🟢 Low | 단순 검증 |
| 2 | 스크립트 없이 검증 | POST /plugin/validation | 400 | P0 | 🟢 Low | 에러 처리 |
| 3 | 새 플러그인 업로드 | POST /plugin/upload | 200 | P0 | 🔴 Critical | **287줄 함수** |
| 4 | 기존 플러그인 덮어쓰기 | POST /plugin/upload | 200 | P1 | 🔴 Critical | 백업/롤백 로직 |
| 5 | Official 플러그인 수정 시도 | POST /plugin/upload | 403 | P0 | 🟢 Low | 권한 체크 |

#### B. 스크립트 업로드 테스트 (Critical)

| No | 테스트 케이스 | HTTP Method | 예상 결과 | 우선순위 | 리팩토링 | 비고 |
|----|------------|-------------|----------|---------|----------|------|
| 6 | 스크립트 파일 업로드 성공 | POST /plugin/upload_scripts | 200 | P0 | 🔴 Critical | **454줄 함수!** |
| 7 | 스크립트 업로드 중 실패 시 롤백 | POST /plugin/upload_scripts | 500 | P1 | 🔴 Critical | 복잡한 롤백 |
| 8 | Official 플러그인 스크립트 업로드 시도 | POST /plugin/upload_scripts | 403 | P0 | 🟢 Low | 권한 체크 |
| 9 | 스테이징 디렉토리 정리 실패 처리 | POST /plugin/upload_scripts | 500 | P1 | 🔴 Critical | 에러 복구 |
| 10 | Reference 폴더 복원 | POST /plugin/upload_scripts | 200 | P1 | 🟡 Medium | 로직 복잡 |

**리팩토링 예상 효과**:
- upload_scripts 함수: 454줄 → 30줄 (93% 감소)
- 테스트 가능성: 20% → 90% (350% 증가)
- 순환 복잡도: 25+ → 5 (80% 감소)

---

### 9.3 워크플로우 API (test_workflows.py) - 우선순위: P0

**파일**: `app/routes/endpoints/workflow.py` (추정)
**복잡도**: 🟢 낮음
**리팩토링 필요도**: 중간

**테스트 시나리오** (8개):

| No | 테스트 케이스 | HTTP Method | 예상 결과 | 우선순위 | 검증 항목 |
|----|------------|-------------|----------|---------|----------|
| 1 | 워크플로우 생성 | POST /workflow | 201 | P0 | DB 저장, drawflow 구조 |
| 2 | 빈 제목으로 생성 | POST /workflow | 400 | P0 | Pydantic 검증 |
| 3 | 잘못된 drawflow 구조 | POST /workflow | 400 | P1 | ValidationService 필요 |
| 4 | 사용자 워크플로우 목록 | GET /workflow | 200 | P0 | 본인 것만 조회 |
| 5 | 타인 워크플로우 조회 | GET /workflow/{id} | 404 | P1 | 권한 체크 |
| 6 | 워크플로우 수정 | PUT /workflow/{id} | 200 | P0 | title, workflow_info 변경 |
| 7 | 타인 워크플로우 수정 | PUT /workflow/{id} | 403 | P1 | 권한 체크 |
| 8 | 워크플로우 삭제 | DELETE /workflow/{id} | 200 | P0 | DB 삭제, Cascade 확인 |

---

### 9.4 태스크 API (test_tasks.py) - 우선순위: P0

**파일**: `app/routes/endpoints/task.py`
**복잡도**: 🟡 중간
**리팩토링 필요도**: 중간

**테스트 시나리오** (10개):

| No | 테스트 케이스 | 우선순위 | 검증 항목 | 리팩토링 포인트 |
|----|------------|---------|----------|---------------|
| 1 | 워크플로우 실행 요청 | P0 | Celery task_id 반환 | TaskService 분리 |
| 2 | 태스크 상태 조회 (PENDING) | P0 | 상태 정확성 | - |
| 3 | 태스크 상태 조회 (RUNNING) | P0 | 진행률 정보 | - |
| 4 | 태스크 상태 조회 (SUCCESS) | P0 | 결과 파일 경로 | - |
| 5 | 태스크 상태 조회 (FAILURE) | P1 | 에러 메시지 | 에러 핸들링 |
| 6 | 태스크 취소 요청 | P1 | Celery revoke 호출 | 취소 로직 개선 |
| 7 | 취소 시 컨테이너 정리 | P1 | Docker 컨테이너 종료 | 컨테이너 관리 |
| 8 | 사용자 태스크 목록 | P0 | 본인 것만 조회 | 페이지네이션 추가 |

---

### 9.5 Celery 태스크 테스트 (test_celery_tasks.py) - 우선순위: P1

**파일**: `app/routes/celery_tasks.py` (200+줄)
**복잡도**: 🔴 높음
**리팩토링 필요도**: 🔴 High

**테스트 시나리오** (8개):

| No | 테스트 케이스 | 우선순위 | 검증 항목 | 리팩토링 포인트 |
|----|------------|---------|----------|---------------|
| 1 | process_data_task 정상 실행 | P0 | Snakemake 호출, 결과 파일 생성 | 로직 분리 |
| 2 | process_data_task 실패 시 컨테이너 정리 | P1 | Docker container 정리 | 컨테이너 관리 개선 |
| 3 | build_plugin_task 성공 | P0 | Docker 이미지 생성 | 빌드 로직 분리 |
| 4 | build_plugin_task 실패 시 로그 | P1 | 빌드 로그 파일 생성 | 로깅 개선 |
| 5 | before_start hook | P1 | DB 태스크 레코드 생성 | Hook 로직 분리 |
| 6 | on_success hook | P0 | DB 상태 업데이트 | - |
| 7 | on_failure hook | P1 | 컨테이너 정리, 상태 업데이트 | - |
| 8 | on_revoke hook | P1 | 컨테이너 강제 종료 | 패턴 매칭 개선 |

---

### 9.6 테스트 우선순위 요약

| API 그룹 | 테스트 수 | P0 테스트 | P1 테스트 | 예상 소요 |
|---------|----------|----------|----------|----------|
| **인증 (auth)** | 9 | 7 | 2 | 2-3일 |
| **플러그인 (plugin)** | 25+ | 15 | 10+ | **10-15일** (Critical) |
| **워크플로우 (workflow)** | 8 | 5 | 3 | 3-4일 |
| **태스크 (task)** | 10 | 6 | 4 | 4-5일 |
| **Celery 태스크** | 8 | 3 | 5 | 5-7일 |
| **파일 (files)** | 7 | 4 | 3 | 2-3일 |
| **합계** | **67+** | **40** | **27+** | **26-40일** |

**Critical Path**: 플러그인 API 테스트 + 리팩토링 (40% 시간 소요 예상)

---

## 10. 코딩 에이전트 작업 체크리스트

### 즉시 착수 항목 (Week 1)
- [ ] `conftest.py` 기본 픽스처 구현
- [ ] `pytest.ini` 설정 파일 작성
- [ ] `app/dependencies.py` DI 컨테이너 생성
- [ ] `app/repositories/base.py` 기본 Repository 구현
- [ ] GitHub Actions 워크플로우 기본 구성

### 테스트 & 리팩토링 병행 작업
```
Week 1-2:
  테스트: conftest.py, pytest.ini, CI 설정
  리팩토링: DI 패턴 도입, Base Repository 생성

Week 3-4:
  테스트: test_auth.py, test_models.py
  리팩토링: AuthService 분리, 중복 코드 제거

Week 5-6:
  테스트: test_workflows.py, test_plugins.py
  리팩토링: WorkflowRepository, ValidationService

Week 7-8:
  테스트: test_files.py, test_tasks.py
  리팩토링: FileService, Celery 태스크 개선

Week 9-12:
  테스트: 통합 테스트, E2E 테스트
  리팩토링: 성능 최적화, 쿼리 최적화
```

---

## 11. 성공 지표 (KPI)

### 테스트 지표
```
✓ 테스트 커버리지: 60% → 75% → 85%
✓ 테스트 실행 시간: < 5분 (전체)
✓ CI/CD 통과율: 95% 이상
✓ P0 테스트: 100% 완료
```

### 리팩토링 지표
```
✓ 평균 함수 복잡도: <10 (McCabe)
✓ 코드 중복률: <5%
✓ Pylint 점수: 8.0 이상
✓ 순환 의존성: 0개
✓ 테스트 가능한 코드 비율: 90% 이상
```

### 품질 개선 지표
```
✓ Repository 패턴 적용: 100% (CRUD 작업)
✓ Service 레이어 적용: 100% (비즈니스 로직)
✓ DI 적용: 100% (API 엔드포인트)
✓ 커스텀 Exception: 모든 에러 케이스 커버
```

---

## 12. 위험 요소 & 대응 방안

| 위험 요소 | 영향도 | 대응 방안 |
|----------|--------|----------|
| 대규모 리팩토링으로 기존 기능 영향 | 높음 | 단계적 리팩토링, 기존 테스트 유지 |
| 리팩토링 시간 부족 | 중간 | 우선순위 명확화, P0만 먼저 처리 |
| 팀원 간 리팩토링 방향 불일치 | 중간 | 아키텍처 문서화, 코드 리뷰 강화 |
| 테스트 작성 시간 예상 초과 | 낮음 | 병렬 작업, 페어 프로그래밍 |

---

## 13. 참고 문서

```
공식 문서:
- FastAPI Testing: https://fastapi.tiangolo.com/tutorial/testing/
- pytest Documentation: https://docs.pytest.org/
- Clean Architecture: https://blog.cleancoder.com/

내부 문서 (작성 예정):
- tests/README.md - 테스트 실행 가이드
- tests/CONTRIBUTING.md - 테스트 작성 가이드라인
- docs/architecture.md - 리팩토링 후 아키텍처
- docs/refactoring-guide.md - 리팩토링 가이드
- docs/refactoring-backlog.md - 리팩토링 백로그
```

---

## 14. 최종 아키텍처 비전

```
리팩토링 완료 후 목표 구조:

app/
├── api/                    # API 라우터 (얇은 계층)
│   └── endpoints/
│       ├── auth.py        # 라우팅만 담당
│       ├── workflows.py
│       └── tasks.py
│
├── services/              # 비즈니스 로직
│   ├── auth_service.py
│   ├── workflow_service.py
│   └── validation_service.py
│
├── repositories/          # 데이터 접근
│   ├── base.py
│   ├── workflow_repository.py
│   └── user_repository.py
│
├── models/                # DB 모델
├── schemas/               # Pydantic 스키마
├── dependencies.py        # DI 컨테이너
└── exceptions.py          # 커스텀 Exception

레이어 의존성:
API → Service → Repository → Model
(단방향, 역의존성 없음)
```

---

**다음 단계**:
1. 코딩 에이전트에게 `conftest.py` + `dependencies.py` 동시 작성 할당
2. 첫 테스트 작성과 동시에 DI 패턴 적용 시작
3. 매주 리팩토링 백로그 리뷰 및 우선순위 조정

---

## 15. 실행 계획 요약 (Executive Summary)

### 프로젝트 현황

**코드베이스 규모**:
- 총 라인 수: ~2,500줄 (주요 백엔드 코드)
- 복잡도 High 파일: 2개 (plugin.py 1576줄, celery_tasks.py 200+줄)
- 현재 테스트 커버리지: **0%**

**주요 문제점**:
1. **🔴 Critical**: plugin.py (1576줄) - 비즈니스 로직이 라우터에 집중
   - `upload_scripts()`: 454줄 단일 함수
   - `upload_plugin()`: 287줄 단일 함수
   - 테스트 불가능한 구조

2. **🔴 High**: 파일 처리 로직 중복
   - 백업/롤백 패턴이 3곳에서 반복
   - 에러 핸들링 일관성 부족

3. **🟡 Medium**: CRUD 레이어의 한계
   - 단순 쿼리 래퍼만 제공
   - Repository 패턴 미적용

### 목표

**테스트 커버리지 목표**:
- Phase 1 (Week 1-4): 0% → 30% (인증, 기본 CRUD)
- Phase 2 (Week 5-8): 30% → 60% (핵심 API, 리팩토링 병행)
- Phase 3 (Week 9-12): 60% → 75% (통합 테스트, 성능 최적화)
- Phase 4 (Week 13-16): 75% → 85% (엣지 케이스, E2E)

**코드 품질 목표**:
- 평균 함수 길이: 150줄 → 30줄
- 순환 복잡도: 15-20 → 5-10
- 코드 중복률: 20% → <5%
- 테스트 가능 코드 비율: 20% → 90%

### 핵심 작업 항목

#### Week 1-2: 인프라 구축
- [ ] `conftest.py` 픽스처 구현 (2일)
- [ ] `pytest.ini` 설정 (1일)
- [ ] `app/dependencies.py` DI 컨테이너 (2일)
- [ ] `app/exceptions.py` 커스텀 예외 (1일)
- [ ] GitHub Actions CI 기본 설정 (1일)

#### Week 3-4: 인증 & 기본 테스트
- [ ] `test_auth.py` 작성 (3일)
- [ ] `test_models.py` 작성 (2일)
- [ ] `test_schemas.py` 작성 (1일)
- [ ] AuthService 분리 (선택, 1일)

#### Week 5-8: 플러그인 API 리팩토링 (Critical)
- [ ] `FileTransactionManager` 유틸리티 작성 (3일)
- [ ] `PluginService` 생성 및 로직 분리 (7일)
- [ ] `PluginRepository` 구현 (3일)
- [ ] `test_plugins.py` 작성 (7일)
- [ ] upload_scripts, upload_plugin 리팩토링 (5일)

**예상 효과**:
- plugin.py: 1576줄 → ~500줄 (68% 감소)
- 테스트 커버리지: plugin API 85%+

#### Week 9-12: 통합 테스트 & 성능
- [ ] `test_workflow_execution.py` E2E (5일)
- [ ] `test_celery_tasks.py` (4일)
- [ ] 데이터베이스 쿼리 최적화 (N+1 해결) (3일)
- [ ] 페이지네이션 추가 (2일)
- [ ] 성능 테스트 추가 (2일)

### 리스크 & 대응

| 리스크 | 확률 | 영향도 | 대응 방안 |
|-------|-----|-------|----------|
| 플러그인 리팩토링 시간 초과 | 높음 | 높음 | Phase별 검증, 단계적 배포 |
| 기존 기능 영향 | 중간 | 높음 | 통합 테스트 우선 작성, Rollback 계획 |
| 테스트 DB 환경 차이 | 중간 | 중간 | PostgreSQL TestContainer 고려 |
| 팀원 리소스 부족 | 중간 | 중간 | 우선순위 명확화, P0만 집중 |

### 성공 지표 (KPI)

**정량 지표**:
- ✅ 테스트 커버리지: 0% → 60% (Week 8) → 85% (Week 16)
- ✅ CI/CD 통과율: 95% 이상
- ✅ 평균 함수 복잡도: <10 (McCabe)
- ✅ Pylint 점수: 8.0 이상
- ✅ P0 테스트: 100% 완료 (40개)

**정성 지표**:
- ✅ 모든 API 엔드포인트에 테스트 존재
- ✅ Repository 패턴 100% 적용
- ✅ Service 레이어 100% 분리
- ✅ 커스텀 Exception 모든 에러 케이스 커버
- ✅ 문서화 완료 (테스트 실행 가이드, 아키텍처 문서)

### 예상 일정 (간트 차트)

```
Week  1-2  | 인프라 구축               | ████████                     |
Week  3-4  | 인증 & 기본 테스트         | ████████                     |
Week  5-8  | 플러그인 리팩토링 (Critical)| ████████████████             |
Week  9-12 | 통합 테스트 & 성능         | ████████████                 |
Week 13-16 | E2E & 최종 마무리          | ████████                     |
           |___________________________|______________________________|
             인프라  인증  플러그인  통합  최종
```

**Critical Path**: Week 5-8 플러그인 리팩토링 (전체 40% 시간 소요)

---

## 부록 A. 참고 자료

### 공식 문서
- FastAPI Testing: https://fastapi.tiangolo.com/tutorial/testing/
- pytest Documentation: https://docs.pytest.org/
- SQLAlchemy Testing: https://docs.sqlalchemy.org/en/14/core/testing.html
- Celery Testing: https://docs.celeryproject.org/en/stable/userguide/testing.html

### 내부 문서 (작성 예정)
- `backend/tests/README.md` - 테스트 실행 가이드
- `backend/tests/CONTRIBUTING.md` - 테스트 작성 가이드라인
- `backend/docs/architecture.md` - 리팩토링 후 아키텍처
- `backend/docs/refactoring-guide.md` - 리팩토링 가이드
- `backend/docs/refactoring-backlog.md` - 리팩토링 백로그

### 코드 예제
- 실제 픽스처: `backend/tests/conftest.py`
- 테스트 예제: `backend/tests/api/test_auth.py`
- 서비스 레이어: `backend/app/services/plugin_service.py` (신규)
- Repository: `backend/app/repositories/plugin_repository.py` (신규)

---

## 부록 B. 용어 사전

| 용어 | 설명 |
|------|------|
| **P0** | 최우선 순위 테스트 (Critical Path) |
| **P1** | 높은 순위 테스트 (Important) |
| **P2** | 낮은 순위 테스트 (Nice to Have) |
| **AAA 패턴** | Arrange-Act-Assert 테스트 패턴 |
| **DI** | Dependency Injection (의존성 주입) |
| **Repository 패턴** | 데이터 접근 로직을 추상화하는 디자인 패턴 |
| **Service 레이어** | 비즈니스 로직을 캡슐화하는 계층 |
| **Celery EAGER 모드** | Celery 태스크를 동기적으로 실행하는 테스트 모드 |
| **순환 복잡도** | 코드의 복잡도를 측정하는 지표 (McCabe Complexity) |
| **N+1 쿼리 문제** | ORM에서 발생하는 비효율적인 쿼리 패턴 |
| **Drawflow** | 시각적 워크플로우 편집 라이브러리 |
| **Snakemake** | Python 기반 워크플로우 엔진 |

---

## 문서 변경 이력

- **2025-01-XX**: 초안 작성
- **2025-01-XX**: 실제 코드베이스 검토 반영 (v2.0)
  - 실제 프로젝트 구조 분석 추가 (섹션 2)
  - 코드 품질 이슈 상세 분석 (섹션 3)
  - 실제 모델 기반 픽스처 구체화 (섹션 4)
  - 67개 실제 테스트 케이스 정의 (섹션 9)
  - 리팩토링 예상 효과 추가
- **2025-01-XX**: 테스트 환경 구축 진행 상황 반영 (v3.0)
  - 🎯 테스트 환경 구축 진행 상황 섹션 추가
  - Section 1.1: environment-test.yml 실제 구현 반영
  - Section 4.1: PostgreSQL test-db + drop/create 격리 전략 반영
  - 주요 의사결정 및 변경 사항 문서화:
    * SQLite → PostgreSQL 전환 (JSONB 호환성)
    * Transaction rollback → Drop/create tables (격리 보장)
    * environment.yml → environment-test.yml (의존성 격리)
    * httpx 버전 호환성 해결 (0.24.1)
  - 현재 테스트 결과 및 다음 단계 정리
  - Docker Compose test-db 설정 문서화
---

## 16. 코드베이스 상세 분석 결과 (Codex 분석)

> 이 섹션은 OpenAI Codex를 활용한 자동 코드베이스 분석 결과입니다.
> 생성 일시: 2025-01-27
> 분석 도구: codex (gpt-5-codex, medium reasoning)

### 16.1 app/common/ 모듈 상세 분석

#### 복잡도 분류 기준
- **Simple**: LOC < 100, Cyclomatic Complexity < 20
- **Moderate**: LOC 100-500, CC 20-100
- **Complex**: LOC > 500, CC > 100

#### Simple 복잡도 모듈 (우선순위 높음 - 먼저 테스트)

| 모듈 | LOC | CC | 주요 기능 | 의존성 | 테스트 가능성 | 우선순위 |
|------|-----|----|-----------| -------|-------------|---------|
| `config.py` | 54 | 5 | 설정 관리 (Celery, DB, CORS) | 없음 | ✅ Easy | **P0** |
| `security.py` | 23 | 2 | JWT 토큰, 비밀번호 해싱 | config | ✅ Easy | **P0** |
| `enums.py` | 4 | 1 | Enum 타입 정의 | 없음 | ✅ Easy | P1 |
| `utils/datatable_utils.py` | 42 | 5 | 데이터테이블 Pydantic 모델 | 없음 | ✅ Easy | P2 |
| `utils/celery_utils.py` | 74 | 5 | Celery 앱 생성 | config | ⚠️ Moderate | P1 |

**테스트 전략**:
- config.py: 환경변수 로딩, 설정값 검증 테스트
- security.py: JWT 토큰 생성/검증, 비밀번호 해싱 테스트
- enums.py: Enum 값 검증 테스트

#### Moderate 복잡도 모듈 (중간 우선순위)

| 모듈 | LOC | CC | 주요 기능 | 의존성 | 테스트 가능성 | 우선순위 |
|------|-----|----|-----------| -------|-------------|---------|
| `utils/cache_utils.py` | 381 | 46 | 시각화 결과 캐싱 (7일 만료) | 없음 | ⚠️ Moderate | P1 |
| `utils/docker_utils.py` | 212 | 44 | Docker 컨테이너 관리 | 없음 | 🔴 Hard | P1 |
| `utils/error_utils.py` | 230 | 9 | 에러 핸들링 유틸리티 | 없음 | ✅ Easy | **P0** |
| `utils/github_registry_client.py` | 161 | 14 | GitHub Container Registry 클라이언트 | 없음 | ⚠️ Moderate | P1 |
| `utils/h5ad_utils.py` | 162 | 26 | H5AD 파일 처리 | 없음 | ⚠️ Moderate | P1 |
| `utils/plugin_sync_manager.py` | 167 | 17 | 플러그인 동기화 관리 | 없음 | ⚠️ Moderate | P1 |
| `utils/plugin_version_validator.py` | 173 | 19 | 플러그인 버전 검증 | github_registry_client, plugin_sync_manager | ⚠️ Moderate | P1 |
| `utils/snakemake_utils.py` | 342 | 61 | Snakemake 워크플로우 실행 | docker_utils, github_registry_client | 🔴 Hard | P1 |

**테스트 전략**:
- error_utils.py: 커스텀 예외 클래스, HTTPException 래퍼 테스트
- cache_utils.py: 캐싱 로직, 만료 정책, 심볼릭 링크 관리 테스트
- docker_utils.py: 모킹 필수 (Docker SDK 의존성)
- github_registry_client.py: HTTP 요청 모킹, 캐시 검증 테스트

#### Complex 복잡도 모듈 (후순위 - 리팩토링 필요)

| 모듈 | LOC | CC | 주요 기능 | 의존성 | 테스트 가능성 | 우선순위 | 리팩토링 필요도 |
|------|-----|----|-----------|--------|-------------|---------|----------------|
| `utils/plugin_utils.py` | 1727 | 257 | 플러그인 파일 처리 (업로드, 검증) | 없음 | 🔴 Hard | **P0** | 🔴 **Critical** |
| `utils/snakefile_dag_parser.py` | 1595 | 323 | Snakefile DAG 파싱 | 없음 | 🔴 Hard | P1 | 🔴 **Critical** |
| `utils/workflow_utils.py` | 794 | 162 | 워크플로우 검증 및 Snakefile 생성 | error_utils | 🔴 Hard | **P0** | 🟡 Medium |
| `utils/snakemake_native_parser.py` | 565 | 96 | Snakemake 네이티브 파서 | 없음 | 🔴 Hard | P1 | 🟡 Medium |

**리팩토링 권장사항**:
- `plugin_utils.py` (1727줄, CC 257): 
  - 🔴 **Critical 리팩토링 필요**
  - 단일 파일이 너무 크고 복잡함
  - 제안: PluginUploadService, PluginValidationService, PluginBuildService로 분리
  - 테스트 전략: 리팩토링 후 서비스 레이어 단위 테스트 작성

- `snakefile_dag_parser.py` (1595줄, CC 323):
  - 🔴 **Critical 리팩토링 필요**
  - DAG 파싱 로직이 매우 복잡함
  - 제안: Parser, Validator, Transformer 클래스로 분리
  - 테스트 전략: 파싱 테스트 케이스 작성 후 리팩토링

- `workflow_utils.py` (794줄, CC 162):
  - 🟡 Medium 리팩토링 필요
  - 워크플로우 검증 및 Snakefile 생성 로직 분리
  - 제안: WorkflowValidator, SnakefileGenerator 클래스 도입
  - 테스트 전략: 검증 로직과 생성 로직 분리 테스트

#### 의존성 그래프 (app/common/)

```
config.py (독립)
  ↓
security.py
  ↓
celery_utils.py

error_utils.py (독립)
  ↓
workflow_utils.py

github_registry_client.py (독립)
  ↓
plugin_sync_manager.py
  ↓
plugin_version_validator.py

docker_utils.py (독립)
  ↓
snakemake_utils.py

독립 모듈:
- enums.py
- datatable_utils.py
- cache_utils.py
- h5ad_utils.py
- plugin_utils.py
- snakefile_dag_parser.py
- snakemake_native_parser.py
```

**테스트 작성 순서 (의존성 기반)**:
1. **Wave 1** (독립 모듈, Simple): config, enums, error_utils, datatable_utils
2. **Wave 2** (Wave 1 의존): security, celery_utils, github_registry_client
3. **Wave 3** (Wave 2 의존): plugin_sync_manager, plugin_version_validator, workflow_utils
4. **Wave 4** (독립, Moderate): cache_utils, h5ad_utils, docker_utils
5. **Wave 5** (Complex, 리팩토링 후): plugin_utils, snakefile_dag_parser, snakemake_native_parser, snakemake_utils

---

### 16.2 app/routes/ 엔드포인트 상세 분석

#### 복잡도 분류

| 파일 | LOC | CC | 엔드포인트 수 | 주요 의존성 | 복잡도 | 우선순위 |
|------|-----|----|-------------|-------------|--------|---------|
| `api.py` | 10 | 1 | 0 (라우터 통합만) | endpoints.* | Simple | P2 |
| `dep.py` | 44 | 5 | 0 (의존성 주입) | config, security | Simple | **P0** |
| `celery_tasks.py` | 265 | 37 | 0 (Celery 태스크) | snakemake_utils, docker_utils, plugin_utils, cache_utils | Moderate | **P0** |
| `endpoints/auth.py` | 73 | 5 | **4개** | config, security, dep | Simple | **P0** |
| `endpoints/datatable.py` | 18 | 2 | **1개** | datatable_utils, workflow_utils, dep | Simple | P2 |
| `endpoints/files.py` | 304 | 31 | **17개** | h5ad_utils, workflow_utils, dep | Moderate | P1 |
| `endpoints/admin.py` | 533 | 62 | **23개** | plugin_sync_manager, plugin_version_validator, dep | Moderate | P1 |
| `endpoints/workflow.py` | 693 | 90 | **13개** | cache_utils, celery_utils, error_utils, plugin_utils, snakemake_utils, workflow_utils, celery_tasks, dep | Moderate | **P0** |
| `endpoints/task.py` | 1005 | 165 | **14개** | enums, celery_utils, docker_utils, snakefile_dag_parser, snakemake_native_parser, dep | **Complex** | **P0** |
| `endpoints/plugin.py` | 1213 | 210 | **24개** | github_registry_client, plugin_utils, celery_tasks, dep | **Complex** | **P0** |

#### P0 우선순위 (핵심 기능) - 총 5개 파일

**이유**: 핵심 비즈니스 로직, 높은 복잡도, 많은 의존성

| 파일 | 테스트 케이스 수 (예상) | 커버리지 목표 | 리팩토링 필요도 |
|------|----------------------|-------------|---------------|
| `dep.py` | 5개 | 90% | 🟢 Low |
| `celery_tasks.py` | 8개 | 75% | 🟡 Medium |
| `endpoints/auth.py` | 13개 (기존 11개 + 2개) | 95% | 🟢 Low |
| `endpoints/workflow.py` | 20개 | 70% | 🟡 Medium |
| `endpoints/task.py` | 25개 | 65% | 🔴 High |
| `endpoints/plugin.py` | 30개 | 60% | 🔴 **Critical** |

**세부 테스트 계획**:

**1. dep.py (의존성 주입)** - 5개 테스트
- get_db() 세션 생성/종료 테스트
- get_current_user() JWT 검증 테스트
- 유효/만료/없음 토큰 시나리오

**2. celery_tasks.py (Celery 태스크)** - 8개 테스트
- process_data_task 정상 실행 테스트
- process_data_task 실패 시 컨테이너 정리 테스트
- build_plugin_task 성공/실패 테스트
- MyTask lifecycle hooks 테스트 (before_start, on_success, on_failure, on_revoke)

**3. endpoints/auth.py** - 13개 테스트 (backend-test-docs.md Section 9.1 참조)
- 기존 11개 + 추가 2개 (토큰 갱신, 플러그인 연결/해제)

**4. endpoints/workflow.py** - 20개 테스트
- POST /compile: 워크플로우 컴파일 (정상, 잘못된 drawflow, 파일 누락) - 4개
- POST /visualization: 시각화 생성 (정상, 플러그인 없음, 파라미터 누락) - 3개
- POST /save: 워크플로우 저장 (정상, 중복 제목, 권한 없음) - 3개
- GET /me: 사용자 워크플로우 목록 - 2개
- POST /find: 워크플로우 검색 - 2개
- POST /delete: 워크플로우 삭제 - 2개
- POST /node/*: 노드 저장/읽기/삭제 - 4개

**5. endpoints/task.py** - 25개 테스트
- GET /info/{task_id}: 태스크 정보 조회 - 3개
- GET /monitoring: 태스크 모니터링 - 2개
- DELETE /revoke/{task_id}: 태스크 취소 - 3개
- DELETE /delete/{task_id}: 태스크 삭제 - 2개
- GET /logs/{task_id}: 로그 조회 - 3개
- GET /{task_id}/dag-structure: DAG 구조 - 2개
- GET /{task_id}/rule-status: Rule 상태 - 2개
- GET /{task_id}/enhanced-progress: 진행률 - 2개
- GET /{task_id}/rule-logs/{rule_name}: Rule 로그 - 2개
- GET /cache/stats, DELETE /cache/clear: 캐시 관리 - 2개
- GET /{task_id}/execution-manifest: 실행 manifest - 2개

**6. endpoints/plugin.py** - 30개 테스트 (backend-test-docs.md Section 9.2 참조)
- 플러그인 생명주기 (검증, 업로드, 빌드) - 10개
- 스크립트 업로드 및 롤백 - 5개
- 버전 관리 (버전 조회, 업데이트) - 5개
- 플러그인 목록 및 정보 - 5개
- 빌드 태스크 관리 - 5개

#### P1 우선순위 (중요 기능) - 총 2개 파일

| 파일 | 테스트 케이스 수 (예상) | 커버리지 목표 | 리팩토링 필요도 |
|------|----------------------|-------------|---------------|
| `endpoints/files.py` | 20개 | 70% | 🟢 Low |
| `endpoints/admin.py` | 25개 | 65% | 🟡 Medium |

**세부 테스트 계획**:

**7. endpoints/files.py** - 20개 테스트
- POST /upload: 파일 업로드 (정상, 크기 초과, 중복) - 4개
- GET /me: 파일 목록 - 2개
- POST /find, POST /folder, POST /delete: 파일 관리 - 6개
- POST /convert: H5AD 변환 - 2개
- GET /check/{file_name}: 파일 존재 확인 - 2개
- POST /columns, POST /clusters: 메타데이터 조회 - 4개

**8. endpoints/admin.py** - 25개 테스트
- GET /users, /users_count: 사용자 관리 - 4개
- GET /files, /files_count: 파일 관리 - 4개
- GET /workflows, /workflows_count: 워크플로우 관리 - 4개
- GET /tasks, /tasks_count: 태스크 관리 - 4개
- GET /plugins, /plugins_count: 플러그인 관리 - 4개
- PUT /users/{user_id}, DELETE /users/{user_id}: 사용자 수정/삭제 - 3개
- POST /plugins/sync, GET /plugins/consistency: 플러그인 동기화 - 2개

#### P2 우선순위 (부가 기능) - 총 2개 파일

| 파일 | 테스트 케이스 수 (예상) | 커버리지 목표 | 리팩토링 필요도 |
|------|----------------------|-------------|---------------|
| `api.py` | 2개 | 80% | 🟢 Low |
| `endpoints/datatable.py` | 3개 | 80% | 🟢 Low |

**세부 테스트 계획**:

**9. api.py** - 2개 테스트
- 라우터 통합 확인 테스트
- 엔드포인트 prefix/tags 검증 테스트

**10. endpoints/datatable.py** - 3개 테스트
- POST /load_data: 데이터테이블 로드 (정상, 파일 없음, 잘못된 형식)

---

### 16.3 테스트 작성 우선순위 매트릭스

#### 위험도 평가 (Risk Assessment)

| 모듈/엔드포인트 | 복잡도 | 의존성 | 비즈니스 중요도 | 위험 점수 (0-10) | 우선순위 |
|----------------|--------|--------|---------------|-----------------|---------|
| `endpoints/plugin.py` | 10 | 8 | 10 | **9.3** | **P0** |
| `endpoints/task.py` | 10 | 7 | 10 | **9.0** | **P0** |
| `utils/plugin_utils.py` | 10 | 3 | 10 | **7.7** | **P0** |
| `endpoints/workflow.py` | 7 | 9 | 10 | **8.7** | **P0** |
| `celery_tasks.py` | 5 | 8 | 10 | **7.7** | **P0** |
| `utils/workflow_utils.py` | 9 | 4 | 9 | **7.3** | **P0** |
| `endpoints/auth.py` | 2 | 5 | 10 | **5.7** | **P0** |
| `dep.py` | 2 | 3 | 10 | **5.0** | **P0** |
| `config.py` | 2 | 0 | 9 | **3.7** | **P0** |
| `security.py` | 1 | 2 | 10 | **4.3** | **P0** |
| `error_utils.py` | 3 | 0 | 8 | **3.7** | **P0** |
| `utils/snakefile_dag_parser.py` | 10 | 0 | 7 | **5.7** | P1 |
| `utils/docker_utils.py` | 6 | 0 | 8 | **4.7** | P1 |
| `utils/snakemake_utils.py` | 7 | 5 | 8 | **6.7** | P1 |
| `endpoints/files.py` | 5 | 4 | 7 | **5.3** | P1 |
| `endpoints/admin.py` | 7 | 5 | 6 | **6.0** | P1 |

**위험 점수 계산식**: (복잡도 × 0.4) + (의존성 × 0.3) + (비즈니스 중요도 × 0.3)

**위험도 기준**:
- 🔴 **Critical** (8.0-10.0): 즉시 테스트 필요, 리팩토링 고려
- 🟡 **High** (6.0-7.9): 우선 테스트 필요
- 🟢 **Medium** (4.0-5.9): 중간 우선순위
- ⚪ **Low** (0-3.9): 낮은 우선순위

#### 테스트 커버리지 목표 (단계별)

**Phase 1** (Week 1-2): 기본 인프라 + Simple 모듈
- 목표 커버리지: **30%**
- 모듈: config, security, enums, error_utils, dep, auth.py, datatable.py
- 테스트 수: ~35개

**Phase 2** (Week 3-4): Moderate 모듈 + 핵심 엔드포인트
- 목표 커버리지: **50%**
- 모듈: celery_tasks, workflow.py, files.py
- 테스트 수: ~45개

**Phase 3** (Week 5-6): Complex 엔드포인트
- 목표 커버리지: **65%**
- 모듈: task.py, plugin.py, admin.py
- 테스트 수: ~80개

**Phase 4** (Week 7-8): Utils 모듈
- 목표 커버리지: **75%**
- 모듈: cache_utils, h5ad_utils, docker_utils, github_registry_client, plugin_sync_manager
- 테스트 수: ~40개

**Phase 5** (Week 9-12): Complex Utils (리팩토링 병행)
- 목표 커버리지: **85%**
- 모듈: plugin_utils, snakefile_dag_parser, workflow_utils, snakemake_utils
- 테스트 수: ~50개

---

### 16.4 테스트 작성 로드맵 (상세 계획)

#### Week 1-2: 기본 인프라 + Simple 모듈 (30% 커버리지)

**목표**: 테스트 환경 안정화 및 기본 모듈 테스트 작성

**작업 항목**:
- [x] 테스트 환경 구축 (PostgreSQL test-db, Conda environment)
- [x] conftest.py 픽스처 구현
- [x] pytest.ini 설정
- [ ] **app/common/config.py** 테스트 작성 (5개)
  - 환경변수 로딩 테스트
  - Settings 클래스 검증
  - route_task 함수 테스트
- [ ] **app/common/security.py** 테스트 작성 (5개)
  - JWT 토큰 생성/검증 테스트
  - 비밀번호 해싱 테스트
  - 만료 토큰 테스트
- [ ] **app/common/enums.py** 테스트 작성 (2개)
- [ ] **app/common/utils/error_utils.py** 테스트 작성 (10개)
  - 커스텀 예외 클래스 테스트
  - HTTPException 래퍼 테스트
- [ ] **app/routes/dep.py** 테스트 작성 (5개)
  - get_db() 세션 관리 테스트
  - get_current_user() 인증 테스트
- [x] **app/routes/endpoints/auth.py** 테스트 작성 (13개) - 기존 11개 완료, 2개 추가
  - 플러그인 연결/해제 테스트 추가
- [ ] **app/routes/endpoints/datatable.py** 테스트 작성 (3개)

**예상 소요**: 10-12일
**산출물**: ~43개 테스트, 30% 커버리지

#### Week 3-4: Moderate 모듈 + 핵심 엔드포인트 (50% 커버리지)

**목표**: 핵심 비즈니스 로직 테스트 및 Celery 태스크 테스트

**작업 항목**:
- [ ] **app/routes/celery_tasks.py** 테스트 작성 (8개)
  - process_data_task 테스트
  - build_plugin_task 테스트
  - MyTask lifecycle hooks 테스트
- [ ] **app/routes/endpoints/workflow.py** 테스트 작성 (20개)
  - 워크플로우 컴파일 테스트
  - 시각화 생성 테스트
  - 워크플로우 CRUD 테스트
  - 노드 관리 테스트
- [ ] **app/routes/endpoints/files.py** 테스트 작성 (20개)
  - 파일 업로드 테스트
  - H5AD 변환 테스트
  - 메타데이터 조회 테스트
- [ ] **app/common/utils/cache_utils.py** 테스트 작성 (8개)
  - 캐싱 로직 테스트
  - 만료 정책 테스트
  - 심볼릭 링크 관리 테스트

**예상 소요**: 10-12일
**산출물**: ~56개 테스트, 50% 커버리지 (누적)

#### Week 5-6: Complex 엔드포인트 (65% 커버리지)

**목표**: 복잡한 엔드포인트 테스트 작성

**작업 항목**:
- [ ] **app/routes/endpoints/task.py** 테스트 작성 (25개)
  - 태스크 정보 조회 테스트
  - 태스크 취소/삭제 테스트
  - 로그 조회 테스트
  - DAG 구조 테스트
  - 진행률 테스트
- [ ] **app/routes/endpoints/plugin.py** 테스트 작성 (30개)
  - 플러그인 검증 테스트
  - 플러그인 업로드 테스트
  - 스크립트 업로드 및 롤백 테스트
  - 빌드 태스크 테스트
  - 버전 관리 테스트
- [ ] **app/routes/endpoints/admin.py** 테스트 작성 (25개)
  - 관리자 CRUD 테스트
  - 플러그인 동기화 테스트
  - 시스템 통계 테스트

**예상 소요**: 12-14일
**산출물**: ~80개 테스트, 65% 커버리지 (누적)

#### Week 7-8: Utils 모듈 (75% 커버리지)

**목표**: 유틸리티 모듈 테스트 작성

**작업 항목**:
- [ ] **app/common/utils/h5ad_utils.py** 테스트 작성 (10개)
- [ ] **app/common/utils/docker_utils.py** 테스트 작성 (12개) - 모킹 필수
- [ ] **app/common/utils/github_registry_client.py** 테스트 작성 (8개)
- [ ] **app/common/utils/plugin_sync_manager.py** 테스트 작성 (10개)
- [ ] **app/common/utils/celery_utils.py** 테스트 작성 (5개)
- [ ] **app/common/utils/datatable_utils.py** 테스트 작성 (3개)

**예상 소요**: 10-12일
**산출물**: ~48개 테스트, 75% 커버리지 (누적)

#### Week 9-12: Complex Utils (리팩토링 병행) (85% 커버리지)

**목표**: 복잡한 유틸리티 모듈 리팩토링 및 테스트 작성

**작업 항목**:
- [ ] **app/common/utils/plugin_utils.py** 리팩토링 + 테스트 (20개)
  - 리팩토링: PluginUploadService, PluginValidationService, PluginBuildService 분리
  - 테스트: 서비스별 단위 테스트 작성
- [ ] **app/common/utils/workflow_utils.py** 리팩토링 + 테스트 (15개)
  - 리팩토링: WorkflowValidator, SnakefileGenerator 분리
  - 테스트: 검증 로직 및 생성 로직 분리 테스트
- [ ] **app/common/utils/snakefile_dag_parser.py** 리팩토링 + 테스트 (15개)
  - 리팩토링: Parser, Validator, Transformer 클래스 분리
  - 테스트: 파싱 테스트 케이스 작성
- [ ] **app/common/utils/snakemake_utils.py** 테스트 작성 (12개)
- [ ] **app/common/utils/snakemake_native_parser.py** 테스트 작성 (10개)
- [ ] **app/common/utils/plugin_version_validator.py** 테스트 작성 (8개)

**예상 소요**: 20-25일
**산출물**: ~80개 테스트, 85% 커버리지 (누적)

---

### 16.5 테스트 작성 지침 (추가)

#### 모킹 전략

**Docker 모킹**:
```python
@pytest.fixture
def mock_docker(mocker):
    """Docker 클라이언트 모킹"""
    mock_client = mocker.patch("docker.from_env")
    mock_container = mocker.MagicMock()
    mock_client.return_value.containers.run.return_value = mock_container
    return mock_client
```

**Celery 모킹**:
```python
@pytest.fixture(autouse=True)
def mock_celery(monkeypatch):
    """Celery EAGER 모드 설정"""
    monkeypatch.setenv("CELERY_ALWAYS_EAGER", "True")
    monkeypatch.setenv("CELERY_EAGER_PROPAGATES", "True")
```

**HTTP 요청 모킹** (GitHub Registry):
```python
@pytest.fixture
def mock_requests(mocker):
    """HTTP 요청 모킹"""
    mock_get = mocker.patch("requests.get")
    mock_get.return_value.status_code = 200
    mock_get.return_value.json.return_value = {"data": "test"}
    return mock_get
```

#### 복잡한 모듈 테스트 전략

**plugin_utils.py 테스트**:
- 리팩토링 전: 통합 테스트 위주 (E2E 시나리오)
- 리팩토링 후: 서비스별 단위 테스트 + 통합 테스트

**snakefile_dag_parser.py 테스트**:
- Snakefile 샘플 파일 준비 (fixtures/sample_snakefile.smk)
- 파싱 결과 JSON 비교 테스트
- 에러 케이스 테스트 (잘못된 구문, 순환 의존성)

**workflow_utils.py 테스트**:
- Drawflow JSON 샘플 준비 (fixtures/sample_drawflow.json)
- Snakefile 생성 결과 검증
- 워크플로우 검증 로직 테스트

---

### 16.6 리팩토링 백로그 (우선순위별)

#### 🔴 Critical (즉시 처리)

1. **app/common/utils/plugin_utils.py** (1727줄, CC 257)
   - 현재 문제: 단일 파일에 플러그인 업로드, 검증, 빌드 로직이 모두 포함
   - 리팩토링 계획:
     ```python
     # services/plugin_upload_service.py
     class PluginUploadService:
         def upload_plugin(self, files, metadata): ...
         def upload_scripts(self, plugin_name, files): ...
     
     # services/plugin_validation_service.py
     class PluginValidationService:
         def validate_metadata(self, metadata): ...
         def validate_scripts(self, scripts): ...
     
     # services/plugin_build_service.py
     class PluginBuildService:
         def build_docker_image(self, plugin_name): ...
         def check_image_exists(self, plugin_name): ...
     ```
   - 예상 효과: 테스트 가능성 300% 증가, CC 257 → 50-70
   - 소요 시간: 5-7일

2. **app/common/utils/snakefile_dag_parser.py** (1595줄, CC 323)
   - 현재 문제: DAG 파싱 로직이 매우 복잡하고 가독성 낮음
   - 리팩토링 계획:
     ```python
     # parsers/snakefile_parser.py
     class SnakefileParser:
         def parse(self, snakefile_path): ...
     
     # parsers/snakefile_validator.py
     class SnakefileValidator:
         def validate_dag(self, dag): ...
         def check_circular_dependencies(self, dag): ...
     
     # parsers/snakefile_transformer.py
     class SnakefileTransformer:
         def transform_to_json(self, dag): ...
     ```
   - 예상 효과: CC 323 → 80-100, 테스트 커버리지 0% → 60%
   - 소요 시간: 7-10일

#### 🟡 Medium (Phase 4-5에 처리)

3. **app/common/utils/workflow_utils.py** (794줄, CC 162)
   - 리팩토링 계획: WorkflowValidator, SnakefileGenerator 분리
   - 소요 시간: 3-5일

4. **app/routes/endpoints/plugin.py** (1213줄, CC 210)
   - 리팩토링 계획: 비즈니스 로직을 PluginService로 이동
   - 소요 시간: 5-7일

5. **app/routes/endpoints/task.py** (1005줄, CC 165)
   - 리팩토링 계획: TaskService, LogService 분리
   - 소요 시간: 4-6일

#### 🟢 Low (Phase 5 이후 또는 선택)

6. **app/common/utils/snakemake_native_parser.py** (565줄, CC 96)
7. **app/routes/endpoints/admin.py** (533줄, CC 62)
8. **app/common/utils/snakemake_utils.py** (342줄, CC 61)

---

### 16.7 예상 일정 및 리소스

**총 기간**: 9-12주 (리팩토링 포함 14-16주)

**리소스 투입**:
- 개발자 1명 full-time 기준
- 주당 30-35시간 테스트 작성 가능
- 리팩토링 병행 시 주당 20-25시간 테스트 작성

**마일스톤**:
- Week 2: 30% 커버리지 달성
- Week 4: 50% 커버리지 달성
- Week 6: 65% 커버리지 달성
- Week 8: 75% 커버리지 달성
- Week 12: 85% 커버리지 달성 (리팩토링 포함 Week 16)

**위험 요소**:
- plugin_utils.py, snakefile_dag_parser.py 리팩토링 시간 초과 가능성
- Celery, Docker 모킹 복잡도로 인한 테스트 작성 지연
- 기존 코드 변경 시 회귀 테스트 필요

**대응 방안**:
- 리팩토링과 테스트 작성을 병행하여 시간 절약
- 모킹 라이브러리 (pytest-mock) 적극 활용
- 기존 테스트 자동화 (CI/CD)로 회귀 테스트 부담 감소

---

