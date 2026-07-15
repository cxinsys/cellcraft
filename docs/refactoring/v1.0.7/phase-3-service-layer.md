# Phase 3 — 서비스 계층 추출 (PR-5 ~ PR-8)

> 예상: 3~5일 | 위험도: 중간 | 선행: PR-4 머지
> 이 Phase가 리팩토링의 핵심 가치. 도메인별 1 PR, **PR-5 → 8 순서 권장이나 도메인이 겹치지 않으므로 병렬 진행 가능.**

## 공통 규칙 (모든 PR 적용)

### 서비스 함수 형태 (Dispatch식)
```python
# app/plugin/service.py
def upload_plugin(*, db: Session, user: User, payload: PluginUploadRequest) -> Plugin:
    """비즈니스 로직만. HTTP 관심사 없음."""
```
- 키워드 전용 인자(`*`), 첫 인자로 `db: Session`
- 반환은 모델 또는 도메인 객체 — `JSONResponse`/`HTTPException` 금지
- 실패는 `core/exceptions.py`의 도메인 예외로 raise

### 라우터의 역할 (이것만 남긴다)
1. 요청 파싱/검증 (pydantic 스키마)
2. 인증/인가 dependency
3. service 호출
4. 응답 변환 (+ 도메인 예외 → HTTPException 매핑은 전역 핸들러가 담당)

### 추출 방법 (동작 보존)
1. 엔드포인트 본문을 그대로 service 함수로 복사-이동 (copy-paste refactoring)
2. HTTP 관심사(HTTPException, status_code, Response 객체)를 경계 밖으로 밀어내기
3. characterization 테스트(PR-1)로 전후 동일성 확인
4. service 단위 테스트 신규 작성 — **이 단계에서 커버리지를 만든다** (도메인 service 80% 목표)
5. 한 엔드포인트씩 커밋 — 커밋당 리뷰 가능 단위 유지

### 응답 형태 절대 불변
JSON 키, 상태코드, 에러 메시지 문자열까지 유지 (OpenAPI diff 0 + characterization 테스트로 검증).

---

## PR-5: plugin 서비스 (`refactor/v107-05-plugin-service`)

대상: `plugin/router.py` 1,616줄 → 목표 ~400줄

### 추출 대상 (엔드포인트 → 서비스 함수)
| 현재 엔드포인트 로직 | 서비스 함수 | 비고 |
|---|---|---|
| `upload_plugin()` 250줄+ | `service.upload_plugin()` — 내부를 단계 함수로 분해: `_save_plugin_files()`, `_build_image()`(→ `builder.py` 위임), `_upsert_plugin()`, `_rollback_upload()` | 파일 I/O·Docker 빌드·DB upsert·Celery 디스패치·롤백이 한 함수에 있음 |
| 플러그인 목록/상세/삭제 | `service.list_plugins()`, `get_plugin()`, `delete_plugin()` | CRUD 위임 위주 |
| 빌드 상태/로그 조회 | `service.get_build_status()` | |
| 공식 플러그인 sync 관련 | 기존 `plugin/sync.py` 재사용, router에서 직접 로직 금지 |

### 체크리스트
- [x] Docker 빌드 로직을 `plugin/builder.py`로 위임 (utils.py에서 빌드 관련 호출부만 — utils 내부 분할은 PR-9)
- [x] 롤백 로직을 명시적 함수로 분리하고 단위 테스트 (실패 주입 테스트: 빌드 실패 시 파일/DB 원복)
- [x] `plugin/router.py` 400줄 이하 확인 (301줄)
- [x] service 단위 테스트 (Docker 모킹)

---

## PR-6: workflow 서비스 (`refactor/v107-06-workflow-service`)

대상: `workflow/router.py` 988줄

### 추출 대상
| 현재 | 서비스 함수 |
|---|---|
| `compileWorkflow()` 200줄+ (DAG 파싱→경로 해석→Snakefile 생성→캐싱 인라인) | `service.compile_workflow()` — 파이프라인 단계 함수: `_parse_dag()`, `_resolve_paths()`, `_generate_snakefile()`, `_cache_result()` (compiler/ 모듈 위임) |
| 워크플로우 CRUD 엔드포인트 | `service.create/get/list/update/delete_workflow()` |

### 체크리스트
- [x] compile 파이프라인의 각 단계가 `workflow/compiler/*`를 호출하는 구조로 정리 (파서 내부 수정은 PR-10)
- [x] 캐시 적중/미적중 경로 characterization 테스트 통과
- [x] service 단위 테스트 (파서는 실제, 파일시스템은 tmp_path)

---

## PR-7: task 서비스 + SSE 분리 (`refactor/v107-07-task-service`)

대상: `task/router.py` 1,447줄

### 추출 대상
| 현재 | 이동 후 |
|---|---|
| `get_task_status()` SSE 루프 80줄+ (수동 타임아웃·리소스 정리) | `task/sse.py` — 제너레이터/스트리밍 유틸로 분리, router는 `EventSourceResponse(sse.stream_task_status(...))` 형태 |
| 태스크 생성·Celery 디스패치 (15+ kwargs 인라인) | `service.submit_task()` + **kwargs를 pydantic 스키마로 타입화**: `task/schemas.py`에 `TaskSubmission(BaseModel)` 정의, `.dict()`로 디스패치 |
| 태스크 취소/재시도/결과 조회 | `service.cancel_task()`, `get_task_result()` 등 |

### Celery 디스패치 타입화 (이 PR의 핵심 개선)
```python
# task/schemas.py (pydantic v1 문법)
class TaskSubmission(BaseModel):
    user_id: int
    workflow_id: int
    algorithm_id: Optional[int]
    plugin_name: str
    cache_key: Optional[str]
    # ... 현재 celery_tasks.py 태스크 시그니처의 kwargs 전부 명시

# task/service.py
def submit_task(*, db: Session, user: User, payload: ...) -> Task:
    submission = TaskSubmission(...)
    celery_task = run_pipeline.apply_async(kwargs=submission.dict(), ...)
```
- [~] `worker/tasks.py` 쪽에서도 수신 시 같은 스키마로 파싱 → 런타임 검증
      (SKIPPED: `process_data_task` 시그니처는 `**kwargs`가 아닌 명시 인자라 워커 변경 없음. `TaskSubmission`은 task/schemas.py에 정의 — 디스패치측 타입화 전용)
- [x] **와이어 포맷(kwargs 키 이름) 불변** — 구버전 워커/신버전 웹 혼재 배포 호환
      (`TaskSubmission.dict(exclude_none=True)` == 기존 compile/visualization kwargs, 왕복 테스트로 검증)

### 체크리스트
- [x] SSE 이벤트 포맷·종료 조건 characterization 테스트 통과 (patch 경로만 `app.task.sse`/`app.task.service`로 갱신)
- [x] 타임아웃/정리(finally) 로직이 sse.py로 이동 후에도 동일 동작 (`stream_task_status` 단위 테스트: 종료조건/타임아웃/예외/cleanup)
- [x] service + sse 단위 테스트 (tests/unit/test_task_service.py, 28건)

---

## PR-8: 나머지 도메인 + 전역 예외 처리 (`refactor/v107-08-remaining-services`)

### 1. 전역 예외 핸들러 (이 PR에서 도입 — 이후 전 도메인에 적용됨)
- [x] `core/exceptions.py`에 도메인 예외 계층 정의 (기존 error_utils 유지, 계층만 추가):
```python
class CellcraftError(Exception):        # status_code=500, .detail 보유, status_code= 오버라이드 지원
class NotFoundError(CellcraftError):        # 404
class PermissionDeniedError(CellcraftError):  # 403
class ValidationFailedError(CellcraftError):  # 400
class ExternalServiceError(CellcraftError):   # 502 (Docker/GitHub/Redis), 500/503 오버라이드
```
- [x] `main.py`(create_app)에 `app.add_exception_handler(CellcraftError, ...)` 등록 — 핸들러가 `JSONResponse(status_code=exc.status_code, content={"detail": exc.detail})`를 반환하여 FastAPI 기본 HTTPException 응답과 **바이트 단위 동일**. 422 RequestValidationError 등 기본 동작은 미변경 (CellcraftError 서브클래스에만 발동)
- [x] PR-5~7 서비스 HTTPException → 도메인 예외 전환 (응답 불변 검증: characterization+contract 통과):
  - **전환**: workflow(list/find/delete_workflow, get_results, save/read/delete_node_data, get_visualization_result, check_visualization_result), task(get_resources 503, get_task_status_simple 400). 모두 감싸는 swallow 가드가 없는 sentinel raise만 전환하고 관련 unit 테스트 예외 타입 갱신
  - **유지(보고)**: `compile_workflow`/`_resolve_plugin_paths`(외곽 `except Exception→HTTPException(400,str(e))`가 상세를 재-문자열화 → 전환 시 본문 변형), `visualize_data`(CellCraftHTTPException 구조화 본문/`create_error_response`), plugin/service.py 전체(모든 raise가 `except HTTPException→raise / except Exception→500` 가드 체인 안 + 기존 unit이 HTTPException(500) 롤백을 단언), task/service.py의 로그/DAG/manifest 함수(`except HTTPException: raise` 가드 체인). 모두 "무리한 전환 금지"에 따라 HTTPException 유지 — 응답 불변 보장

### 2. file 서비스
- [x] `file/router.py` 598줄 → 178줄 로직을 `file/service.py`로 (업로드/다운로드/삭제/컬럼/클러스터/setup/shared/tutorial). router 195줄 (≤300 목표 달성)
- [x] `file/security.py`(validate_file_upload/validate_file_path) 호출을 service 경유로 일원화. path-traversal 거부는 security.py의 HTTPException 그대로 (FastAPI 기본 핸들러 처리, 응답 불변)

### 3. admin 서비스
- [x] `admin/router.py` 760줄 → 314줄. 통계/사용자·파일·워크플로우·태스크·플러그인 관리, 시스템 통계(Docker), 플러그인 sync, 컨테이너 매니저 로직을 `admin/service.py`로. 반복되던 `is_superuser` 가드는 `_require_admin()`(PermissionDeniedError)로 일원화 (router 400 목표 달성)
- [x] 타 도메인 데이터 접근 **보고**: `admin/crud.py`가 users/files/workflows/tasks/plugins 테이블을 직접 조회 (대시보드용). 이는 PR-8 이전부터 존재하며 도메인 경계 강화는 범위 밖 — 그대로 이동만 함 (계획의 "일단 그대로 옮기고 보고" 지침대로)

### 4. auth·user·datatable 서비스
- [x] `auth/service.py`: 로그인/토큰 발급. router는 OAuth2 폼 처리 + response_model만
- [x] `user/service.py`: 가입/조회/수정. **보고**: register/me/plugins 3개 엔드포인트는 물리적으로 `auth/router.py`에 잔류 (다른 router로 이동 시 마운트 경로가 바뀌어 OpenAPI 스냅샷 계약 위반). 로직만 user/service.py로 추출
- [x] `datatable/service.py`: h5ad 조회·캐싱 오케스트레이션. router 29줄 → 13줄 (utils 내부 미변경, service 경유로만 전환)

### 체크리스트
- [x] 신규 5개 도메인(file/admin/auth/user/datatable) router의 `try/except + HTTPException` 반복 패턴 제거 (도메인 예외 + 전역 핸들러로 대체). PR-5~7 서비스는 위 §1 "유지" 항목대로 일부 잔존
- [x] 에러 응답 본문·상태코드 불변 (characterization file/auth·task 통과 + contract OpenAPI diff 0)
- [~] import-linter 계약 `router는 crud를 직접 import하지 않는다`: 신규 5개 router는 이미 crud 미import(service 경유). 단, plugin/workflow router는 아직 `crud`를 `# noqa: F401` 재-export(기존 테스트 patch 경로 호환용)로 보유 — 계약 활성화는 전 도메인 동시 충족이 필요하므로 후속 PR로 이관 (현 계약 2종 KEPT 유지, 신규 계약 미추가)

---

## Phase 3 완료 기준
- [x] 8개 router 파일 각각 500줄 이하 (file 195 · admin 314 · auth 53 · datatable 13 · plugin 301 · task 245 · workflow 177 · version 15)
- [x] 모든 도메인에 `service.py` 존재 (plugin/workflow/task/file/admin/auth/user/datatable), 신규 service 단위 테스트 신규 작성 (test_{file,admin,auth,datatable,user}_service.py, 44건). 도메인 예외 매핑·검증 브랜치 커버
- [x] router에 남은 것: 파싱/deps/service 호출/응답 변환뿐 (admin의 반복 superuser 가드는 service `_require_admin`로 이동)
- [x] 공통 머지 게이트: 전체 테스트 실패 집합 baseline과 **완전 동일**(68 failed, 신규 regression 0·우발 수정 0), 신규 테스트 통과로 570 passed(=524+46), contract OpenAPI diff 0, import-linter 2 contracts KEPT
