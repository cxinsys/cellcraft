# FastAPI 백엔드 리팩토링 마스터 플랜 (release/v1.0.7)

작성일: 2026-07-11
대상: `backend/app` (약 20,000 LOC, 60개 Python 파일, 8개 라우터 / 99개 엔드포인트)
참고 아키텍처: [fastapi-best-architecture](https://github.com/fastapi-practices/fastapi-best-architecture), [Netflix Dispatch](https://github.com/Netflix/dispatch), [Onyx](https://github.com/onyx-dot-app/onyx)

---

## 1. 배경과 문제 정의

### 현재 구조의 문제
| # | 문제 | 근거 |
|---|---|---|
| 1 | **서비스 계층 부재** — 비즈니스 로직이 엔드포인트에 인라인 | `plugin.py` 1,616줄(upload_plugin 한 함수 250줄+: 파일 I/O·Docker 빌드·DB·Celery·롤백), `task.py` 1,447줄, `workflow.py` 988줄 |
| 2 | **utils 잡동사니** — 도메인 경계 없는 26개 파일 11.5K LOC | `plugin_utils.py` 2,220줄, `snakefile_dag_parser.py` 2,227줄, utils 간 교차 의존 11건+, 함수 내부 지연 import 다수 |
| 3 | **main.py 과부하 (337줄)** — import 시점 `run_migrations()` 실행, Celery 앱 생성, 시그널 핸들러, 플러그인 sync, Docker pull 혼재 | `app/main.py:84` |
| 4 | **Celery 워커가 웹 앱 전체를 import** — 워커 부팅마다 마이그레이션 실행 | `Dockerfile.celery`: `celery -A app.main.celery worker` |
| 5 | **애매한 위치의 코드** — 루트 `app/model.py`(auth 스키마), `routes/dep.py`(DB 세션 + 인증 혼재) | |
| 6 | **전역 가변 싱글턴** — `docker_utils.container_manager` (스레드 안전성 없음) | |

### 기술 스택 제약 (1.0.7에서 변경하지 않음)
| 항목 | 버전 | 제약 |
|---|---|---|
| FastAPI | 0.78.0 | lifespan 파라미터 미지원 → `@app.on_event` 유지 |
| Pydantic | 1.9.1 | v2 문법(`model_validate` 등) 사용 금지 |
| SQLAlchemy | 1.4.36 | 동기 ORM 유지 |
| Celery | 5.2.7 | 태스크 이름 하위호환 필수 |
| Snakemake | 7.14.0 | conda 환경 고정 |

### 유지할 자산 (건드리지 않음)
- Alembic 단일 스키마 관리 + 3-상태 감지 부팅 로직
- pytest 구조 (`tests/unit`, `tests/integration`, `tests/benchmarks`, conftest 격리)
- CRUD 계층 분리 (`crud/base.py` 패턴)
- 최근 보안 커밋 (CORS 제한, 환경 인지 앱 설정)

---

## 2. 아키텍처 결정 (ADR)

### 결정: Dispatch식 도메인 패키지 + FBA식 계층 규율 + Onyx식 워커 분리

| 레퍼런스 | 채택 | 비채택 | 이유 |
|---|---|---|---|
| **Netflix Dispatch** ★ | 도메인 플랫 패키지(도메인 폴더에 router/service/crud/models/schemas 동거), 서비스=평범한 함수 | ― | 동기 SQLAlchemy + 중간 규모 + 명확한 도메인 — cellcraft와 조건 최유사 |
| **fastapi-best-architecture** | `router → service → crud` 단방향 계층 규율 | SQLAlchemy 2.0 async, Pydantic v2, Casbin RBAC, Snowflake ID | 스택 불일치·과설계 |
| **Onyx** | 웹 계층 / Celery `worker/` 대분리 | 초대형 규모용 server/db 완전 분리 | 워커의 `app.main` import 문제 직접 해결 |

### 설계 규칙 (모든 PR에 적용)
1. **의존 방향**: `router → service → crud/모듈` 단방향. 신규 코드에서 router의 crud 직접 호출 금지.
2. **service는 평범한 함수**: `def upload_plugin(*, db: Session, user: User, payload: PluginUpload) -> Plugin`. 클래스/DI 컨테이너 불필요.
3. **HTTPException은 router에서만**. service는 `core/exceptions.py`의 도메인 예외를 raise.
4. **도메인 간 참조는 상대 도메인의 service 경유** (crud 직접 접근 금지).
5. **모델 분리 후에도 단일 `Base` metadata 유지**. relationship은 문자열 참조(`"User"`)로 순환 import 방지.
6. **worker는 웹 계층(FastAPI) 미의존**.
7. 파일당 800줄 이하, 함수당 50줄 이하 (신규/수정 코드 기준).

### 목표 디렉터리 구조
```
backend/app/
├── main.py                  # create_app() 팩토리 + 라우터 마운트 (~60줄)
├── api.py                   # api_router 집계
├── core/                    # config, security, exceptions, logging, startup, enums
├── db/                      # session(engine/SessionLocal/get_db), base(Base/CRUDBase)
├── auth/                    # router, service, schemas, deps(get_current_user)
├── user/                    # models, crud, schemas, service
├── plugin/                  # router, service, crud, models, schemas, builder, registry, sync, versioning, cache
├── workflow/                # router, service, crud, models, schemas, compiler/{dag_parser,native_parser,snakefile}
├── task/                    # router, service, crud, models, schemas, sse
├── file/                    # router, service, crud, models, schemas, storage, security
├── datatable/               # router, service, schemas, h5ad, cache
├── admin/                   # router, service, crud, schemas
├── version/                 # router
├── shared/                  # docker, redis, cache, resources, archive (도메인 무관 인프라)
└── worker/                  # app(create_celery), tasks (FastAPI 미의존)
```

---

## 3. PR 로드맵 (총 10개 PR, 순차 머지)

모든 PR의 base 브랜치는 **`release/v1.0.7`**. 브랜치 네이밍: `refactor/v107-<번호>-<이름>`.
Phase 간 의존이 있으므로 **머지 순서 고정**, 병렬 진행은 Phase 3의 PR-5~8만 제한적으로 허용(충돌 범위가 도메인별로 분리되므로).

| PR | 브랜치 | 내용 | 예상 규모 | 위험도 | 문서 |
|---|---|---|---|---|---|
| PR-1 | `refactor/v107-01-safety-net` | Phase 0: OpenAPI 스냅샷 + baseline 테스트 | +300줄 | 낮음 | [phase-0](./phase-0-safety-net.md) |
| PR-2 | `refactor/v107-02-app-assembly` | Phase 1: main.py 분해, worker 분리 | ±500줄 | 낮음 | [phase-1](./phase-1-app-assembly.md) |
| PR-3 | `refactor/v107-03-infra-skeleton` | Phase 2a: core/db/shared/worker 이동 | 이동 위주 | 낮음 | [phase-2](./phase-2-domain-skeleton.md) |
| PR-4 | `refactor/v107-04-domain-skeleton` | Phase 2b: 도메인 이동 + models 분리 + shim 제거 | 이동 위주 | 중간 | [phase-2](./phase-2-domain-skeleton.md) |
| PR-5 | `refactor/v107-05-plugin-service` | Phase 3a: plugin 서비스 추출 | ±1,600줄 | 중간 | [phase-3](./phase-3-service-layer.md) |
| PR-6 | `refactor/v107-06-workflow-service` | Phase 3b: workflow 서비스 추출 | ±1,000줄 | 중간 | [phase-3](./phase-3-service-layer.md) |
| PR-7 | `refactor/v107-07-task-service` | Phase 3c: task 서비스 + SSE 분리 + Celery kwargs 타입화 | ±1,400줄 | 중간 | [phase-3](./phase-3-service-layer.md) |
| PR-8 | `refactor/v107-08-remaining-services` | Phase 3d: files/admin/auth/datatable + 전역 예외 핸들러 | ±1,500줄 | 중간 | [phase-3](./phase-3-service-layer.md) |
| PR-9 | `refactor/v107-09-split-plugin-utils` | Phase 4a: plugin_utils.py(2,220줄) 분할 | ±2,200줄 | 중간 | [phase-4](./phase-4-monolith-split.md) |
| PR-10 | `refactor/v107-10-split-dag-parser` | Phase 4b: snakefile_dag_parser.py(2,227줄) 분할 + 지연 import 제거 | ±2,300줄 | 중간 | [phase-4](./phase-4-monolith-split.md) |

Phase 5(프레임워크 업그레이드 등)는 1.0.7 범위 밖 → [phase-5-backlog.md](./phase-5-backlog.md)

### 진행 흐름
```
release/v1.0.7
  ← PR-1 (안전망)          ─ 이후 모든 PR의 게이트 기준 확립
  ← PR-2 (앱 조립부)        ─ 워커 진입점 변경 포함 (배포 주의)
  ← PR-3 (인프라 골격)      ─┐ 기계적 이동
  ← PR-4 (도메인 골격)      ─┘
  ← PR-5 ~ PR-8 (서비스 추출) ─ 도메인별 독립, 5→8 순 권장 (병렬 가능)
  ← PR-9 ~ PR-10 (모놀리스 분할)
```

---

## 4. 공통 머지 게이트 (모든 PR 필수)

머지 전 아래 4개를 전부 통과해야 한다. PR 본문에 결과 첨부.

1. **테스트**: `pytest backend/tests` 전체 통과 (baseline 대비 실패 증가 0)
2. **API 계약 동결**: OpenAPI 스냅샷 diff == 0 (PR-1에서 구축하는 테스트)
   - 예외: 스키마 이름 변경 등 의도된 diff는 PR 본문에 명시하고 스냅샷 갱신 커밋 분리
3. **DB 스키마 동결**: `alembic revision --autogenerate` 결과가 빈 마이그레이션 (모델 이동이 스키마를 바꾸지 않았음을 증명, 확인 후 파일 삭제)
4. **스모크**: 웹 + Celery 워커 컨테이너 기동 확인 (`docker compose up` 후 `/routes/version` 응답, 워커 ready 로그)

### 롤백 정책
- PR 단위 `git revert`로 롤백 가능해야 함 → PR 간 커밋 섞임 금지, 순차 머지 유지
- PR-2(워커 진입점 변경)만 배포 절차 필요: 큐 드레인 후 배포 (문서 내 상세)

## 5. 비목표 (1.0.7에서 하지 않는 것)
- API 경로·요청/응답 스키마 변경 (프론트엔드 호환 — OpenAPI diff 0)
- DB 스키마 변경 (Alembic revision 생성 없음)
- FastAPI/Pydantic/SQLAlchemy 업그레이드, async 전환 → Phase 5 백로그
- 인증 개편 (refresh token 등) → Phase 5 백로그
- 성능 최적화 (구조 개선과 분리)

## 6. 일정 요약
| Phase | PR | 예상 |
|---|---|---|
| 0 | PR-1 | 0.5~1일 |
| 1 | PR-2 | 1일 |
| 2 | PR-3, PR-4 | 1~2일 |
| 3 | PR-5~8 | 3~5일 |
| 4 | PR-9, PR-10 | 2~3일 |
| **합계** | 10 PR | **8~12 작업일** |
