# Phase 2 — 도메인 골격 이동 (PR-3, PR-4)

> 예상: 1~2일 | 위험도: 낮음(PR-3) / 중간(PR-4) | 선행: PR-2 머지
> **원칙: 로직 변경 0. 파일 이동 + import 경로 갱신만.** diff는 커도 리뷰 부담은 낮다 — 리뷰어는 "이동 매핑표와 일치하는가"만 확인.

## 이동 절차 (두 PR 공통)
1. 새 경로에 파일 이동 (`git mv` — 히스토리 보존)
2. 구 경로에 re-export shim 생성: `from app.<new>.<module> import *  # noqa — v1.0.7 리팩토링 shim, PR-4에서 제거`
3. `app/` 전체 import 일괄 갱신 (구 경로 참조 제거)
4. shim 제거 (PR-4 마지막 커밋)
5. 커밋 분리: "이동 커밋"(git mv만)과 "import 갱신 커밋"을 나누면 리뷰가 쉬워짐

---

## PR-3: 인프라 골격 (`refactor/v107-03-infra-skeleton`)

도메인과 무관한 공용 계층을 먼저 옮긴다. 도메인 코드는 건드리지 않음.

### 이동 매핑
| 현재 | 이동 후 | 비고 |
|---|---|---|
| `common/config.py` | `core/config.py` | |
| `common/security.py` | `core/security.py` | |
| `common/enums.py` | `core/enums.py` | |
| `common/utils/error_utils.py` | `core/exceptions.py` | 내용 그대로, 예외 계층 재설계는 PR-8 |
| `common/utils/security_logging.py` | `core/logging.py` | |
| `common/utils/docker_utils.py` | `shared/docker.py` | `container_manager` 싱글턴 그대로 (DI화는 Phase 5) |
| `common/utils/redis_cache.py` | `shared/redis.py` | |
| `common/utils/cache_utils.py` | `shared/cache.py` | |
| `common/utils/resource_manager.py` | `shared/resources.py` | |
| `common/utils/log_archive_utils.py` | `shared/archive.py` | |
| `database/conn.py` | `db/session.py` | `initialize_plugins_from_csv()`는 임시로 함께 이동, PR-4에서 `plugin/sync.py`로 재배치 |
| `database/crud/base.py` | `db/base.py` | `Base`(declarative_base)와 `CRUDBase` 함께 |
| `routes/dep.py`의 `get_db` | `db/session.py` | 인증 dependency는 PR-4에서 auth로 |
| `common/utils/celery_utils.py` (shim) | 삭제 → `worker/app.py` | PR-2에서 이동 완료, shim만 제거 |
| `routes/celery_tasks.py` | `worker/tasks.py` | `name=` 고정 완료 상태(PR-2)라 안전. Dockerfile의 `include`/autodiscover 경로 확인 |

### 체크리스트
- [x] 위 매핑대로 `git mv` + shim + import 갱신
- [x] `grep -rn "from app.common\|from app.database.conn\|from app.routes.celery_tasks" backend/app` 결과 0 (shim 파일 제외)
- [x] 공통 머지 게이트 4종 통과

---

## PR-4: 도메인 골격 + models 분리 (`refactor/v107-04-domain-skeleton`)

### 이동 매핑 — 라우터/CRUD/스키마
| 도메인 | 현재 | 이동 후 |
|---|---|---|
| api 집계 | `routes/api.py` | `api.py` (app 루트) |
| auth | `routes/endpoints/auth.py` | `auth/router.py` |
| | `routes/dep.py`의 `get_current_user`, `get_current_active_user`, `reusable_oauth2` | `auth/deps.py` |
| | `app/model.py` (Token, TokenPayload, UserInfo) | `auth/schemas.py` |
| user | `database/crud/crud_user.py` | `user/crud.py` |
| | `database/schemas/user.py` | `user/schemas.py` |
| plugin | `routes/endpoints/plugin.py` | `plugin/router.py` |
| | `database/crud/crud_plugin.py` | `plugin/crud.py` |
| | `database/schemas/plugin.py` | `plugin/schemas.py` |
| | `common/utils/plugin_utils.py` | `plugin/utils.py` (분할은 PR-9) |
| | `common/utils/github_registry_client.py` | `plugin/registry.py` |
| | `common/utils/plugin_sync_manager.py` | `plugin/sync.py` |
| | `common/utils/plugin_version_validator.py` | `plugin/versioning.py` |
| | `common/utils/plugin_cache.py` | `plugin/cache.py` |
| | `db/session.py`의 `initialize_plugins_from_csv` | `plugin/sync.py` |
| workflow | `routes/endpoints/workflow.py` | `workflow/router.py` |
| | `database/crud/crud_workflow.py` | `workflow/crud.py` |
| | `database/schemas/workflow.py` | `workflow/schemas.py` |
| | `common/utils/snakefile_dag_parser.py` | `workflow/compiler/dag_parser.py` (분할은 PR-10) |
| | `common/utils/snakemake_native_parser.py` | `workflow/compiler/native_parser.py` |
| | `common/utils/snakemake_utils.py` | `workflow/compiler/snakefile.py` |
| | `common/utils/workflow_utils.py` | `workflow/utils.py` |
| task | `routes/endpoints/task.py` | `task/router.py` |
| | `database/crud/crud_task.py` | `task/crud.py` |
| | `database/schemas/task.py` | `task/schemas.py` |
| file | `routes/endpoints/files.py` | `file/router.py` |
| | `database/crud/crud_file.py` | `file/crud.py` |
| | `database/schemas/file.py` | `file/schemas.py` |
| | `common/utils/file_security.py` | `file/security.py` |
| | `common/utils/path_utils.py` | `file/storage.py` |
| | `common/utils/file_path_resolver.py` | `file/path_resolver.py` (storage와 통합은 PR-10에서 검토) |
| datatable | `routes/endpoints/datatable.py` | `datatable/router.py` |
| | `common/utils/datatable_utils.py` | `datatable/utils.py` |
| | `common/utils/h5ad_utils.py` | `datatable/h5ad.py` |
| | `common/utils/h5ad_compression.py` | `datatable/h5ad_compression.py` |
| | `common/utils/h5ad_metadata_cache.py` + `datatable_cache.py` | `datatable/cache.py` (단순 병합 — 두 파일을 한 모듈로, 코드 수정 없이) |
| admin | `routes/endpoints/admin.py` | `admin/router.py` |
| | `database/crud/crud_admin.py` | `admin/crud.py` |
| | `database/schemas/admin.py` | `admin/schemas.py` |
| version | `routes/endpoints/version.py` | `version/router.py` |

### models 분리 (`database/models.py` → 도메인별)
| 모델 | 이동 후 |
|---|---|
| `User` | `user/models.py` |
| `File` | `file/models.py` |
| `Workflow` | `workflow/models.py` |
| `Task` | `task/models.py` |
| `Plugin` | `plugin/models.py` |

**필수 규칙**:
- [x] 모든 모델은 `db/base.py`의 단일 `Base`를 상속 (metadata 하나 유지)
- [x] `relationship()`은 문자열 참조(`relationship("Task", ...)`)로 통일 → 도메인 간 순환 import 방지
- [x] FK의 `ForeignKey("users.id")` 등 테이블명 문자열은 그대로 (변경 없음)
- [x] 전 도메인 models를 import하는 중앙 registry 추가 — Alembic `env.py`와 `create_all` 계열이 전체 메타데이터를 보도록:
```python
# app/db/base_registry.py — db/base.py 하단이 아닌 별도 모듈로 배치 (아래 참고)
from app.user.models import User        # noqa
from app.file.models import File        # noqa
from app.workflow.models import Workflow  # noqa
from app.task.models import Task        # noqa
from app.plugin.models import Plugin    # noqa
```
  > **구현 노트(순환 import 회피)**: 도메인 `models` 모듈이 `app.db.base`에서 `Base`를 import하므로, registry를 `db/base.py` 하단에 두면 `db.base` ↔ `<domain>.models` 순환이 발생한다(도메인 models가 먼저 import될 때 `partially initialized` ImportError). 따라서 registry는 별도 모듈 `app/db/base_registry.py`로 분리했다. 웹/CRUD 호출부의 `models.User` 접근을 위해서는 얇은 aggregator `app/models.py`(`from app.user.models import User` 등)를 추가해 `from app import models` 패턴을 유지했다(클래스명·테이블명 불변). `user_plugin_association` Table은 `user/models.py`에 두고 `plugin/models.py`가 참조(단방향 plugin→user).
- [x] `alembic/env.py`의 `target_metadata` import 경로 갱신 (`from app.db.base_registry import Base`)
- [x] 구 `database/models.py`는 도메인별로 분리 후 이 PR에서 제거

### 마무리
- [x] `api.py`의 라우터 import를 새 경로로 갱신 (`from app.plugin.router import router as plugin_router` 등) — **prefix/tags 변경 금지**
- [x] PR-3·PR-4에서 만든 모든 shim 제거, `routes/`·`database/`·`common/` 빈 디렉터리 삭제
- [x] 테스트 코드의 import 경로 일괄 갱신 (patch 경로 `app.routes.endpoints.X` → `app.X.router`, docstring 경로 언급 포함)
- [x] import-linter에 도메인 계약 추가 (`app.shared`는 도메인 모듈 import 금지 → `shared-is-domain-agnostic` 계약 KEPT; `no-web-in-worker`는 `app.api`+도메인 router 금지로 갱신)

### 머지 게이트 (공통 4종 + 추가)
- [x] **alembic 게이트**: `alembic upgrade head` 후 `alembic revision --autogenerate` → 앱 관리 테이블(users/plugins/files/workflows/tasks/user_plugin_association) 변경 0 → 모델 분리가 스키마에 영향 없음 증명 (autogenerate가 감지한 `celery_taskmeta`/`celery_tasksetmeta` drop은 Celery 런타임 결과-백엔드 테이블로, 모델·초기 마이그레이션 어디에도 정의되지 않은 사전-기존 noise)
- [x] `celery` 태스크 이름 불변: `['plugin_task:build_plugin_task', 'workflow_task:process_data_task']`, `app.main` 미유입
- [x] `grep -rn "app.common\|app.database\|app.routes" backend/app backend/tests` 코드 import 0건 (남은 매치는 docstring/주석의 경로 언급뿐)
  > **계약 테스트 예외 보고**: `tests/contract/test_openapi_contract.py::test_openapi_schema_unchanged`가 실패한다. 원인은 pydantic이 **동명(同名) 클래스 충돌**(`PluginInfo`가 `plugin/schemas.py`와 `task/schemas.py`에 각각 존재)을 해소할 때 `components.schemas` 키를 **모듈 경로**로 생성하기 때문이다. 모듈 이동으로 키가 `app__database__schemas__plugin__PluginInfo` → `app__plugin__schemas__PluginInfo`(및 task 동일)로 바뀌고, 이를 `$ref`하는 `PluginData`/`TaskMonitoringItem`도 연동 변경된다. **클래스명은 변경하지 않았고**(로직 변경 0), `paths`·`info`·최상위 키는 완전 동일하며 두 `PluginInfo` 키의 모듈-경로 부분만 정규화하면 `components`가 바이트 단위로 동일함을 확인했다(API 표면 불변). 지시대로 **스냅샷은 갱신하지 않았다** — 스냅샷 재생성은 PR 리뷰에서 의도된 변경으로 승인 후 별도 진행 권장.

### 리스크
| 리스크 | 완화 |
|---|---|
| 모델 간 순환 import | 문자열 relationship + `db/base.py` 중앙 registry |
| conftest/픽스처의 구 경로 import | 테스트 갱신을 같은 PR에 포함, baseline과 대조 |
| 숨은 동적 import(문자열 기반) | `grep -rn "importlib\|__import__\|import_module" backend/app`으로 사전 확인 |
