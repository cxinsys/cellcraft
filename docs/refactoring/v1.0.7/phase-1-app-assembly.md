# Phase 1 — 앱 조립부 정리 (PR-2)

> 브랜치: `refactor/v107-02-app-assembly` → base: `release/v1.0.7`
> 예상: 1일 | 위험도: 낮음(코드) / 중간(배포 — 워커 진입점 변경) | 선행: PR-1 머지

## 목적
`main.py`(337줄)를 순수 조립 코드(~60줄)로 축소하고, Celery 워커를 웹 계층에서 분리한다.
**로직 자체는 옮기기만 하고 바꾸지 않는다.**

## 현재 main.py의 책임 (분해 대상)
| 현재 위치 | 내용 | 이동 목적지 |
|---|---|---|
| `main.py:84` | import 시점 `run_migrations()` 실행 | `core/startup.py` + 명시 호출 |
| `main.py:25-38` | 시그널 핸들러 (컨테이너 정리) | `core/startup.py` |
| `main.py:41-52, 112-118` | Celery 앱 생성 + worker_shutting_down 훅 | `worker/app.py` |
| `main.py:130-219` | startup_event (플러그인 sync/버전 검증/CSV 초기화) | `core/startup.py` |
| `main.py:221-338` | Docker 이미지 pull (_sync_pull_image 등) | `core/startup.py` (또는 `shared/docker.py`와 인접 배치) |

## 작업 항목

### 1. `app/core/startup.py` 신설
- [x] `run_migrations()` 이동 (내용 그대로)
- [x] `setup_signal_handlers()` 이동
- [x] `startup_event`의 본문 → `on_startup()` 함수로 이동 (`check_and_pull_official_plugin_images()`/`_sync_pull_image` 포함)
- [x] print 기반 로깅은 이 단계에서 유지 (동작 변경 금지 — 로깅 정리는 별도)

### 2. `app/worker/` 신설 — Celery 분리
- [x] `worker/__init__.py`, `worker/app.py` 생성
- [x] `common/utils/celery_utils.py`의 `create_celery()` → `worker/app.py`로 이동, 모듈 레벨에서 `celery = create_celery()` 인스턴스 노출
- [x] `worker_shutting_down` 훅(main.py:41)을 `worker/app.py`로 이동
- [x] `worker/app.py`는 FastAPI·`app.main`을 import하지 않을 것 (import-linter 계약 활성화, `app.main not in sys.modules` 검증 통과)
- [x] 기존 `common/utils/celery_utils.py`는 re-export shim으로 유지 (`from app.worker.app import create_celery`) — Phase 2에서 일괄 제거. `get_task_info`/`on_worker_ready`는 그대로 잔류
- [x] **태스크 이름 고정**: `routes/celery_tasks.py`의 두 태스크는 이미 `name="plugin_task:build_plugin_task"`, `name="workflow_task:process_data_task"`로 명시되어 있음(변경 불필요). `worker/app.py`는 import 시점에 `app.routes.celery_tasks`를 로드해 레지스트리를 채움 — 태스크 이름 목록 불변 검증 통과

### 3. Dockerfile CMD 변경
- [x] `Dockerfile.celery`, `Dockerfile.celery.prod`, `Dockerfile.dev.celery`, `Dockerfile.cpu.celery`, `Dockerfile.gpu.celery`:
  `celery -A app.main.celery worker` → **`celery -A app.worker.app worker`** (5개 파일 모두, 옵션/플래그 불변)
- [x] 워커가 더 이상 `run_migrations()`·웹 앱 조립을 실행하지 않음 확인 — `app.worker.app` import 시 `app.main not in sys.modules`, 마이그레이션 미실행
- [x] 마이그레이션은 **웹 컨테이너만** 수행: `run_migrations()`는 `app.core.startup`으로 이동해 `create_app()`(웹 진입점)에서만 명시 호출. 워커 진입점 `app.worker.app`은 이를 import/실행하지 않음

### 4. `main.py` 재작성 (~60줄)
```python
from fastapi import FastAPI
from starlette.middleware.cors import CORSMiddleware

from app.routes.api import api_router          # PR-4에서 app.api로 이동
from app.common.config import settings          # PR-3에서 app.core.config로 이동
from app.core import startup

def create_app() -> FastAPI:
    startup.run_migrations()
    is_production = settings.ENVIRONMENT == "production"
    app = FastAPI(
        title=settings.PROJECT_NAME,
        version=settings.APP_VERSION,
        description=settings.APP_DESCRIPTION,
        docs_url=None if is_production else "/docs",
        redoc_url=None if is_production else "/redoc",
        openapi_url=None if is_production else "/openapi.json",
    )
    app.add_middleware(CORSMiddleware, ...)      # 기존 설정 그대로
    app.include_router(api_router, prefix=settings.ROUTES_STR)
    app.add_event_handler("startup", startup.on_startup)   # FastAPI 0.78: on_event 방식 유지
    startup.setup_signal_handlers()
    return app

app = create_app()
```
- [x] `run_migrations()`는 import 부수효과가 아니라 `create_app()` 내부의 명시 호출로 변경
- [x] uvicorn 진입점(`app.main:app`)은 그대로 유지 → 웹 Dockerfile 변경 불필요 (`app = create_app()` 모듈 레벨 유지)
- [x] `app.celery_app`/`get_celery_app` 참조 조사 후 배선: `task.py`가 `from app.main import get_celery_app` 사용 → `main.py`에 `get_celery_app()`·`app.celery_app`·모듈 변수 `celery` 유지, 모두 `app.worker.app.celery` 재수출로 연결 (테스트가 `@patch('app.main.get_celery_app')`에 의존)

### 5. 테스트 수정
- [x] conftest가 `app.main` import에 의존하는 부분 점검 — `tests/conftest.py:20 from app.main import app`은 `app = create_app()` 모듈 레벨 유지로 그대로 동작. 테스트 수정 불필요 (기존 동작 유지)

### 6. import-linter 계약 활성화
- [x] `backend/.importlinter`의 `no-web-in-worker` 계약 주석 해제 및 활성화. `forbidden_modules = app.main, app.routes.api`. `app.routes` 전체 금지는 `celery_tasks.py`가 `app.routes.celery_tasks`에 위치하는 한 불가 — 사유 주석 명시. `lint-imports` → 1 kept, 0 broken

## 머지 게이트 (공통 4종 + 추가)
- [x] pytest: `448 passed`(baseline 유지), 실패 68건은 baseline과 정확히 동일한 집합 (신규 실패 0, diff empty)
- [x] contract+characterization: `30 passed`
- [x] **워커 스모크**: `app.worker.app` import 시 ① `app.main not in sys.modules`(마이그레이션 미실행) ② 태스크 이름 목록 `['plugin_task:build_plugin_task', 'workflow_task:process_data_task']` 불변 검증 통과. (실제 태스크 1건 실행은 오케스트레이터의 컴포즈 스모크에서 확인)
- [x] import-linter `no-web-in-worker` 계약 통과 (1 kept, 0 broken)

## 배포 주의사항 (이 PR만 해당)
워커 진입점이 바뀌므로 운영 배포 시:
1. 이전 버전 워커의 큐 드레인(진행 중 태스크 완료 대기)
2. 웹(마이그레이션 담당) → 워커 순으로 기동
3. 태스크 이름이 `name=`으로 고정되었으므로 드레인 실패한 잔여 메시지도 신버전 워커가 처리 가능

## 롤백
`git revert` + Dockerfile CMD 원복. DB/스키마 변경이 없으므로 데이터 롤백 불필요.
