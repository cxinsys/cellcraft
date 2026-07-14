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
- [ ] `run_migrations()` 이동 (내용 그대로)
- [ ] `setup_signal_handlers()` 이동
- [ ] `startup_event`의 본문 → `initialize_plugin_system(db)` / `check_and_pull_official_plugin_images()` 함수로 이동
- [ ] print 기반 로깅은 이 단계에서 유지 (동작 변경 금지 — 로깅 정리는 별도)

### 2. `app/worker/` 신설 — Celery 분리
- [ ] `worker/__init__.py`, `worker/app.py` 생성
- [ ] `common/utils/celery_utils.py`의 `create_celery()` → `worker/app.py`로 이동, 모듈 레벨에서 `celery = create_celery()` 인스턴스 노출
- [ ] `worker_shutting_down` 훅(main.py:41)을 `worker/app.py`로 이동
- [ ] `worker/app.py`는 FastAPI·`app.main`을 import하지 않을 것 (PR-1의 import-linter 계약 활성화)
- [ ] 기존 `common/utils/celery_utils.py`는 re-export shim으로 유지 (`from app.worker.app import *`) — Phase 2에서 일괄 제거
- [ ] **태스크 이름 고정**: `routes/celery_tasks.py`의 모든 `@shared_task`/`@celery.task`에 현재 등록 이름을 `name="..."`으로 명시 (`celery -A ... inspect registered`로 현재 이름 확인 후). 이후 파일 이동에도 브로커 메시지 호환 유지

### 3. Dockerfile CMD 변경
- [ ] `Dockerfile.celery`, `Dockerfile.celery.prod`, `Dockerfile.dev.celery`, `Dockerfile.cpu.celery`, `Dockerfile.gpu.celery`:
  `celery -A app.main.celery worker` → **`celery -A app.worker.app worker`**
- [ ] 워커가 더 이상 `run_migrations()`·웹 앱 조립을 실행하지 않음을 기동 로그로 확인
- [ ] 워커에 필요한 startup 로직이 있는지 검토(현재 워커는 main.py import로 마이그레이션까지 얻고 있었음) — 마이그레이션은 **웹 컨테이너만** 수행하는 것으로 명시. 웹보다 워커가 먼저 뜨는 경우를 대비해 compose의 `depends_on` 확인

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
- [ ] `run_migrations()`는 import 부수효과가 아니라 `create_app()` 내부의 명시 호출로 변경
- [ ] uvicorn 진입점(`app.main:app`)은 그대로 유지 → 웹 Dockerfile 변경 불필요
- [ ] `app.celery_app` 속성에 의존하는 코드가 있는지 검색(`grep -rn "app.celery_app\|get_celery_app"`) 후 `app.worker.app.celery` 참조로 교체

### 5. 테스트 수정
- [ ] conftest가 `app.main` import에 의존하는 부분 점검 — `create_app()` 팩토리 덕에 테스트에서 마이그레이션 시점 제어 가능해짐 (기존 동작 유지가 우선, 개선은 선택)

## 머지 게이트 (공통 4종 + 추가)
- [ ] 공통 게이트: pytest / OpenAPI diff 0 / alembic autogenerate 빈 결과 / 스모크
- [ ] **워커 스모크 강화**: `docker compose up` 후 ① 워커 로그에 마이그레이션 출력이 없어야 함 ② `celery inspect registered`의 태스크 이름 목록이 변경 전과 동일 ③ 실제 태스크 1건 실행 확인
- [ ] import-linter `no-web-in-worker` 계약 통과

## 배포 주의사항 (이 PR만 해당)
워커 진입점이 바뀌므로 운영 배포 시:
1. 이전 버전 워커의 큐 드레인(진행 중 태스크 완료 대기)
2. 웹(마이그레이션 담당) → 워커 순으로 기동
3. 태스크 이름이 `name=`으로 고정되었으므로 드레인 실패한 잔여 메시지도 신버전 워커가 처리 가능

## 롤백
`git revert` + Dockerfile CMD 원복. DB/스키마 변경이 없으므로 데이터 롤백 불필요.
