# Phase 0 — 안전망 구축 (PR-1)

> 브랜치: `refactor/v107-01-safety-net` → base: `release/v1.0.7`
> 예상: 0.5~1일 | 위험도: 낮음 | 선행 조건: 없음 (모든 후속 PR의 게이트 기준)

## 목적
리팩토링 전 과정에서 "동작이 바뀌지 않았음"을 기계적으로 증명할 수단을 만든다.
이 PR이 머지되기 전에는 어떤 리팩토링 PR도 머지하지 않는다.

## 작업 항목

### 1. OpenAPI 스냅샷 테스트
API 계약(경로·메서드·요청/응답 스키마) 동결을 검증하는 테스트.

- [x] `backend/tests/contract/__init__.py` 생성
- [x] `backend/tests/contract/openapi_snapshot.json` — 생성 완료 (paths 104개, schemas 40개)
- [x] `backend/tests/contract/test_openapi_contract.py` 작성:

```python
import json
from pathlib import Path

SNAPSHOT = Path(__file__).parent / "openapi_snapshot.json"

def test_openapi_schema_unchanged(client):  # conftest의 기존 client 픽스처 재사용
    current = client.app.openapi()
    expected = json.loads(SNAPSHOT.read_text())
    assert current == expected, (
        "OpenAPI 계약이 변경되었습니다. 의도된 변경이면 스냅샷을 갱신하고 "
        "PR 본문에 변경 사유를 명시하세요: "
        "pytest tests/contract --update-snapshot"
    )
```

- [x] 스냅샷 갱신 헬퍼(`scripts/update_openapi_snapshot.py`) 추가 — TESTING=1 설정 후 app import, `json.dump(sort_keys=True, indent=2)`
- [x] **주의**: `app.main`은 import 시점에 `run_migrations()`를 실행하므로(main.py:84) 스냅샷 생성은 반드시 기존 테스트 환경(conftest의 DB 셋업) 안에서 수행한다. PR-2에서 이 부수효과가 제거되면 단순화 가능. (스크립트/테스트가 이 조건을 문서화 및 준수)

### 2. pytest baseline 고정
- [x] 전체 실행 결과를 `docs/refactoring/v1.0.7/baseline-test-report.txt`로 저장 — **448 passed / 68 failed / 3 skipped / 3 xpassed**
- [x] 실패 68건 목록화 및 원인 분류(A: 상태코드 드리프트, B: 테스트 자체 결함, C: asyncio 루프) — 리팩토링 이전부터 존재하던 실패 (백엔드 CI 부재로 드리프트). 게이트 판정 기준: "목록 외 새로운 실패 = 리그레션"
- [ ] 커버리지 수치 기록 (Phase 3 착수 시 측정)

### 3. 핵심 경로 characterization 테스트 보강
리팩토링 대상 중 테스트가 얇은 핵심 경로에 대해, **현재 동작 그대로**를 고정하는 테스트를 추가한다 (동작이 이상해 보여도 수정하지 않고 그대로 고정 — 수정은 Phase 3에서).

| 대상 | 고정할 동작 |
|---|---|
| `POST /plugin` upload 경로 | 성공 응답 형태, 잘못된 입력 시 상태코드/에러 메시지, 롤백 후 DB 상태 |
| `POST /workflow` compile 경로 | 컴파일 성공/실패 응답, 캐시 적중 시 응답 동일성 |
| `GET /task` 상태 조회(SSE) | 이벤트 포맷, 종료 조건 (Docker/Celery는 conftest 패턴대로 모킹) |
| `POST /files` 업로드 | 경로 검증 거부 케이스(file_security), 성공 시 DB 레코드 |
| `POST /auth/login/access-token` | 토큰 발급 형태, 비활성 유저 거부 |

- [x] 각 테스트는 `tests/integration/` 기존 구조와 픽스처를 따른다 (5개 파일: auth/plugin_upload/workflow_compile/task_status/file_upload)
- [x] Docker·GitHub 외부 의존은 모킹 (`unittest.mock.patch`로 `get_plugin_path`, `process_data_task.apply_async`, `plugin_utils.*`, `get_task_info` 등 모킹)

### 4. import 계약 린트 (선택이지만 권장)
- [x] `import-linter` 도입, `backend/.importlinter` 계약 정의 (Phase 2 이후 활성화할 계약을 미리 작성 — `no-web-in-worker` 계약을 주석으로 정의하고 PR-2 활성화 시점 명시):

```ini
[importlinter]
root_package = app

[importlinter:contract:no-web-in-worker]
name = worker는 웹 계층에 의존하지 않는다
type = forbidden
source_modules = app.worker
forbidden_modules = app.main, app.api
```

- 계약은 PR-2~4 진행에 맞춰 단계적으로 추가 (도메인 단방향 의존 등)

### 5. 테스트 실행 방법 (compose 오버라이드)
- [x] `docker-compose.test.yml` 신규 작성 — 격리된 `test-network` + `test-db`(postgres:15, `cellcraft_test`/`test_user`/`test_pass`, 호스트 포트 `${TEST_DB_PORT:-5433}`) + `backend-test`(image `cellcraft-backend:refactor-test`, `./backend` 마운트, `TESTING=1`/`TEST_DATABASE_URI` 설정, `command: sleep infinity` — 오케스트레이터가 `exec`으로 pytest 실행)
- [x] `backend/tests/conftest.py`의 테스트 DB URI를 `os.environ.get("TEST_DATABASE_URI", <기존 기본값>)`으로 환경변수화 (기존 환경 호환)
- 상세 실행/스냅샷 갱신 방법: `backend/tests/contract/README.md` 참조.

```bash
# 이미지 빌드 → 스택 기동 → 스냅샷 생성 → 테스트
docker build -f backend/Dockerfile.dev -t cellcraft-backend:refactor-test ./backend
docker compose -f docker-compose.test.yml up -d
docker compose -f docker-compose.test.yml exec backend-test python scripts/update_openapi_snapshot.py
docker compose -f docker-compose.test.yml exec backend-test pytest
docker compose -f docker-compose.test.yml down
# 호스트 5433 점유 시: TEST_DB_PORT=5544 docker compose -f docker-compose.test.yml up -d
```

## 머지 게이트
- [ ] `pytest backend/tests` 통과 (신규 테스트 포함)
- [ ] 스냅샷 테스트가 의도적 변경 시 명확한 에러 메시지와 갱신 방법을 출력
- [ ] baseline 리포트 커밋됨

## 산출물
- [x] `tests/contract/` (계약 테스트 + `README.md`; 스냅샷 JSON은 컨테이너에서 생성)
- [x] `backend/scripts/update_openapi_snapshot.py` (스냅샷 생성/갱신 헬퍼)
- [x] `tests/integration/` characterization 테스트 5건 (auth/plugin_upload/workflow_compile/task_status/file_upload)
- [x] `docker-compose.test.yml` (격리 테스트 스택)
- [x] `backend/.importlinter` (PR-2 활성화 예정 계약 정의)
- [x] `docs/refactoring/v1.0.7/baseline-test-report.txt`
