# Phase 4 — 모놀리스 utils 내부 분할 (PR-9, PR-10)

> 예상: 2~3일 | 위험도: 중간 | 선행: PR-5(plugin), PR-6(workflow) 머지 (각 해당 PR만)
> **원칙: 함수 시그니처·동작 불변. 파일 내부를 책임별 모듈로 쪼개는 것이 전부.** 파일당 800줄 이하 목표.

---

## PR-9: plugin_utils.py 분할 (`refactor/v107-09-split-plugin-utils`)

대상: `plugin/utils.py` (구 `common/utils/plugin_utils.py`, 2,220줄)

### 절차
1. 파일 내 함수/클래스를 책임별로 그룹핑 (착수 시 실제 코드 기준으로 확정 — 아래는 예상 분류)
2. 그룹별로 새 모듈로 이동:

| 예상 그룹 | 이동 후 |
|---|---|
| Docker 이미지 빌드/태깅/푸시 | `plugin/builder.py` (PR-5에서 생성된 파일에 통합) |
| 플러그인 파일 저장/검증/압축 해제 | `plugin/files.py` |
| 플러그인 메타데이터 파싱(config/schema) | `plugin/metadata.py` |
| 플러그인 실행 준비(컨테이너 인자 조립 등) | `plugin/runtime.py` |
| 나머지 순수 헬퍼 | `plugin/utils.py`에 잔류 (300줄 이하 목표) |

3. `plugin/service.py`·`worker/tasks.py` 등 호출부 import 갱신
4. 각 모듈 단위 테스트 보강 (분할 과정에서 발견되는 죽은 코드는 **삭제하지 말고 목록만 기록** — 삭제는 별도 PR)

### 체크리스트
- [x] 모든 신규 모듈 800줄 이하 — files.py 456 / metadata.py 657 / runtime.py 727 / builder.py 515 / utils.py 63 (잔류 facade)
- [x] `plugin/` 내부 지연 import(함수 내부 import) 제거 — 새 모듈 함수 내부 import 0건 (create_plugin_folder의 `import tempfile`를 files.py 상단으로 승격). builder↔worker 순환은 utils facade가 builder 심볼을 re-export하지 않도록 경계 조정(builder 심볼은 `app.plugin.builder`에서 직접 import)하여 해소
- [x] 공통 머지 게이트 4종 + plugin characterization 테스트 통과 — 전체 570 passed / 68 failed(baseline 집합 동일), lint-imports 2 kept/0 broken, WORKER OK(app.main 미로드)

---

## PR-10: snakefile_dag_parser 분할 + 지연 import 정리 (`refactor/v107-10-split-dag-parser`)

### 1. `workflow/compiler/dag_parser.py` (2,227줄) 분할
| 예상 그룹 | 이동 후 |
|---|---|
| Snakefile 토크나이징/파싱 코어 | `workflow/compiler/parsing.py` |
| DAG 구성/노드·엣지 모델 | `workflow/compiler/dag.py` |
| 규칙(rule) 해석/와일드카드 처리 | `workflow/compiler/rules.py` |
| 출력 직렬화(프론트용 그래프 변환) | `workflow/compiler/serialize.py` |
| 공개 API (기존 함수 시그니처 유지) | `workflow/compiler/dag_parser.py` — 파사드로 잔류, 내부 모듈 위임 |

- **파사드 유지**: 외부 호출부는 `dag_parser`의 기존 함수명을 그대로 사용 → 호출부 수정 최소화

### 2. `workflow/utils.py` (1,215줄) 검토
- [ ] compile 파이프라인 관련은 `compiler/`로, 나머지는 service로 흡수 가능한지 분류
- [ ] 800줄 초과 시 분할, 이하로 정리되면 유지

### 3. `file/storage.py` + `file/path_resolver.py` 통합 검토 (PR-4에서 보류한 항목)
- [ ] 중복 로직 확인 후 `file/storage.py`로 통합 또는 역할 명확화

### 4. 지연 import 제거 (전역)
PR-2~9를 거치며 대부분 해소되지만 잔여분 일괄 정리:
- [ ] `grep -rn "^\s\+import \|^\s\+from " backend/app --include="*.py"`로 함수 내부 import 전수 조사
- [ ] 각 건에 대해: (a) 순환 의존이 원인이면 모듈 경계 수정, (b) 무거운 모듈 지연 로딩이 목적이면 주석으로 사유 명시 후 유지
- [ ] import-linter 최종 계약 활성화:
  - 도메인 간 직접 crud 접근 금지
  - `shared/` → 도메인 import 금지
  - `worker/` → `main`/`api` import 금지

### 체크리스트
- [ ] workflow compile characterization 테스트 통과 (파싱 결과 스냅샷 비교 권장 — 대표 Snakefile 3~5개에 대한 파서 출력 고정)
- [ ] 공통 머지 게이트 4종

---

## Phase 4 완료 = v1.0.7 리팩토링 완료 기준
- [ ] `backend/app` 내 800줄 초과 파일 0개
- [ ] 함수 내부 지연 import 0건 (사유 명시된 예외 제외)
- [ ] import-linter 전 계약 통과
- [ ] OpenAPI 스냅샷: Phase 0 대비 diff 0
- [ ] DB 스키마: Phase 0 대비 revision 0개
- [ ] 최종 회귀: 전체 pytest + docker compose 스모크 + 실제 워크플로우 1건 end-to-end 실행
