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

### 1. `workflow/compiler/dag_parser.py` (2,227줄) 분할 ✅

**경계 조정 사유**: 실제 코드는 문서의 예상 4분할(parsing/dag/rules/serialize)과
맞지 않았다. `SnakemakeDAGParser`(약 510줄) 단일 클래스가 파싱/DAG 엣지/직렬화를
이미 하나의 응집 단위로 처리하므로 클래스를 4개 파일로 쪼갤 수 없었고, 파일의 대부분
(약 1,455줄)은 문서 그룹에 없던 `SnakemakeRuleStatusTracker`(런타임 상태 추적)였다.
이 클래스 하나가 800줄을 초과하여, 시그니처/동작을 보존하면서 믹스인으로 재분할했다.

| 실제 그룹 | 이동 후 | 줄수 |
|---|---|---|
| 캐시 (`DAGCache`/`LogFileCache` + 캐시 관리 함수) | `workflow/compiler/cache.py` | 258 |
| Snakefile 파싱/DAG/직렬화 코어 (`DAGParsingError`, `SnakemakeDAGParser`) | `workflow/compiler/parsing.py` | 536 |
| 상태 추적기 본체 (오케스트레이션 + 로그 IO + 후처리) | `workflow/compiler/status.py` | 608 |
| 개별 rule 판정 믹스인 (`_RuleStatusAnalysisMixin`) | `workflow/compiler/status_analysis.py` | 538 |
| 진행률/타이밍/병목 리포트 믹스인 (`_ProgressReportMixin`) | `workflow/compiler/status_progress.py` | 364 |
| 공개 API 파사드 (기존 심볼 re-export) | `workflow/compiler/dag_parser.py` | 35 |

- **파사드 유지**: 외부 호출부(`task/service.py`)와 테스트 patch 경로
  (`app.workflow.compiler.dag_parser.SnakemakeDAGParser` 등)는 그대로 동작.
- `SnakemakeRuleStatusTracker(_RuleStatusAnalysisMixin, _ProgressReportMixin)` —
  믹스인 메서드/속성은 합성된 인스턴스의 MRO로 해석되어 동작 불변.
- 지연 `import re` 1건(`_check_workflow_completion`)은 `status.py` 상단으로 승격.

### 2. `workflow/utils.py` (1,215줄) 검토 ✅
- [x] compile 파이프라인(Snakefile rule 추출/생성) → `compiler/visualization.py`로 이동.
      나머지(데이터 로딩/그래프 탐색/파일 검증)는 service 흡수 후보로만 보고, 이동 최소화.
- [x] 800줄 초과 해소: utils.py 1,216 → 685줄 (facade re-export 포함). 이동 함수
      7개(`extract_rule_block`, `generate_visualization_snakefile`,
      `validate_visualization_rule`, `get_rule_metadata`, `_get_available_rules`,
      `_apply_replacements`, `_extract_filename_from_pattern`) → `compiler/visualization.py` (549줄).
- **순환 방지**: 공유 예외 계층을 `workflow/errors.py`(44줄)로 분리 —
      utils(잔류 함수)와 visualization(이동 함수)이 예외를 공유하는데 utils가 visualization을
      re-export 하므로, 예외를 leaf 모듈로 빼 순환을 제거. utils가 예외/함수를 모두 re-export.

### 3. `file/storage.py` + `file/path_resolver.py` 통합 검토 (PR-4에서 보류한 항목) ✅
- [x] **중복 아님 → 현상 유지**. `storage.py`(273줄)는 로그 경로 구성·파일명 sanitize
      (`construct_logs_path`, `is_safe_path`, `sanitize_filename`, `validate_export_filename`),
      `path_resolver.py`(31줄)는 데이터 파일 경로 해석(`resolve_data_file_path`, user/tutorials)로
      책임이 다르다. 둘 다 800줄 이하이므로 통합/shim 불필요.

### 4. 지연 import 정리 (전역) ✅
- [x] 함수 내부 import 전수 조사: 총 81건 (아래 보고 표 참조).
- [x] 제거 1건: `dag_parser`의 `import re`(→ `status.py` 상단 승격, 재작성분에 흡수).
      나머지 80건은 (a) 순환 방지, (b) 무거운 모듈 지연 로딩(scanpy/docker/GPUtil/psutil/
      anndata/scipy/filelock), (c) `app.main` 지연(worker no-web 계약 + 테스트 patch 호환),
      (d) config late-binding, (e) 로컬 stdlib 로 정당한 사유가 있어 유지.
      `task/service.py`의 `from app.main import get_celery_app`에는 사유 주석을 추가.
- [x] import-linter 계약 통과 (`lint-imports`): no-web-in-worker, shared-is-domain-agnostic.
      (crud 접근 금지 계약은 `.importlinter`에 아직 미정의 — 별도 백로그.)

### 체크리스트
- [x] workflow compile characterization 테스트 통과 + 파서 출력 스냅샷 고정
      (`tests/unit/test_dag_parser_snapshot.py` 신규 — 대표 3-rule Snakefile 파서 출력 pin, 6 tests)
- [x] 공통 머지 게이트 4종 (subset 회귀 동일, 전체 pytest 집합 동일, lint-imports, WORKER OK)

---

## Phase 4 완료 = v1.0.7 리팩토링 완료 기준
- [ ] `backend/app` 내 800줄 초과 파일 0개 — **미충족 (PR-10 범위 밖)**. PR-10으로
      compiler/utils 계열은 전부 800줄 이하로 정리됐으나, Phase 3 서비스 추출로 생성된
      service.py 3개가 초과 상태다: `plugin/service.py` 1,531줄 / `task/service.py` 1,295줄 /
      `workflow/service.py` 1,014줄 (실측 2026-07). **서비스 계층 분할은 별도 백로그**.
- [x] 함수 내부 지연 import 0건 (사유 명시된 예외 제외) — 정당 사유 없는 지연 import는
      제거(dag_parser의 `import re`). 잔여 80건은 순환 방지/무거운 모듈/`app.main`/config
      late-binding 등 사유가 명시된 예외 (PR-10 §4 표).
- [x] import-linter 전 계약 통과 — `.importlinter`에 정의된 2개 계약(no-web-in-worker,
      shared-is-domain-agnostic) 통과. (crud 접근 금지 계약은 미정의 상태 = 백로그.)
- [x] OpenAPI 스냅샷: Phase 0 대비 의도된 diff 1건만 존재 (PR-4에서 문서화 — 동명
      클래스 PluginInfo×2의 자동 생성 컴포넌트 키 개명, 응답 페이로드 불변). 그 외 diff 0.
- [x] DB 스키마: 신규 revision 0개. (0001 마이그레이션의 fresh-install 버그 수정과
      모델-스키마 정합화 2건은 스키마 변경이 아닌 선언 수정 — PR-1·PR-4에 기록)
- [~] 최종 회귀: 전체 pytest ✅ (576 passed / 68 failed = baseline 집합 동일) +
      격리 스택 스모크 ✅. **실제 워크플로우 1건 end-to-end 실행은 rabbitmq/워커 포함
      풀스택 필요 — 릴리즈 전 dev 환경에서 수동 확인 항목으로 이관.**
