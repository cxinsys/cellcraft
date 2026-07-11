# 리팩토링 PR 본문 템플릿

모든 `refactor/v107-*` PR은 아래 템플릿을 사용한다. base 브랜치: `release/v1.0.7`.

```markdown
## refactor(v1.0.7): <PR 제목>

**Phase**: <0~4> | **PR 번호**: PR-<n> | **계획 문서**: docs/refactoring/v1.0.7/phase-<n>-*.md
**선행 PR**: PR-<n-1> (머지 완료 여부: ✅/⬜)

### 범위
- <이 PR이 하는 것>

### 범위 아님
- <이 PR이 의도적으로 하지 않는 것 — 후속 PR 번호 명시>

### 변경 유형
- [ ] 파일 이동만 (로직 변경 0)
- [ ] 서비스 추출 (동작 보존)
- [ ] 모듈 분할 (시그니처 보존)
- [ ] 신규 테스트/도구

### 머지 게이트 결과 (필수 첨부)
| 게이트 | 결과 |
|---|---|
| 1. pytest 전체 (baseline 대비 실패 증가 0) | ⬜ 통과 / 로그: |
| 2. OpenAPI 스냅샷 diff == 0 | ⬜ 통과 (의도된 diff 있으면 사유:) |
| 3. alembic autogenerate == 빈 마이그레이션 | ⬜ 통과 |
| 4. docker compose 스모크 (웹 + 워커 기동) | ⬜ 통과 |
| (해당 시) celery inspect registered 이름 불변 | ⬜ 통과 |
| (해당 시) import-linter 계약 | ⬜ 통과 |

### 리뷰 가이드
- 이동 커밋과 로직 커밋이 분리되어 있음: <커밋 해시 안내>
- 이동 매핑표: <phase 문서 링크 또는 표>

### 롤백
- `git revert <merge commit>` 로 롤백 가능 (DB 스키마 변경 없음)
- <PR-2만: Dockerfile CMD 원복 절차>
```

## 커밋 규칙
- 타입: `refactor:` (도구/테스트 추가는 `test:`, `chore:`)
- 이동 커밋은 `git mv`만 포함, 메시지에 `[move-only]` 표기
- 예: `refactor: [move-only] common/utils 플러그인 모듈을 plugin/ 도메인으로 이동`
