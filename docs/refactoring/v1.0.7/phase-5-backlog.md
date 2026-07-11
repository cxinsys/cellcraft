# Phase 5 — 백로그 (v1.0.7 범위 제외)

v1.0.7 리팩토링(Phase 0~4)이 완료된 후, 별도 릴리즈에서 진행할 항목.
Phase 0~4가 만든 구조(서비스 계층, 도메인 경계, 계약 테스트) 덕분에 아래 작업들의 위험이 크게 줄어든다.

## 1. 프레임워크 업그레이드 (최우선 백로그, 일괄 진행 필요)

| 항목 | 현재 | 목표 | 비고 |
|---|---|---|---|
| FastAPI | 0.78.0 (2022) | 최신 | Pydantic v2와 함께 올려야 함 |
| Pydantic | 1.9.1 | 2.x | 전 스키마 마이그레이션 (`bump-pydantic` 도구 활용). Phase 2~3에서 스키마가 도메인별로 분리되어 있어 도메인 단위 진행 가능 |
| SQLAlchemy | 1.4.36 | 2.0 | 1.4는 2.0 스타일 과도기 지원 → 먼저 1.4에서 `future=True`로 2.0 스타일 쿼리로 전환 후 버전 업 |
| Starlette/uvicorn | 0.19/0.17 | FastAPI에 종속 | |

**순서 권고**: SQLAlchemy 2.0 스타일 전환(1.4 유지) → Pydantic v2 + FastAPI 동시 업그레이드 → `@app.on_event` → lifespan 전환

**주의**: conda `environment.yml`(snakemake 7.14 고정)과의 의존성 충돌 사전 검증 필요.

## 2. 비동기 전환 (선택적·부분적)
- 파일 업로드, SSE, 외부 API(GitHub registry) 호출 등 I/O 바운드 경로부터 부분 적용
- DB 비동기(asyncpg + SQLAlchemy async)는 전면 전환 비용이 크므로 실측(응답 지연·동시성 병목) 후 결정 — "Measure First"
- 현재 sync Celery 구조와의 일관성 고려

## 3. 인증 강화
- [ ] Refresh token 도입 (현재 32시간 단일 액세스 토큰)
- [ ] 토큰 폐기(revocation) — Redis 블랙리스트
- [ ] 로그인 rate limiting
- [ ] `SECRET_KEY` 로테이션 절차 문서화

## 4. 전역 상태 제거
- [ ] `shared/docker.py`의 `container_manager` 싱글턴 → 명시적 수명주기 관리 (app.state 또는 dependency 주입, 스레드 안전성 확보)
- [ ] Redis 커넥션 풀 lazy 생성 → 앱 기동 시 초기화 + 헬스체크

## 5. 관측성(Observability)
- [ ] `core/startup.py`의 print 로깅 → 구조화 로깅(logging + JSON formatter) 통일
- [ ] 요청 ID 미들웨어 + Celery 태스크 상관관계 ID
- [ ] Sentry 또는 동급 에러 트래킹 (uptime-kuma 모니터링과 연계)

## 6. Celery 파이프라인 개선
- [ ] 태스크 멱등성 보장 (재시도 안전)
- [ ] 엔드포인트 트랜잭션과 태스크 디스패치 사이 원자성 (transactional outbox 또는 커밋 후 디스패치 패턴)
- [ ] 죽은 코드/미사용 태스크 정리 (PR-9~10에서 기록한 죽은 코드 목록 처리)

## 7. 테스트 인프라
- [ ] Docker/GitHub 통합 테스트용 공식 모킹 레이어 (현재 개별 테스트에서 임시 모킹)
- [ ] E2E: 대표 워크플로우 실행 시나리오 자동화 (CI nightly)
- [ ] 커버리지 게이트 CI 통합 (80% 기준)
