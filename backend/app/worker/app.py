"""Celery worker entry point.

This module owns the Celery application instance and worker-lifecycle hooks.
It must NOT import the FastAPI web layer (``app.main`` / ``app.routes.api``);
the ``no-web-in-worker`` import-linter contract enforces this.

Worker Dockerfiles boot via ``celery -A app.worker.app worker``.
"""
from celery import current_app as current_celery_app
from celery.signals import worker_shutting_down

from app.common.config import settings
from app.common.utils.docker_utils import container_manager


def create_celery():
    celery_app = current_celery_app
    celery_app.config_from_object(settings, namespace='CELERY')
    celery_app.conf.update(task_track_started=True)
    celery_app.conf.update(task_acks_late=False)
    celery_app.conf.update(task_serializer='json')
    celery_app.conf.update(result_serializer='json')
    celery_app.conf.update(accept_content=['json'])
    celery_app.conf.update(enable_unsafe_serializers=False)
    celery_app.conf.update(result_expires=200)
    celery_app.conf.update(result_persistent=True)
    celery_app.conf.update(worker_send_task_events=True)
    celery_app.conf.update(worker_prefetch_multiplier=1)
    celery_app.conf.update(worker_max_tasks_per_child=5)  # 메모리 누수 방지: 5개 태스크 처리 후 워커 프로세스 재시작 (24건/일 ÷ 36워커 = ~7.5일/워커)

    # broker_transport_options 설정
    # - visibility_timeout: acks_late=False이므로 실질적 영향 없음. 방어적으로 유지.
    # - confirm_publish: 메시지 발행 확인 활성화
    # - confirm_timeout: 발행 확인 타임아웃
    celery_app.conf.update(broker_transport_options={
        'visibility_timeout': 259200,
        'confirm_publish': True,
        'confirm_timeout': 10.0
    })

    # Task Time Limits 비활성화 — GRN inference 작업은 데이터셋에 따라 수일~수주 소요 가능
    # None = 무제한. 작업 중단은 사용자가 UI에서 revoke로 수행.
    celery_app.conf.update(task_time_limit=None)
    celery_app.conf.update(task_soft_time_limit=None)

    # Celery 설정에서 ampq 연결 끊김 방지를 위한 연결 관련 옵션 조정
    celery_app.conf.update(
        broker_heartbeat=0,  # 브로커 하트비트 비활성화
        broker_connection_timeout=60,  # 연결 시도 제한 시간 (초 단위)
        broker_connection_retry=True,  # 연결 실패 시 재시도 활성화
        broker_connection_max_retries=10,  # 최대 재시도 횟수
        broker_connection_retry_delay=1,  # 재시도 간격 (초 단위)
        broker_connection_retry_jitter=False,  # 재시도 간격 랜덤화 비활성화
    )

    # 작업자(worker) 동시성 제한 설정
    celery_app.conf.update(
        task_reject_on_worker_lost=True  # 워커 손실 시 작업 거부
    )

    # 태스크 등록 보장: 워커가 app.routes.celery_tasks의 태스크를 인식하도록 include에 명시.
    # (celery_tasks.py는 웹 계층을 import하지 않으므로 워커가 FastAPI를 끌어오지 않음)
    celery_app.conf.update(include=['app.routes.celery_tasks'])

    return celery_app


celery = create_celery()

# 태스크 모듈을 import 시점에 로드하여 태스크 레지스트리를 채운다.
# (celery_tasks.py는 웹 계층을 import하지 않으므로 워커가 FastAPI를 끌어오지 않음)
import app.routes.celery_tasks  # noqa: E402,F401


# Celery 이벤트 핸들러를 정의합니다
def on_worker_shut_down(sender=None, conf=None, **kwargs):
    print("Worker is shutting down. Cleaning up containers and trying to reconnect...")
    try:
        # 컨테이너 정리
        container_manager.cleanup_all_task_containers()

        # 재연결 시도
        with celery.connection() as connection:
            connection.ensure_connection(max_retries=3)
            print("재연결에 성공했습니다.")
    except Exception as e:
        print(f"Worker shutdown 처리 중 오류 발생: {e}")


# Celery 이벤트 핸들러를 Celery 애플리케이션 인스턴스에 연결합니다
worker_shutting_down.connect(on_worker_shut_down, sender=celery)

# worker_ready 훅(on_worker_ready)은 celery_utils에 잔류(Phase-1 계획). 구 워커는
# app.main 경유로 celery_utils를 import하며 자동 등록됐으나, 신 워커 진입점(app.worker.app)은
# celery_utils를 import하지 않는다. 이 import로 @worker_ready.connect 핸들러를 명시 등록하여
# 워커 기동 시 리소스 총량 초기화/스테일 정리가 정상 수행되게 한다.
# (하단 배치 이유: celery_utils가 create_celery를 역참조하므로 celery 정의 이후에 import.)
import app.common.utils.celery_utils  # noqa: E402,F401
