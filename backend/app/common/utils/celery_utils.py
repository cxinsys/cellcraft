from celery import current_app as current_celery_app
from celery.result import AsyncResult
from celery.signals import worker_ready
import logging

from app.common.config import settings

logger = logging.getLogger(__name__)


@worker_ready.connect
def on_worker_ready(sender, **kwargs):
    """Initialize resource totals and clean up stale allocations on worker startup."""
    from app.common.utils.resource_manager import (
        initialize_resource_totals, cleanup_stale_resources
    )
    result = initialize_resource_totals()
    cleanup_stale_resources()
    logger.info(
        f"Worker ready — resources initialized: "
        f"CPU={result['cpu_total']}, GPU={result['gpu_total']}, "
        f"concurrency={result['worker_concurrency']}"
    )


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

    return celery_app

def get_task_info(task_id):
    """
    return task info for the given task_id
    """
    task_result = AsyncResult(task_id)
    result = {
        "task_id": task_id,
        "task_status": task_result.status,
        "task_result": task_result.result,
    }
    return result
