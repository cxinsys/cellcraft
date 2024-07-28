from celery import current_app as current_celery_app
from celery.result import AsyncResult

from app.common.config import settings

def create_celery():
    celery_app = current_celery_app
    celery_app.config_from_object(settings, namespace='CELERY')
    celery_app.conf.update(task_track_started=True)
    celery_app.conf.update(task_acks_late=True)
    celery_app.conf.update(task_serializer='json')
    celery_app.conf.update(result_serializer='json')
    celery_app.conf.update(accept_content=['json'])
    celery_app.conf.update(enable_unsafe_serializers=False)
    celery_app.conf.update(result_expires=200)
    celery_app.conf.update(result_persistent=True)
    celery_app.conf.update(worker_send_task_events=False)
    celery_app.conf.update(worker_prefetch_multiplier=1)
    # 긴 작업을 위한 타임아웃 설정 추가
    celery_app.conf.update(broker_transport_options={
        'visibility_timeout': 400000,  # 긴 작업을 위한 visibility timeout 설정 (초 단위)
        'confirm_publish': True,  # 메시지 발행 확인 활성화
        'confirm_timeout': 60.0  # 메시지 발행 확인을 위한 타임아웃 설정 (초 단위)
    })
    # Celery 설정에서 ampq 연결 끊김 방지를 위한 연결 관련 옵션 조정
    celery_app.conf.update(
        broker_heartbeat=0,  # 브로커 하트비트 비활성화
        broker_connection_timeout=60,  # 연결 시도 제한 시간 (초 단위)
        broker_connection_retry=True,  # 연결 실패 시 재시도 활성화
        broker_connection_max_retries=10,  # 최대 재시도 횟수
        broker_connection_retry_delay=1,  # 재시도 간격 (초 단위)
        broker_connection_retry_jitter=False,  # 재시도 간격 랜덤화 비활성화
    )
    celery_app.conf.broker_transport_options = {'confirm_publish': True, 'confirm_timeout': 10.0}


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