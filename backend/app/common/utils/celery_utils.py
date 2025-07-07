from celery import current_app as current_celery_app
from celery.result import AsyncResult
from celery.signals import celeryd_after_setup
from kombu import Queue, Exchange
import threading
import GPUtil
import psutil
import time

from app.common.config import settings

# 자원 사용량 제한 설정
CPU_USAGE_LIMIT = 85  # CPU 사용량 제한
MEMORY_USAGE_LIMIT = 90  # 메모리 사용량 제한
GPU_USAGE_LIMIT = 90  # GPU 사용량 제한

def check_system_usage():
    """ CPU, 메모리 및 GPU 사용량 확인 후 Celery Task 큐 컨트롤 """
    while True:
        cpu_usage = psutil.cpu_percent(interval=1)
        memory_usage = psutil.virtual_memory().percent
        gpus = GPUtil.getGPUs()
        gpu_usage = max((gpu.load * 100 for gpu in gpus), default=0)

        celery_app = current_celery_app

        print(f"💡 시스템 리소스 체크: CPU {cpu_usage}%, 메모리 {memory_usage}%, GPU {gpu_usage}%")

        if cpu_usage >= CPU_USAGE_LIMIT or memory_usage >= MEMORY_USAGE_LIMIT or gpu_usage >= GPU_USAGE_LIMIT:
            print(f"⚠️ 리소스 초과! (CPU: {cpu_usage}%, Memory: {memory_usage}%, GPU: {gpu_usage}%) → 큐 대기 중...")
            # celery_app.control.pause_consumer("workflow_task")  # 특정 Task Queue 중지
            celery_app.control.cancel_consumer("workflow_task")  # 특정 Task Queue 중지
            celery_app.control.cancel_consumer("plugin_task")  # 플러그인 Task Queue 중지
        else:
            print(f"✅ 정상 상태 (CPU: {cpu_usage}%, Memory: {memory_usage}%, GPU: {gpu_usage}%) → 큐 실행 중...")
            # celery_app.control.resume_consumer("workflow_task")  # 특정 Task Queue 다시 실행
            celery_app.control.add_consumer("workflow_task")  # 특정 Task Queue 다시 실행
            celery_app.control.add_consumer("plugin_task")  # 플러그인 Task Queue 다시 실행

        time.sleep(5)  # 5초마다 상태 체크

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

    # 긴 작업을 위한 타임아웃 설정 추가
    celery_app.conf.update(broker_transport_options={
        'visibility_timeout': 400000,  # 긴 작업을 위한 visibility timeout 설정 (초 단위)
        'confirm_publish': True,  # 메시지 발행 확인 활성화
        'confirm_timeout': 60.0  # 메시지 발행 확인을 위한 타임아웃 설정 (초 단위)
    })

    # Task Time Limits 설정 (24시간 제한, 소프트 타임아웃 23시간 56분)
    celery_app.conf.update(task_time_limit=86400)  # 작업 시간 제한
    celery_app.conf.update(task_soft_time_limit=86200)  # 소프트 타임아웃 설정

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

    # 라우팅 설정 활성화 (큐별 라우팅을 위해 필요)
    # celery_app.conf.task_routes = None  # 기본 설정으로 돌아감

    # 작업자(worker) 동시성 제한 설정 수정
    celery_app.conf.update(
        task_reject_on_worker_lost=True  # 워커 손실 시 작업 거부
    )

    # monitoring_thread = threading.Thread(target=check_system_usage, daemon=True)
    # monitoring_thread.start()

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