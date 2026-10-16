"""
Task status SSE (Server-Sent Events) streaming.

Extracted from ``task/router.py`` in PR-7 (Phase 3c). The ``event_generator``
that previously lived inline inside ``get_task_status`` now lives here as the
reusable async streaming utility ``stream_task_status``. The router endpoint is
reduced to ``EventSourceResponse(sse.stream_task_status(task_id))``.

Behavior is preserved exactly: a 1-hour (3600s) timeout guard, a 5-second poll
interval, termination on the terminal Celery states (SUCCESS / FAILURE /
REVOKED), cancellation handling, and the ``finally`` cleanup log.

NOTE on test patching: ``get_task_info`` is bound into this module's namespace
at import time, so tests exercising the stream must patch
``app.task.sse.get_task_info`` (not ``app.worker.utils.get_task_info``);
otherwise the loop polls the real Celery backend.
"""
import asyncio
import time

from app.worker.utils import get_task_info


async def stream_task_status(task_id: str):
    """Yield task status frames for an SSE connection until a terminal state.

    Emits the current Celery ``task_status`` every 5 seconds, terminates on
    SUCCESS / FAILURE / REVOKED, and emits ``TIMEOUT`` after 3600s to prevent
    memory leaks from abandoned connections. Cancellation and unexpected errors
    are swallowed (logged via ``print``) and the ``finally`` block logs cleanup —
    identical to the original inline generator.
    """
    timeout = 3600  # 1시간 타임아웃 (메모리 누수 방지)
    start_time = time.time()
    try:
        while True:
            # 타임아웃 체크
            if time.time() - start_time > timeout:
                print(f"SSE timeout reached for task {task_id}")
                yield f"TIMEOUT"
                break

            if task_id:
                task = get_task_info(task_id)
                task_status = task.get('task_status', 'UNKNOWN')

                if task_status in ['SUCCESS', 'FAILURE', 'REVOKED']:
                    yield f"{task_status}"
                    break

                print(f"Task {task_id} status: {task_status}")
                yield f"{task_status}"
                await asyncio.sleep(5)
            else:
                break
    except asyncio.CancelledError:
        # 클라이언트 연결 끊김 처리
        print(f"SSE connection cancelled for task {task_id}")
    except Exception as e:
        print(f"SSE error for task {task_id}: {e}")
    finally:
        # 정리 로직 - 리소스 해제
        print(f"SSE generator cleanup for task {task_id}")
