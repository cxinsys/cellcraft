import docker
import threading
import signal
import os
import fnmatch
from typing import Dict, Optional, Set
from datetime import datetime

class ContainerManager:
    """작업별 Docker 컨테이너 추적 및 관리"""
    
    def __init__(self):
        self._task_containers: Dict[str, str] = {}  # task_id -> container_id
        self._container_tasks: Dict[str, str] = {}  # container_id -> task_id
        self._cleanup_in_progress: Set[str] = set()  # 정리 중인 컨테이너 ID
        self._lock = threading.Lock()
        self._docker_client = None
        self._setup_signal_handlers()
    
    @property
    def docker_client(self):
        """Docker 클라이언트 lazy initialization"""
        if self._docker_client is None:
            self._docker_client = docker.from_env()
        return self._docker_client
    
    def register_container(self, task_id: str, container_id: str):
        """작업 ID와 컨테이너 ID 매핑 등록"""
        with self._lock:
            self._task_containers[task_id] = container_id
            self._container_tasks[container_id] = task_id
            print(f"Container registered: Task {task_id} -> Container {container_id}")
    
    def unregister_container(self, task_id: str):
        """작업 완료 시 매핑 해제"""
        with self._lock:
            if task_id in self._task_containers:
                container_id = self._task_containers.pop(task_id)
                self._container_tasks.pop(container_id, None)
                # 정리 중 상태도 제거
                self._cleanup_in_progress.discard(container_id)
                print(f"Container unregistered: Task {task_id}")
    
    def get_container_id(self, task_id: str) -> Optional[str]:
        """작업 ID로 컨테이너 ID 조회"""
        with self._lock:
            return self._task_containers.get(task_id)
    
    def get_task_id(self, container_id: str) -> Optional[str]:
        """컨테이너 ID로 작업 ID 조회"""
        with self._lock:
            return self._container_tasks.get(container_id)
    
    def _is_cleanup_in_progress(self, container_id: str) -> bool:
        """컨테이너 정리가 진행 중인지 확인"""
        with self._lock:
            return container_id in self._cleanup_in_progress
    
    def _mark_cleanup_in_progress(self, container_id: str) -> bool:
        """컨테이너 정리 진행 상태로 마킹. 이미 진행 중이면 False 반환"""
        with self._lock:
            if container_id in self._cleanup_in_progress:
                return False
            self._cleanup_in_progress.add(container_id)
            return True
    
    def _unmark_cleanup_in_progress(self, container_id: str):
        """컨테이너 정리 완료 상태로 변경"""
        with self._lock:
            self._cleanup_in_progress.discard(container_id)
    
    def stop_task_container(self, task_id: str, timeout: int = 10) -> bool:
        """특정 작업의 컨테이너 강제 중단"""
        container_id = self.get_container_id(task_id)
        if not container_id:
            print(f"No container found for task {task_id}")
            # 매핑에 없어도 라벨을 통해 검색 시도
            return self._stop_container_by_task_label(task_id, timeout)
        
        # Race condition 방지: 이미 정리 중인지 확인
        if not self._mark_cleanup_in_progress(container_id):
            print(f"Container {container_id} cleanup already in progress")
            return True  # 다른 곳에서 정리 중이므로 성공으로 간주
        
        try:
            container = self.docker_client.containers.get(container_id)
            print(f"Stopping container {container_id} for task {task_id}")
            
            # 컨테이너 강제 중단
            container.kill(signal='SIGTERM')
            container.wait(timeout=timeout)
            container.remove(force=True)
            
            print(f"Container {container_id} stopped and removed successfully")
            self.unregister_container(task_id)
            return True
            
        except docker.errors.NotFound:
            print(f"Container {container_id} not found (already stopped)")
            self.unregister_container(task_id)
            return True
        except docker.errors.APIError as e:
            if "already in progress" in str(e).lower() or "409" in str(e):
                print(f"Container {container_id} removal already in progress")
                self.unregister_container(task_id)
                return True  # 이미 정리 중이므로 성공으로 간주
            else:
                print(f"API error stopping container {container_id}: {e}")
                return False
        except Exception as e:
            print(f"Error stopping container {container_id}: {e}")
            return False
        finally:
            self._unmark_cleanup_in_progress(container_id)
    
    def _stop_container_by_task_label(self, task_id: str, timeout: int = 10) -> bool:
        """작업 ID 라벨을 통한 컨테이너 검색 및 중단"""
        try:
            # 1. 정확한 task_id 라벨로 검색
            containers = self.docker_client.containers.list(
                filters={"label": f"celery.task_id={task_id}"}
            )
            
            # 2. task_id가 unknown인 경우를 위한 추가 검색
            if not containers:
                # 컨테이너 환경변수나 이름으로 검색
                all_containers = self.docker_client.containers.list(
                    filters={"label": "container.type=plugin-execution"}
                )
                
                for container in all_containers:
                    # 환경변수 확인
                    env_vars = container.attrs.get('Config', {}).get('Env', [])
                    for env in env_vars:
                        if env.startswith('CELERY_TASK_ID=') and task_id in env:
                            containers.append(container)
                            break
                    
                    # 컨테이너 이름 확인
                    if task_id[:8] in container.name:
                        containers.append(container)
            
            stopped_any = False
            for container in containers:
                print(f"Found container {container.id} ({container.name}) for task {task_id} via label/search")
                try:
                    container.kill(signal='SIGTERM')
                    container.wait(timeout=timeout)
                    container.remove(force=True)
                    print(f"Container {container.id} stopped successfully via label")
                    stopped_any = True
                except Exception as e:
                    print(f"Error stopping container {container.id} via label: {e}")
            
            return stopped_any
            
        except Exception as e:
            print(f"Error searching containers by task label {task_id}: {e}")
            return False
    
    def stop_container_by_name(self, container_name_pattern: str, timeout: int = 10) -> bool:
        """컨테이너 이름 패턴으로 컨테이너 강제 중단"""
        try:
            # 모든 컨테이너 목록 조회
            containers = self.docker_client.containers.list()
            
            matched_containers = []
            for container in containers:
                # fnmatch를 사용한 패턴 매칭
                if fnmatch.fnmatch(container.name, container_name_pattern):
                    matched_containers.append(container)
            
            if not matched_containers:
                print(f"No containers found matching pattern: {container_name_pattern}")
                return False
            
            for container in matched_containers:
                print(f"Stopping container by name pattern: {container.name}")
                try:
                    container.kill(signal='SIGTERM')
                    container.wait(timeout=timeout)
                    container.remove(force=True)
                    
                    # 매핑에서도 제거
                    task_id = self.get_task_id(container.id)
                    if task_id:
                        self.unregister_container(task_id)
                    
                    print(f"Container {container.id} ({container.name}) stopped successfully")
                except Exception as e:
                    print(f"Error stopping container {container.name}: {e}")
            
            return len(matched_containers) > 0
            
        except Exception as e:
            print(f"Error stopping container by name pattern {container_name_pattern}: {e}")
            return False
    
    def cleanup_plugin_containers(self) -> int:
        """모든 플러그인 실행 컨테이너 정리"""
        try:
            containers = self.docker_client.containers.list(
                filters={"label": "container.type=plugin-execution"}
            )
            
            cleaned_count = 0
            for container in containers:
                try:
                    print(f"Cleaning up plugin container: {container.name}")
                    container.kill(signal='SIGTERM')
                    container.wait(timeout=5)
                    container.remove(force=True)
                    cleaned_count += 1
                    
                    # 매핑에서도 제거
                    task_id = self.get_task_id(container.id)
                    if task_id:
                        self.unregister_container(task_id)
                        
                except Exception as e:
                    print(f"Error cleaning up container {container.name}: {e}")
            
            print(f"Cleaned up {cleaned_count} plugin containers")
            return cleaned_count
            
        except Exception as e:
            print(f"Error during plugin container cleanup: {e}")
            return 0
    
    def cleanup_all_task_containers(self):
        """모든 등록된 작업 컨테이너 및 플러그인 컨테이너 정리"""
        with self._lock:
            task_ids = list(self._task_containers.keys())
        
        print(f"Cleaning up {len(task_ids)} tracked task containers...")
        for task_id in task_ids:
            self.stop_task_container(task_id)
        
        # 추가로 라벨 기반 정리
        plugin_count = self.cleanup_plugin_containers()
        print(f"Total cleanup completed. Plugin containers cleaned: {plugin_count}")
    
    def _setup_signal_handlers(self):
        """시그널 핸들러 설정"""
        def signal_handler(signum, frame):
            print(f"ContainerManager received signal {signum}. Cleaning up containers...")
            self.cleanup_all_task_containers()
            # 원래 시그널 핸들러 호출
            if signum == signal.SIGTERM:
                os._exit(1)
        
        # SIGTERM과 SIGINT 핸들러 등록
        try:
            signal.signal(signal.SIGTERM, signal_handler)
            signal.signal(signal.SIGINT, signal_handler)
        except ValueError:
            # 메인 스레드가 아닌 경우 시그널 핸들러 등록 불가
            print("Warning: Could not register signal handlers (not in main thread)")
    
    def get_container_name_for_task(self, task_id: str, plugin_name: str) -> str:
        """작업 ID와 플러그인 이름으로 컨테이너 이름 생성"""
        timestamp = int(datetime.now().timestamp())
        return f"plugin-{plugin_name.lower()}-task-{task_id[:8]}-{timestamp}"
    
    def get_status(self) -> dict:
        """컨테이너 매니저 상태 조회"""
        with self._lock:
            return {
                "tracked_tasks": len(self._task_containers),
                "task_container_mapping": dict(self._task_containers),
                "container_task_mapping": dict(self._container_tasks)
            }

# 전역 컨테이너 매니저 인스턴스
container_manager = ContainerManager() 