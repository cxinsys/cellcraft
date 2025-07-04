from fastapi import APIRouter, Depends, HTTPException, Query
from typing import List,Optional,Any
from sqlalchemy.orm import Session
import docker
from sqlalchemy import asc, desc, or_

from app.routes import dep
from app.database.schemas.admin import Conditions
from app.database.crud import crud_admin
from app.database import models

router = APIRouter()

@router.get("/users", response_model=Any)
def get_filtered_users(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
):
    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    users, total_count = crud_admin.get_filtered_users(db, conditions)

    if not users:
        raise HTTPException(status_code=404, detail="Users not found")
    return users

@router.get("/users_count", response_model=Any)
def get_users_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    users_num = crud_admin.get_users_count(db)
    return users_num

@router.get("/files", response_model=Any)
def get_filtered_files(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):
    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    files, total_count = crud_admin.get_filtered_files(db, conditions)

    if not files:
        raise HTTPException(status_code=404, detail="Files not found")
    
    # 결과 포맷팅
    formatted_files = []
    for file, username in files:
        formatted_file = {
            'id': file.id,
            'file_name': file.file_name,
            'file_path': file.file_path,
            'file_size': file.file_size,
            'folder': file.folder,
            'username': username,
            'user_id': file.user_id
        }
        formatted_files.append(formatted_file)
    
    return {
        'data': formatted_files,
        'total_count': total_count
    }

@router.get("/files_count", response_model=Any)
def get_files_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    files_num = crud_admin.get_files_count(db)
    return files_num

@router.get("/workflows", response_model=Any)
def get_filtered_workflows(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
):
    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    workflows, total_count = crud_admin.get_filtered_workflows(db, conditions)

    if not workflows:
        raise HTTPException(status_code=404, detail="Workflows not found")
    
    # 결과 포맷팅
    formatted_workflows = []
    for workflow, username in workflows:
        formatted_workflow = {
            'id': workflow.id,
            'title': workflow.title,
            'username': username,
            'updated_at': workflow.updated_at,
            'user_id': workflow.user_id
        }
        formatted_workflows.append(formatted_workflow)
    
    return {
        'data': formatted_workflows,
        'total_count': total_count
    }

@router.get("/workflows_count", response_model=Any)
def get_workflows_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):

    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    workflows_num = crud_admin.get_workflows_count(db)
    return workflows_num

@router.get("/tasks", response_model=Any)
def get_filtered_tasks(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):

    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    tasks, total_count = crud_admin.get_filtered_tasks(db, conditions)

    if not tasks:
        raise HTTPException(status_code=404, detail="Tasks not found")
    
    # 결과 포맷팅
    formatted_tasks = []
    for task, username, workflow_title in tasks:
        formatted_task = {
            'id': task.id,
            'user_id': task.user_id,
            'workflow_id': task.workflow_id,
            'username': username,
            'workflow_title': workflow_title,
            'status': task.status,
            'start_time': task.start_time
        }
        formatted_tasks.append(formatted_task)
    
    return {
        'data': formatted_tasks,
        'total_count': total_count
    }

@router.get("/tasks_count", response_model=Any)
def get_tasks_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    tasks_num = crud_admin.get_tasks_count(db)
    return tasks_num

@router.get("/plugins", response_model=Any)
def get_filtered_plugins(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):
    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    plugins, total_count = crud_admin.get_filtered_plugins(db, conditions)

    if not plugins:
        raise HTTPException(status_code=404, detail="Plugins not found")
    return plugins

@router.get("/plugins_count", response_model=Any)
def get_plugins_count(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    # 관리자 권한 확인
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")

    plugins_num = crud_admin.get_plugins_count(db)
    return plugins_num

@router.get("/system/stats", response_model=Any)
def get_system_stats():
    # Docker 클라이언트 연결
    client = docker.DockerClient(base_url='unix://var/run/docker.sock')

    try:
        containers = client.containers.list()
        print("현재 실행 중인 컨테이너 목록:")
        
        # 모든 컨테이너의 통계 정보를 저장할 리스트
        container_stats = []
        
        for container in containers:
            container_name = container.name
            print(container_name)
            
            # 컨테이너 필터링 로직
            if "cellcraft-rabbitmq" in container_name:
                continue  # rabbitmq 컨테이너는 건너뛰기
            
            if "cellcraft" in container_name and "celery" not in container_name:
                continue  # cellcraft가 포함되어 있지만 celery가 아닌 컨테이너는 건너뛰기
            
            try:
                # 컨테이너 상세 정보 가져오기
                stats = container.stats(stream=False)
                container_info = container.attrs

                # CPU 상세 정보 계산
                cpu_stats = stats['cpu_stats']
                precpu_stats = stats['precpu_stats']
                cpu_usage = cpu_stats['cpu_usage']['total_usage'] - precpu_stats['cpu_usage']['total_usage']
                system_cpu_usage = cpu_stats['system_cpu_usage'] - precpu_stats['system_cpu_usage']
                num_cpus = len(cpu_stats['cpu_usage'].get('percpu_usage', []))
                cpu_percent = (cpu_usage / system_cpu_usage) * 100 * num_cpus if system_cpu_usage > 0 else 0

                # 메모리 상세 정보 계산
                memory_stats = stats['memory_stats']
                memory_usage = memory_stats.get('usage', 0)
                memory_limit = memory_stats.get('limit', 0)
                memory_percent = (memory_usage / memory_limit) * 100 if memory_limit > 0 else 0
                memory_stats_detailed = {
                    'total_bytes': memory_limit,
                    'used_bytes': memory_usage,
                    'available_bytes': memory_limit - memory_usage,
                    'percent': memory_percent,
                    'stats': memory_stats.get('stats', {})
                }

                # GPU 정보 가져오기
                gpu_stats = None
                try:
                    exec_result = container.exec_run(
                        'nvidia-smi --query-gpu=index,name,temperature.gpu,utilization.gpu,memory.total,memory.used,memory.free,power.draw,power.limit --format=csv,noheader,nounits'
                    )
                    if exec_result.exit_code == 0:
                        gpu_stats = []
                        gpu_output = exec_result.output.decode('utf-8')
                        for line in gpu_output.splitlines():
                            index, name, temp, util, mem_total, mem_used, mem_free, power_draw, power_limit = line.split(',')
                            gpu_stats.append({
                                'id': int(index.strip()),
                                'name': name.strip(),
                                'temperature_c': float(temp.strip()),
                                'utilization_percent': float(util.strip()),
                                'memory': {
                                    'total_bytes': int(mem_total.strip()) * 1024 * 1024,
                                    'used_bytes': int(mem_used.strip()) * 1024 * 1024,
                                    'free_bytes': int(mem_free.strip()) * 1024 * 1024,
                                    'utilization_percent': (int(mem_used.strip()) / int(mem_total.strip())) * 100
                                },
                                'power': {
                                    'draw_watts': float(power_draw.strip()),
                                    'limit_watts': float(power_limit.strip())
                                }
                            })
                except Exception as e:
                    print(f"GPU 정보 조회 실패: {e}")

                # 네트워크 정보
                network_stats = stats.get('networks', {})
                
                # 컨테이너별 통계 정보 저장
                container_stats.append({
                    'container_info': {
                        'id': container.id,
                        'name': container.name,
                        'status': container.status,
                        'created': container_info['Created'],
                        'image': container_info['Config']['Image'],
                        'command': container_info['Config']['Cmd'],
                        'labels': container_info['Config'].get('Labels', {})
                    },
                    'cpu': {
                        'usage_percent': cpu_percent,
                        'num_cpus': num_cpus,
                        'total_usage': cpu_usage,
                        'system_usage': system_cpu_usage,
                        'per_cpu_usage': cpu_stats['cpu_usage'].get('percpu_usage', [])
                    },
                    'memory': memory_stats_detailed,
                    'gpu': gpu_stats,
                    'network': {
                        interface: {
                            'rx_bytes': data['rx_bytes'],
                            'tx_bytes': data['tx_bytes'],
                            'rx_packets': data['rx_packets'],
                            'tx_packets': data['tx_packets'],
                            'rx_errors': data['rx_errors'],
                            'tx_errors': data['tx_errors']
                        } for interface, data in network_stats.items()
                    }
                })
            except Exception as e:
                print(f"컨테이너 {container.name} 통계 정보 조회 실패: {e}")
                continue

        # 전체 시스템 통계 계산
        total_cpu_usage = sum(stat['cpu']['usage_percent'] for stat in container_stats)
        total_memory_usage = sum(stat['memory']['used_bytes'] for stat in container_stats)
        total_memory_limit = sum(stat['memory']['total_bytes'] for stat in container_stats)
        
        system_stats = {
            'total_containers': len(container_stats),
            'total_cpu_usage_percent': total_cpu_usage,
            'total_memory_usage_bytes': total_memory_usage,
            'total_memory_limit_bytes': total_memory_limit,
            'total_memory_usage_percent': (total_memory_usage / total_memory_limit * 100) if total_memory_limit > 0 else 0,
            'containers': container_stats
        }
        
        return system_stats

    except docker.errors.NotFound:
        return {'message': "Docker 데몬에 연결할 수 없습니다."}
    except Exception as e:
        return {'message': f"시스템 정보 조회 중 오류 발생: {str(e)}"}

@router.put("/users/{user_id}", response_model=Any)
def update_user(
    user_id: int,
    user_data: dict,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")
    
    user = crud_admin.update_user(db, user_id, user_data)
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    return user

@router.delete("/users/{user_id}")
def delete_user(
    user_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")
    
    success = crud_admin.delete_user(db, user_id)
    if not success:
        raise HTTPException(status_code=404, detail="User not found")
    return {"message": "User deleted successfully"}

@router.delete("/files/{file_id}")
def delete_file(
    file_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")
    
    success = crud_admin.delete_file(db, file_id)
    if not success:
        raise HTTPException(status_code=404, detail="File not found")
    return {"message": "File deleted successfully"}

@router.delete("/workflows/{workflow_id}")
def delete_workflow(
    workflow_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")
    
    success = crud_admin.delete_workflow(db, workflow_id)
    if not success:
        raise HTTPException(status_code=404, detail="Workflow not found")
    return {"message": "Workflow deleted successfully"}

@router.post("/tasks/{task_id}/cancel")
def cancel_task(
    task_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")
    
    success = crud_admin.cancel_task(db, task_id)
    if not success:
        raise HTTPException(status_code=404, detail="Task not found")
    return {"message": "Task cancelled successfully"}

@router.post("/plugins/{plugin_id}/install-dependencies")
def install_plugin_dependencies(
    plugin_id: int,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
):
    if not current_user.is_superuser:
        raise HTTPException(status_code=403, detail="Access denied: Admins only")
    
    success = crud_admin.install_plugin_dependencies(db, plugin_id)
    if not success:
        raise HTTPException(status_code=404, detail="Plugin not found")
    return {"message": "Plugin dependencies installed successfully"}