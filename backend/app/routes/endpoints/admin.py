from fastapi import APIRouter, Depends, HTTPException
from typing import List,Optional,Any
from sqlalchemy.orm import Session
import psutil
import GPUtil
import docker
import subprocess
import json

from app.routes import dep
from app.database.schemas.user import User
from app.database.schemas.admin import Conditions
from app.database.crud import crud_user, crud_admin

from app.database import models

router = APIRouter()

@router.get("/users", response_model=Any)
def get_filtered_users(
    *,
    db: Session = Depends(dep.get_db),
    # conditions: Conditions, # 이거 안됨
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    users = crud_user.get_filtered_users(db, conditions)

    if users is None:
        raise HTTPException(status_code=404, detail="Users not found")
    return users

@router.get("/users_count", response_model=Any)
def get_users_count(
    *,
    db: Session = Depends(dep.get_db),
):
    users_num = crud_user.get_users_count(db)
    return users_num

@router.get("/files", response_model=Any)
def get_filtered_files(
    *,
    db: Session = Depends(dep.get_db),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    files = crud_admin.get_filtered_files(db, conditions)

    if files is None:
        raise HTTPException(status_code=404, detail="Files not found")
    return files

@router.get("/files_count", response_model=Any)
def get_files_count(
    *,
    db: Session = Depends(dep.get_db),
):
    files_num = crud_admin.get_files_count(db)
    return files_num

@router.get("/workflows", response_model=Any)
def get_filtered_workflows(
    *,
    db: Session = Depends(dep.get_db),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    workflows = crud_admin.get_filtered_workflows(db, conditions)

    if workflows is None:
        raise HTTPException(status_code=404, detail="Workflows not found")
    return workflows

@router.get("/workflows_count", response_model=Any)
def get_workflows_count(
    *,
    db: Session = Depends(dep.get_db),
):
    workflows_num = crud_admin.get_workflows_count(db)
    return workflows_num

@router.get("/tasks", response_model=Any)
def get_filtered_tasks(
    *,
    db: Session = Depends(dep.get_db),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    tasks = crud_admin.get_filtered_tasks(db, conditions)

    if tasks is None:
        raise HTTPException(status_code=404, detail="Tasks not found")
    return tasks

@router.get("/tasks_count", response_model=Any)
def get_tasks_count(
    *,
    db: Session = Depends(dep.get_db),
):
    tasks_num = crud_admin.get_tasks_count(db)
    return tasks_num

@router.get("/plugins", response_model=Any)
def get_filtered_plugins(
    *,
    db: Session = Depends(dep.get_db),
    amount: int,
    page_num: int,
    sort: str,
    order: str,
    searchTerm: str,
    ):

    conditions = Conditions(
        amount=amount,
        page_num=page_num,
        sort=sort,
        order=order,
        searchTerm=searchTerm
    )

    plugins = crud_admin.get_filtered_plugins(db, conditions)

    if plugins is None:
        raise HTTPException(status_code=404, detail="Plugins not found")
    return plugins

@router.get("/plugins_count", response_model=Any)
def get_plugins_count(
    *,
    db: Session = Depends(dep.get_db),
):
    plugins_num = crud_admin.get_plugins_count(db)
    return plugins_num

@router.get("/system/stats", response_model=Any)
def get_celery_stats():
    # Docker 클라이언트 연결
    client = docker.DockerClient(base_url='unix://var/run/docker.sock')

    try:
        containers = client.containers.list()
        print("현재 실행 중인 컨테이너 목록:")
        
        celery_container = None
        for container in containers:
            print(container.name)
            if 'celery' in container.name:
                celery_container = container
                break

        if not celery_container:
            return {'message': "Celery 컨테이너를 찾을 수 없습니다."}
        
        # 컨테이너 상세 정보 가져오기
        stats = celery_container.stats(stream=False)
        container_info = celery_container.attrs

        # CPU 상세 정보 계산
        cpu_stats = stats['cpu_stats']
        precpu_stats = stats['precpu_stats']
        cpu_usage = cpu_stats['cpu_usage']['total_usage'] - precpu_stats['cpu_usage']['total_usage']
        system_cpu_usage = cpu_stats['system_cpu_usage'] - precpu_stats['system_cpu_usage']
        num_cpus = len(cpu_stats['cpu_usage'].get('percpu_usage', []))
        cpu_percent = (cpu_usage / system_cpu_usage) * 100 * num_cpus

        # 메모리 상세 정보 계산
        memory_stats = stats['memory_stats']
        memory_usage = memory_stats['usage']
        memory_limit = memory_stats['limit']
        memory_percent = (memory_usage / memory_limit) * 100
        memory_stats_detailed = {
            'total_bytes': memory_limit,
            'used_bytes': memory_usage,
            'available_bytes': memory_limit - memory_usage,
            'percent': memory_percent,
            'stats': memory_stats.get('stats', {})
        }

        # GPU 정보 가져오기
        try:
            gpu_stats = []
            exec_result = celery_container.exec_run(
                'nvidia-smi --query-gpu=index,name,temperature.gpu,utilization.gpu,memory.total,memory.used,memory.free,power.draw,power.limit --format=csv,noheader,nounits'
            )
            if exec_result.exit_code == 0:
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
            gpu_stats = None
            print(f"GPU 정보 조회 실패: {e}")

        # 네트워크 정보
        network_stats = stats.get('networks', {})
        
        return {
            'container_info': {
                'id': celery_container.id,
                'name': celery_container.name,
                'status': celery_container.status,
                'created': container_info['Created'],
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
        }

    except docker.errors.NotFound:
        return {'message': "Celery 컨테이너를 찾을 수 없습니다."}
    except Exception as e:
        return {'message': f"시스템 정보 조회 중 오류 발생: {str(e)}"}

# @router.get("/system/stats", response_model=Any)
# def get_celery_stats():
#     # Docker 클라이언트 연결
#     client = docker.DockerClient(base_url='unix://var/run/docker.sock')

#     try:
#         containers = client.containers.list()
#         print("현재 실행 중인 컨테이너 목록:")

#         celery_container = None
#         for container in containers:
#             print(container.name)
#             if 'celery' in container.name:
#                 celery_container = container
#                 break

#         if not celery_container:
#             return {'message': "Celery 컨테이너를 찾을 수 없습니다."}

#         # 컨테이너 상세 정보 가져오기
#         stats = celery_container.stats(stream=False)
#         container_info = celery_container.attrs

#         # CPU 상세 정보 계산
#         cpu_stats = stats['cpu_stats']
#         precpu_stats = stats['precpu_stats']
#         cpu_usage = cpu_stats['cpu_usage']['total_usage'] - precpu_stats['cpu_usage']['total_usage']
#         system_cpu_usage = cpu_stats['system_cpu_usage'] - precpu_stats['system_cpu_usage']
#         num_cpus = len(cpu_stats['cpu_usage'].get('percpu_usage', []))
#         cpu_percent = (cpu_usage / system_cpu_usage) * 100 * num_cpus

#         # 메모리 상세 정보 계산
#         memory_stats = stats['memory_stats']
#         memory_usage = memory_stats['usage']
#         memory_limit = memory_stats['limit']
#         memory_percent = (memory_usage / memory_limit) * 100
#         memory_stats_detailed = {
#             'total_bytes': memory_limit,
#             'used_bytes': memory_usage,
#             'available_bytes': memory_limit - memory_usage,
#             'percent': memory_percent,
#             'stats': memory_stats.get('stats', {})
#         }

#         # **Celery 컨테이너 내부에서 GPU 정보 가져오기 (`pynvml` 실행)**
#         try:
#             python_script = """
# import json
# from pynvml import nvmlInit, nvmlDeviceGetCount, nvmlDeviceGetHandleByIndex, \
#     nvmlDeviceGetName, nvmlDeviceGetTemperature, nvmlDeviceGetMemoryInfo, \
#     nvmlDeviceGetPowerUsage, nvmlDeviceGetUtilizationRates, nvmlShutdown

# nvmlInit()
# gpu_stats = []
# device_count = nvmlDeviceGetCount()

# for i in range(device_count):
#     handle = nvmlDeviceGetHandleByIndex(i)
#     mem_info = nvmlDeviceGetMemoryInfo(handle)
#     utilization = nvmlDeviceGetUtilizationRates(handle)

#     gpu_stats.append({
#         'id': i,
#         'name': nvmlDeviceGetName(handle),  # decode() 제거
#         'temperature_c': nvmlDeviceGetTemperature(handle, 0),
#         'utilization_percent': utilization.gpu,
#         'memory': {
#             'total_bytes': mem_info.total,
#             'used_bytes': mem_info.used,
#             'free_bytes': mem_info.free,
#             'utilization_percent': (mem_info.used / mem_info.total) * 100 if mem_info.total > 0 else 0
#         },
#         'power': {
#             'draw_watts': nvmlDeviceGetPowerUsage(handle) / 1000
#         }
#     })

# nvmlShutdown()
# print(json.dumps(gpu_stats))
# """

#             exec_result = celery_container.exec_run(f'/opt/conda/envs/snakemake/bin/python -c "{python_script}"')

#             if exec_result.exit_code == 0:
#                 gpu_stats = json.loads(exec_result.output.decode('utf-8'))
#             else:
#                 gpu_stats = None
#                 print(f"GPU 정보 조회 실패: {exec_result.output.decode('utf-8')}")

#         except Exception as e:
#             gpu_stats = None
#             print(f"GPU 정보 조회 실패: {e}")

#         # 네트워크 정보
#         network_stats = stats.get('networks', {})

#         return {
#             'container_info': {
#                 'id': celery_container.id,
#                 'name': celery_container.name,
#                 'status': celery_container.status,
#                 'created': container_info['Created'],
#             },
#             'cpu': {
#                 'usage_percent': cpu_percent,
#                 'num_cpus': num_cpus,
#                 'total_usage': cpu_usage,
#                 'system_usage': system_cpu_usage,
#                 'per_cpu_usage': cpu_stats['cpu_usage'].get('percpu_usage', [])
#             },
#             'memory': memory_stats_detailed,
#             'gpu': gpu_stats,
#             'network': {
#                 interface: {
#                     'rx_bytes': data['rx_bytes'],
#                     'tx_bytes': data['tx_bytes'],
#                     'rx_packets': data['rx_packets'],
#                     'tx_packets': data['tx_packets'],
#                     'rx_errors': data['rx_errors'],
#                     'tx_errors': data['tx_errors']
#                 } for interface, data in network_stats.items()
#             }
#         }

#     except docker.errors.NotFound:
#         return {'message': "Celery 컨테이너를 찾을 수 없습니다."}
#     except Exception as e:
#         return {'message': f"시스템 정보 조회 중 오류 발생: {str(e)}"}