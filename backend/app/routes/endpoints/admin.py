from fastapi import APIRouter, Depends, HTTPException
from typing import List,Optional,Any
from sqlalchemy.orm import Session
import psutil
import GPUtil
import docker
import subprocess

from app.routes import dep  # get_db를 정의한 python 파일을 import
from app.database.schemas.user import User  # User를 정의한 python 파일을 import
from app.database.schemas.admin import Conditions  # Conditions를 정의한 python 파일을 import
from app.database.crud import crud_user  # get_filtered_users를 정의한 python 파일을 import

from app.database import models

router = APIRouter()

@router.get("/users", response_model=List[User])
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

# def get_gpu_stats(container_name):
#     try:
#         # 'nvidia-smi' 명령어를 Docker 컨테이너 내에서 실행하여 GPU 정보를 얻음
#         gpu_stats_output = subprocess.check_output(
#             ['docker', 'exec', container_name, 'nvidia-smi', '--query-gpu=index,name,utilization.gpu,memory.total,memory.used,memory.free', '--format=csv,noheader,nounits']
#         ).decode('utf-8')

#         # 여러 GPU가 있을 경우를 대비하여 리스트로 처리
#         gpu_stats = []
#         for line in gpu_stats_output.splitlines():
#             gpu_index, gpu_name, gpu_utilization, gpu_mem_total, gpu_mem_used, gpu_mem_free = line.split(',')
#             gpu_stats.append({
#                 'id': int(gpu_index.strip()),
#                 'name': gpu_name.strip(),
#                 'load_percent': float(gpu_utilization.strip()),  # GPU 사용량 (%)
#                 'total_memory_bytes': int(gpu_mem_total.strip()) * 1024 * 1024,  # MB를 바이트로 변환
#                 'used_memory_bytes': int(gpu_mem_used.strip()) * 1024 * 1024,
#                 'available_memory_bytes': int(gpu_mem_free.strip()) * 1024 * 1024
#             })

#         # 만약 GPU가 한 개라면 첫 번째 GPU만 반환
#         return gpu_stats[0] if len(gpu_stats) > 0 else None

#     except subprocess.CalledProcessError as e:
#         print(f"GPU stats retrieval failed: {e}")
#         return None

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

# def get_system_stats():
#     # CPU 사용량
#     cpu_usage = psutil.cpu_percent(interval=1)
    
#     # 메모리 사용량
#     memory_info = psutil.virtual_memory()
#     memory_usage = memory_info.percent
#     total_memory = memory_info.total
#     available_memory = memory_info.available
#     used_memory = memory_info.used

#     # GPU 사용량 (GPUtil을 이용하여 GPU 정보 조회)
#     gpus = GPUtil.getGPUs()
#     gpu_info = []
    
#     for gpu in gpus:
#         gpu_info.append({
#             "id": gpu.id,
#             "name": gpu.name,
#             "load_percent": gpu.load * 100,  # GPU 사용량 (0.0~1.0) -> 퍼센트로 변환
#             "memory_used_percent": gpu.memoryUtil * 100,  # GPU 메모리 사용량 (0.0~1.0) -> 퍼센트로 변환
#             "total_memory_bytes": gpu.memoryTotal * 1024 * 1024,  # 메모리 총량 (MB 단위)
#             "used_memory_bytes": gpu.memoryUsed * 1024 * 1024,  # 사용 중인 메모리 (MB 단위)
#             "available_memory_bytes": gpu.memoryFree * 1024 * 1024  # 사용 가능한 메모리 (MB 단위)
#         })

#     print(gpu_info)

#     # 반환할 데이터
#     system_stats = {
#         "cpu_usage_percent": cpu_usage,
#         "memory_usage_percent": memory_usage,
#         "total_memory_bytes": total_memory,
#         "available_memory_bytes": available_memory,
#         "used_memory_bytes": used_memory,
#         "gpu_stats": gpu_info
#     }

#     return system_stats
