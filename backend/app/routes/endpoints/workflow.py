from fastapi import APIRouter, Depends, HTTPException
from typing import Any
from sqlalchemy.orm import Session
import os
import shutil
import json
from fastapi.responses import FileResponse, JSONResponse

from app.common.utils.plugin_utils import verify_dependencies
from app.common.utils.celery_utils import get_task_info
from app.common.utils.snakemake_utils import change_snakefile_parameter
from app.common.utils.workflow_utils import extract_rule_block, extract_all_algorithms, extract_algorithm_data, extract_visualization_data, extract_target_data, generate_user_input, generate_plugin_params, generate_visualization_params
from app.database.crud import crud_workflow
from app.database.schemas.workflow import WorkflowDelete, WorkflowCreate, WorkflowUpdate, WorkflowResult, WorkflowFind, WorkflowNodeFileCreate, WorkflowNodeFileDelete, WorkflowNodeFileRead
from app.routes import dep
from app.database import models
from app.routes.celery_tasks import process_data_task

router = APIRouter()                                

#compile workflow
@router.post("/compile")
def compileWorkflow(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowCreate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    try:
        user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow.id)
        user_path = f"./user/{current_user.username}/"
        task_ids = []  # 여러 개의 태스크 ID를 저장할 리스트

        if user_workflow:
            crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.thumbnail, workflow.workflow_info)
            algorithms = extract_all_algorithms(user_workflow.workflow_info['drawflow']['Home']['data'])

            if not algorithms:
                raise HTTPException(status_code=400, detail="No algorithm nodes found in workflow.")

            for algorithm in algorithms:
                selected_plugin = algorithm['selectedPlugin']['name']
                user_workflow_task_path = f"{user_path}workflow_{workflow.id}/algorithm_{algorithm['id']}"

                # 플러그인 폴더 내의 dependency 폴더 경로 설정
                plugin_dependency_path = f"./plugin/{selected_plugin}/dependency"

                # 입력 데이터 및 파라미터 추출
                user_input = generate_user_input(algorithm['selectedPluginInputOutput'])
                plugin_params = generate_plugin_params(algorithm['selectedPluginRules'])
                target_list = extract_target_data(algorithm['selectedPluginInputOutput'], user_workflow_task_path)

                additional_data = {
                    "user_name": current_user.username,
                    "workflow_id": str(workflow.id),
                    "algorithm_id": str(algorithm['id']),
                    "plugin_name": selected_plugin,
                }
                user_input.update(additional_data)

                # 작업 폴더 생성
                if not os.path.exists(user_workflow_task_path):
                    os.makedirs(user_workflow_task_path)

                # Snakefile 생성
                plugin_snakefile_path = f"./plugin/{selected_plugin}/Snakefile"
                user_snakefile_path = change_snakefile_parameter(plugin_snakefile_path, user_workflow_task_path + "/Snakefile", user_input, plugin_params)

                # Celery 작업 실행
                process_task = process_data_task.apply_async(
                    args=[current_user.username, user_snakefile_path, selected_plugin, target_list],
                    kwargs={'user_id': current_user.id, 'workflow_id': workflow.id, 'algorithm_id': algorithm['id'], 'plugin_name': selected_plugin, 'task_type': 'compile'},
                    ignore_result=False
                )

                # 실행된 태스크 ID 저장
                task_ids.append(process_task.id)

            return {
                "message": "Multiple tasks added to queue",
                "task_ids": task_ids,
                "results": [get_task_info(task_id) for task_id in task_ids]
            }

    except Exception as e:
        raise HTTPException(
                status_code=400,
                detail=str(e),
        )

# visualize compile
@router.post("/visualization")
def visualizeData(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowUpdate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    try:
        user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow.id)
        user_path = f"./user/{current_user.username}/"
        if user_workflow:
            crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.thumbnail, workflow.workflow_info)
            extract_algorithm = extract_algorithm_data(user_workflow.workflow_info['drawflow']['Home']['data'], workflow.algorithm_id)
            extract_visualization = extract_visualization_data(user_workflow.workflow_info['drawflow']['Home']['data'], workflow.current_node_id)

            user_task_path = f"{user_path}workflow_{workflow.id}/algorithm_{extract_algorithm['id']}"

            selected_plugin = extract_algorithm['selectedPlugin']['name']
            # 플러그인 폴더 내의 dependency 폴더에 있는 파일 리스트를 순회하면서 검증
            plugin_dependency_path = "None"

            # print("extract_algorithm:", extract_algorithm)
            # print("extract_visualization:", extract_visualization)

            visualization_params, visualization_inputs, visualization_outputs = generate_visualization_params(extract_visualization['selectedVisualizationParams'])

            print("visualization_inputs:", visualization_inputs)
            print("visualization_outputs:", visualization_outputs)
            print("visualization_params:", visualization_params)

            user_workflow_visualization_result_directory = f"{user_task_path}/results"
            user_workflow_visualization_result_name = ""

            for key in visualization_inputs:
                user_workflow_visualization_result_name += f"{visualization_inputs[key]}_"

            for key in visualization_params:
                user_workflow_visualization_result_name += f"{visualization_params[key]}_"

            additional_data = {
                "user_name": current_user.username,
                "workflow_id": str(workflow.id),
                "algorithm_id": str(extract_algorithm['id']),
                "plugin_name": selected_plugin,
                "visualization_name": extract_visualization['selectedVisualizationTitle'],
                "visualization_result_path": user_workflow_visualization_result_name
            }

            # visualization_outputs의 첫 번째 값을 사용
            output_key = list(visualization_outputs.keys())[0]
            user_workflow_visualization_result_path = user_workflow_visualization_result_directory + "/" + user_workflow_visualization_result_name + f"{visualization_outputs[output_key]}"

            # user_workflow_visualization_result_path 있는 지 확인 후, 있으면 해당 파일 반환 없으면 코드 계속 실행
            if os.path.exists(user_workflow_visualization_result_path):
                print("user_workflow_visualization_result_path:", user_workflow_visualization_result_path)
                # message로 파일 있다고 반환
                message = "Visualization result already exists"
                return {
                    "message": message,
                    "result_path": user_workflow_visualization_result_name + f"{visualization_outputs[output_key]}"
                }

            else:
                visualization_snakefile_path = f"./plugin/{extract_algorithm['selectedPlugin']['name']}/visualization_Snakefile"
                
                # rule 추출 추가
                rule_content, rule_path = extract_rule_block(visualization_snakefile_path, extract_visualization['selectedVisualizationTitle'], user_task_path + "/visualization_Snakefile")
                print("rule_content:", rule_content)
                visualization_inputs.update(additional_data)

                user_visualization_snakefile_path = change_snakefile_parameter(rule_path, user_task_path + "/visualization_Snakefile", visualization_inputs, visualization_params)
                # user_workflow_visualization_result_path를 target_list에 추가
                target_list = [user_workflow_visualization_result_path]

                process_task = process_data_task.apply_async(
                    (current_user.username, user_visualization_snakefile_path, selected_plugin, target_list),
                    kwargs={'user_id': current_user.id, 'workflow_id': workflow.id, 'algorithm_id': extract_algorithm['id'], 'plugin_name': selected_plugin, 'task_type': 'visualization'}
                )

                message = "Tasks added to queue"
                task_id = process_task.id
                result = get_task_info(process_task.id)

                return {
                    "message": message,
                    "task_id": task_id,
                    "result": result,
                    "result_path": user_workflow_visualization_result_path
                }

    except Exception as e:
        raise HTTPException(
                status_code=400,
                detail=str(e),
                )

# get visualization result
@router.post("/visualize/result")
def getVisualizationResult(
    WorkflowResult: WorkflowResult, 
    current_user: models.User = Depends(dep.get_current_active_user)
):
    PATH_VISUALIZAE_RESULT = WorkflowResult.filename
    print("PATH_VISUALIZAE_RESULT:", PATH_VISUALIZAE_RESULT)

    try:
        with open(PATH_VISUALIZAE_RESULT, "r") as f:
            plotly_data = json.load(f)
        return JSONResponse(content=plotly_data)
        
    except Exception as e:
        raise HTTPException(
                status_code=400,
                detail=str(e),
                )

#save workflow data
@router.post("/save")
def update_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    workflow: WorkflowCreate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflow.id)
    if user_workflow:
        # workflow 수정
        return crud_workflow.update_workflow(db, current_user.id, workflow.id, workflow.title, workflow.thumbnail, workflow.workflow_info)
    else :
        # workflow 생성
        return crud_workflow.create_workflow(db, workflow.title, workflow.thumbnail, workflow.workflow_info, current_user.id)

#User workflow delete
@router.post("/delete")
def delete_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    workflowInfo: WorkflowDelete,
    ) -> Any:
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowInfo.id)
    if user_workflow:
        delete_workflow = crud_workflow.delete_user_workflow(db, current_user.id, workflowInfo.id)
        user_path = f"./user/{current_user.username}/"
        workflow_path = f"{user_path}workflow{workflowInfo.id}"
        if os.path.exists(workflow_path):
            shutil.rmtree(workflow_path)
        return delete_workflow
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#response workflow data
@router.get("/me")
def get_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    ) -> Any:
    user_workflows = crud_workflow.get_user_workflows(db, current_user.id)
    if user_workflows:
        res = []
        for item in user_workflows:
            workflow_res = {
                'id': item.id, 
                'title': item.title, 
                'thumbnail': item.thumbnail,
                'updated_at': item.updated_at, 
                'user_id': item.user_id
            }
            res.append(workflow_res)
        # print(res)
        return res
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#find user workflow
@router.post("/find", response_model=WorkflowCreate)
def find_user_workflow(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    workflowInfo: WorkflowFind,
    ) -> Any:
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowInfo.id)
    if user_workflow:
        return {
            'title': user_workflow.title,
            'thumbnail': user_workflow.thumbnail,
            'workflow_info': user_workflow.workflow_info,
        }
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#response workflow results
@router.post("/results")
def getResults(
    WorkflowResult: WorkflowResult, 
    current_user: models.User = Depends(dep.get_current_active_user)
):
    PATH_COMPILE_RESULT = f'./user/{current_user.username}/workflow_{WorkflowResult.id}/algorithm_{WorkflowResult.algorithm_id}/results'
    
    # 파일 리스트와 각 파일의 부가 정보 가져오기
    file_info_list = []
    if os.path.exists(PATH_COMPILE_RESULT):
        for file_name in os.listdir(PATH_COMPILE_RESULT):
            file_path = os.path.join(PATH_COMPILE_RESULT, file_name)
            if os.path.isfile(file_path):
                file_info = {
                    "name": file_name,
                    "size": os.path.getsize(file_path),  # 파일 크기
                    "modified_time": os.path.getmtime(file_path)  # 마지막 수정 시간
                }
                file_info_list.append(file_info)
    else:
        raise HTTPException(status_code=404, detail="Results directory not found.")
    
    return file_info_list

#response workflow result
@router.post("/result")
def checkResult(WorkflowResult: WorkflowResult, current_user: models.User = Depends(dep.get_current_active_user)):
    PATH_COMPILE_RESULT = f'./user/{current_user.username}/workflow_{WorkflowResult.id}/algorithm_{WorkflowResult.algorithm_id}/results'
    file_list = os.listdir(PATH_COMPILE_RESULT)
    # print(file_list)
    FILE_NAME = WorkflowResult.filename
    
    for item_file in file_list:
        if FILE_NAME in item_file:
            FILE_NAME = item_file
    # print(FILE_NAME)
    FILE_PATH = os.path.join(PATH_COMPILE_RESULT, FILE_NAME)
    # print(FILE_PATH)

    return FileResponse(FILE_PATH)

#response workflow visualization result
@router.post("/visualization/result")
def checkVisualizationResult(WorkflowResult: WorkflowResult, current_user: models.User = Depends(dep.get_current_active_user)):
    PATH_COMPILE_RESULT = f'./user/{current_user.username}/workflow_{WorkflowResult.id}/algorithm_{WorkflowResult.algorithm_id}/results'
    file_list = os.listdir(PATH_COMPILE_RESULT)
    # print(file_list)
    FILE_NAME = WorkflowResult.filename
    
    for item_file in file_list:
        if FILE_NAME in item_file:
            FILE_NAME = item_file
    # print(FILE_NAME)
    FILE_PATH = os.path.join(PATH_COMPILE_RESULT, FILE_NAME)
    # print(FILE_PATH)

    try:
        with open(FILE_PATH, "r") as f:
            plotly_data = json.load(f)
        return JSONResponse(content=plotly_data)
        
    except Exception as e:
        raise HTTPException(
                status_code=400,
                detail=str(e),
                )

#save workflow node modal data
@router.post("/node/save")
def saveNodeData(
    *,
    db: Session = Depends(dep.get_db),
    workflowNodeFileInfo: WorkflowNodeFileCreate, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowNodeFileInfo.id)
    if user_workflow:
        # workflowNodeFileInfo.id로 사용자 폴더에 workflow 폴더 생성
        user_path = f"./user/{current_user.username}/"
        workflow_path = f"{user_path}workflow_{workflowNodeFileInfo.id}"
        # user 폴더에 workflow 폴더 존재하지 않으면 생성
        if not os.path.exists(workflow_path):
            os.makedirs(workflow_path)
        # workflowNodeFileInfo.node_name이랑 workflowNodeFileInfo.node_id로 파일 이름 생성
        file_name = f"{workflowNodeFileInfo.node_name}_{workflowNodeFileInfo.node_id}.{workflowNodeFileInfo.file_extension}"
        result_file_path = f"{workflow_path}/{file_name}"
        # user 폴더에 파일 생성
        with open(result_file_path, "w") as f:
            # workflowNodeFileInfo.file_content가 List이면 json.dump으로 파일 생성
            if workflowNodeFileInfo.file_extension == "json":
                json.dump(workflowNodeFileInfo.file_content, f)
            if workflowNodeFileInfo.file_extension == "txt" or workflowNodeFileInfo.file_extension == "tsv" or workflowNodeFileInfo.file_extension == "csv":
                f.write(workflowNodeFileInfo.file_content)
        return {
            'id': workflowNodeFileInfo.id,
            'node_id': workflowNodeFileInfo.node_id,
            'node_name': workflowNodeFileInfo.node_name,
            'file_extension': workflowNodeFileInfo.file_extension,
            'file_path': result_file_path
        }
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#read workflow node modal data
@router.post("/node/read")
def readNodeData(
    *,
    db: Session = Depends(dep.get_db),
    workflowNodeFileInfo: WorkflowNodeFileRead, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowNodeFileInfo.id)
    if user_workflow:
        # user 폴더에 workflow 폴더 존재하지 않으면 에러 발생
        user_path = f"./user/{current_user.username}/"
        workflow_path = f"{user_path}workflow_{workflowNodeFileInfo.id}"
        print(workflow_path)
        if not os.path.exists(workflow_path):
            raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )
        # workflowNodeFileInfo.node_name이랑 workflowNodeFileInfo.node_id로 파일 읽어오기
        file_name = f"{workflowNodeFileInfo.node_name}_{workflowNodeFileInfo.node_id}.{workflowNodeFileInfo.file_extension}"
        with open(f"{workflow_path}/{file_name}", "r") as f:
            if workflowNodeFileInfo.file_extension == "json":
                file_content = json.load(f)
            else:
                file_content = f.read()

        return {
            'id': workflowNodeFileInfo.id,
            'node_id': workflowNodeFileInfo.node_id,
            'node_name': workflowNodeFileInfo.node_name,
            'file_content': file_content,
            'file_extension': workflowNodeFileInfo.file_extension
        }
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )

#delete workflow node modal data
@router.post("/node/delete")
def deleteNodeData(
    *,
    db: Session = Depends(dep.get_db),
    workflowNodeFileInfo: WorkflowNodeFileDelete, 
    current_user: models.User = Depends(dep.get_current_active_user)
    ):
    user_workflow = crud_workflow.get_user_workflow(db, current_user.id, workflowNodeFileInfo.id)
    if user_workflow:
        # user 폴더에 workflow 폴더 존재하지 않으면 에러 발생
        user_path = f"./user/{current_user.username}/"
        workflow_path = f"{user_path}workflow_{workflowNodeFileInfo.id}"
        if not os.path.exists(workflow_path):
            raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )
        # workflowNodeFileInfo.node_name이랑 workflowNodeFileInfo.node_id로 파일 삭제
        file_name = f"{workflowNodeFileInfo.node_name}_{workflowNodeFileInfo.node_id}.{workflowNodeFileInfo.file_extension}"
        os.remove(f"{workflow_path}/{file_name}")
        return {
            'id': workflowNodeFileInfo.id,
            'node_id': workflowNodeFileInfo.node_id,
            'node_name': workflowNodeFileInfo.node_name,
            'file_extension': workflowNodeFileInfo.file_extension
        }
    else:
        raise HTTPException(
                status_code=400,
                detail="this workflow not exists in your workflows",
                )
