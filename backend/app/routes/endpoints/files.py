from typing import Any
from venv import create
from fastapi import APIRouter, Body, Depends, HTTPException, Request, UploadFile, File
from fastapi.responses import HTMLResponse, FileResponse
from typing import List, Union
from sqlalchemy.orm import Session
from datetime import datetime
import os
import pandas as pd
import scanpy as sc
import json
import gc
import asyncio

from app.routes import dep
from app.database.crud import crud_file
from app.database import models
from app.database.schemas.file import FileCreate, FileDelete, FileUpdate, FileFind, FolderFind, FileGet, FileResultFind
from app.common.utils.h5ad_utils import organize_column_dtypes, get_annotation_columns, get_pseudotime_columns
from app.common.utils.h5ad_metadata_cache import get_cached_columns, set_cached_columns, get_cached_clusters, set_cached_clusters, invalidate_file
from app.common.utils.workflow_utils import load_tab_file

router = APIRouter()

#h5ad file-upload
@router.post("/upload")
async def fileUpload(
    *,
    db: Session = Depends(dep.get_db),
    files: List[UploadFile] = File(),
    current_user: models.User = Depends(dep.get_current_active_user),
) -> Any:
    from app.common.config import MAX_FILES_PER_UPLOAD
    from app.common.utils.file_security import validate_file_upload

    # Validation 1: Check batch upload limit
    if len(files) > MAX_FILES_PER_UPLOAD:
        raise HTTPException(
            status_code=400,
            detail=f"Too many files. Maximum {MAX_FILES_PER_UPLOAD} files per upload."
        )

    # Setup upload directory
    UPLOAD_DIRECTORY = f'./user/{current_user.username}/data'
    os.makedirs(UPLOAD_DIRECTORY, exist_ok=True)

    uploaded_files = []

    for item_file in files:
        # Validation 2: Comprehensive file validation
        safe_filename = validate_file_upload(item_file)

        # Parse folder and filename (format: "folder_filename.ext")
        folder_file = safe_filename.split('_', 1)
        if len(folder_file) != 2:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid filename format. Expected 'folder_filename.ext', got '{safe_filename}'"
            )

        folder, actual_filename = folder_file
        final_filename = actual_filename

        # Check for duplicate files
        user_file = crud_file.get_user_file(db, current_user.id, final_filename)
        if user_file:
            raise HTTPException(
                status_code=400,
                detail=f"File '{final_filename}' already exists in your files"
            )

        # Stream file to disk in chunks to prevent memory exhaustion
        from app.common.config import UPLOAD_CHUNK_SIZE, MAX_UPLOAD_SIZE

        file_path = os.path.join(UPLOAD_DIRECTORY, final_filename)
        file_size = 0

        try:
            with open(file_path, "wb") as f:
                while chunk := await item_file.read(UPLOAD_CHUNK_SIZE):
                    f.write(chunk)
                    file_size += len(chunk)

                    # Enforce size limit during streaming to fail fast
                    if file_size > MAX_UPLOAD_SIZE:
                        # Clean up partial file
                        f.close()
                        if os.path.exists(file_path):
                            os.remove(file_path)
                        raise HTTPException(
                            status_code=413,
                            detail=f"File too large. Maximum size is {MAX_UPLOAD_SIZE / (1024*1024):.0f}MB"
                        )
        except HTTPException:
            raise
        except Exception as e:
            # Clean up partial file on error
            if os.path.exists(file_path):
                os.remove(file_path)
            raise HTTPException(
                status_code=500,
                detail=f"File upload failed: {str(e)}"
            )

        # Create database record
        created_file = crud_file.create_file(
            db,
            final_filename,
            file_size,
            UPLOAD_DIRECTORY,
            folder,
            current_user.id
        )
        uploaded_files.append(created_file)

    # Return last uploaded file (maintain backward compatibility)
    # or return list if multiple files
    if len(uploaded_files) == 1:
        return uploaded_files[0]
    else:
        return {"files": uploaded_files, "count": len(uploaded_files)}

#User Files get
@router.get("/me")
def read_user_folder(
    current_user: models.User = Depends(dep.get_current_active_user)
    ) -> Any:
    USER_FOLDER = f'./user/{current_user.username}/data'
    res = []
    for (dir_path, dir_names, file_names) in os.walk(USER_FOLDER):
        folder = dir_path.replace(USER_FOLDER , 'data')
        res.extend({folder : file_names}.items())
        # print(f"Directories: {dir_path}, Files: {file_names}")
    # print(res)

    # print(user_files)
    return res

#User File find
@router.post("/find")
def find_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        return user_file
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
                )

#User find File of folder
@router.post("/folder")
def find_user_folder(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    folder: FolderFind,
    ) -> Any:
    user_file = crud_file.get_user_folder(db, current_user.id, folder.folder_name)
    if user_file:
        return user_file
    else:
        raise HTTPException(
                status_code=400,
                detail="this folder not exists in your folders",
                )

#User File delete
@router.post("/delete")
def delete_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileDelete,
    ) -> Any:
    from app.common.utils.file_security import validate_file_path

    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        # 실제 파일 경로 구성
        user_folder = f'./user/{current_user.username}/data'
        file_path = os.path.join(user_folder, fileInfo.file_name)

        # 보안: Path traversal 방지 (centralized validation)
        validate_file_path(user_folder, file_path)

        # 파일 존재 여부 확인 및 삭제
        if os.path.exists(file_path) and os.path.isfile(file_path):
            try:
                from app.common.utils.datatable_cache import invalidate_datatable
                invalidate_file(file_path)
                invalidate_datatable(file_path)
                os.remove(file_path)
            except Exception as e:
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to delete file: {str(e)}"
                )

        # DB에서 파일 정보 삭제
        delete_file = crud_file.delete_user_file(db, current_user.id, fileInfo.file_name)
        return delete_file
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
                )

#User File Update
@router.post("/update")
def update_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileUpdate,
    ) -> Any:
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        update_file = crud_file.update_user_file(db, current_user.id, fileInfo)
        # print(update_file)
        return update_file
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
                )
    
@router.post("/convert")
def user_file_convert(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        # 해당 유저의 폴더에 있는 파일을 가져와서 csv로 변환
        # 변환된 파일을 해당 유저의 폴더에 저장
        # 변환된 파일의 정보를 db에 저장
        # 유저 폴더 내 파일 경로 "user/{username}/data/{filename}.h5ad"
        # 변환된 파일 경로 "user/{username}/result/{filename}.csv"
        folder_path = './user' + '/' + current_user.username
        input_filename = fileInfo.file_name
        output_filename = fileInfo.file_name.replace('.h5ad', '') + '_obs_umap.csv'
        input_filepath = f"{folder_path}/data/{input_filename}"
        output_filepath = f"{folder_path}/result/{output_filename}"

        # Check if directory exists, if not, create it
        if not os.path.exists(folder_path + '/result'):
            os.makedirs(folder_path + '/result')

        # Check if file exists, if not, continue
        if not os.path.isfile(output_filepath):
            # 메모리 누수 방지를 위해 adata를 명시적으로 해제
            adata = None
            try:
                # Read the h5ad file
                adata = sc.read_h5ad(input_filepath)

                # Process and combine data to form the desired dataframe structure
                df = pd.concat([
                    adata.obs.copy(),
                    pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'], index=adata.obs.index)
                ], axis=1)

                # Save the dataframe to a specific folder
                df.to_csv(output_filepath, index=False)
                return {'file_name': output_filename}
            finally:
                if adata is not None:
                    del adata
                gc.collect()
        else:
            return {'file_name': output_filename}
            
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
                )

@router.get("/check/{file_name}")
def check_user_file_convert(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    file_name: str,
    ) -> Any:
    folder_path = './user' + '/' + current_user.username
    output_filename = file_name.replace('.h5ad', '') + '_obs_umap.csv'
    output_filepath = f"{folder_path}/result/{output_filename}"
    if os.path.isfile(output_filepath):
        return { 'file_name': output_filename }
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
                )


@router.post("/columns")
def h5ad_columns (
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    from app.common.utils.file_path_resolver import resolve_data_file_path

    # 공용 파일은 DB 소유권 검사 스킵
    if fileInfo.source != "shared":
        user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
        if not user_file:
            raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
            )

    input_filepath = resolve_data_file_path(fileInfo.file_name, current_user.username, fileInfo.source)

    if not os.path.isfile(input_filepath):
        raise HTTPException(status_code=404, detail="File not found")

    cached = get_cached_columns(input_filepath)
    if cached is not None:
        return cached

    adata = None
    try:
        adata = sc.read_h5ad(input_filepath)
        obs = organize_column_dtypes(adata.obs)
        anno_columns = get_annotation_columns(obs, organized=True)
        pseudo_columns = get_pseudotime_columns(obs, organized=True)
        result = {'anno_columns': anno_columns, 'pseudo_columns': pseudo_columns}
        set_cached_columns(input_filepath, result)
        return result
    finally:
        if adata is not None:
            del adata
        gc.collect()

@router.post("/clusters")
def h5ad_cluster (
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    from app.common.utils.file_path_resolver import resolve_data_file_path

    # 공용 파일은 DB 소유권 검사 스킵
    if fileInfo.source != "shared":
        user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
        if not user_file:
            raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
            )

    input_filepath = resolve_data_file_path(fileInfo.file_name, current_user.username, fileInfo.source)

    if not os.path.isfile(input_filepath):
        raise HTTPException(status_code=404, detail="File not found")

    cached = get_cached_clusters(input_filepath, fileInfo.anno_column)
    if cached is not None:
        return {'clusters': cached}

    adata = None
    try:
        adata = sc.read_h5ad(input_filepath)
        adata.obs = organize_column_dtypes(adata.obs)
        clusters = list(map(str, adata.obs[fileInfo.anno_column].value_counts().index))
        set_cached_clusters(input_filepath, fileInfo.anno_column, clusters)
        return {'clusters': clusters}
    finally:
        if adata is not None:
            del adata
        gc.collect()
    
@router.get("/setup/check")
def algorithm_setup_check (
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    ) -> Any:
    folder_path = './user' + '/' + current_user.username
    input_filepath = f"{folder_path}/data/"

    # input_filepath 안에 .option.json 파일이 있는지 확인 후, 있으면 모두 return
    files = os.listdir(input_filepath)
    option_files = []
    for file in files:
        if file.endswith('.json'):
            option_files.append(file)
    # print(option_files)
    return {'option_files': option_files}

@router.get("/setup/{file_name}")
def read_user_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    file_name: str,
    ) -> Any:
    from app.common.utils.file_security import validate_file_path

    folder_path = f'./user/{current_user.username}/data'
    input_filepath = os.path.join(folder_path, file_name)

    # Security: Prevent path traversal
    validate_file_path(folder_path, input_filepath)

    if not os.path.exists(input_filepath):
        raise HTTPException(status_code=404, detail="File not found")

    # file_name에 해당하는 파일이 .json 확장자면 json으로 로드해서 return
    # 아니면 파일을 읽어서 return
    if file_name.endswith('.json'):
        with open(input_filepath) as json_file:
            data = json.load(json_file)
            return data
    else:
        with open(input_filepath, 'rb') as file:
            contents = file.read()
            return contents
        
@router.get("/shared")
def get_shared_files(
    current_user: models.User = Depends(dep.get_current_active_user),
) -> Any:
    """tutorials/ 디렉토리의 데이터 파일 목록 반환 (html 제외)"""
    SHARED_DIR = './tutorials'
    DATA_EXTENSIONS = {'.h5ad', '.csv', '.txt', '.json'}
    shared_files = []

    if not os.path.isdir(SHARED_DIR):
        return shared_files

    for f in os.listdir(SHARED_DIR):
        ext = os.path.splitext(f)[1].lower()
        if ext in DATA_EXTENSIONS:
            file_path = os.path.join(SHARED_DIR, f)
            shared_files.append({
                "file_name": f,
                "file_size": str(os.path.getsize(file_path)),
                "file_path": SHARED_DIR,
                "folder": "tutorials",
                "source": "shared",
            })
    return shared_files


@router.post("/shared/find")
def find_shared_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
) -> Any:
    """공용 파일 단건 조회 (InputFile.vue에서 사용)"""
    from app.common.utils.file_security import validate_file_path

    SHARED_DIR = './tutorials'
    file_path = os.path.join(SHARED_DIR, fileInfo.file_name)
    validate_file_path(SHARED_DIR, file_path)

    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="Shared file not found")

    return {
        "file_name": fileInfo.file_name,
        "file_size": str(os.path.getsize(file_path)),
        "file_path": SHARED_DIR,
        "folder": "tutorials",
        "source": "shared",
    }


@router.get("/html/{filename}", response_class=HTMLResponse)
async def read_html_file(
    filename: str,
) -> Any:
    from app.common.utils.file_security import validate_file_path

    folder_path = './tutorials'
    file_path = os.path.join(folder_path, f'{filename}.html')

    validate_file_path(folder_path, file_path)

    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found")

    with open(file_path, "r") as f:
        html = f.read()
    return HTMLResponse(content=html)

@router.post("/result")
async def read_result_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileResultFind,
    ) -> Any:
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        folder_path = './user' + '/' + current_user.username + "/result/"
        filename = fileInfo.file_name
        option_filename = fileInfo.option_file_name
        ## folder_path 안에 파일들의 이름에 option_filename, filename이 포함되어 있는지 확인
        ## 둘 다 포함되어 있으면 파일들을 읽어서 return
        ## 파일들을 리스트 형식으로 return
        files = os.listdir(folder_path)
        result_files = []
        for file in files:
            if option_filename in file and filename in file:
                result_files.append(file)
        # print(result_files)
        return {'result_files': result_files}
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
        )

@router.get("/result/{filename}")
async def download_result_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    filename: str,
    ) -> Any:
    from app.common.utils.file_security import validate_file_path

    folder_path = f'./user/{current_user.username}/result'
    file_path = os.path.join(folder_path, filename)

    # Security: Prevent path traversal
    validate_file_path(folder_path, file_path)

    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(file_path, filename=filename)

@router.get("/data/{filename}")
async def download_data_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    filename: str,
    source: str = None,
    ) -> Any:
    from app.common.utils.file_path_resolver import resolve_data_file_path

    file_path = resolve_data_file_path(filename, current_user.username, source)

    # 파일이 존재하는지 확인
    if not os.path.isfile(file_path):
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
        )

    if filename.endswith('.h5ad'):
        # asyncio.to_thread()를 사용하여 blocking I/O를 별도 스레드에서 실행
        # 이벤트 루프 블로킹 방지
        df = await asyncio.to_thread(load_tab_file, file_path)
        return df.to_dict(orient="records")
    return FileResponse(file_path, filename=filename)

@router.get("/tutorials/{filename}")
async def download_tutorial_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    filename: str,
    ) -> Any:
    from app.common.utils.file_security import validate_file_path

    folder_path = './tutorials'
    file_path = os.path.join(folder_path, filename)

    # Security: Prevent path traversal
    validate_file_path(folder_path, file_path)

    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(file_path, filename=filename)