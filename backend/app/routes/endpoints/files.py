from typing import Any
from venv import create
from fastapi import APIRouter, Body, Depends, HTTPException, Request, UploadFile, File
from typing import List, Union
from sqlalchemy.orm import Session
import os
import pandas as pd
import scanpy as sc
import json

from app.routes import dep
from app.database.crud import crud_file
from app.database import models
from app.database.schemas.file import FileCreate, FileDelete, FileUpdate, FileFind, FolderFind, FileSetup
from app.routes.h5ad_utils import organize_column_dtypes, get_annotation_columns, get_pseudotime_columns

router = APIRouter()

#workflow file-upload
@router.post("/upload")
async def fileUpload(
    *,
    db: Session = Depends(dep.get_db),
    files: List[UploadFile] = File(),
    current_user: models.User = Depends(dep.get_current_active_user),
    ) -> Any:
    #upload to directory
    UPLOAD_DIRECTORY = f'./user/{current_user.username}/data'
    for item_file in files:
        contents = await item_file.read()
        print(item_file.filename.split("_"))
        folder_file = item_file.filename.split('_')
        with open(os.path.join(UPLOAD_DIRECTORY, folder_file[1]), "wb") as f:
            f.write(contents)
        user_file = crud_file.get_user_file(db, current_user.id, folder_file[1])
        if user_file:
            raise HTTPException(
                status_code=400,
                detail="this file already exists in your files",
                )
        print(len(contents))
        create_file = crud_file.create_file(db, folder_file[1], len(contents), UPLOAD_DIRECTORY, folder_file[0], current_user.id)
        
    return create_file

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
        print(f"Directories: {dir_path}, Files: {file_names}")
    print(res)

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
def find_user_file(
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
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        delete_file = crud_file.delete_user_file(db, current_user.id, fileInfo.file_name)
        print(delete_file)
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
        print(update_file)
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

        adata = sc.read_h5ad(input_filepath)

        # Process and combine data to form the desired dataframe structure
        df = pd.concat([adata.obs, pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'], index=adata.obs.index)], axis=1)

        # Save the dataframe to a specific folder
        df.to_csv(output_filepath, index=False)
        return {'file_name': output_filename}
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
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        # 해당 유저의 폴더에 있는 파일을 가져와서 csv로 변환
        # 변환된 파일을 해당 유저의 폴더에 저장
        # 변환된 파일의 정보를 db에 저장
        # 유저 폴더 내 파일 경로 "user/{username}/data/{filename}.h5ad"
        # 변환된 파일 경로 "user/{username}/result/{filename}.csv"
        folder_path = './user' + '/' + current_user.username
        input_filename = fileInfo.file_name
        input_filepath = f"{folder_path}/data/{input_filename}"

        adata = sc.read_h5ad(input_filepath)
        adata.obs = organize_column_dtypes(adata.obs)
        anno_columns = get_annotation_columns(adata.obs)
        pseudo_columns = get_pseudotime_columns(adata.obs)
        return {'anno_columns': anno_columns, 'pseudo_columns': pseudo_columns}
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
        )
    
@router.post("/clusters")
def h5ad_cluster (
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
        folder_path = './user' + '/' + current_user.username
        input_filename = fileInfo.file_name
        input_filepath = f"{folder_path}/data/{input_filename}"

        adata = sc.read_h5ad(input_filepath)
        adata.obs = organize_column_dtypes(adata.obs)
        clusters = map(str, adata.obs[fileInfo.anno_column].value_counts().index)
        print(clusters)
        return {'clusters': list(clusters)}
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
        )
    
@router.post("/setup")
def algorithm_setup (
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    options: FileSetup,
    ) -> Any:
    user_file = crud_file.get_user_file(db, current_user.id, options.file_name)
    if user_file:
        folder_path = './user' + '/' + current_user.username
        input_filename = options.file_name
        input_filepath = f"{folder_path}/data/{input_filename}"

        # folder_path에 input_filename + '_option.json' 파일을 만들어서 저장
        # option 파일에는 anno_of_interest, pseudo_of_interest, clusters_of_interest, make_binary, device, device_ids, batch_size 정보를 저장
        # option 파일의 경로 "user/{username}/data/{filename}_option.json"

        options = {
            "anno_of_interest": options.anno_of_interest,
            "pseudo_of_interest": options.pseudo_of_interest,
            "clusters_of_interest": options.clusters_of_interest,
            "make_binary": options.make_binary,
            "device": options.device,
            "device_ids": options.device_ids,
            "batch_size": options.batch_size,
            "kp": options.kp,
            "percentile": options.percentile,
            "win_length": options.win_length,
            "polyorder": options.polyorder
        }

        with open(input_filepath + '_option.json', 'w', encoding='utf-8') as f:
            json.dump(options, f, ensure_ascii=False, indent=4)

        # 파일 생성 되었는지 확인
        if os.path.isfile(input_filepath + '_option.json'):
            return {'file_name': input_filename + '_option.json'}
        else:
            raise HTTPException(
                status_code=400,
                detail="option file not created",
            )
    else:
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
        )