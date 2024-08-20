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

from app.routes import dep
from app.database.crud import crud_file
from app.database import models
from app.database.schemas.file import FileCreate, FileDelete, FileUpdate, FileFind, FolderFind, FileSetup, FileGet, FileResultFind
from app.common.utils.h5ad_utils import organize_column_dtypes, get_annotation_columns, get_pseudotime_columns
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
    #upload to directory
    UPLOAD_DIRECTORY = f'./user/{current_user.username}/data'
    for item_file in files:
        contents = await item_file.read()
        # print(item_file.filename.split("_", 1))
        # item_file.filename에 h5ad가 포함되어 있으면 다음 과정 진행
        # 포함되어 있지 않으면 DB에 저장하지 않고 파일 생성만 진행
        if 'h5ad' not in item_file.filename:
            # 폴더 안에 파일이 있는지 확인
            # 있으면 파일 이름을 그대로 파일 내용만 변경해서 저장
            # 없으면 그냥 저장
            with open(os.path.join(UPLOAD_DIRECTORY, item_file.filename), "wb") as f:
                f.write(contents)
            return {'file_name': item_file.filename}

        folder_file = item_file.filename.split('_', 1)
        with open(os.path.join(UPLOAD_DIRECTORY, folder_file[1]), "wb") as f:
            f.write(contents)
        user_file = crud_file.get_user_file(db, current_user.id, folder_file[1])
        if user_file:
            raise HTTPException(
                status_code=400,
                detail="this file already exists in your files",
                )
        # print(len(contents))
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

        # Check if file exists, if not, continue
        if not os.path.isfile(output_filepath):
            # Read the h5ad file
            adata = sc.read_h5ad(input_filepath)

            # Process and combine data to form the desired dataframe structure
            df = pd.concat([adata.obs, pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'], index=adata.obs.index)], axis=1)

            # Save the dataframe to a specific folder
            df.to_csv(output_filepath, index=False)
            return {'file_name': output_filename}
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
    user_file = crud_file.get_user_file(db, current_user.id, fileInfo.file_name)
    if user_file:
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
        # print(clusters)
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
    current_time = datetime.now().strftime("%Y%m%d%H%M")
    user_file = crud_file.get_user_file(db, current_user.id, options.file_name)
    if user_file:
        folder_path = './user' + '/' + current_user.username + "/data/"
        input_filename = f"{current_time}_{options.option_name}_{options.file_name}"
        input_filepath = f"{folder_path}{input_filename}"

        gene_list_filename = options.gene_list
        gene_list_filepath = f"{folder_path}{gene_list_filename}"
        # gene_list_filename이 존재하면 해당 파일을 읽어서 gene_list에 할당 (csv 파일)
        # 존재하지 않으면 gene_list에 None 할당
        if os.path.isfile(gene_list_filepath):
            # pd.read_csv(gene_list_filepath, header=None) : csv 파일을 읽어서 dataframe으로 만듦
            # .iloc[:, 0] : dataframe의 첫번째 열을 가져옴
            # .tolist() : 첫번째 열을 list로 변환
            gene_list = pd.read_csv(gene_list_filepath, header=None).iloc[:, 0].tolist()
        else:
            gene_list = []

        setOptions = {
            "algorithm": options.algorithm,
            "anno_of_interest": options.anno_of_interest,
            "pseudo_of_interest": options.pseudo_of_interest,
            "clusters_of_interest": options.clusters_of_interest,
            "num_of_threads": options.num_of_threads,
            "history_length": options.history_length,
            "species": options.species,
            "gene_list_file": options.gene_list,
            "gene_list": gene_list,
            "cutoff_for_fdr": options.cutoff_for_fdr,
            "num_of_links": options.num_of_links,
            "trimming_indirect_edges": options.trimming_indirect_edges,
        }

        if len(options.clusters_of_interest) == 0:
            setOptions["selected_indices"] = options.selected_indices

        # print(setOptions)

        with open(input_filepath + '_option.json', 'w', encoding='utf-8') as f:
            json.dump(setOptions, f, ensure_ascii=False, indent=4)

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
    folder_path = './user' + '/' + current_user.username
    input_filepath = f"{folder_path}/data/{file_name}"
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
        
@router.get("/html/{filename}", response_class=HTMLResponse)
async def read_html_file(filename: str):
    with open(f'./tutorials/{filename}.html', "r") as f:
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
    folder_path = './user' + '/' + current_user.username + "/result/"
    return FileResponse(folder_path + filename ,filename=filename)

@router.get("/data/{filename}")
async def download_data_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    filename: str,
    ) -> Any:
    PATH_DATA_FILE = './user' + '/' + current_user.username + "/data/" + filename
    # 파일이 존재하는지 확인
    if not os.path.isfile(PATH_DATA_FILE):
        raise HTTPException(
                status_code=400,
                detail="this file not exists in your files",
        )

    if filename.endswith('.h5ad'):
        df = load_tab_file(PATH_DATA_FILE)
        return df.to_dict(orient="records")
    return FileResponse(PATH_DATA_FILE, filename=filename)