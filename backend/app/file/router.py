from typing import Any, List
from fastapi import APIRouter, Depends, UploadFile, File
from fastapi.responses import HTMLResponse
from sqlalchemy.orm import Session

from app.auth import deps as dep
from app import models
from app.file.schemas import FileDelete, FileUpdate, FileFind, FolderFind, FileResultFind
from app.file import service

router = APIRouter()


# h5ad file-upload
@router.post("/upload")
async def fileUpload(
    *,
    db: Session = Depends(dep.get_db),
    files: List[UploadFile] = File(),
    current_user: models.User = Depends(dep.get_current_active_user),
) -> Any:
    return await service.save_upload(db=db, files=files, current_user=current_user)


# User Files get
@router.get("/me")
def read_user_folder(
    current_user: models.User = Depends(dep.get_current_active_user)
    ) -> Any:
    return service.read_user_folder(current_user=current_user)


# User File find
@router.post("/find")
def find_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    return service.find_user_file(db=db, current_user=current_user, file_info=fileInfo)


# User find File of folder
@router.post("/folder")
def find_user_folder(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    folder: FolderFind,
    ) -> Any:
    return service.find_user_folder(db=db, current_user=current_user, folder=folder)


# User File delete
@router.post("/delete")
def delete_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileDelete,
    ) -> Any:
    return service.delete_file(db=db, current_user=current_user, file_info=fileInfo)


# User File Update
@router.post("/update")
def update_user_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileUpdate,
    ) -> Any:
    return service.update_file(db=db, current_user=current_user, file_info=fileInfo)


@router.post("/convert")
def user_file_convert(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    return service.convert_file(db=db, current_user=current_user, file_info=fileInfo)


@router.get("/check/{file_name}")
def check_user_file_convert(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    file_name: str,
    ) -> Any:
    return service.check_convert(current_user=current_user, file_name=file_name)


@router.post("/columns")
def h5ad_columns(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    return service.get_h5ad_columns(db=db, current_user=current_user, file_info=fileInfo)


@router.post("/clusters")
def h5ad_cluster(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
    ) -> Any:
    return service.get_h5ad_clusters(db=db, current_user=current_user, file_info=fileInfo)


@router.get("/setup/check")
def algorithm_setup_check(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    ) -> Any:
    return service.list_option_files(current_user=current_user)


@router.get("/setup/{file_name}")
def read_user_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    file_name: str,
    ) -> Any:
    return service.read_setup_file(current_user=current_user, file_name=file_name)


@router.get("/shared")
def get_shared_files(
    current_user: models.User = Depends(dep.get_current_active_user),
) -> Any:
    """tutorials/ 디렉토리의 데이터 파일 목록 반환 (html 제외)"""
    return service.get_shared_files()


@router.post("/shared/find")
def find_shared_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileFind,
) -> Any:
    """공용 파일 단건 조회 (InputFile.vue에서 사용)"""
    return service.find_shared_file(file_info=fileInfo)


@router.get("/html/{filename}", response_class=HTMLResponse)
async def read_html_file(
    filename: str,
) -> Any:
    return service.read_html_file(filename=filename)


@router.post("/result")
async def read_result_file(
    *,
    db: Session = Depends(dep.get_db),
    current_user: models.User = Depends(dep.get_current_active_user),
    fileInfo: FileResultFind,
    ) -> Any:
    return service.read_result_file(db=db, current_user=current_user, file_info=fileInfo)


@router.get("/result/{filename}")
async def download_result_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    filename: str,
    ) -> Any:
    return service.download_result_file(current_user=current_user, filename=filename)


@router.get("/data/{filename}")
async def download_data_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    filename: str,
    source: str = None,
    ) -> Any:
    return await service.download_data_file(
        current_user=current_user, filename=filename, source=source
    )


@router.get("/tutorials/{filename}")
async def download_tutorial_file(
    *,
    current_user: models.User = Depends(dep.get_current_active_user),
    filename: str,
    ) -> Any:
    return service.download_tutorial_file(current_user=current_user, filename=filename)
