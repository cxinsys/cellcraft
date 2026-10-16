"""
File domain service layer.

Extracted from ``file/router.py`` in PR-8 (Phase 3d). This module holds the
business logic that previously lived inline in the endpoints: streaming uploads
(with size/extension validation via ``file/security.py``), user file CRUD, h5ad
column/cluster extraction, setup/result/shared/tutorial file access, and the
download handlers. The router now only parses requests, resolves dependencies,
calls a service function, and returns its result.

Exception policy (PR-8): business-rule failures raise the domain exceptions from
``app.core.exceptions`` (``NotFoundError`` -> 404, ``ValidationFailedError`` ->
400, ``CellcraftError`` -> 500). The global handler in ``app.main`` maps these
onto the exact ``{"detail": ...}`` wire format FastAPI produces for
``HTTPException``, so responses (status code + detail string) are unchanged.
Path-traversal rejections still come from ``file/security.py`` as ``HTTPException``
and are handled by FastAPI's default handler (also unchanged).
"""
import os
import gc
import json
import asyncio
from typing import Any, List, Optional

import pandas as pd
import scanpy as sc
from fastapi import HTTPException, UploadFile
from fastapi.responses import HTMLResponse, FileResponse
from sqlalchemy.orm import Session

from app import models
from app.core.exceptions import NotFoundError, ValidationFailedError, CellcraftError
from app.file import crud as crud_file
from app.file.schemas import FileFind, FileDelete, FileUpdate, FolderFind, FileResultFind
from app.datatable.h5ad import (
    organize_column_dtypes, get_annotation_columns, get_pseudotime_columns
)
from app.datatable.cache import (
    get_cached_columns, set_cached_columns,
    get_cached_clusters, set_cached_clusters, invalidate_file
)
from app.workflow.utils import load_tab_file

SHARED_DIR = './tutorials'
DATA_EXTENSIONS = {'.h5ad', '.csv', '.txt', '.json'}


async def save_upload(*, db: Session, files: List[UploadFile],
                      current_user: models.User) -> Any:
    """Stream-validate-persist one or more uploaded data files for a user."""
    from app.core.config import MAX_FILES_PER_UPLOAD
    from app.file.security import validate_file_upload

    # Validation 1: Check batch upload limit
    if len(files) > MAX_FILES_PER_UPLOAD:
        raise ValidationFailedError(
            f"Too many files. Maximum {MAX_FILES_PER_UPLOAD} files per upload."
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
            raise ValidationFailedError(
                f"Invalid filename format. Expected 'folder_filename.ext', got '{safe_filename}'"
            )

        folder, actual_filename = folder_file
        final_filename = actual_filename

        # Check for duplicate files
        user_file = crud_file.get_user_file(db, current_user.id, final_filename)
        if user_file:
            raise ValidationFailedError(
                f"File '{final_filename}' already exists in your files"
            )

        # Stream file to disk in chunks to prevent memory exhaustion
        from app.core.config import UPLOAD_CHUNK_SIZE, MAX_UPLOAD_SIZE

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
                        raise ValidationFailedError(
                            f"File too large. Maximum size is {MAX_UPLOAD_SIZE / (1024*1024):.0f}MB",
                            status_code=413,
                        )
        except (HTTPException, CellcraftError):
            raise
        except Exception as e:
            # Clean up partial file on error
            if os.path.exists(file_path):
                os.remove(file_path)
            raise CellcraftError(f"File upload failed: {str(e)}")

        # Compress H5AD files to reduce storage
        if final_filename.lower().endswith('.h5ad'):
            from app.datatable.h5ad_compression import compress_h5ad_file
            from app.core.config import H5AD_COMPRESSION_ENABLED, H5AD_COMPRESSION_MIN_SIZE

            if H5AD_COMPRESSION_ENABLED:
                _compressed, file_size = await asyncio.to_thread(
                    compress_h5ad_file,
                    file_path,
                    min_size_bytes=H5AD_COMPRESSION_MIN_SIZE,
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


def read_user_folder(*, current_user: models.User) -> Any:
    """Walk the user's data folder and return {folder: [filenames]} entries."""
    USER_FOLDER = f'./user/{current_user.username}/data'
    res = []
    for (dir_path, dir_names, file_names) in os.walk(USER_FOLDER):
        folder = dir_path.replace(USER_FOLDER, 'data')
        res.extend({folder: file_names}.items())
    return res


def find_user_file(*, db: Session, current_user: models.User, file_info: FileFind) -> Any:
    """Return a user's file record by name (400 if the user has no such file)."""
    user_file = crud_file.get_user_file(db, current_user.id, file_info.file_name)
    if user_file:
        return user_file
    raise ValidationFailedError("this file not exists in your files")


def find_user_folder(*, db: Session, current_user: models.User, folder: FolderFind) -> Any:
    """Return a user's files within a folder (400 if the folder is unknown)."""
    user_file = crud_file.get_user_folder(db, current_user.id, folder.folder_name)
    if user_file:
        return user_file
    raise ValidationFailedError("this folder not exists in your folders")


def delete_file(*, db: Session, current_user: models.User, file_info: FileDelete) -> Any:
    """Delete a user's file from disk (path-traversal guarded) and its DB row."""
    from app.file.security import validate_file_path

    user_file = crud_file.get_user_file(db, current_user.id, file_info.file_name)
    if user_file:
        # 실제 파일 경로 구성
        user_folder = f'./user/{current_user.username}/data'
        file_path = os.path.join(user_folder, file_info.file_name)

        # 보안: Path traversal 방지 (centralized validation)
        validate_file_path(user_folder, file_path)

        # 파일 존재 여부 확인 및 삭제
        if os.path.exists(file_path) and os.path.isfile(file_path):
            try:
                from app.datatable.cache import invalidate_datatable
                invalidate_file(file_path)
                invalidate_datatable(file_path)
                os.remove(file_path)
            except Exception as e:
                raise CellcraftError(f"Failed to delete file: {str(e)}")

        # DB에서 파일 정보 삭제
        return crud_file.delete_user_file(db, current_user.id, file_info.file_name)
    raise ValidationFailedError("this file not exists in your files")


def update_file(*, db: Session, current_user: models.User, file_info: FileUpdate) -> Any:
    """Update a user's file metadata (400 if the user has no such file)."""
    user_file = crud_file.get_user_file(db, current_user.id, file_info.file_name)
    if user_file:
        return crud_file.update_user_file(db, current_user.id, file_info)
    raise ValidationFailedError("this file not exists in your files")


def convert_file(*, db: Session, current_user: models.User, file_info: FileFind) -> Any:
    """Convert a user's h5ad file to a `_obs_umap.csv` result (idempotent)."""
    user_file = crud_file.get_user_file(db, current_user.id, file_info.file_name)
    if not user_file:
        raise ValidationFailedError("this file not exists in your files")

    folder_path = './user' + '/' + current_user.username
    input_filename = file_info.file_name
    output_filename = file_info.file_name.replace('.h5ad', '') + '_obs_umap.csv'
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


def check_convert(*, current_user: models.User, file_name: str) -> Any:
    """Return the converted csv filename if it already exists (400 otherwise)."""
    folder_path = './user' + '/' + current_user.username
    output_filename = file_name.replace('.h5ad', '') + '_obs_umap.csv'
    output_filepath = f"{folder_path}/result/{output_filename}"
    if os.path.isfile(output_filepath):
        return {'file_name': output_filename}
    raise ValidationFailedError("this file not exists in your files")


def get_h5ad_columns(*, db: Session, current_user: models.User, file_info: FileFind) -> Any:
    """Return organized annotation/pseudotime columns for an h5ad file (cached)."""
    from app.file.path_resolver import resolve_data_file_path

    # 공용 파일은 DB 소유권 검사 스킵
    if file_info.source != "shared":
        user_file = crud_file.get_user_file(db, current_user.id, file_info.file_name)
        if not user_file:
            raise ValidationFailedError("this file not exists in your files")

    input_filepath = resolve_data_file_path(file_info.file_name, current_user.username, file_info.source)

    if not os.path.isfile(input_filepath):
        raise NotFoundError("File not found")

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


def get_h5ad_clusters(*, db: Session, current_user: models.User, file_info: FileFind) -> Any:
    """Return the distinct cluster labels for an h5ad annotation column (cached)."""
    from app.file.path_resolver import resolve_data_file_path

    # 공용 파일은 DB 소유권 검사 스킵
    if file_info.source != "shared":
        user_file = crud_file.get_user_file(db, current_user.id, file_info.file_name)
        if not user_file:
            raise ValidationFailedError("this file not exists in your files")

    input_filepath = resolve_data_file_path(file_info.file_name, current_user.username, file_info.source)

    if not os.path.isfile(input_filepath):
        raise NotFoundError("File not found")

    cached = get_cached_clusters(input_filepath, file_info.anno_column)
    if cached is not None:
        return {'clusters': cached}

    adata = None
    try:
        adata = sc.read_h5ad(input_filepath)
        adata.obs = organize_column_dtypes(adata.obs)
        clusters = list(map(str, adata.obs[file_info.anno_column].value_counts().index))
        set_cached_clusters(input_filepath, file_info.anno_column, clusters)
        return {'clusters': clusters}
    finally:
        if adata is not None:
            del adata
        gc.collect()


def list_option_files(*, current_user: models.User) -> Any:
    """List `.json` option files present in the user's data folder."""
    folder_path = './user' + '/' + current_user.username
    input_filepath = f"{folder_path}/data/"

    files = os.listdir(input_filepath)
    option_files = []
    for file in files:
        if file.endswith('.json'):
            option_files.append(file)
    return {'option_files': option_files}


def read_setup_file(*, current_user: models.User, file_name: str) -> Any:
    """Read a setup file from the user's data folder (json parsed, else raw bytes)."""
    from app.file.security import validate_file_path

    folder_path = f'./user/{current_user.username}/data'
    input_filepath = os.path.join(folder_path, file_name)

    # Security: Prevent path traversal
    validate_file_path(folder_path, input_filepath)

    if not os.path.exists(input_filepath):
        raise NotFoundError("File not found")

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


def get_shared_files() -> list:
    """List data files (html excluded) available under the shared tutorials dir."""
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


def find_shared_file(*, file_info: FileFind) -> Any:
    """Return a single shared file's metadata (path-traversal guarded)."""
    from app.file.security import validate_file_path

    file_path = os.path.join(SHARED_DIR, file_info.file_name)
    validate_file_path(SHARED_DIR, file_path)

    if not os.path.exists(file_path):
        raise NotFoundError("Shared file not found")

    return {
        "file_name": file_info.file_name,
        "file_size": str(os.path.getsize(file_path)),
        "file_path": SHARED_DIR,
        "folder": "tutorials",
        "source": "shared",
    }


def read_html_file(*, filename: str) -> HTMLResponse:
    """Return a tutorial HTML file's content (path-traversal guarded)."""
    from app.file.security import validate_file_path

    folder_path = './tutorials'
    file_path = os.path.join(folder_path, f'{filename}.html')

    validate_file_path(folder_path, file_path)

    if not os.path.exists(file_path):
        raise NotFoundError("File not found")

    with open(file_path, "r") as f:
        html = f.read()
    return HTMLResponse(content=html)


def read_result_file(*, db: Session, current_user: models.User,
                     file_info: FileResultFind) -> Any:
    """List result files whose names contain both the file and option filenames."""
    user_file = crud_file.get_user_file(db, current_user.id, file_info.file_name)
    if user_file:
        folder_path = './user' + '/' + current_user.username + "/result/"
        filename = file_info.file_name
        option_filename = file_info.option_file_name
        files = os.listdir(folder_path)
        result_files = []
        for file in files:
            if option_filename in file and filename in file:
                result_files.append(file)
        return {'result_files': result_files}
    raise ValidationFailedError("this file not exists in your files")


def download_result_file(*, current_user: models.User, filename: str) -> FileResponse:
    """Return a FileResponse for a file in the user's result folder (traversal guarded)."""
    from app.file.security import validate_file_path

    folder_path = f'./user/{current_user.username}/result'
    file_path = os.path.join(folder_path, filename)

    # Security: Prevent path traversal
    validate_file_path(folder_path, file_path)

    if not os.path.exists(file_path):
        raise NotFoundError("File not found")

    return FileResponse(file_path, filename=filename)


async def download_data_file(*, current_user: models.User, filename: str,
                             source: Optional[str] = None) -> Any:
    """Download a user/shared data file; h5ad is loaded and returned as records."""
    from app.file.path_resolver import resolve_data_file_path

    file_path = resolve_data_file_path(filename, current_user.username, source)

    # 파일이 존재하는지 확인
    if not os.path.isfile(file_path):
        raise ValidationFailedError("this file not exists in your files")

    if filename.endswith('.h5ad'):
        # asyncio.to_thread()를 사용하여 blocking I/O를 별도 스레드에서 실행
        # 이벤트 루프 블로킹 방지
        df = await asyncio.to_thread(load_tab_file, file_path)
        return df.to_dict(orient="records")
    return FileResponse(file_path, filename=filename)


def download_tutorial_file(*, current_user: models.User, filename: str) -> FileResponse:
    """Return a FileResponse for a tutorial file (path-traversal guarded)."""
    from app.file.security import validate_file_path

    folder_path = './tutorials'
    file_path = os.path.join(folder_path, filename)

    # Security: Prevent path traversal
    validate_file_path(folder_path, file_path)

    if not os.path.exists(file_path):
        raise NotFoundError("File not found")

    return FileResponse(file_path, filename=filename)
