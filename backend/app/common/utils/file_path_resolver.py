import os
from fastapi import HTTPException
from app.common.utils.file_security import validate_file_path

SHARED_DIR = './tutorials'


def resolve_data_file_path(filename: str, username: str, source: str = None) -> str:
    """사용자 파일 또는 공용(tutorials) 파일의 경로를 반환한다.

    Args:
        filename: 파일명
        username: 현재 사용자 이름
        source: "user" | "shared" | None (None은 "user"로 취급)

    Returns:
        검증된 파일 경로 문자열

    Raises:
        HTTPException: source 값이 유효하지 않거나 경로 순회 공격이 감지된 경우
    """
    if source == "shared":
        folder_path = SHARED_DIR
    elif source == "user" or source is None:
        folder_path = f"./user/{username}/data"
    else:
        raise HTTPException(status_code=400, detail=f"Invalid file source: {source}")

    file_path = os.path.join(folder_path, filename)
    validate_file_path(folder_path, file_path)
    return file_path
