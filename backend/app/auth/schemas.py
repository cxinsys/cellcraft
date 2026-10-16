from typing import Optional
from pydantic.main import BaseModel

class UserInfo(BaseModel):
    is_superuser: bool
    email: str
    is_active: bool
    must_change_password: bool = False


class PasswordChange(BaseModel):
    current_password: str
    new_password: str

class Token(BaseModel):
    access_token: str
    token_type: str
    user_info: UserInfo


class TokenPayload(BaseModel):
    sub: Optional[int] = None