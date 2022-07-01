from typing import List, Optional
from typing import Optional
from pydantic import BaseModel, EmailStr

# Shared properties
class UserBase(BaseModel):
    email: Optional[EmailStr] = None
    is_active: Optional[bool] = True
    is_superuser: bool = False
    full_name: Optional[str] = None


class UserCreate(UserBase):
    email: EmailStr
    password: str


class User(UserBase):
    id: Optional[int] = None
    is_active: bool

    class Config:
        orm_mode = True

# # Properties to receive via API on update
# class UserUpdate(UserBase):
#     password: Optional[str] = None
