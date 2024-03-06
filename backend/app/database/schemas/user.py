from typing import Optional
from pydantic import BaseModel, EmailStr

# Shared properties
class UserBase(BaseModel):
    email: Optional[EmailStr] = None
    is_active: Optional[bool] = True
    is_superuser: Optional[bool] = False
    username: Optional[str] = None


class UserCreate(UserBase):
    email: EmailStr
    password: str
    username: Optional[str] = None

class UserProfile(UserBase):
    email: EmailStr
    username: str
    is_superuser: bool

class User(UserBase):
    id: Optional[int] = None
    is_active: bool
    hashed_password: str

    class Config:
        orm_mode = True

