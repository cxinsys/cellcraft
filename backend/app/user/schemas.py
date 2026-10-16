from typing import Optional
from pydantic import BaseModel, EmailStr, Field, validator

# Shared properties
class UserBase(BaseModel):
    email: Optional[EmailStr] = None
    is_active: Optional[bool] = True
    is_superuser: Optional[bool] = False
    username: Optional[str] = None


class UserCreate(UserBase):
    email: EmailStr = Field(...)
    password: str = Field(..., min_length=1)
    username: str = Field(..., min_length=1)

    @validator('username', 'password', pre=True)
    def prevent_empty_strings(cls, v, field):
        """Prevent empty strings for username and password fields.

        Fixes BUG-002: Empty string validation gap.
        Empty strings should be rejected with ValidationError.
        """
        if v is not None and isinstance(v, str) and len(v.strip()) == 0:
            raise ValueError(f'{field.name} cannot be empty')
        return v

class UserUpdate(UserBase):
    password: Optional[str] = None
    plugins: Optional[list] = None

class UserProfile(UserBase):
    id: int
    email: EmailStr
    username: str
    is_superuser: bool

    class Config:
        orm_mode = True

class User(UserBase):
    id: Optional[int] = None
    is_active: bool
    hashed_password: str

    class Config:
        orm_mode = True

