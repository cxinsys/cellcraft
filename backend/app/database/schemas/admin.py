from typing import Optional
from pydantic import BaseModel, EmailStr

from app.database.schemas.user import User

class Conditions(BaseModel):
    amount: int = 20, 
    page_num: int = 1, 
    sort: str = "id", 
    order: str = "asc", 
    searchTerm: str = None

class Users(BaseModel):
    users: list[User] = []
