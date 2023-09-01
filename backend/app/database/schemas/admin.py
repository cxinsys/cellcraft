from typing import Optional
from pydantic import BaseModel, EmailStr

class Conditions(BaseModel):
    amount: int 
    page_num: int
    sort: str
    order: str
    searchTerm: Optional[str] = None
