from fastapi import APIRouter
from fastapi import Depends

from app.auth import deps as dep
from app.datatable.utils import DataTableEvent
from app.datatable import service

router = APIRouter()


@router.post("/load_data")
def load_data(vgt_info: DataTableEvent, current_user: str = Depends(dep.get_current_active_user)):
    return service.load_data(vgt_info=vgt_info, current_user=current_user)
