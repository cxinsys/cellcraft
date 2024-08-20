from fastapi import APIRouter, HTTPException
from fastapi import Depends

from app.routes import dep
from app.database import models
from app.common.utils.workflow_utils import load_tab_file, transform_df_to_vgt_format
from app.common.utils.datatable_utils import RemoteDataTable, DataTableEvent

router = APIRouter()

@router.post("/load_data")
def load_data(vgt_info: DataTableEvent, current_user: str = Depends(dep.get_current_active_user)):
    try:
        PATH_DATATABLE_FILE = f'./user/{current_user.username}/data/{vgt_info.file_name}'
        df = load_tab_file(PATH_DATATABLE_FILE)
        remote_table = RemoteDataTable(df, vgt_info)
        data = remote_table.get_filtered_sorted_paginated_data()
        vgt_data = transform_df_to_vgt_format(data["df"])
        return {"status": "success", "columns": vgt_data['columns'], "rows": vgt_data["rows"], "totalRecords": data["totalRecords"]}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
