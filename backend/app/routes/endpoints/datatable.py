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
        from app.common.utils.file_path_resolver import resolve_data_file_path
        PATH_DATATABLE_FILE = resolve_data_file_path(vgt_info.file_name, current_user.username, vgt_info.source)
        df = load_tab_file(PATH_DATATABLE_FILE)
        remote_table = RemoteDataTable(df, vgt_info)
        data = remote_table.get_filtered_sorted_paginated_data()
        vgt_data = transform_df_to_vgt_format(data["df"])
        return {"status": "success", "columns": vgt_data['columns'], "rows": vgt_data["rows"], "totalRecords": data["totalRecords"]}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
