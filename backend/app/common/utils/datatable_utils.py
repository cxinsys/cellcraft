import pandas as pd
from typing import Dict, List
from pydantic import BaseModel

class DataTableBase(BaseModel):
    file_name: str

class Sort(BaseModel):
    field: str
    type: str

class DataTableEvent(DataTableBase):
    page: int
    perPage: int
    columnFilters: Dict[str, str]
    sort: List[Sort]

class DataTable(DataTableBase):
    data: List[Dict[str, str]]

class RemoteDataTable:
    def __init__(self, df, vgt_info : DataTableEvent):
        self.df = df
        self.server_params = vgt_info
        self.total_records = len(df)
        self.filtered_sorted_df = df

    def apply_filters(self):
        df = self.df
        for column, value in self.server_params.columnFilters.items():
            df = df[df[column].astype(str).str.contains(value)]
        self.filtered_sorted_df = df

    def apply_sorting(self):
        if self.server_params.sort:
            sort_field = self.server_params.sort[0].field
            sort_type = self.server_params.sort[0].type
            if sort_field and sort_type:
                ascending = sort_type == 'asc'
                self.filtered_sorted_df = self.filtered_sorted_df.sort_values(by=sort_field, ascending=ascending)

    def get_paginated_data(self):
        start = (self.server_params.page - 1) * self.server_params.perPage
        end = start + self.server_params.perPage
        paginated_data = self.filtered_sorted_df.iloc[start:end]
        return {"totalRecords": self.total_records, "df": paginated_data}

    def get_filtered_sorted_paginated_data(self):
        self.apply_filters()
        self.apply_sorting()
        return self.get_paginated_data()
