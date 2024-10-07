import os
import pandas as pd
import scanpy as sc

# def load_tab_file(file_path: str, max_rows: int = 5000):
def load_tab_file(file_path: str):
    # Check if the file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    # Extract the file extension
    file_extension = os.path.splitext(file_path)[1].lower()
    
    # Check if the file extension is either csv or tsv
    if file_extension not in ['.csv', '.tsv', '.h5ad']:
        raise ValueError(f"The file {file_path} must have a .csv or .tsv or .h5ad extension.")
    
    if file_extension == '.h5ad':
        adata = sc.read_h5ad(file_path)
        # Process and combine data to form the desired dataframe structure
        df = pd.concat([adata.obs, pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'], index=adata.obs.index)], axis=1)
        # slice the dataframe to the maximum number of rows
        # df_trimmed = df.head(max_rows)
        # return df_trimmed
        return df
    else :
        # Determine the separator based on the file extension
        sep = ',' if file_extension == '.csv' else '\t'
        
        # Load the file with pandas
        # df = pd.read_csv(file_path, sep=sep, nrows=max_rows)
        df = pd.read_csv(file_path, sep=sep)
        
        # Return the trimmed dataframe
        return df

def transform_df_to_vgt_format(df):
    columns = [
        {
            "label": col,
            "field": col,
            "type": "number" if df[col].dtype in [int, float] else "string"
        }
        for col in df.columns
    ]

    rows = df.to_dict(orient='records')
    
    return {"columns": columns, "rows": rows}

def extract_algorithm_data(workflow_info):
    # 탐색하여 class가 "Algorithm"인 객체를 찾음
    for key, value in workflow_info.items():
        if value.get("class") == "Algorithm":
            algorithm_data = value.get("data", {})
            algorithm_id = value.get("id")
            
            # 필드 추출
            files = algorithm_data.get("files")
            selected_plugin = algorithm_data.get("selectedPlugin")
            selected_plugin_input_output = algorithm_data.get("selectedPluginInputOutput")
            selected_plugin_rules = algorithm_data.get("selectedPluginRules")
            
            # 결과 반환
            return {
                "id": algorithm_id,
                "files": files,
                "selectedPlugin": selected_plugin,
                "selectedPluginInputOutput": selected_plugin_input_output,
                "selectedPluginRules": selected_plugin_rules
            }
    
    # "Algorithm" 클래스를 가진 객체가 없을 경우 빈 딕셔너리 반환
    return {}

def generate_user_input(selectedPluginInputOutput):
    user_input = {}

    # selectedPluginInputOutput에서 type이 "input"인 파라미터 추출 및 사용자 입력에 추가
    for parameter in selectedPluginInputOutput:
        if parameter.get("type") == "inputFile":
            # 파라미터 이름 추출
            parameter_key = parameter.get("defaultValue")
            
            # 사용자 입력에 추가
            user_input[parameter_key] = parameter.get("file_name")
    
    return user_input

def generate_plugin_params(selectedPluginRules):
    plugin_params = {}

    # selectedPluginRules에서 파라미터 추출 및 사용자 입력에 추가
    for rule in selectedPluginRules:
        for parameter in rule.get("parameters", []):
            # 파라미터 이름 추출
            parameter_name = parameter.get("name")
            
            plugin_params[parameter_name] = parameter.get("defaultValue")

    return plugin_params

def extract_target_data(selectedPluginInputOutput, user_workflow_task_path):
    data_list = []

    user_workflow_task_path = os.path.join(user_workflow_task_path, "results")

    # type이 "output"인 데이터들 중, activated가 True인 데이터들을 찾아서 리스트에 추가
    for data in selectedPluginInputOutput:
        if data.get("type") == "outputFile" and data.get("activate"):
            data_list.append(os.path.join(user_workflow_task_path, data.get("defaultValue")))

    # 데이터 리스트 반환
    return data_list
