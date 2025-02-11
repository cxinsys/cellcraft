import os
import pandas as pd
import scanpy as sc
import numpy as np

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
        # 기본 데이터프레임으로 adata.obs 사용
        df = adata.obs.copy()
        
        # X_umap이 있는 경우에만 추가
        if 'X_umap' in adata.obsm:
            umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'], index=adata.obs.index)
            df = pd.concat([df, umap_df], axis=1)
        # X_pca가 있는 경우 대체 사용
        elif 'X_pca' in adata.obsm:
            pca_df = pd.DataFrame(adata.obsm['X_pca'][:, :2], columns=['X', 'Y'], index=adata.obs.index)
            df = pd.concat([df, pca_df], axis=1)
        # 둘 다 없는 경우 raw 데이터에서 처음 두 컬럼 사용
        else:
            # 데이터가 1차원인 경우 처리
            if adata.X.shape[1] == 1:
                raw_df = pd.DataFrame({
                    'X': adata.X[:, 0],
                    'Y': np.zeros(adata.X.shape[0])  # Y값을 0으로 설정
                }, index=adata.obs.index)
            else:
                raw_df = pd.DataFrame(adata.X[:, :2], columns=['X', 'Y'], index=adata.obs.index)
            df = pd.concat([df, raw_df], axis=1)
            
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

def extract_all_algorithms(workflow_info):
    algorithms = []
    
    # 탐색하여 class가 "Algorithm"인 모든 객체를 찾음
    for key, value in workflow_info.items():
        if value.get("class") == "Algorithm":
            algorithm_data = value.get("data", {})
            algorithm_id = value.get("id")

            algorithms.append({
                "id": algorithm_id,
                "files": algorithm_data.get("files"),
                "selectedPlugin": algorithm_data.get("selectedPlugin"),
                "selectedPluginInputOutput": algorithm_data.get("selectedPluginInputOutput"),
                "selectedPluginRules": algorithm_data.get("selectedPluginRules")
            })

    return algorithms

def extract_visualization_data(workflow_info, node_id):
    # 탐색하여 class가 "Visualization"이고 id가 node_id인 객체를 찾음
    for key, value in workflow_info.items():
        print(value.get("class"))
        print(value.get("id"))
        if value.get("class") == "Visualization" and str(value.get("id")) == node_id:
            visualization_data = value.get("data", {})
            
            # 필드 추출
            selected_visualization_params= visualization_data.get("selectedVisualizationParams")
            selected_visualization_title = visualization_data.get("selectedVisualizationTitle")

            # 결과 반환
            return {
                "selectedVisualizationParams": selected_visualization_params,
                "selectedVisualizationTitle": selected_visualization_title
            }
    
    # "Visualization" 클래스를 가진 객체가 없을 경우 빈 딕셔너리 반환
    return {}

def generate_user_input(selectedPluginInputOutput):
    user_input = {}

    # selectedPluginInputOutput에서 type이 "input"인 파라미터 추출 및 사용자 입력에 추가
    for parameter in selectedPluginInputOutput:
        if parameter.get("type") == "inputFile" or parameter.get("type") == "optionalInputFile":
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

def generate_visualization_params(selectedVisualizationParams):
    visualization_params = {}
    visualization_inputs = {}
    visualization_outputs = {}

    # 모든 파라미터를 적절한 딕셔너리에 분류
    for parameter in selectedVisualizationParams:
        param_type = parameter.get("type")
        param_name = parameter.get("name")
        param_value = parameter.get("defaultValue")
        
        # inputFile과 optionalInputFile은 visualization_inputs에 할당
        if param_type in ["inputFile", "optionalInputFile"]:
            if "target" in param_name:
                param_name = param_name + parameter.get("fileExtension")
            visualization_inputs[param_name] = param_value
            
        # outputFile은 visualization_outputs에 할당
        elif param_type == "outputFile":
            visualization_outputs[param_name] = param_value
            
        # 나머지는 visualization_params에 할당
        else:
            visualization_params[param_name] = param_value

    return visualization_params, visualization_inputs, visualization_outputs

def extract_target_data(selectedPluginInputOutput, user_workflow_task_path):
    data_list = []

    user_workflow_task_path = os.path.join(user_workflow_task_path, "results")

    # type이 "output"인 데이터들 중, activated가 True인 데이터들을 찾아서 리스트에 추가
    for data in selectedPluginInputOutput:
        if data.get("type") == "outputFile" and data.get("activate"):
            data_list.append(os.path.join(user_workflow_task_path, data.get("defaultValue")))

    # 데이터 리스트 반환
    return data_list

def extract_rule_block(snakefile_path: str, visualization_title: str, result_path: str) -> str:
    # rule 블록들을 분리하기 위한 시작 패턴
    rule_pattern = f"rule {visualization_title}:"
    
    # 파일 내용을 줄 단위로 분리
    with open(snakefile_path, 'r') as file:
        lines = file.readlines()
    
    result = []
    is_target_rule = False
    current_indent = 0
    
    for line in lines:
        # 대상 rule 블록의 시작을 찾음
        if line.startswith(rule_pattern):
            is_target_rule = True
            current_indent = len(line) - len(line.lstrip())
            result.append(line)
            continue
            
        # rule 블록 내부의 내용을 수집
        if is_target_rule:
            # 현재 줄의 들여쓰기 수준 확인
            if line.strip() == "":
                continue
                
            line_indent = len(line) - len(line.lstrip())
            
            # 들여쓰기가 같거나 큰 경우에만 해당 rule에 속한 것으로 간주
            if line_indent > current_indent:
                result.append(line)
            else:
                break
                
    # 빈 줄 제거하고 결과 문자열 생성
    result = [line.rstrip() for line in result if line.strip()]
    extracted_content = 'import os\n' + '\n'.join(result)
    
    # 추출한 내용을 파일에 쓰기
    with open(result_path, 'w') as file:
        file.write(extracted_content)

    return extracted_content, result_path