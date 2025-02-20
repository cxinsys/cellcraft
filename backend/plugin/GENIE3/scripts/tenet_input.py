import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np
import sys
import json

def organize_column_dtypes(data_frame):
    df = data_frame.copy()
    for column in df.columns:
        try:
            if all(df[column].astype(float) % 1 == 0):
                df[column] = df[column].astype(int)
                check_if_counts = len(df[column].value_counts())   
                if column not in ['nCount_RNA','nFeature_RNA','n_genes','n_genes_by_counts', 'total_counts','total_counts_mt']:
                    df[column] = df[column].astype('category')
            else:
                df[column] = df[column].astype(float)
        except TypeError as e:
            df[column] = df[column].astype(str)
        except ValueError as e:
            df[column] = df[column].astype('category')
    return df

def check_raw_data(data):
    if isinstance(data, csr_matrix):
        if np.all(data.data % 1 != 0):  # 어떤 원소라도 정수가 아니라면
            return False
        else:
            return True
    elif isinstance(data, np.ndarray):
        if np.all(data % 1 != 0):  # 어떤 원소라도 정수가 아니라면
            return False
        else:
            return True
    else:
        return False

def make_cell_select(adata, categories, annotation_column = 'seurat_annotation'):
    categories = [int(cat) if cat.isdigit() else cat for cat in categories]
    adata.obs['cell_select'] = np.where(adata.obs[annotation_column].isin(categories),'1','0')

def make_cell_select_lasso(adata, selected_indices): # selected_indices: list of indices
    cell_selected_bools = np.array(['0'] * adata.obs.shape[0])
    cell_selected_bools[selected_indices] = '1'
    adata.obs['cell_select'] = cell_selected_bools

def select_cells(adata, pseudotime_column = 'pseudotime'):
    doi = adata[adata.obs['cell_select'] == '1']
    # print(doi.shape)
    if any(doi.obs[pseudotime_column] == float('inf')) or any(doi.obs[pseudotime_column].isna()):
        doi = doi[doi.obs[pseudotime_column] != float('inf')]
        doi = doi[~doi.obs[pseudotime_column].isna()]
    return doi

def synch_raw_counts(adata):
    try:
        if hasattr(adata, 'X') and check_raw_data(adata.X):
            print("Your adata object contains raw counts adata.X")
        if hasattr(adata, 'raw') and hasattr(adata.raw, 'X') and check_raw_data(adata.raw.X):
            print("Your adata object contains raw counts adata.raw.X")

        if 'counts' in adata.layers and check_raw_data(adata.layers['counts']):
            print("Your adata object already contains raw counts")
            return
        if hasattr(adata, 'X') and check_raw_data(adata.X):
            adata.layers['counts'] = adata.X
            return
        if hasattr(adata, 'raw') and hasattr(adata.raw, 'X') and check_raw_data(adata.raw.X):
            adata.layers['counts'] = adata.raw.X
            return
    except AttributeError as e:
        print("Your adata object does not contain raw counts", e)
    except KeyError as e:
        print("Your adata object does not contain raw counts", e)

def tenet_input_make(adata, path_exp, path_pseudo, path_cell_select, pseudotime_column = 'dpt_pseudotime', gene_list = []):
    try:
        if 'counts' not in adata.layers:
            print("adata.layers에 'counts' 키가 없습니다. 먼저 raw count 데이터를 설정해주세요.")
            return

        input_data = sc.get.obs_df(adata, keys=[pseudotime_column,'cell_select'])
        
        # gene_list가 비어있으면 모든 유전자 사용
        if not gene_list:
            gene_list = list(adata.var_names)
        else:
            # gene_list와 adata.var_names의 교집합 찾기
            available_genes = [gene for gene in gene_list if gene in adata.var_names]
            if not available_genes:
                print("입력된 gene_list의 유전자들이 데이터에 존재하지 않습니다.")
                return
            if len(available_genes) < len(gene_list):
                missing_genes = set(gene_list) - set(available_genes)
                print(f"다음 유전자들은 데이터에서 찾을 수 없습니다: {', '.join(missing_genes)}")
            gene_list = available_genes
        
        # 선택된 유전자들의 인덱스 찾기
        gene_indices = [adata.var_names.get_loc(gene) for gene in gene_list]
        
        # 매트릭스 데이터 추출
        if isinstance(adata.layers['counts'], csr_matrix):
            matrix_data = adata.layers['counts'][:, gene_indices].todense()
        else:
            matrix_data = adata.layers['counts'][:, gene_indices]
        
        # 데이터 저장
        pd.DataFrame(data=matrix_data,
                    index=adata.obs_names, 
                    columns=gene_list).to_csv(path_exp)
        
        #save pseudotime
        input_data[pseudotime_column].to_csv(path_pseudo, 
                                        index=False, header=False)
        
        #save cell sellect
        input_data['cell_select'].to_csv(path_cell_select, 
                                    index=False, header=False)
    except Exception as e:
        print(f"Error: {e}")

h5ad = sys.argv[1]
gene_list_path = None if sys.argv[2].lower() == "none" else sys.argv[2]
selected_indices_json_path = None if sys.argv[3].lower() == "none" else sys.argv[3]
expMatrix = sys.argv[4]
pseudotime = sys.argv[5]
cellSelect = sys.argv[6]
anno_of_interest = sys.argv[7]
pseudo_of_interest = sys.argv[8]
clusters_of_interest_string = sys.argv[9]

# gene_list.txt 파일에서 유전자 목록 읽기
gene_list = []
if gene_list_path:  # gene_list 파일 경로가 제공된 경우
    try:
        with open(gene_list_path, 'r') as f:
            gene_list = [line.strip() for line in f if line.strip()]  # 빈 줄 제거
    except FileNotFoundError:
        print(f"gene_list 파일을 찾을 수 없습니다: {gene_list_path}")
    except Exception as e:
        print(f"gene_list 파일을 읽는 중 오류가 발생했습니다: {e}")

clusters_of_interest = clusters_of_interest_string.split(';')
clusters_of_interest = [cluster for cluster in clusters_of_interest if cluster]

orig_data = sc.read_h5ad(h5ad)
adata = orig_data.copy()
adata.obs = organize_column_dtypes(adata.obs)

make_cell_select(adata, annotation_column = anno_of_interest, categories = clusters_of_interest)

# JSON 파일에서 선택된 인덱스 읽기 (선택적)
if selected_indices_json_path:  # JSON 파일 경로가 제공된 경우
    try:
        with open(selected_indices_json_path, 'r') as f:
            selected_indices = json.load(f)
        # make_cell_select_lasso 함수 호출하여 adata 업데이트
        make_cell_select_lasso(adata, selected_indices)
    except FileNotFoundError:
        print(f"선택된 인덱스 JSON 파일을 찾을 수 없습니다: {selected_indices_json_path}")
    except Exception as e:
        print(f"JSON 파일을 읽는 중 오류가 발생했습니다: {e}")

# cell_select column을 기반으로 cell 선택 및
# pseudotime column의 inf와 NA값 제거
adata = select_cells(adata, pseudotime_column = pseudo_of_interest)
    
# 뽑힌 cell이 없을시 프로그램 종료
if adata.shape[0] == 0:
    print('no cell selected')
    sys.exit(1)

# cell X gene matrix를 만들기 위해 raw counts가 
# adata.layers['counts']에 존재하도록 데이터 구조 동기화
synch_raw_counts(adata)
    
# tenet inputfile 생성
tenet_input_make(adata, expMatrix, pseudotime, cellSelect, pseudotime_column = pseudo_of_interest, gene_list = gene_list)