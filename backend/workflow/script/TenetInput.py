import matplotlib.pyplot as plt
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np
import json

def organize_column_dtypes(data_frame):
    '''
    dataframe의 각 column의 dtype을 정리합니다. 범주형 변수의 dtype은
    object, count data의 dtype은 int, 그 외 연속형 변수의 dtype은 float로
    맞춰줍니다. 
    '''
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
    '''
    data에 raw data(정수값)이 들어있는지 확인합니다.
    '''
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
    '''
    make cell select함수는 유저가 관심있는 annotation column을 드롭다운 메뉴로 
    선택 -> UMAP(scatterplot)에 annotation을 기반으로 coloring이 이루어지면 
    유저가 마우스를 통해 클러스터들을 선택하는 형태를 구현하면 좋을 것 같습니다.
    categories(cluster / celltype / cell state 등등)을 체크박스의 형태로 선택
    하도록 해도 좋을 것 같습니다. adata.obs에 'cell_select' column을 추가합니다.
    '''
    categories = [int(cat) for cat in categories]
    adata.obs['cell_select'] = np.where(adata.obs[annotation_column].isin(categories),'1','0')

def select_cells(adata, pseudotime_column = 'pseudotime'):
    '''
    cell_select column을 기반으로 관심있는 cell만을 추출한 adata object를
    반환합니다. tenet을 돌릴 때 pseudotime 값이 누락되거나 inf인 값이 있으면
    안됩니다.
    '''
    doi = adata[adata.obs['cell_select'] == '1']
    # print(doi.shape)
    if any(doi.obs[pseudotime_column] == float('inf')) or any(doi.obs[pseudotime_column].isna()):
        doi = doi[doi.obs[pseudotime_column] != float('inf')]
        doi = doi[~doi.obs[pseudotime_column].isna()]
    return doi

def synch_raw_counts(adata):
    '''
    TENET을 위해 cell X gene matrix 파일을 만들기 위한 전처리 과정입니다. 
    cell X gene matrix를 일관되게추출할 수 있도록 raw count expression을 
    adata.layers['counts']에 저장하도록 합니다. adata.layers['counts']가
    이미 있을 경우 adata.layers['counts'] 내에 있는 숫자들이 정수인지
    확인합니다. raw count expression은 정수 값이기 때문에 정수가 아니라면
    counts로 normalized data 대신 raw data를 넣어달라고 요구합니다.
    '''
    try:
        print(check_raw_data(adata.X), check_raw_data(adata.raw.X))
        # check_data_type(adata)
        if 'counts' in adata.layers and check_raw_data(adata.layers['counts']):
            print("Your adata object already contains raw counts")
            return
        if check_raw_data(adata.X):
            print("Your adata object contains raw counts adata.X", adata.X)
            adata.layers['counts'] = adata.X
            return
        if hasattr(adata, 'raw') and check_raw_data(adata.raw.X):
            print("Your adata object contains raw counts adata.raw.X", adata.raw.X)
            adata.layers['counts'] = adata.raw.X
            return
    except AttributeError as e:
        if check_raw_data(adata.layers['counts']):
            print("Your adata object contains raw counts adata.layers['counts']", adata.layers['counts'])
            return 
        else:
            print("Your adata object does not contain raw counts", e)
    except KeyError as e:
        print("Your adata object does not contain raw counts", e)

def tenet_input_make(adata, path_result, path_exp, path_pseudo, path_cell_select, pseudotime_column = 'dpt_pseudotime', gene_list = []):
    '''
    전처리가 완료된 adata를 매개변수를 받아서 3개의 tenet input file을 생성합니다.
    사용자가 만든 adata object에 있는 gene들과 사용자가 실제로 tenet 분석을 할 때
    관심이 있는 유전자 set이 다를 수 있어 사용자로 하여금 csv 혹은 txt파일의 형태로
    gene list를 따로 업로드 하게끔 할 수도 있을 것 같습니다.. 이 부분은 아직
    논의가 필요합니다.
    
    만약 업로드 한 gene
    list가 없을 경우 adata object에 있는 gene들로 cell X gene matrix를 만들어서
    분석을 진행하도록 합니다.
    '''
    try:
        if 'counts' not in adata.layers:
            print("adata.layers에 'counts' 키가 없습니다. 먼저 raw count 데이터를 설정해주세요.")
            return

        input_data = sc.get.obs_df(adata, keys=[pseudotime_column,'cell_select'])
        if not gene_list:
            gene_list = list(adata.var_names)

        if isinstance(adata.layers['counts'], csr_matrix):
            matrix_data = adata.layers['counts'][:, [gene_list.index(gene) for gene in gene_list]].todense()
        else:
            matrix_data = adata.layers['counts'][:, [gene_list.index(gene) for gene in gene_list]]
        
        
        #exp matrix
        pd.DataFrame(data=matrix_data,
                index=adata.obs_names, 
                columns=adata.var_names).to_csv(path_exp)
        
        #save pseudotime
        input_data[pseudotime_column].to_csv(path_pseudo, 
                                        index=False, header=False)
        
        #save cell sellect
        input_data['cell_select'].to_csv(path_cell_select, 
                                    index=False, header=False)
        with open(path_result, "w") as f:
            f.write("success")
    except:
        raise
    

options = json.load(open(snakemake.input[1]))

orig_data = sc.read_h5ad(snakemake.input[0])
adata = orig_data.copy()
adata.obs = organize_column_dtypes(adata.obs)

anno_of_interest = options['anno_of_interest']
pseudo_of_interest = options['pseudo_of_interest']
clusters_of_interest = options['clusters_of_interest']

# cell_select column 생성 및 '1' 혹은 '0'의 값 할당
make_cell_select(adata, annotation_column = anno_of_interest, categories = clusters_of_interest)

# cell_select column을 기반으로 cell 선택 및
# pseudotime column의 inf와 NA값 제거
adata = select_cells(adata, pseudotime_column = pseudo_of_interest)
    
# 뽑힌 cell이 없을시 프로그램 종료
if adata.shape[0] == 0:
    print('no cell selected')

# cell X gene matrix를 만들기 위해 raw counts가 
# adata.layers['counts']에 존재하도록 데이터 구조 동기화
synch_raw_counts(adata)
    
# tenet inputfile 생성
tenet_input_make(adata, snakemake.output[0], snakemake.output[1], snakemake.output[2], snakemake.output[3], pseudotime_column = pseudo_of_interest)