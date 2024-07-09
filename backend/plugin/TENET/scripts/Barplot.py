import sys
import pandas as pd
import plotly.express as px
import json

# 명령줄 인수에서 파일 경로 받기
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

# 데이터 로드
gene_list = pd.read_csv(input_file_path, sep='\t', header=None, names=['Source', 'Weight', 'Target'])
gene_list = pd.DataFrame(gene_list)

# 상위 10개 Source 계산
source_counts = gene_list['Source'].value_counts()
top_sources = source_counts.head(10).index

# Plotly로 시각화 및 JSON 저장
gene_list_top_sources = gene_list[gene_list['Source'].isin(top_sources)]

# Source별로 빈도를 계산하여 새로운 데이터프레임 생성
source_freq = gene_list_top_sources['Source'].value_counts().reset_index()
source_freq.columns = ['Source', 'Frequency']

# 내림차순으로 정렬
source_freq = source_freq.sort_values(by='Frequency', ascending=False)

# 사용자 정의 색상 팔레트
color_palette = px.colors.qualitative.Set3

fig = px.bar(
    source_freq,
    x='Frequency',
    y='Source',
    color='Source',
    orientation='h',
    title='Top 10 Source Frequency Barplot',
    color_discrete_sequence=color_palette
)

fig.update_layout(yaxis={'categoryorder':'total ascending'})

fig_json = fig.to_json()
with open(output_file_path, "w") as json_file:
    json_file.write(fig_json)
