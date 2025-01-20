import sys
import pandas as pd
import plotly.express as px
import json

# 명령줄 인수에서 파일 경로 받기
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]
top_n = int(sys.argv[3])

# 데이터 로드 및 형식 확인
df = pd.read_csv(input_file_path, sep=None, engine='python', header=None)

if len(df.columns) == 3:  # fdrTrim.sif 형식 (3열)
    df.columns = ['Source', 'Weight', 'Target']
    # Source별로 Target 개수 계산
    source_freq = df['Source'].value_counts().reset_index()
    source_freq.columns = ['Genes', 'Outbound']
elif len(df.columns) == 2:  # numlinksOutdegree.txt 형식 (2열)
    df.columns = ['Genes', 'Outbound']
    source_freq = df
else:
    raise ValueError("입력 파일은 2열 또는 3열 형식이어야 합니다.")

# 상위 N개 선택
source_freq = source_freq.head(top_n)

# 내림차순으로 정렬
source_freq = source_freq.sort_values(by='Outbound', ascending=False)

# 사용자 정의 색상 팔레트
color_palette = px.colors.qualitative.Set3

# Plotly로 시각화
fig = px.bar(
    source_freq,
    x='Outbound',
    y='Genes',
    color='Genes',
    orientation='h',
    title=f"Top {top_n} Genes Outbound Barplot",
    color_discrete_sequence=color_palette
)

fig.update_layout(yaxis={'categoryorder':'total ascending'})

# JSON 저장
fig_json = fig.to_json()
with open(output_file_path, "w") as json_file:
    json_file.write(fig_json)
