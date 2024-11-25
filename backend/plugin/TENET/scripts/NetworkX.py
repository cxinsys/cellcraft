import networkx as nx
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import random
from plotly.colors import qualitative
import sys

def load_graph(file_path):
    G = nx.DiGraph()
    source_nodes = set()
    target_nodes = set()

    try:
        with open(file_path, "r") as file:
            lines = file.readlines()
            # 첫 번째 줄로 파일 형식 확인
            first_line = lines[0].strip().split("\t")
            
            if len(first_line) == 3:  # 3열 형식 (source, weight, target)
                for line in lines:
                    parts = line.strip().split("\t")
                    source, weight, target = parts
                    G.add_edge(source, target, weight=round(float(weight), 4))
                    source_nodes.add(source)
                    target_nodes.add(target)
            
            elif len(first_line) == 2:  # 2열 형식 (gene, weight)
                for line in lines:
                    parts = line.strip().split("\t")
                    gene, weight = parts
                    # 노드 추가 (self-loop로 처리)
                    G.add_node(gene, weight=round(float(weight), 4))
                    source_nodes.add(gene)
            else:
                raise ValueError("입력 파일은 2열 또는 3열 형식이어야 합니다.")
                
    except FileNotFoundError:
        print(f"File {file_path} not found.")
        return None, None, None
    except ValueError as e:
        print(str(e))
        return None, None, None
    
    return G, source_nodes, target_nodes

def sample_graph_by_degree(G, top_n=10):
    degree_dict = dict(G.degree())
    if not degree_dict:  # 그래프가 비어있는 경우
        return G
    # 실제 그래프 크기와 요청된 top_n 중 작은 값 사용
    actual_n = min(top_n, len(degree_dict))
    top_nodes = sorted(degree_dict, key=degree_dict.get, reverse=True)[:actual_n]
    subgraph = G.subgraph(top_nodes).copy()
    return subgraph

def calculate_node_size(G, method="degree", scale=20, min_size=10):
    if len(G.nodes()) == 0:  # 빈 그래프 체크
        return []
    if method == "degree":
        centrality = dict(G.degree())
    elif method == "betweenness":
        centrality = nx.betweenness_centrality(G)
    elif method == "pagerank":
        centrality = nx.pagerank(G)
    else:
        raise ValueError(f"Unknown centrality method: {method}")

    max_value = max(centrality.values()) or 1  # Avoid division by zero
    min_value = min(centrality.values())

    # Normalize node sizes between min_size and scale
    return [
        ((value - min_value) / (max_value - min_value)) * (scale - min_size) + min_size
        for value in centrality.values()
    ]

def plot_graph_plotly(G, source_nodes, target_nodes, node_sizes):
    pos = nx.spring_layout(G, seed=42)

    edge_colors = random.choices(qualitative.Plotly, k=len(G.edges))
    node_colors = random.choices(qualitative.Set3, k=len(G.nodes))

    edge_trace = []
    for i, edge in enumerate(G.edges(data=True)):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        weight = edge[2].get('weight', 1)
        edge_trace.append(go.Scatter(
            x=[x0, x1, None],
            y=[y0, y1, None],
            line=dict(width=weight, color=edge_colors[i]),
            hoverinfo='none',
            mode='lines'
        ))

    node_x = []
    node_y = []
    node_text = []
    for node in G.nodes:
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(node)

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        text=node_text,
        mode='markers+text',
        textposition='top center',
        hoverinfo='text',
        marker=dict(
            showscale=False,
            colorscale='Viridis',
            color=node_colors,
            size=node_sizes,
            line_width=1
        )
    )

    fig = go.Figure(data=edge_trace + [node_trace],
                    layout=go.Layout(
                        title='Gene Regulatory Network',
                        titlefont=dict(size=16),
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        annotations=[dict(
                            text="",
                            showarrow=False,
                            xref="paper", yref="paper"
                        )],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                    ))

    return fig

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]
top_n = int(sys.argv[3])
centrality_method = sys.argv[4] if len(sys.argv) > 4 else "degree"

G, source_nodes, target_nodes = load_graph(input_file_path)

if G is not None:
    sampled_G = sample_graph_by_degree(G, top_n=top_n)
    node_sizes = calculate_node_size(sampled_G, method=centrality_method, scale=80, min_size=10)

    fig = plot_graph_plotly(sampled_G, source_nodes, target_nodes, node_sizes)
    fig_json = fig.to_json()
    with open(output_file_path, "w") as json_file:
        json_file.write(fig_json)
