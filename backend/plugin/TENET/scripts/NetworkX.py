import networkx as nx
import matplotlib.pyplot as plt
import math
import mpld3
import plotly.graph_objs as go
import json
import sys

def load_graph(file_path):
    G = nx.DiGraph()
    source_nodes = set()
    target_nodes = set()

    try:
        with open(file_path, "r") as file:
            for line in file:
                parts = line.strip().split("\t")
                if len(parts) == 3:
                    source, weight, target = parts
                    G.add_edge(source, target, weight=round(float(weight), 4))
                    source_nodes.add(source)
                    target_nodes.add(target)
    except FileNotFoundError:
        print(f"File {file_path} not found.")
        return None, None, None
    
    return G, source_nodes, target_nodes

def plot_graph_plotly(G, source_nodes, target_nodes, parameter_node_size=30):
    pos = nx.shell_layout(G)
    edge_weights = [data['weight'] for (source, target, data) in G.edges(data=True)]
    node_sizes = [max(len(list(G.neighbors(node))), 1) * parameter_node_size for node in G.nodes]

    edge_trace = []
    for edge in G.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        weight = edge[2]['weight']
        edge_trace.append(go.Scatter(
            x=[x0, x1, None],
            y=[y0, y1, None],
            line=dict(width=weight * 2, color='grey'),
            hoverinfo='none',
            mode='lines'
        ))

    node_x = []
    node_y = []
    node_text = []
    node_color = []

    for node in G.nodes:
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(node)
        if node in source_nodes:
            node_color.append('yellow')
        else:
            node_color.append('skyblue')

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        text=node_text,
        mode='markers+text',
        textposition='top center',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='Viridis',
            color=node_color,
            size=node_sizes,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2
        )
    )

    fig = go.Figure(data=edge_trace + [node_trace],
                    layout=go.Layout(
                        title='Network Graph',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        annotations=[dict(
                            text="Plotly Network Graph",
                            showarrow=False,
                            xref="paper", yref="paper"
                        )],
                        xaxis=dict(showgrid=False, zeroline=False),
                        yaxis=dict(showgrid=False, zeroline=False)
                    ))

    return fig

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

G, source_nodes, target_nodes = load_graph(input_file_path)

if G is not None:
    # Plotly 그래프 생성 및 JSON 파일로 저장
    fig = plot_graph_plotly(G, source_nodes, target_nodes)
    fig_json = fig.to_json()
    with open(output_file_path, "w") as json_file:
        json_file.write(fig_json)

