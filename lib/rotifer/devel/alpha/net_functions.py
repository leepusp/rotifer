from operator import pos
from dash import Dash, dcc, html, Input, Output, ctx, State
import dash_cytoscape as cyto
import networkx as nx
import numpy as np
import pandas as pd
import warnings
from rotifer.devel.alpha import gian_func as gf

# Enable SVG export
cyto.load_extra_layouts()

def network_annotation(networkdf,
                     node_color=None,
                     node_shape='ellipse',
                     second_shape=None,
                     second_shape_list=None,
                     node_outline=None,
                     node_outline_color=None,
                     node_outline_list=None,
                     edge_color=True,
                     unidirected_graph=True,
                     split_types=False):
    """
    networkdf: get_relations output (source, target, edge_type, edge_count, edge_total).
               Without edge_total, it is estimated by summing edge_count over the pair's
               edge types, which counts a block once per type it occurs in.
    unidirected_graph: if False, a MultiDiGraph with one edge per (source, target, edge_type).
    split_types: undirected only. If True, a MultiGraph with one edge per edge type
                 (keyed by edge_type); if False, a Graph with one edge per pair carrying
                 its most frequent type. In both, mirrored relations (A->B Cterminal and
                 B->A Nterminal) are merged, and edge_type is read from the node in the
                 edge's type_from attribute (the alphabetically first of the pair) to the
                 other. edge_total is the pair's block count either way.
    """

    import networkx as nx
    import pandas as pd

    second_shape = {} if second_shape is None else second_shape
    second_shape_list = [] if second_shape_list is None else second_shape_list
    node_outline_color = {} if node_outline_color is None else node_outline_color
    node_outline_list = [] if node_outline_list is None else node_outline_list

    networkdf = networkdf.copy()
    if 'edge_total' not in networkdf.columns:
        warnings.warn("networkdf has no edge_total column; summing edge_count over edge types, "
                      "which counts a block once per type it occurs in")
        networkdf['edge_total'] = networkdf.groupby(['source', 'target'])['edge_count'].transform('sum')

    edge_attr = ['edge_type', 'edge_count', 'edge_total']
    if unidirected_graph:
        # One row per undirected relation: orient every row source <= target, flipping its type
        flip = {'Cterminal': 'Nterminal', 'Nterminal': 'Cterminal', 'Downstream': 'Upstream', 'Upstream': 'Downstream'}
        swap = networkdf.source > networkdf.target
        networkdf.loc[swap, ['source', 'target']] = networkdf.loc[swap, ['target', 'source']].values
        networkdf.loc[swap, 'edge_type'] = networkdf.loc[swap, 'edge_type'].replace(flip)
        # A self relation is its own mirror
        self_rel = networkdf.source == networkdf.target
        networkdf.loc[self_rel, 'edge_type'] = networkdf.loc[self_rel, 'edge_type'].replace({'Nterminal': 'Cterminal', 'Upstream': 'Downstream'})
        networkdf = networkdf.groupby(['source', 'target', 'edge_type'], as_index=False)[['edge_count', 'edge_total']].max()
        # An undirected edge does not keep its source: record which end edge_type refers to
        networkdf['type_from'] = networkdf.source
        edge_attr = edge_attr + ['type_from']
        if split_types:
            G = nx.from_pandas_edgelist(networkdf, edge_attr=edge_attr, edge_key='edge_type', create_using=nx.MultiGraph())
        else:
            networkdf = networkdf.sort_values(['edge_count', 'edge_type'], ascending=[False, True]).drop_duplicates(['source', 'target'])
            G = nx.from_pandas_edgelist(networkdf, edge_attr=edge_attr)
    else:
        G = nx.from_pandas_edgelist(networkdf, edge_attr=edge_attr, create_using=nx.MultiDiGraph())


    default_color = {item: "#D3D3D3" for item in G.nodes()}
    if node_color:
        missing = [node for node in G.nodes() if node not in node_color]
        if missing:
            warnings.warn(f"nodes missing from node_color, drawn grey: {sorted(missing)}")
        node_color = {**default_color, **node_color}
    else:
        node_color = default_color
    nx.set_node_attributes(G, node_color, name="node_color") 

    #Getting the size of nodes
    # In the future we can add more ways of getting the node size
    node_degree = dict(_simple_graph(G).degree())
    nx.set_node_attributes(G, node_degree, name="degree") 
    

    #Getting the shape of nodes
    #Creating a dictionary to shape nodes
    if second_shape == {}:
        second_shape={item: "hexagon" for item in second_shape_list}
    elif isinstance(second_shape, str):
        second_shape={item: second_shape for item in second_shape_list}
    else:
        pass

    missing = [node for node in G.nodes() if node in second_shape_list and node not in second_shape]
    if missing:
        warnings.warn(f"nodes missing from second_shape, drawn as {node_shape}: {sorted(missing)}")
    node_shape_dict  = {node: second_shape.get(node, node_shape) if node in second_shape_list else node_shape for node in G.nodes()}
    nx.set_node_attributes(G, node_shape_dict, name="node_shape")


    # Getting the outline of nodes:
    if node_outline_color == {}:
        node_outline_color={item: "#000000" for item in node_outline_list}
    elif isinstance(node_outline_color, str):
        node_outline_color={item: node_outline_color for item in node_outline_list}
    else:
        pass
    if node_outline:
        missing = [node for node in G.nodes() if node in node_outline_list and node not in node_outline_color]
        if missing:
            warnings.warn(f"nodes missing from node_outline_color, outlined in their node colour: {sorted(missing)}")
        node_outline_dict = {node: node_outline_color.get(node, node_color[node]) if node in node_outline_list else node_color[node] for node in G.nodes()}
    else:
        # The viewer always draws a border from node_outline: blend it into the node
        node_outline_dict = node_color
    nx.set_node_attributes(G, node_outline_dict, name="node_outline") 
    
    if edge_color:
        edge_colors=[]
        edge_alpha=[]
        edge_styles=[]
    #Getting the color of edges():
        if G.is_multigraph():
            edges_obj = G.edges(keys=True)
            for source,target, key in edges_obj:
                if node_color[source] == node_color[target]:
                    edge_colors.append(node_color[source])
                    edge_styles.append('solid')
                    edge_alpha.append(0.7)
                else:
                    edge_colors.append('gray')
                    edge_styles.append('dotted')
                    edge_alpha.append(0.2)    
                   
            dict_edge_colors = { (u, v, key): color for (u, v, key), color in zip(edges_obj, edge_colors)}
            dict_edge_style = { (u, v, key): style for (u, v, key), style in zip(edges_obj, edge_styles)}
            dict_edge_alpha = { (u, v, key): alpha for (u, v, key), alpha in zip(edges_obj, edge_alpha)}

        else:    
            edges_obj = G.edges() 
            for source,target in edges_obj:
                if node_color[source] == node_color[target]:
                    edge_colors.append(node_color[source])
                    edge_styles.append('solid')
                    edge_alpha.append(0.7)
                else:
                    edge_colors.append('gray')
                    edge_styles.append('dotted')
                    edge_alpha.append(0.2)    
            dict_edge_colors = {k: v for k, v in zip(G.edges(), edge_colors)}
            dict_edge_style = {k: v for k, v in zip(G.edges(), edge_styles)}
            dict_edge_alpha = {k: v for k, v in zip(G.edges(), edge_alpha)}

        nx.set_edge_attributes(G, dict_edge_style, 'line_style')
        nx.set_edge_attributes(G, dict_edge_colors, 'edge_color')
        nx.set_edge_attributes(G, dict_edge_alpha, 'opacity')        

    return G

def _simple_graph(G):
    # Parallel edges of an undirected MultiGraph (split_types) share one pair: collapse them
    # so degree, layouts and communities match the one-edge-per-pair graph
    return nx.Graph(G) if G.is_multigraph() and not G.is_directed() else G

def calculate_positions(G, iterations):
    # Recalculate positions as needed
    G = _simple_graph(G)
    pos1 = nx.spring_layout(G, iterations=iterations, seed=1)
    pos2 = nx.kamada_kawai_layout(G, pos=pos1)
    pos3 = nx.kamada_kawai_layout(G)
    pos4 = nx.circular_layout(G)
    return {"Spring Layout": pos1,
            "Kamada-Kawai from Spring layout": pos2,
            "Kamada-Kawai": pos3,
            "Circular": pos4}

def convert_to_cytoscape(G, positions, layout_key, scaling_factor=900, x_scale=1, y_scale=1):
    cy = nx.cytoscape_data(G)
    # Modify node positions and style
    for node in cy["elements"]["nodes"]:
        node_id = node["data"]["id"]
        if node_id in positions[layout_key]:
            p = positions[layout_key][node_id]
            node["position"] = {
                "x": int(p[1] * scaling_factor * x_scale),
                "y": int(p[0] * scaling_factor * y_scale),
            }
    return cy["elements"]["nodes"] + cy["elements"]["edges"]


def network_dash(G,pos_dict, xx, function, color_schemes, **run_kwargs):
    """
    Extra keyword arguments (e.g. host, port, debug, jupyter_mode) are passed to app.run;
    defaults: jupyter_mode="external", debug=True.
    """
    from dash import Dash, dcc, html, Input, Output, ctx, State
    import dash
    import dash_cytoscape as cyto
    import networkx as nx
    from networkx.readwrite import json_graph
    import numpy as np
    import pandas as pd
    from rotifer.devel.alpha import net_functions
    from dash.exceptions import PreventUpdate


    # Serializing G to save it od dash  as dcc.store, that way we could manipulate in any way that we need
    G_serialized = json_graph.node_link_data(G.copy())

    app = dash.Dash(__name__)
    default_stylesheet = [
        {
            'selector': 'node',
            'style': {
                'content': 'data(name)',
                'background-color': 'data(node_color)',
                'border-color': 'data(node_outline)',
                "border-width": 3,
                'shape': 'data(node_shape)',
                'width': 'data(size)',
                'height': 'data(size)',
                'text-valign': 'bottom',
                'text-halign': 'center',
                'font-size': '12px', 
                'color': 'black'
            }
        },
        {
            'selector': 'edge',
            'style': {
                'line-color': 'data(edge_color)',
                'line-style': 'data(line_style)',
                'opacity': 'data(opacity)',
                'curve-style': 'bezier'  # draw parallel edges apart
            }
        }
    ]
    
    
    app.layout = html.Div([
        dcc.Dropdown(
            id='layout-selector',
            options=[{'label': key, 'value': key} for key, value in pos_dict.items()],
            value=list(pos_dict.keys())[0],
            clearable=False,
            style={'width': '30%', 'margin-right': '10p'}
        ),
        html.Label("Spring Layout iterations:", style={'margin-right': '5px'}),
        dcc.Input(
                id='iterations-input',
                type='number',
                placeholder='Iterations',
                value=10,
                min=1,
                step=1,
                style={'width': '30px'}
            ),
        html.Label("Minimmun connection:", style={'margin-right': '5px'}),
        dcc.Input(
                id='connections-input',
                type='number',
                placeholder='min_connections',
                value=0,
                min=0,
                step=1,
                style={'width': '30px'}
            ),
        html.Label("x scale:", style={'margin-right': '5px'}),
        dcc.Input(
                id='x_scale',
                type='number',
                placeholder='x_scale',
                value=1,
                min=1,
                step=0.1,
                style={'width': '30px'}
            ),
        html.Label("y scale:", style={'margin-right': '5px'}),
        dcc.Input(
                id='y_scale',
                type='number',
                placeholder='y scale',
                value=1,
                min=1,
                step=0.1,
                style={'width': '30px'}
            ),
        html.Button("Update layouts", id="up_lay"),
        dcc.Input(id='new-group-key', type='text', placeholder='Group name'),
        dcc.Input(id='new-group-color', type='text', placeholder='Color hex'),
        html.Button('Add group', id='add-group-button', n_clicks=0),
    html.Div([
        html.Div([
            # Source option with radio and input
            html.Div([
                dcc.RadioItems(
                    id="select-mode",
                    options=[
                        {"label": "Source", "value": "source"},
                        {"label": "Target", "value": "target"}
                    ],
                    value="source",
                    labelStyle={'display': 'inline-block', 'margin-right': '15px'}
                )
            ], style={"margin-bottom": "10px"}),
            html.Div([
                html.Label("Source Node", htmlFor="select_source", style={"margin-right": "5px"}),
                dcc.Input(id="select_source", type="text", placeholder="Source Node", style={"width": "150px"})
            ], style={"display": "inline-block", "margin-right": "30px"}),
            html.Div([
                html.Label("Target Node", htmlFor="select_target", style={"margin-right": "5px"}),
                dcc.Input(
                    id="select_target",
                    type="text",
                    placeholder="Target Node",
                    style={"width": "150px", "margin-right": "10px"}
                    ),
                html.Button("Add Edge", id="add-edge-button", n_clicks=0, style={"margin-right": "5px"}),
                html.Button("Remove Edge", id="remove-edge-button", n_clicks=0),
                html.Button("Update Communities", id="update-communities-btn", n_clicks=0)
                ], style={"display": "inline-block"}),
            ### Feedback from Add/Remove Edge
            html.Div(id="edge-edit-status", style={"display": "inline-block", "margin-left": "10px"})
        ])
    ]),

    html.Div(id="node-input-box"),

         html.Div([
             dcc.Dropdown(
                 id="colorscheme",
                 options=[],   # preenchido no callback
                 value=None,   # será definido no callback também
                 placeholder="Select color scheme",
                 style={"flex": "0 0 10%", "margin-right": "10px"}
                 ),
             dcc.Dropdown(
                 id='group_selection',
                 options=[],
                 value=[],
                 multi=True,
                 placeholder="Select groups",
                 style={"flex": "1"}
                 )
             ], style={"display": "flex", "align-items": "center", "margin-bottom": "10px"}),

        cyto.Cytoscape(
            id='cytoscape',
            layout={'name': 'preset'},
            style={'height': '85vh', 'width': '100%'},
            stylesheet=default_stylesheet,
            elements=[]  # Initially empty
        ),
        html.Label("Font Size:"),
        dcc.Input(
            id='font-size-input',
            type='number',
            value=12,  # Default font size
            min=0,
            max=24,
            step=0.5,
            style={'width': '30px'}
        ),
        html.Label("Edge width:"),
        dcc.Input(
            id='edge-size',
            type='number',
            value=1.5,  # Default font size
            min=0,
            max=10,
            step=0.1,
            style={'width': '30px'}
        ),
        html.Label("Outline width:"),
        dcc.Input(
            id='outline-size-input',
            type='number',
            value=1,  # Default font size
            min=0,
            max=20,
            step=0.1,
            style={'width': '30px'}
        ),
        html.Label("Node size:", style={'margin-right': '5px'}),
        dcc.Input(
            id='Node_scale',
            type='number',
            value=3,  # Default font size
            min=0,
            max=50,
            step=0.1,
            style={'width': '30px'}
        ),
        html.Label("Edge distances:",style={'margin-right': '5px'}),
        dcc.Input(
            id='scale-size-input',
            type='range',
            value=900,  # Default font size
            min=100,
            max=3000,
            step=100,
            style={'width': '90px'}
        ),
        dcc.Input(
            id='Saved_position_text',
            value='Pos_name',  # Default font size
            type='text',
            style={'width': '70px'}
        ),
        html.Button("Save Position", id="save-button", n_clicks=0),
        html.Button("Print Stored Data", id="print-button"),  # Button to trigger action
        html.Label("Save network as:"),
        html.Button("JPG", id="btn-get-jpg" ),
        html.Button("PNG", id="btn-get-png"),
        html.Button("SVG", id="btn-get-svg"),
        #### Store  to do the double clik
        dcc.Store(id="edge-color-cycle", data=[]),
        dcc.Store(id="edge-style-overrides", data={}),
        dcc.Store(id="last-edge-click", data={"id": None, "time": 0}),
        dcc.Store(id="last-node-click", data={"id": None, "time": 0}),
        dcc.Store(id="node-color-cycle", data=[]),
        ### Adding stroe for groups
        dcc.Store(id="xx-store", data=xx),
        dcc.Store(id="function-store", data=function),
        dcc.Store(id="color-schemes-store", data=color_schemes),

        ### Saving position
        dcc.Store(id="pos-dict-store", data=pos_dict),
        #### Saving the networkX G object 
        dcc.Store(id="graph-store", data=G_serialized),

        ### Dummy store to update add and remove edges on fly
        dcc.Store(id="edge-modification-trigger", data=0),




        dcc.Store(id="saved-pos-store", data={}, clear_data=True, storage_type='memory'),
        html.Div(id="output"),  # Output for display
        
    ])
   


    @app.callback(
        Output('cytoscape', 'stylesheet'),
        [Input('font-size-input', 'value'),  # Font size input
         Input('outline-size-input', 'value'),
         Input('edge-size', 'value')]       # Edge size input
    )
    def update_stylesheet(font_size,outline_size, edge_size):
        return [
            {
                'selector': 'node',
                'style': {
                    'content': 'data(name)',
                    'background-color': 'data(node_color)',
                    'border-color': 'data(node_outline)',
                    'shape': 'data(node_shape)',
                    "border-width": f'{outline_size}px',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'text-valign': 'bottom',
                    'text-halign': 'center',
                    'font-size': f'{font_size}px',  # Dynamic font size
                    'color': 'black'
                }
            },
            {
                'selector': 'edge',
                'style': {
                    'line-color': 'data(edge_color)',
                    'line-style': 'data(line_style)',
                    'opacity': 'data(opacity)',
                    'curve-style': 'bezier',
                    'width': f'{edge_size}px'  # Dynamic edge width
                }
            }
        ]
    """    
    @app.callback(
        Output('cytoscape', 'elements'),
        [Input('layout-selector', 'value'),
        Input('group_selection', 'value'),
        Input('connections-input', 'value'),
        Input('Node_scale', 'value'),
        Input('x_scale', 'value'),
        Input('y_scale', 'value'),
        Input('scale-size-input', 'value')]
    )
    def update_graph(selected_option, groups,min_connections,node_scale,x_scale,y_scale, scale_value):
        scale_value = int(scale_value)
       # y_scale = int(scale_value)
       # x_scale = int(scale_value)
        # Dynamically generate elements based on selected_option
        #### Node filter based on groups:
        filter_nodes = list({key: value for key, value in function.items() if value in groups}.keys())
        G_filtered = G.subgraph(filter_nodes)
        ### Edge filtered based in total Edge count
        filter_edges = [(u, v) for u, v, attrs in G_filtered.edges(data=True) if attrs.get('edge_total', 0) >= min_connections]
        subgraph = G_filtered.edge_subgraph(filter_edges).copy()
        subgraph.remove_nodes_from(list(nx.isolates(subgraph)))
        for edge in subgraph.edges(data=True):
            if function[edge[0]] == function[edge[1]]:
                edge[2]['line_style']='solid'
                edge[2]['opacity']=1
                edge[2]['edge_color']= xx[function[edge[0]]]
            else:
                edge[2]['line_style']='dotted'
                edge[2]['opacity']=0.3
                edge[2]['edge_color']= 'gray'

        ## Scale the node size:
        for node, attrs in subgraph.nodes(data=True):
            attrs['size'] = np.sqrt(attrs['degree']) * node_scale
            attrs['node_color'] = xx[function[node]]
            
        element = net_functions.convert_to_cytoscape(subgraph, pos_dict, selected_option,x_scale=x_scale, y_scale=y_scale, scaling_factor=scale_value)
        return element
    """ 

    @app.callback(
        Output("select_source", "value"),
        Output("select_target", "value"),
        Input("cytoscape", "tapNodeData"),
        State("select-mode", "value"),
    )
    def update_input_boxes(node_data,toggle):
        if node_data and "id" in node_data:
            node_id = node_data["id"]
            if  toggle == "source":
                return node_id, dash.no_update
            elif toggle == "target":
                return dash.no_update, node_id
        return dash.no_update, dash.no_update


    @app.callback(
        Output("cytoscape", "elements"),
        Output("last-edge-click", "data"),
        Output("edge-color-cycle", "data"),
        Output("edge-style-overrides", "data"),
        Output("last-node-click", "data"),
        Output("node-color-cycle", "data"),
        Output("function-store", "data"),
        Input("cytoscape", "tapEdgeData"),
        Input("cytoscape", "tapNodeData"),
        Input("layout-selector", "value"),
        Input("group_selection", "value"),
        Input("connections-input", "value"),
        Input("Node_scale", "value"),
        Input("x_scale", "value"),
        Input("y_scale", "value"),
        Input("scale-size-input", "value"),
        Input("pos-dict-store", "data"),
        Input("edge-modification-trigger", "data"),
        State("cytoscape", "elements"),
        State("last-edge-click", "data"),
        State("edge-color-cycle", "data"),
        State("edge-style-overrides", "data"),
        State("last-node-click", "data"),
        State("node-color-cycle", "data"),
        State('xx-store', 'data'),
        State("function-store", "data"),
        State("pos-dict-store", "data"),
        State("graph-store", "data"),
    )
    def update_elements_or_cycle(
        edge_data, node_data,
        selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value,pos_dict_input,dummy_call,
        elements, last_click, color_cycle, style_overrides, last_node_click, node_color_cycle, func_color_dict, function_data, pos_dict_data, serial_graph
    ):

        from dash import Input, Output, State, ctx, no_update
        import time
        import numpy as np
        triggered_id = ctx.triggered_id

        # Edge style overrides are kept per session, tagged with the colour scheme they were made under:
        # a scheme change discards them, together with the style cycle built under that scheme
        style_overrides = style_overrides or {}
        if any(func_color_dict.get(g) != c for g, c in style_overrides.get("xx", {}).items()):
            style_overrides, color_cycle = {}, []
        style_overrides = {"xx": dict(func_color_dict), "overrides": style_overrides.get("overrides", {})}

        # Duplo clique em aresta
        if triggered_id == "cytoscape" and ctx.triggered[0]["prop_id"] == "cytoscape.tapEdgeData":
            if not edge_data or "id" not in edge_data:
                return no_update, last_click, color_cycle, style_overrides, last_node_click, node_color_cycle, function_data

            edge_id = f"{edge_data['source']}|{edge_data['target']}"
            now = time.time()
            double_click = (last_click["id"] == edge_id and (now - last_click["time"]) < 0.3)

            if not double_click:
                return no_update, {"id": edge_id, "time": now}, color_cycle, style_overrides, last_node_click, node_color_cycle, function_data

            # Captura estilos únicos de arestas visíveis
            if not color_cycle:
                seen = set()
                color_cycle = []
                for el in elements:
                    if "source" not in el.get("data", {}):
                        continue
                    color = el["data"].get("edge_color", "#ccc")
                    style = el["data"].get("line_style", "solid")
                    opacity = round(float(el["data"].get("opacity", 1)), 2)
                    key = (color, style, opacity)
                    if key not in seen:
                        seen.add(key)
                        color_cycle.append(key)

            # Aplica estilo seguinte com contador e salva override
            new_elements = []
            for el in elements:
                data = el.get("data", {})
                if f"{data.get('source')}|{data.get('target')}" == edge_id:
                    cycle_index = (data.get("cycle_index", -1) + 1) % len(color_cycle)
                    data["cycle_index"] = cycle_index
                    next_style = color_cycle[cycle_index]
                    data["edge_color"], data["line_style"], data["opacity"] = next_style

                    el["style"] = el.get("style", {})
                    el["style"].update({
                        "line-color": next_style[0],
                        "target-arrow-color": next_style[0],
                        "line-style": next_style[1],
                        "opacity": next_style[2],
                    })

                    style_overrides["overrides"][edge_id] = list(next_style)

                new_elements.append(el)

            return new_elements, {"id": edge_id, "time": now}, color_cycle, style_overrides, last_node_click, node_color_cycle, function_data

        # Duplo clique em nó
        if triggered_id == "cytoscape" and ctx.triggered[0]["prop_id"] == "cytoscape.tapNodeData":
            if not node_data or "id" not in node_data:
                return no_update, last_click, color_cycle, style_overrides, last_node_click, node_color_cycle, function_data

            node_id = node_data["id"]
            now = time.time()
            double_click = (last_node_click["id"] == node_id and (now - last_node_click["time"]) < 0.3)

            if not double_click:
                return no_update, last_click, color_cycle, style_overrides, {"id": node_id, "time": now}, node_color_cycle, function_data
            """
            if not node_color_cycle:
                node_color_cycle = list(func_color_dict.values())
            """
            node_color_cycle = list(func_color_dict.values())

            new_elements = []
            for el in elements:
                data = el.get("data", {})
                if data.get("id") == node_id:
                    current_color = data.get("node_color", "#ccc")
                    try:
                        next_color = node_color_cycle[(node_color_cycle.index(current_color) + 1) % len(node_color_cycle)]
                    except ValueError:
                        next_color = node_color_cycle[0]

                    data["node_color"] = next_color
                    el["style"] = el.get("style", {})
                    el["style"]["background-color"] = next_color
                    """
                    for func_name, color in func_color_dict.items():
                        if color == next_color:
                            function_data[node_id] = func_name
                            break
                    """        
                    # Procura o nome da função correspondente à cor
                    matched_group = next((func_name for func_name, color in func_color_dict.items() if color == next_color), None)

                    if matched_group is not None:
                        function_data[node_id] = matched_group
                    else:
                        # Cor não está no dicionário original, tenta buscar no xx-store atualizado
                        for func_name, color in xx.items():
                            if color == next_color:
                                function_data[node_id] = func_name
                                break

                new_elements.append(el)

            return new_elements, last_click, color_cycle, style_overrides, {"id": node_id, "time": now}, node_color_cycle, function_data

        # Atualização normal do grafo
        scale_value = int(scale_value)
        filter_nodes = list({key: value for key, value in function_data.items() if value in groups}.keys())
        stored_G = json_graph.node_link_graph(serial_graph)
        G_filtered = stored_G.subgraph(filter_nodes)

        # Multigraph edges need their key; manually added edges are always kept
        edges = G_filtered.edges(keys=True, data=True) if G_filtered.is_multigraph() else G_filtered.edges(data=True)
        filter_edges = [e[:-1] for e in edges
                        if e[-1].get('edge_total', 0) >= min_connections or e[-1].get('edge_type') == 'manual']
        subgraph = G_filtered.edge_subgraph(filter_edges).copy()
        subgraph.remove_nodes_from(list(nx.isolates(subgraph)))

        for edge in subgraph.edges(data=True):
            uid = f"{edge[0]}|{edge[1]}"
            if uid in style_overrides["overrides"]:
                color, style, opacity = style_overrides["overrides"][uid]
                edge[2]['line_style'] = style
                edge[2]['opacity'] = opacity
                edge[2]['edge_color'] = color
            elif function_data[edge[0]] == function_data[edge[1]]:
                edge[2]['line_style'] = 'solid'
                edge[2]['opacity'] = 1
                edge[2]['edge_color'] = func_color_dict[function_data[edge[0]]]
            else:
                edge[2]['line_style'] = 'dotted'
                edge[2]['opacity'] = 0.3
                edge[2]['edge_color'] = 'gray'

        # Size from the degree in the displayed subgraph; colour always from the stored group
        for node, attrs in subgraph.nodes(data=True):
            attrs['size'] = np.sqrt(net_functions._simple_graph(subgraph).degree(node)) * node_scale
            attrs['node_color'] = func_color_dict[function_data[node]]

        new_elements = convert_to_cytoscape(
            subgraph, pos_dict_data, selected_option,
            x_scale=x_scale, y_scale=y_scale, scaling_factor=scale_value
        )
        return new_elements, last_click, color_cycle, style_overrides, last_node_click, node_color_cycle, function_data

    @app.callback(
        [Output("saved-pos-store", "data"),
        Output("layout-selector", "options"),
        Output("pos-dict-store","data"),
        Output("layout-selector", "value")],
        [Input("save-button", "n_clicks"),
        State("Saved_position_text", "value"),
        State('x_scale', 'value'),
        State('y_scale', 'value'),
        State('scale-size-input', 'value'),
        State("cytoscape", "elements"),
        State("pos-dict-store", "data"),],
        prevent_initial_call=True,
    )
    def save_positions(n_clicks, text,x_scale,y_scale, scale, n_elements, current_pos_dict):
        # Extract positions from the Cytoscape graph
        positions = {
            node["data"]["id"]: node.get("position", {})
            for node in n_elements if "position" in node
        }
        # A type='range' dcc.Input sends its value as a string once the slider is moved
        scale = int(scale)
        positions = {key: np.array([value['y']/scale/y_scale, value['x']/scale/x_scale]) for key, value in positions.items()}
        #global p
        current_pos_dict[text] = positions
        to_selector = [{'label': key, 'value': key} for key, value in current_pos_dict.items()]
        return positions, to_selector, current_pos_dict, text
    
    
    @app.callback(
        Output("output", "children"),
        Input("print-button", "n_clicks"),
        State("saved-pos-store", "data"),
        prevent_initial_call=True,  # Prevent execution on page load
    )
    def print_stored_data(n_clicks, toprint):
        
        print("Stored data:", toprint)  # Print to console
        return f"Printed data: {toprint}"  # Show feedback in the app
    
    
    @app.callback(
        Output("cytoscape", "generateImage"),
        [
            Input("btn-get-jpg", "n_clicks"),
            Input("btn-get-png", "n_clicks"),
            Input("btn-get-svg", "n_clicks"),
        ],
        prevent_initial_call=True
    )
    def get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks):
        triggered_id = ctx.triggered_id
        file_type_map = {
            'btn-get-jpg': 'jpg',
            'btn-get-png': 'png',
            'btn-get-svg': 'svg'
        }
        ftype = file_type_map.get(triggered_id, 'jpg')
        action = 'download' if ftype == 'svg' else 'both'
        return {'type': ftype, 'action': action}
  
    ##### Adding dinamically groups to dropdown menu:
    @app.callback(
        Output("group_selection", "options"),
        Output("group_selection", "value"),
        Output("function-store", "data", allow_duplicate=True),  # <-- adicionado
        Output("xx-store", "data", allow_duplicate=True),  # <-- adicionado
        Input("colorscheme", "value"),
        State("color-schemes-store", "data"),
        prevent_initial_call=True
    )
    def update_group_selection(selected_scheme, color_schemes_data):
        if not selected_scheme or selected_scheme not in color_schemes_data:
            return [], [], dash.no_update, dash.no_update

        function_dict, color_dict = color_schemes_data[selected_scheme]

        # Transforma o dicionário de funções para string-string (caso necessário)
        function_dict = {str(k): str(v) for k, v in function_dict.items()}

        options = [{"label": g, "value": g} for g in color_dict.keys()]
        return options, list(color_dict.keys()), function_dict, color_dict

    
    @app.callback(
        Output('xx-store', 'data'),
        Input('add-group-button', 'n_clicks'),
        State('new-group-key', 'value'),
        State('new-group-color', 'value'),
        State('xx-store', 'data'),
        prevent_initial_call=True
    )
    def add_group(n_clicks, new_key, new_color, current_xx):
        if not new_key or not new_color:
            return dash.no_update
        current_xx[new_key] = new_color
        return current_xx



    @app.callback(
        Output("layout-selector", "options", allow_duplicate=True),
        Output("pos-dict-store","data", allow_duplicate=True),
        Input("up_lay", "n_clicks"),
        State('group_selection', 'value'),
        State('iterations-input', 'value'),
        State('connections-input', 'value'),
        State("pos-dict-store", "data"),
        State("graph-store", "data"),
        State("function-store", "data"),
        prevent_initial_call=True,  # Prevent execution on page load
    )
    def update_layouts(n_clicks, groups,Spring_lay_iteractions, min_connections, pos_dict, serial_graph, function_data):
        stored_G = json_graph.node_link_graph(serial_graph)
        filter_nodes = list({key: value for key, value in function_data.items() if value in groups}.keys())
        G_filtered = stored_G.subgraph(filter_nodes)
        ### Edge filtered based in total Edge count (same filter as the graph render)
        edges = G_filtered.edges(keys=True, data=True) if G_filtered.is_multigraph() else G_filtered.edges(data=True)
        filter_edges = [e[:-1] for e in edges
                        if e[-1].get('edge_total', 0) >= min_connections or e[-1].get('edge_type') == 'manual']
        subgraph = G_filtered.edge_subgraph(filter_edges).copy()
        subgraph.remove_nodes_from(list(nx.isolates(subgraph)))
        newp = net_functions.calculate_positions(subgraph, Spring_lay_iteractions)
        pos_dict['updated Spring'] = newp['Spring Layout']
        pos_dict['updated Kamada from spring'] = newp['Kamada-Kawai from Spring layout']
        pos_dict['updated Kamada'] = newp['Kamada-Kawai']
        pos_dict['updated Circular'] = newp['Circular']
        to_selector = [{'label': key, 'value': key} for key, value in pos_dict.items()]
        return to_selector, pos_dict


    @app.callback(
        Output("graph-store", "data", allow_duplicate=True),
        Output("edge-modification-trigger", "data", allow_duplicate=True),
        Output("edge-edit-status", "children"),
        Input("add-edge-button", "n_clicks"),
        Input("remove-edge-button", "n_clicks"),
        State("select_source", "value"),
        State("select_target", "value"),
        State("graph-store", "data"),
        State("edge-modification-trigger", "data"),
        prevent_initial_call=True
    )
    def modify_edge(add_clicks, remove_clicks, source, target, graph_data, trigger_count):
        triggered_id = ctx.triggered_id
        if not source or not target:
            raise PreventUpdate

        G = json_graph.node_link_graph(graph_data)

        # The node boxes are free text: reject names that are not in the graph
        missing = [n for n in (source, target) if n not in G]
        if missing:
            return dash.no_update, dash.no_update, f"Unknown node(s): {', '.join(missing)}"

        if triggered_id == "add-edge-button":
            # Do not overwrite the attributes (edge_total etc.) of an existing edge
            if G.has_edge(source, target):
                return dash.no_update, dash.no_update, f"Edge {source} - {target} already exists"
            G.add_edge(source, target, edge_type="manual", edge_count=1, edge_total=1)
            message = f"Added edge {source} - {target}"
        elif triggered_id == "remove-edge-button":
            if not G.has_edge(source, target):
                return dash.no_update, dash.no_update, f"No edge {source} - {target} to remove"
            try:
                if G.is_multigraph():
                    # Only source/target are chosen, so remove all parallel edges between them
                    G.remove_edges_from([(source, target, k) for k in list(G[source][target])])
                else:
                    G.remove_edge(source, target)
            except nx.NetworkXError as err:
                return dash.no_update, dash.no_update, f"Could not remove edge {source} - {target}: {err}"
            message = f"Removed edge {source} - {target}"
        else:
            raise PreventUpdate

        return json_graph.node_link_data(G), trigger_count + 1, message

    @app.callback(
        Output("colorscheme", "options"),
        Output("colorscheme", "value"),
        Input("color-schemes-store", "data")
    )
    def populate_colorscheme_dropdown(color_data):
        if not color_data:
            return [], None
        keys = list(color_data.keys())
        options = [{"label": k, "value": k} for k in keys]
        return options, keys[0]  # seleciona o primeiro por padrão




    run_options = {"jupyter_mode": "external", "debug": True}
    run_options.update(run_kwargs)
    app.run(**run_options)
 

def _cpm_resolution(cpm_resolution, resolution_parameter):
    # resolution_parameter is the old name of cpm_resolution; it wins when given
    if resolution_parameter is None:
        return cpm_resolution
    warnings.warn("resolution_parameter is deprecated, use cpm_resolution "
                  "(it only sets the CPM resolution of the Leidein_2 partition)",
                  DeprecationWarning, stacklevel=3)
    return resolution_parameter

def _undirected_max(G, weight):
    # A pair's directed or parallel edges become one undirected edge with the largest weight
    # (net_functions graphs give every edge of a pair the same edge_total; summing would double it)
    H = nx.Graph()
    H.add_nodes_from(G.nodes(data=True))
    for u, v, d in G.edges(data=True):
        if not H.has_edge(u, v) or (weight is not None and d.get(weight, 1) > H[u][v].get(weight, 1)):
            H.add_edge(u, v)
            H[u][v].clear()
            H[u][v].update(d)
    return H

def _community_scheme(nodes, communities, prefix):
    # [node -> label, label -> colour]: labels numbered by size, largest first (ties by smallest node name);
    # pastel has 10 colours and cycles, so more communities get evenly spaced hues
    import seaborn as sns
    ordered = sorted(communities, key=lambda c: (-len(c), min(str(n) for n in c)))
    label = {n: f"{prefix}{i}" for i, c in enumerate(ordered) for n in c}
    palette = sns.color_palette('pastel' if len(ordered) <= 10 else 'husl', len(ordered)).as_hex()
    return [{n: label[n] for n in nodes}, {f"{prefix}{i}": colour for i, colour in enumerate(palette)}]

def get_network_community(G,
                 community_to_color='leidein', weight='edge_total', cpm_resolution=0.05, seed=0,
                 *, resolution_parameter=None):
    """
    Detect communities with Louvain, Leiden and Leiden-CPM, as network_dash colour schemes.

    Partitions ``G`` three ways: Louvain [1]_ and Leiden [2]_ maximising standard modularity [3]_
    (resolution 1), and Leiden with the Constant Potts Model (CPM) [4]_ at resolution
    ``cpm_resolution``. Each partition is returned as a colour scheme ready for ``network_dash``.

    Parameters
    ----------
    G : networkx.Graph
        Graph to partition. A directed graph is partitioned as an undirected copy in which a pair's
        edges keep the largest weight; the parallel edges of a MultiGraph are collapsed to one edge
        per pair.
    community_to_color : str, default 'leidein'
        Unused, kept so positional calls still work. Any other value raises a DeprecationWarning.
    weight : str or None, default 'edge_total'
        Edge attribute used as weight, or None for unweighted. Edges lacking the attribute count as 1.
    cpm_resolution : float, default 0.05
        CPM resolution of ``Leidein_2``: a group of nodes forms a community only if its edges weigh
        more than this on average per pair of nodes, so higher values give more, smaller
        communities. Louvain and ``Leidein`` do not use it.
    seed : int or None, default 0
        Random seed for all three algorithms, so results are reproducible; None for unseeded runs.
    resolution_parameter : float, optional
        Deprecated name of ``cpm_resolution``. When given, it is used and raises a
        DeprecationWarning.

    Returns
    -------
    dict
        ``{"Louvain" | "Leidein" | "Leidein_2": [node_to_label, label_to_colour]}``. Labels are
        ``Louvain_<i>``, ``Leiden_<i>`` and ``LeidenCPM_<i>``, numbered by community size
        (0 = largest, ties broken by the smallest node name). Colours are pastel up to 10
        communities and evenly spaced hues beyond. An empty graph gives empty dicts.

    Warns
    -----
    UserWarning
        If ``G`` is empty or directed, or if some or all edges lack ``weight``.
    DeprecationWarning
        If ``community_to_color`` or ``resolution_parameter`` is passed.

    See Also
    --------
    community_stability : How much these partitions depend on ``seed``.
    clique_jaccard_community : Communities of overlapping maximal cliques.

    Notes
    -----
    A fixed seed gives a reproducible partition, not necessarily a stable one: modularity has many
    structurally different partitions with nearly the same score [5]_, and on weak community
    structure another seed can give another partition. Use ``community_stability`` to measure it.

    References
    ----------
    .. [1] Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. "Fast unfolding of
       communities in large networks." J. Stat. Mech. P10008 (2008).
       https://doi.org/10.1088/1742-5468/2008/10/P10008
    .. [2] Traag, V. A., Waltman, L. & van Eck, N. J. "From Louvain to Leiden: guaranteeing
       well-connected communities." Sci. Rep. 9, 5233 (2019).
       https://doi.org/10.1038/s41598-019-41695-z
    .. [3] Newman, M. E. J. & Girvan, M. "Finding and evaluating community structure in networks."
       Phys. Rev. E 69, 026113 (2004). https://doi.org/10.1103/PhysRevE.69.026113
    .. [4] Traag, V. A., Van Dooren, P. & Nesterov, Y. "Narrow scope for resolution-limit-free
       community detection." Phys. Rev. E 84, 016114 (2011).
       https://doi.org/10.1103/PhysRevE.84.016114
    .. [5] Good, B. H., de Montjoye, Y.-A. & Clauset, A. "Performance of modularity maximization in
       practical contexts." Phys. Rev. E 81, 046106 (2010).
       https://doi.org/10.1103/PhysRevE.81.046106
    """
    import leidenalg as la
    import igraph as ig

    cpm_resolution = _cpm_resolution(cpm_resolution, resolution_parameter)
    if community_to_color != 'leidein':
        warnings.warn("community_to_color is unused and deprecated; all three partitions are always returned",
                      DeprecationWarning, stacklevel=2)
    if G.number_of_nodes() == 0:
        warnings.warn("G has no nodes; returning empty community schemes.")
        return {name: [{}, {}] for name in ["Louvain", "Leidein", "Leidein_2"]}
    if G.is_directed():
        warnings.warn("G is directed; detecting communities on an undirected copy (a pair's edges keep the largest weight).")
        G = _undirected_max(G, weight)
    else:
        G = _simple_graph(G)

    #Edge weights: weight=None means unweighted; edges lacking the attribute count as 1
    g = ig.Graph.from_networkx(G)
    leiden_weight = weight
    if weight is not None:
        n_missing = sum(1 for u, v, d in G.edges(data=True) if weight not in d)
        if G.number_of_edges() == 0:
            weight = leiden_weight = None
        elif n_missing == G.number_of_edges():
            warnings.warn(f"No edge has the '{weight}' attribute; running community detection unweighted.")
            weight = leiden_weight = None
        elif n_missing:
            warnings.warn(f"{n_missing} of {G.number_of_edges()} edges lack the '{weight}' attribute; using 1 for them.")
            leiden_weight = [1 if w is None else w for w in g.es[weight]]

    #Partitions as sets of node names, never by position
    names = g.vs['_nx_name']
    def leiden(partition_type, **kwargs):
        part = la.find_partition(g, partition_type, weights=leiden_weight, n_iterations=-1, seed=seed, **kwargs)
        return [{names[i] for i in c} for c in part]

    louvain = nx.community.louvain_communities(G, weight=weight, resolution=1, seed=seed)
    return {"Louvain": _community_scheme(G, louvain, 'Louvain_'),
            "Leidein": _community_scheme(G, leiden(la.ModularityVertexPartition), 'Leiden_'),
            #Leiden with CPM, the only one using the resolution:
            "Leidein_2": _community_scheme(G, leiden(la.CPMVertexPartition, resolution_parameter=cpm_resolution), 'LeidenCPM_')}

def community_stability(G, weight='edge_total', cpm_resolution=0.05, n_seeds=50, seed=0, warn_below=0.8):
    """
    Measure how much the partitions of get_network_community depend on the seed.

    Runs ``get_network_community`` with seeds ``seed, seed + 1, ..., seed + n_seeds - 1``. The run
    with ``seed`` is the reference: the partition ``get_network_community(G, seed=seed)`` returns,
    i.e. the one being reported. Modularity has many structurally different partitions with nearly
    the same score [1]_, and Louvain [2]_ and Leiden [3]_ land on one of them depending on the seed;
    this function reports how much that changes the result, for each partition and each node.

    Parameters
    ----------
    G : networkx.Graph
        Graph to partition, as in ``get_network_community``.
    weight : str or None, default 'edge_total'
        Edge attribute used as weight, or None for unweighted.
    cpm_resolution : float, default 0.05
        CPM resolution of ``Leidein_2``.
    n_seeds : int, default 50
        Number of runs, including the reference.
    seed : int, default 0
        Seed of the reference run; the other runs use the following seeds.
    warn_below : float or None, default 0.8
        Warn for each partition whose median ARI is below this; None disables the warning.

    Returns
    -------
    summary : pandas.DataFrame
        One row per partition (``Louvain``, ``Leidein``, ``Leidein_2``) with columns
        ``communities`` (communities in the reference), ``distinct_partitions`` (distinct
        partitions among the runs) and ``ari_median``, ``ari_min`` (adjusted Rand index [4]_ over
        all pairs of runs: 1 = identical, about 0 = chance-level agreement).
    nodes : pandas.DataFrame
        One row per node with, for each partition, ``<name>`` (the node's reference community) and
        ``<name>_stability``: the mean Jaccard index [5]_ between the node's community in each other
        run and its reference community (1 = always the same community; splits and merges both
        lower it).

    Raises
    ------
    ValueError
        If ``n_seeds`` is less than 2.

    Warns
    -----
    UserWarning
        For each partition whose median ARI is below ``warn_below``. Warnings raised by the
        repeated ``get_network_community`` calls are shown once each.

    See Also
    --------
    get_network_community : The partitions being measured.

    Notes
    -----
    The node stability follows the co-assignment idea of consensus clustering [6]_ (how often nodes
    end up together across runs), measured against the reference instead of building a node x node
    matrix, so it costs O(nodes) per run. This measures stability across seeds only: not robustness
    to sampling of the data, nor whether the graph has significant community structure. Runtime is
    about ``n_seeds`` ``get_network_community`` calls.

    References
    ----------
    .. [1] Good, B. H., de Montjoye, Y.-A. & Clauset, A. "Performance of modularity maximization in
       practical contexts." Phys. Rev. E 81, 046106 (2010).
       https://doi.org/10.1103/PhysRevE.81.046106
    .. [2] Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. "Fast unfolding of
       communities in large networks." J. Stat. Mech. P10008 (2008).
       https://doi.org/10.1088/1742-5468/2008/10/P10008
    .. [3] Traag, V. A., Waltman, L. & van Eck, N. J. "From Louvain to Leiden: guaranteeing
       well-connected communities." Sci. Rep. 9, 5233 (2019).
       https://doi.org/10.1038/s41598-019-41695-z
    .. [4] Hubert, L. & Arabie, P. "Comparing partitions." J. Classif. 2, 193-218 (1985).
       https://doi.org/10.1007/BF01908075
    .. [5] Jaccard, P. "The distribution of the flora in the alpine zone." New Phytol. 11, 37-50
       (1912). https://doi.org/10.1111/j.1469-8137.1912.tb05611.x
    .. [6] Lancichinetti, A. & Fortunato, S. "Consensus clustering in complex networks."
       Sci. Rep. 2, 336 (2012). https://doi.org/10.1038/srep00336
    """
    import igraph as ig
    from collections import Counter

    if n_seeds < 2:
        raise ValueError(f"n_seeds must be >= 2 to compare runs, got {n_seeds}")
    # The repeated calls repeat their warnings (directed graph, missing weights): show each once
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        runs = [get_network_community(G, weight=weight, cpm_resolution=cpm_resolution, seed=s)
                for s in range(seed, seed + n_seeds)]
    for message, category in dict.fromkeys((str(w.message), w.category) for w in caught):
        warnings.warn(message, category, stacklevel=2)

    names = ["Louvain", "Leidein", "Leidein_2"]
    nodes = list(runs[0]["Louvain"][0])
    summary = pd.DataFrame(index=pd.Index(names, name='partition'),
                           columns=['communities', 'distinct_partitions', 'ari_median', 'ari_min'], dtype=float)
    node_table = pd.DataFrame(index=pd.Index(nodes, name='node'))
    if not nodes:
        summary[['communities', 'distinct_partitions']] = [0, 1]
        for name in names:
            node_table[name], node_table[f"{name}_stability"] = pd.Series(dtype=object), pd.Series(dtype=float)
        return summary, node_table

    for name in names:
        # Codes by first appearance in node order: relabelled copies of a partition get the same codes
        memberships = [pd.factorize(pd.Series([run[name][0][n] for n in nodes]))[0] for run in runs]
        counts = Counter(m.tobytes() for m in memberships)
        partitions = {m.tobytes(): m for m in memberships}

        # ARI over all pairs of runs, computed once per pair of distinct partitions
        aris, pairs = [], []
        keys = list(counts)
        for i, a in enumerate(keys):
            if counts[a] > 1:
                aris.append(1.0)
                pairs.append(counts[a] * (counts[a] - 1) // 2)
            for b in keys[i + 1:]:
                aris.append(ig.compare_communities(partitions[a].tolist(), partitions[b].tolist(), method='adjusted_rand'))
                pairs.append(counts[a] * counts[b])
        all_aris = np.repeat(aris, pairs)

        # Node stability: Jaccard of its community in each run with its reference community
        ref = memberships[0]
        ref_size = np.bincount(ref)[ref]
        stability = np.zeros(len(nodes))
        for m in memberships[1:]:
            pair = ref * (m.max() + 1) + m
            shared = np.bincount(pair)[pair]
            stability += shared / (ref_size + np.bincount(m)[m] - shared)

        summary.loc[name] = [ref.max() + 1, len(counts), np.median(all_aris), all_aris.min()]
        node_table[name] = [runs[0][name][0][n] for n in nodes]
        node_table[f"{name}_stability"] = stability / (n_seeds - 1)
        if warn_below is not None and np.median(all_aris) < warn_below:
            warnings.warn(f"{name}: median ARI {np.median(all_aris):.2f} across {n_seeds} seeds is below {warn_below}; "
                          f"its communities depend on the seed (the nodes table shows which nodes are unstable).",
                          stacklevel=2)
    summary[['communities', 'distinct_partitions']] = summary[['communities', 'distinct_partitions']].astype(int)
    return summary, node_table

def _clique_pairs(cliques, nodes, min_jaccard_index, warn_checked, warn_kept, block=1000):
    # Pairs of cliques sharing at least one node, from a sparse clique x node matrix times its transpose,
    # one block of cliques at a time so memory stays bounded; each pair once (i < j), no clique with itself.
    # Returns positions i, j (sorted) and their Jaccard index.
    from scipy import sparse
    node_id = {n: k for k, n in enumerate(nodes)}
    size = np.array([len(c) for c in cliques], dtype=np.int64)
    n = len(cliques)
    rows = np.repeat(np.arange(n), size)
    cols = np.fromiter((node_id[v] for c in cliques for v in c), dtype=np.int64, count=size.sum())
    M = sparse.csr_matrix((np.ones(len(cols), dtype=np.float32), (rows, cols)), shape=(n, len(node_id)))
    MT = M.T.tocsr()
    kept_i, kept_j, kept_jaccard = [], [], []
    for start in range(0, n, block):
        part = (M[start:start + block] @ MT).tocoo()
        i = part.row.astype(np.int64) + start
        j = part.col.astype(np.int64)
        upper = j > i
        i, j, shared = i[upper], j[upper], part.data[upper].astype(np.float64)
        if start == 0 and warn_checked is not None:
            # Overlap rate in the first block, extrapolated to all pairs
            block_pairs = sum(n - 1 - r for r in range(min(block, n)))
            estimate = len(i) / block_pairs * n * (n - 1) / 2 if block_pairs else 0
            if estimate > warn_checked:
                warnings.warn(f"about {estimate:,.0f} pairs of the {n:,} cliques share nodes and must be checked "
                              f"(more than warn_checked={warn_checked:,}); this can take minutes. "
                              "Raise min_clique_size or min_jaccard_index to reduce it.", stacklevel=3)
        jaccard = shared / (size[i] + size[j] - shared)
        keep = jaccard >= min_jaccard_index
        order = np.lexsort((j[keep], i[keep]))
        kept_i.append(i[keep][order]); kept_j.append(j[keep][order]); kept_jaccard.append(jaccard[keep][order])
    i, j, jaccard = (np.concatenate(x) if x else np.array([], dtype=np.int64) for x in (kept_i, kept_j, kept_jaccard))
    if warn_kept is not None and len(i) > warn_kept:
        warnings.warn(f"{len(i):,} pairs of cliques pass min_jaccard_index (more than warn_kept={warn_kept:,}); "
                      "the result table will be very large. Raise min_clique_size or min_jaccard_index to reduce it.",
                      stacklevel=3)
    return i, j, jaccard

def clique_jaccard_community(G,
                             min_jaccard_index=0.4,
                             cpm_resolution=0.5,
                             min_clique_size=10,
                             seed=0,
                             *, resolution_parameter=None,
                             warn_checked=100_000_000,
                             warn_kept=5_000_000):
    """
    Detect communities of overlapping maximal cliques.

    Finds the maximal cliques of ``G`` with at least ``min_clique_size`` nodes, connects each pair
    of cliques whose Jaccard index [1]_ (nodes in both over nodes in either) is at least
    ``min_jaccard_index``, and partitions that clique graph with ``get_network_community``, using
    the Jaccard index as edge weight.

    Parameters
    ----------
    G : networkx.Graph
        Graph whose maximal cliques are compared.
    min_jaccard_index : float, default 0.4
        Minimum Jaccard index for two cliques to be connected. Must be greater than 0.
    cpm_resolution : float, default 0.5
        CPM resolution of the ``Leidein_2`` partition (the ``*_leidein_resolution`` columns),
        passed to ``get_network_community``.
    min_clique_size : int, default 10
        Smallest clique, in nodes, that is compared.
    seed : int or None, default 0
        Random seed passed to ``get_network_community``.
    resolution_parameter : float, optional
        Deprecated name of ``cpm_resolution``.
    warn_checked : int or None, default 100_000_000
        Warn when the estimated number of clique pairs sharing a node, which must all be checked,
        exceeds this. Estimated after the first block of 1,000 cliques. None disables it.
    warn_kept : int or None, default 5_000_000
        Warn when more pairs than this pass ``min_jaccard_index``, i.e. the returned table will be
        very large. None disables it.

    Returns
    -------
    pandas.DataFrame
        One row per connected pair of cliques, with ``source < target``. Columns: ``source`` and
        ``target`` (clique numbers, their position in ``networkx.find_cliques`` output),
        ``cliques_x``, ``clique_size_x``, ``cliques_y``, ``clique_size_y`` (each clique's nodes and
        size), ``t`` (constant 1), ``uni`` and ``inter`` (sets of nodes in either and in both
        cliques), ``Jaccard_index``, and the community of each clique in the ``Leidein_2``,
        ``Leidein`` and ``Louvain`` partitions: ``source_leidein_resolution``,
        ``target_leidein_resolution``, ``source_leidein``, ``target_leidein``, ``source_louvain``,
        ``target_louvain``. Empty, with the same columns, when no pair qualifies.

    Raises
    ------
    ValueError
        If ``min_jaccard_index`` is not greater than 0.

    Warns
    -----
    UserWarning
        If the estimated pairs to check exceed ``warn_checked``, or the kept pairs exceed
        ``warn_kept``.

    See Also
    --------
    get_network_community : The partitions computed on the clique graph.

    Notes
    -----
    Pairs come from the number of nodes each pair of cliques shares, computed as a sparse
    clique x node matrix times its transpose in blocks of 1,000 cliques. Pairs sharing no node
    (Jaccard index 0) are never formed, which is why ``min_jaccard_index`` must be positive, and a
    clique that overlaps no other above the threshold does not appear in the result. The cost grows
    with the number of pairs sharing a node, not with all pairs, and memory is bounded by the
    block. When nearly every pair overlaps (many maximal cliques on few nodes) the work is quadratic
    again; raising ``min_clique_size`` or ``min_jaccard_index`` reduces it.

    References
    ----------
    .. [1] Jaccard, P. "The distribution of the flora in the alpine zone." New Phytol. 11, 37-50
       (1912). https://doi.org/10.1111/j.1469-8137.1912.tb05611.x
    """
    cpm_resolution = _cpm_resolution(cpm_resolution, resolution_parameter)
    if not min_jaccard_index > 0:
        raise ValueError(f"min_jaccard_index must be > 0 (pairs sharing no node are never connected), got {min_jaccard_index}")
    cliques = list(nx.find_cliques(G))
    cliques = pd.Series(cliques, dtype=object).rename('cliques').to_frame()
    cliques['clique_size'] = cliques.cliques.apply(lambda x: len(x))
    cliques = cliques.query('clique_size >= @min_clique_size')
    cliques['t']  = 1
    cliques = cliques.reset_index().rename({'index':'clique_n'}, axis=1)
    i, j, jaccard = _clique_pairs(cliques.cliques.tolist(), list(G), min_jaccard_index, warn_checked, warn_kept)
    # Same columns, order and index as the former all-pairs merge of cliques with itself
    left = cliques.iloc[i].reset_index(drop=True).add_suffix('_x').rename({'t_x': 't'}, axis=1)
    right = cliques.iloc[j].reset_index(drop=True).add_suffix('_y').drop(columns='t_y')
    cm = pd.concat([left, right], axis=1)
    cm.index = i * len(cliques) + j
    cm['uni'] = [set(a) | set(b) for a, b in zip(cm.cliques_x, cm.cliques_y)]
    cm['inter'] = [set(a) & set(b) for a, b in zip(cm.cliques_x, cm.cliques_y)]
    cm['Jaccard_index'] = jaccard.astype(float)
    cm.rename({"clique_n_x" :'source', 'clique_n_y' :'target'}, axis=1, inplace=True)
    community_columns = ['source_leidein_resolution', 'target_leidein_resolution', 'source_leidein',
                         'target_leidein', 'source_louvain', 'target_louvain']
    if cm.empty:
        return cm.assign(**{c: pd.Series(dtype=object) for c in community_columns})
    G = nx.from_pandas_edgelist(cm, edge_attr='Jaccard_index')
    z = get_network_community(G, weight='Jaccard_index', cpm_resolution=cpm_resolution, seed=seed)
    cm['source_leidein_resolution'] = cm.source.replace(z['Leidein_2'][0])
    cm['target_leidein_resolution'] = cm.target.replace(z['Leidein_2'][0])
    cm['source_leidein'] = cm.source.replace(z['Leidein'][0])
    cm['target_leidein'] = cm.target.replace(z['Leidein'][0])
    cm['source_louvain'] = cm.source.replace(z['Louvain'][0])
    cm['target_louvain'] = cm.target.replace(z['Louvain'][0])
    return cm

def subgraph_stats(net,
                   subnet,
                   theme_dict,
                   theme_to_stats,
                   alternative = "two-sided",
                   return_df = False):
    from  scipy.stats import fisher_exact
    import pandas as pd
    to_test = pd.Series(list(subnet.nodes),name="node",dtype='object').to_frame()
    to_test['theme'] = to_test.node.map(theme_dict)
    total = pd.Series(list(net.nodes),name="node").to_frame()
    total['theme'] = total.node.map(theme_dict)
    sub_total = to_test.node.nunique()
    sub_theme = to_test.query('theme == @theme_to_stats').node.nunique()
    ttt = to_test.node.tolist()
    total = total[~total.node.isin(ttt)]
    net_total = total.node.nunique()
    net_theme = total.query('theme == @theme_to_stats').node.nunique()
    to_fisher = [[sub_theme, sub_total - sub_theme], [net_theme, net_total - net_theme]] 
    if return_df:
        print(fisher_exact(to_fisher, alternative=alternative))
        return pd.DataFrame(to_fisher, columns=['Theme', 'Total'], index=['Subnetwork', 'net - subnet'] )
    return fisher_exact(to_fisher,alternative=alternative)

def get_relations(domdf,
                  dom='full',
                  min_connections = 5,
                  self_relation = False,
                  to_yaml=False,
                  iteration=1,
                  filter_rename_yaml  = False,
                  *,
                  gene_pos,
                  feature_pos):
    """
    domdf to network, it can create dictionary to create a YAML file or a dataframe.
    dom: list of domains to be used as query in the network, if full all domdf will be used
    min_connections: Minimum number of relationship (Operon or fusion) anode should have to be presente in the network.
    yaml: If decided to transform the network in dictionary to latter be exported as YAML file.
    Iteration: If a list of domains was supplied, the number of iteration in the network to collect other nodes.
               Must be >= 1. For a network of the seed domains only, filter the result:
               edges[edges.source.isin(dom) & edges.target.isin(dom)]
    gene_pos: column ordering genes along the genome (smaller = upstream), e.g. start coordinate.
    feature_pos: column ordering domains along the protein, N->C (smaller = N-terminal),
                 e.g. domain start in amino acid coordinates.
    Dataframe output: one row per (source, target, edge_type); edge_count is the number of
    blocks with that relation, edge_total the number of blocks with the pair in any type.
    Row order does not matter. For gian_func.split_domain output, whose domp
    follows genome order (C->N on - strand), use gene_pos='blockp' and:
        domdf['dom_aa_order'] = np.where(domdf.strand == -1, -domdf.domp, domdf.domp)
    """
    from rotifer.devel.alpha import gian_func as gf
    domdf = domdf.copy()
    domdf.dom = domdf.dom.str.strip()
    if filter_rename_yaml:
        import yaml
        with open(filter_rename_yaml, "r") as file:
            data = yaml.safe_load(file)
            data = pd.DataFrame(data).T
            data.loc[data.Display_name =="", 'Display_name'] = data.loc[data.Display_name ==""].index.tolist()
            rename_dict = data.Display_name.to_dict()
            domdf.dom = domdf.dom.replace(rename_dict)
            to_display = data.query('Include ==1').Display_name.unique().tolist()
            domdf = domdf.query('dom in @to_display')

    if isinstance(dom, str):
        dom = [dom]
    if dom[0] =='full':
        block_id = domdf.block_id.tolist()
    else:
        if iteration < 1:
            raise ValueError(f"iteration must be >= 1, got {iteration}")
        missing = set(dom) - set(domdf.dom)
        if missing:
            warnings.warn(f"seed domains not found in domdf: {sorted(missing)}")
        for x in range(iteration):
            block_id = domdf.query('dom in @dom').block_id.tolist()
            dom = domdf.query('block_id in @block_id').dom.unique().tolist()
    to_net = domdf.query('block_id in @block_id')[['pid', 'block_id', 'dom', gene_pos, feature_pos]]
    to_net.columns = ['pid', 'block_id', 'dom', 'gpos', 'fpos']
    to_net = to_net.merge(to_net, on ='block_id')
    # Edge type from positions, not row order: genes by gene_pos, domains of the same protein by feature_pos
    same_pid = to_net.pid_x == to_net.pid_y
    to_net['edge_type'] = np.select(
            [~same_pid & (to_net.gpos_x < to_net.gpos_y),
             ~same_pid & (to_net.gpos_x > to_net.gpos_y),
             same_pid & (to_net.fpos_x < to_net.fpos_y),
             same_pid & (to_net.fpos_x > to_net.fpos_y)],
            ['Downstream', 'Upstream', 'Cterminal', 'Nterminal'],
            default='')

    to_net= to_net.query('edge_type != ""')
    to_net = to_net[['block_id','dom_x','dom_y','edge_type']]
    to_net.columns = ['block_id','source','target','edge_type']
    # Count each relation once per block, however many times it occurs there
    to_net = to_net.drop_duplicates()
    if not self_relation:
        to_net = to_net[to_net.source != to_net.target]
    if to_yaml:
        to_net['edge_count'] = to_net.groupby(['source', 'target']).block_id.transform('nunique')
        to_net = to_net.query('edge_count >= @ min_connections').drop(['edge_count'], axis=1)
        if to_net.empty:
            return {}
        yamldf = pd.pivot_table(to_net,
                            index='source',
                            columns='edge_type',
                            values = ['target'],
                            aggfunc={
                                'target':lambda x : gf.count_series(x, as_list=True)}
                            )
        # Edge types absent from the filtered relations get no pivot column
        yamldf = yamldf.loc[:,'target'].reindex(columns=['Cterminal', 'Downstream', 'Nterminal', 'Upstream'], fill_value='')
        yamldf = to_net.groupby('source').block_id.nunique().rename('total').to_frame().join(yamldf).sort_values('total')
        yamldf.fillna('', inplace=True)
        yamldf.Nterminal = yamldf.Nterminal.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Cterminal = yamldf.Cterminal.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Downstream = yamldf.Downstream.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Upstream = yamldf.Upstream.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf['L'] = yamldf[['Cterminal', 'Downstream', 'Nterminal', 'Upstream']].apply(
                lambda row: [item for sublist in row for item in sublist if item != ''],
                axis=1
                )
        yamldf.L = yamldf.L.apply(len)
        yamldf = yamldf.query('L >0')
        yamldf = yamldf.drop(['L'], axis=1)

        return(yamldf.to_dict(orient='index'))
    else:
        to_net['edge_total'] = to_net.groupby(['source', 'target']).block_id.transform('nunique')
        to_net = to_net.groupby(['source', 'target', 'edge_type']).agg(
                {'block_id' :'count',
                 'edge_total' : 'first'}).reset_index().sort_values('block_id')
        to_net = to_net.rename({'block_id':'edge_count'},axis=1)
        to_net = to_net.query('edge_total >= @min_connections')
        return(to_net)
