"""
Interactive viewer for domain networks.

:func:`network_dash` opens a graph made by
:func:`rotifer.graph.network.network_annotation` in a Dash app with a
Cytoscape canvas, where nodes can be filtered, recoloured, moved and
connected by hand. :func:`convert_to_cytoscape` turns a graph and a
layout into the elements the canvas draws. The analysis behind the
graph (relations, layouts, communities) is in :mod:`rotifer.graph.network`.

This module needs Dash and dash-cytoscape (the ``view`` extra,
``pip install rotifer[view]``). Importing it loads the extra Cytoscape
layouts of ``dash_cytoscape``, which also enables SVG export in the
viewer. :func:`network_dash` starts a local web server; nothing here
contacts a remote service.
"""
import dash_cytoscape as cyto
import networkx as nx

# Enable SVG export
cyto.load_extra_layouts()

def convert_to_cytoscape(G, positions, layout_key, scaling_factor=900, x_scale=1, y_scale=1):
    """
    Convert a graph and one of its layouts to Dash Cytoscape elements.

    :func:`network_dash` calls it on every redraw, with the layout and the
    scaling chosen in the viewer.

    Parameters
    ----------
    G : networkx.Graph
        Graph to convert. All node and edge attributes go into each
        element's ``data``.
    positions : dict
        Layout name to ``{node: [x, y]}``, as built for
        :func:`network_dash`.
    layout_key : str
        Name of the layout in ``positions`` to use.
    scaling_factor : float, default 900
        Pixels per layout unit. Larger values spread the nodes further
        apart.
    x_scale : float, default 1
        Extra horizontal stretch.
    y_scale : float, default 1
        Extra vertical stretch.

    Returns
    -------
    list of dict
        Node elements followed by edge elements, in the format of
        ``networkx.cytoscape_data``. Node elements found in the layout also
        get a ``position`` with integer ``x`` and ``y`` in pixels.

    See Also
    --------
    network_dash : The viewer that draws these elements.

    Notes
    -----
    The layout's second coordinate becomes the horizontal position and its
    first coordinate the vertical one, so the drawing is the layout
    mirrored along its diagonal. Saved positions use the same convention,
    so they come back unchanged.
    """
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
    Open an interactive viewer for a domain network.

    Starts a Dash app that draws ``G`` with Dash Cytoscape, a wrapper around
    Cytoscape.js [1]_. In the viewer you can switch layouts, show only some
    groups of nodes, hide pairs seen in few blocks, recolour nodes and
    edges, add or remove edges, save node positions as new layouts, and
    export the picture.

    Parameters
    ----------
    G : networkx.Graph
        Graph from :func:`~rotifer.graph.network.network_annotation`. Its ``node_shape`` and
        ``node_outline`` attributes are used as they are; node and edge
        colours are recomputed from the selected groups.
    pos_dict : dict
        Layout name to ``{node: [x, y]}``, with plain lists of floats (the
        data is stored in the browser as JSON). The first layout is shown
        first.
    xx : dict
        Colour of each group, as a CSS colour.
    function : dict
        Group of each node. Every node needs a group, and every group a
        colour in ``xx``.
    color_schemes : dict
        Named alternatives to ``function`` and ``xx``, each
        ``[node_to_group, group_to_colour]``, for example the output of
        :func:`~rotifer.graph.network.get_network_community` or ``{"Theme": [function, xx]}``.
        Choosing one in the viewer replaces the current groups and colours.
    **run_kwargs
        Passed to ``dash.Dash.run``, for example ``port=8051`` or
        ``jupyter_mode="inline"``. Defaults: ``jupyter_mode="external"``
        and ``debug=True``.

    Returns
    -------
    None
        The function returns when the server stops. Outside Jupyter it
        blocks until then.

    See Also
    --------
    rotifer.graph.network.network_annotation : Builds ``G``.
    rotifer.graph.network.calculate_positions : Builds the layouts for ``pos_dict``.
    rotifer.graph.network.get_network_community : Builds colour schemes from communities.

    Notes
    -----
    The graph is copied into the page when the app starts. Edges added or
    removed in the viewer, saved positions and recoloured nodes live only in
    that page and are lost when it reloads; ``G`` itself is never changed.

    What the controls do:

    - "Minimmun connection" hides edges whose ``edge_total`` (number of
      blocks with the pair) is below the value. Edges added by hand have
      ``edge_type='manual'`` and are always shown.
    - Node size is the square root of the node's degree in the shown
      graph, times "Node size".
    - Double-clicking an edge cycles it through the edge styles on screen;
      double-clicking a node moves it to the next group colour, which also
      changes its group.
    - "Update layouts" recomputes :func:`~rotifer.graph.network.calculate_positions` on the shown
      graph and adds the results as layouts named ``updated ...``.
    - The "Update Communities" button has no action yet.

    References
    ----------
    .. [1] Franz, M. et al. "Cytoscape.js: a graph theory library for
       visualisation and analysis." Bioinformatics 32, 309-311 (2016).
       https://doi.org/10.1093/bioinformatics/btv557

    Examples
    --------
    Positions from :func:`~rotifer.graph.network.calculate_positions` hold numpy arrays, which the
    browser store cannot hold, so convert them to lists first:

    >>> from rotifer.graph.network import calculate_positions  # doctest: +SKIP
    >>> layouts = calculate_positions(G, iterations=50)  # doctest: +SKIP
    >>> pos_dict = {name: {n: [float(c) for c in p] for n, p in layout.items()}
    ...             for name, layout in layouts.items()}  # doctest: +SKIP
    >>> network_dash(G, pos_dict, xx, function, {"Theme": [function, xx]})  # doctest: +SKIP
    """
    from dash import Dash, dcc, html, Input, Output, ctx, State
    import dash
    import dash_cytoscape as cyto
    import networkx as nx
    from networkx.readwrite import json_graph
    import numpy as np
    import pandas as pd
    from rotifer.graph import network as net_functions
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
        """Rebuild the stylesheet from the font size, outline width and edge width boxes."""
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
            
        element = convert_to_cytoscape(subgraph, pos_dict, selected_option,x_scale=x_scale, y_scale=y_scale, scaling_factor=scale_value)
        return element
    """ 

    @app.callback(
        Output("select_source", "value"),
        Output("select_target", "value"),
        Input("cytoscape", "tapNodeData"),
        State("select-mode", "value"),
    )
    def update_input_boxes(node_data,toggle):
        """Copy a clicked node's name into the source or target box, as the toggle says."""
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

        """Redraw the graph, or restyle a double-clicked edge or recolour a double-clicked node."""
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
        """Save the nodes' current positions as a new layout named by the text box."""
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
        
        """Print the last saved positions to the console and show them on the page."""
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
        """Export the drawing as JPG, PNG or SVG."""
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
        """Load the groups and colours of the chosen colour scheme."""
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
        """Add a group and its colour to the current colours."""
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
        """Recompute the layouts on the graph currently shown."""
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
        """Add or remove the edge between the source and target nodes, and report the result."""
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
        """List the colour schemes in the dropdown and select the first."""
        if not color_data:
            return [], None
        keys = list(color_data.keys())
        options = [{"label": k, "value": k} for k in keys]
        return options, keys[0]  # seleciona o primeiro por padrão




    run_options = {"jupyter_mode": "external", "debug": True}
    run_options.update(run_kwargs)
    app.run(**run_options)
