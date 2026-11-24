#!/usr/bin/env python3
"""Open rotifer.view.network.network_dash on the toy domain table.

Usage (inside the rotifer env):  python test/data/graph/run_dash_demo.py [--split] [--port N]
  --split  one edge per edge type (network_annotation split_types=True)
Then browse to http://127.0.0.1:8050
"""
import argparse
import os
import pandas as pd
from rotifer.graph import network as nf
from rotifer.view.network import network_dash

parser = argparse.ArgumentParser()
parser.add_argument("--split", action="store_true", help="one edge per edge type")
parser.add_argument("--port", type=int, default=8050)
args = parser.parse_args()

here = os.path.dirname(os.path.abspath(__file__))
domdf = pd.read_csv(os.path.join(here, "toy_domains.tsv"), sep="\t")

# Node -> group, group -> colour
function = {
    **dict.fromkeys(["SBP_bac_5", "BPD_transp_1", "ABC_tran", "oligo_HPY"], "transport"),
    **dict.fromkeys(["Response_reg", "Trans_reg_C", "HAMP", "HisKA", "HATPase_c"], "signalling"),
    **dict.fromkeys(["PhdYeFM_antitox", "ParE_toxin", "MethyltransfD12", "RE_AlwI"], "defence"),
}
xx = {"transport": "#4C9BE8", "signalling": "#E8A84C", "defence": "#D9534F"}

edges = nf.get_relations(domdf, min_connections=2, gene_pos="start", feature_pos="dom_start")
G = nf.network_annotation(edges, node_color={n: xx[function[n]] for n in function}, split_types=args.split)
print(f"{type(G).__name__}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# dcc.Store needs JSON, so turn numpy coordinates into lists
pos_dict = {
    name: {n: [float(c) for c in p] for n, p in layout.items()}
    for name, layout in nf.calculate_positions(G, iterations=50).items()
}

color_schemes = {"Theme": [function, xx]}
try:
    color_schemes.update(nf.get_network_community(G, cpm_resolution=0.05))
except ImportError as e:
    print(f"Skipping community colour schemes ({e})")

network_dash(G, pos_dict, xx, function, color_schemes, port=args.port)
