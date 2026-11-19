#!/usr/bin/env python3
"""Tests for NeighborhoodDF.to_blocks and NeighborhoodDF.to_yaml.

rotifer.pandas is genome-neutral, so the neighborhood defaults (grouping by c80e3,
the shown columns, the statistics, the sort orders) come from the presets in
rotifer.genome.blocks. These tests pin the text those presets produce.
"""

import pandas as pd

from rotifer.genome.data import NeighborhoodDF
from rotifer.genome import blocks as presets
from rotifer.pandas.yaml_blocks import read_yaml


def text(lines):
    return '\n'.join(lines)


def genome_frame():
    """Has every column the to_blocks preset looks for."""
    return pd.DataFrame({
        'pid': ['p1', 'p2', 'p3', 'p4'], 'c80i70': ['i1', 'i1', 'i2', 'i3'],
        'c100i100': ['j1', 'j2', 'j3', 'j4'], 'plen': [100, 20, 300, 40],
        'pfam': ['PF1', None, 'PF2', 'PF3'], 'aravind': ['a', 'b', 'c', 'd'],
        'product': ['prod|1', 'prod2', 'prod3', 'prod4'],
        'organism': ['Org B', 'Org A', 'Org C', 'Org A'],
        'lineage': ['L2', 'L1', 'L1', 'L3'], 'classification': ['k2', 'k1', 'k3', 'k1'],
        'block_id': [1, 1, 2, 3], 'c80e3': ['X', 'X', 'Y', 'Z']})


def yaml_frame():
    """Adds the taxonomy columns the to_yaml preset aggregates and sorts on."""
    df = genome_frame()
    df['superkingdom'] = ['Bacteria', 'Bacteria', 'Archaea', 'Bacteria']
    df['phylum'] = ['P1', 'P2', 'P3', 'P1']
    df['taxid'] = [10, 20, 30, 40]
    return df


GENOME_DEFAULT_BLOCKS = [
    '#',
    '# c80e3=X; lineage=2; block_id=1; c80i70=1; pid=2',
    '  : pid : c80i70 : c100i100 : plen : pfam : aravind : product : organism : lineage : classification',
    '1 : p2  : i1     : j2       : 20   :      : b       : prod2   : Org A    : L1      : k1            ',
    '0 : p1  : i1     : j1       : 100  : PF1  : a       : prod|1  : Org B    : L2      : k2            ',
    '#',
    '# c80e3=Y; lineage=1; block_id=1; c80i70=1; pid=1',
    '  : pid : c80i70 : c100i100 : plen : pfam : aravind : product : organism : lineage : classification',
    '2 : p3  : i2     : j3       : 300  : PF2  : c       : prod3   : Org C    : L1      : k3            ',
    '#',
    '# c80e3=Z; lineage=1; block_id=1; c80i70=1; pid=1',
    '  : pid : c80i70 : c100i100 : plen : pfam : aravind : product : organism : lineage : classification',
    '3 : p4  : i3     : j4       : 40   : PF3  : d       : prod4   : Org A    : L3      : k1            ',
]


def test_to_blocks_genome_defaults():
    assert NeighborhoodDF(genome_frame()).to_blocks() == text(GENOME_DEFAULT_BLOCKS)


def test_to_blocks_keywords_replace_the_preset():
    out = NeighborhoodDF(genome_frame()).to_blocks(columns=['pid'], stats={'pid': 'nunique'},
                                                   sortrows={})
    assert out.splitlines()[:5] == ['#', '# c80e3=X; pid=2', '  : pid', '0 : p1 ', '1 : p2 ']


def test_to_blocks_buf(tmp_path):
    target = tmp_path / 'blocks.txt'
    assert NeighborhoodDF(genome_frame()).to_blocks(buf=target) is None
    assert target.read_text() == text(GENOME_DEFAULT_BLOCKS)


def test_to_yaml_preset_tables_round_trip():
    doc = NeighborhoodDF(yaml_frame()).to_yaml()
    tables = read_yaml(doc, merge=False)
    assert set(tables) == {'annotation', 'stats', 'dist', 'rows'}
    assert list(tables['stats'].columns) == ['c80e3', 'superkingdom', 'lineage', 'phylum',
                                             'taxid', 'c80i70', 'pid']
    # sortblocks: phylum descending, then c80i70 descending. X has 2 phyla (counts
    # of distinct values per block), Y and Z one each.
    assert list(tables['stats']['c80e3']) == ['X', 'Y', 'Z']
    assert list(tables['stats']['pid']) == [2, 1, 1]
    assert list(tables['rows'].columns) == ['c80e3', 'pid', 'c80i70', 'c100i100', 'plen', 'pfam',
                                            'aravind', 'product', 'organism', 'lineage',
                                            'classification']
    assert len(tables['rows']) == 4


def test_to_yaml_keywords_replace_the_preset():
    doc = NeighborhoodDF(yaml_frame()).to_yaml(columns=['pid'], sortblocks={'c80e3': False})
    tables = read_yaml(doc, merge=False)
    assert list(tables['rows'].columns) == ['c80e3', 'pid']
    assert list(tables['stats']['c80e3']) == ['Z', 'Y', 'X']


def test_presets_are_not_shared_between_calls():
    NeighborhoodDF(genome_frame()).to_blocks(columns=['pid'])
    assert presets.NEIGHBORHOOD_BLOCKS['columns'][0] == 'pid'
    assert len(presets.NEIGHBORHOOD_BLOCKS['columns']) == 10
