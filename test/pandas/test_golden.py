#!/usr/bin/env python3
"""Golden-output (characterization) tests for the public API of rotifer.pandas.

These pin what the package emits today, byte for byte, so that the internals can
be reorganised -- module split, a shared block engine, a padding kernel, a new
read_yaml pipeline -- without the output drifting. Everything is reached through
public names only (to_string, to_blocks, default_header, yaml_blocks.to_yaml /
read_yaml / md_table / read_md_table, functions.pad / to_color /
optimize_memory_usage / print_everything / show_display_options).

Expected values are literals captured from the current implementation, never
recomputed by the code under test. A failure here means the output changed: either
a regression, or a deliberate change that must update the literal in the same
commit. Some literals record behaviour that looks accidental (see the QUIRK
comments); they are pinned on purpose and may be changed together with a fix.

The dtype strings assume pandas >= 3, where text columns are 'str' (not 'object').

Run under pytest, or standalone: python test/pandas/test_golden.py
"""

import io

import numpy as np
import pandas as pd
import pytest

import rotifer.pandas as rp
import rotifer.pandas.functions as rpf
from rotifer.pandas import yaml_blocks as yb


def frame():
    """Three groups (B, A, C, in first-appearance order), a NaN float, a missing
    string, a pipe, an embedded newline."""
    return pd.DataFrame({'grp': ['B', 'A', 'B', 'C', 'A'], 'sub': ['x', 'y', 'x', 'z', 'y'],
                         'name': ['alpha|x', 'beta\ngamma', 'delta', 'eps', None],
                         'score': [1.5, np.nan, 2.25, 3.0, 10.125],
                         'n': [3, 10, 200, 4, 5]})


def multilevel_frame():
    df = frame()
    df.columns = pd.MultiIndex.from_tuples(
        [('m', 'grp'), ('m', 'sub'), ('d', 'name'), ('d', 'score'), ('d', 'n')])
    return df


def genome_frame():
    """Has every column the genome-specific defaults of to_blocks look for."""
    return pd.DataFrame({
        'pid': ['p1', 'p2', 'p3', 'p4'], 'c80i70': ['i1', 'i1', 'i2', 'i3'],
        'c100i100': ['j1', 'j2', 'j3', 'j4'], 'plen': [100, 20, 300, 40],
        'pfam': ['PF1', None, 'PF2', 'PF3'], 'aravind': ['a', 'b', 'c', 'd'],
        'product': ['prod|1', 'prod2', 'prod3', 'prod4'],
        'organism': ['Org B', 'Org A', 'Org C', 'Org A'],
        'lineage': ['L2', 'L1', 'L1', 'L3'], 'classification': ['k2', 'k1', 'k3', 'k1'],
        'block_id': [1, 1, 2, 3], 'c80e3': ['X', 'X', 'Y', 'Z']})


def text(lines):
    return '\n'.join(lines)


def snap(result):
    """Frame (or mapping of frames) -> plain literals: columns, dtype strings and
    values, with missing values as None so they compare equal."""
    if isinstance(result, dict):
        return {name: snap(frame) for name, frame in result.items()}
    objects = result.astype(object).where(result.notna(), None)
    return {'columns': list(result.columns), 'dtypes': [str(t) for t in result.dtypes],
            'data': objects.to_dict('list')}


# ------------------------------------------------------------------ to_string
TO_STRING = {
    'default': ({},
     [
        '  : grp : sub : name       : score  : n  ',
        '0 : B   : x   : alpha|x    : 1.5    : 3  ',
        '1 : A   : y   : beta',
        'gamma :        : 10 ',
        '2 : B   : x   : delta      : 2.25   : 200',
        '3 : C   : z   : eps        : 3.0    : 4  ',
        '4 : A   : y   :            : 10.125 : 5  ',
    ]),
    'noindex': ({'index': False},
     [
        'grp : sub : name       : score  : n  ',
        'B   : x   : alpha|x    : 1.5    : 3  ',
        'A   : y   : beta',
        'gamma :        : 10 ',
        'B   : x   : delta      : 2.25   : 200',
        'C   : z   : eps        : 3.0    : 4  ',
        'A   : y   :            : 10.125 : 5  ',
    ]),
    'right': ({'align': 'right'},
     [
        '  : grp : sub :       name :  score :   n',
        '0 :   B :   x :    alpha|x :    1.5 :   3',
        '1 :   A :   y : beta',
        'gamma :        :  10',
        '2 :   B :   x :      delta :   2.25 : 200',
        '3 :   C :   z :        eps :    3.0 :   4',
        '4 :   A :   y :            : 10.125 :   5',
    ]),
    'float_format': ({'float_format': '%.2f'},
     [
        '  : grp : sub : name       : score : n  ',
        '0 : B   : x   : alpha|x    : 1.50  : 3  ',
        '1 : A   : y   : beta',
        'gamma :       : 10 ',
        '2 : B   : x   : delta      : 2.25  : 200',
        '3 : C   : z   : eps        : 3.00  : 4  ',
        '4 : A   : y   :            : 10.12 : 5  ',
    ]),
    'na_rep': ({'na_rep': '-'},
     [
        '  : grp : sub : name       : score  : n  ',
        '0 : B   : x   : alpha|x    : 1.5    : 3  ',
        '1 : A   : y   : beta',
        'gamma : -      : 10 ',
        '2 : B   : x   : delta      : 2.25   : 200',
        '3 : C   : z   : eps        : 3.0    : 4  ',
        '4 : A   : y   : -          : 10.125 : 5  ',
    ]),
    'sep': ({'sep': ' | '},
     [
        '  | grp | sub | name       | score  | n  ',
        '0 | B   | x   | alpha|x    | 1.5    | 3  ',
        '1 | A   | y   | beta',
        'gamma |        | 10 ',
        '2 | B   | x   | delta      | 2.25   | 200',
        '3 | C   | z   | eps        | 3.0    | 4  ',
        '4 | A   | y   |            | 10.125 | 5  ',
    ]),
    'insert': ({'insert': [(0, 'new', ['v', 'w', 'x', 'y', 'z'])]},
     [
        '  : new : grp : sub : name       : score  : n  ',
        '0 : v   : B   : x   : alpha|x    : 1.5    : 3  ',
        '1 : w   : A   : y   : beta',
        'gamma :        : 10 ',
        '2 : x   : B   : x   : delta      : 2.25   : 200',
        '3 : y   : C   : z   : eps        : 3.0    : 4  ',
        '4 : z   : A   : y   :            : 10.125 : 5  ',
    ]),
    'combined': ({'float_format': '%.2f', 'na_rep': '-', 'align': 'right', 'index': False},
     [
        'grp : sub :       name : score :   n',
        '  B :   x :    alpha|x :  1.50 :   3',
        '  A :   y : beta',
        'gamma :     - :  10',
        '  B :   x :      delta :  2.25 : 200',
        '  C :   z :        eps :  3.00 :   4',
        '  A :   y :          - : 10.12 :   5',
    ]),
}


@pytest.mark.parametrize('name', sorted(TO_STRING))
def test_to_string_golden(name):
    kwargs, expected = TO_STRING[name]
    assert rp.to_string(frame(), **kwargs) == text(expected)


# ------------------------------------------------------------------ to_blocks
# changed: 'multilevel_column_labels' used to render the key as a 1-tuple ('B',)
# although ('m', 'grp') is one column label; the key is now the scalar, as in
# default_header called on the column.
TO_BLOCKS = {
    'sum_nunique_unsorted_blocks': ('frame', {'groupby': 'grp', 'stats': {'n': 'sum', 'name': 'nunique'}, 'sortrows': {'n': False}, 'columns': ['name', 'score', 'n']},
     [
        '#',
        '# grp=B; n=203; name=2',
        '  : name    : score : n  ',
        '2 : delta   : 2.25  : 200',
        '0 : alpha|x : 1.5   : 3  ',
        '#',
        '# grp=A; n=15; name=1',
        '  : name       : score  : n ',
        '1 : beta',
        'gamma :        : 10',
        '4 :            : 10.125 : 5 ',
        '#',
        '# grp=C; n=4; name=1',
        '  : name : score : n',
        '3 : eps  : 3.0   : 4',
    ]),
    'sortblocks_and_sortrows': ('frame', {'groupby': 'grp', 'stats': {'n': 'max'}, 'sortrows': {'n': True}, 'sortblocks': {'n': True}, 'columns': ['name', 'score', 'n']},
     [
        '#',
        '# grp=C; n=4',
        '  : name : score : n',
        '3 : eps  : 3.0   : 4',
        '#',
        '# grp=A; n=10',
        '  : name       : score  : n ',
        '4 :            : 10.125 : 5 ',
        '1 : beta',
        'gamma :        : 10',
        '#',
        '# grp=B; n=200',
        '  : name    : score : n  ',
        '0 : alpha|x : 1.5   : 3  ',
        '2 : delta   : 2.25  : 200',
    ]),
    'custom_sep': ('frame', {'groupby': 'grp', 'stats': {'n': 'max'}, 'sortrows': {}, 'sep': '\n=====\n', 'columns': ['name', 'score', 'n']},
     [
        '#',
        '# grp=B; n=200',
        '  : name    : score : n  ',
        '0 : alpha|x : 1.5   : 3  ',
        '2 : delta   : 2.25  : 200',
        '=====',
        '#',
        '# grp=A; n=10',
        '  : name       : score  : n ',
        '1 : beta',
        'gamma :        : 10',
        '4 :            : 10.125 : 5 ',
        '=====',
        '#',
        '# grp=C; n=4',
        '  : name : score : n',
        '3 : eps  : 3.0   : 4',
    ]),
    'two_groupby_columns': ('frame', {'groupby': ['grp', 'sub'], 'stats': {'n': 'max'}, 'sortrows': {}, 'columns': ['name', 'n']},
     [
        '#',
        "# ('grp', 'sub')=('B', 'x'); n=200",
        '  : name    : n  ',
        '0 : alpha|x : 3  ',
        '2 : delta   : 200',
        '#',
        "# ('grp', 'sub')=('A', 'y'); n=10",
        '  : name       : n ',
        '1 : beta',
        'gamma : 10',
        '4 :            : 5 ',
        '#',
        "# ('grp', 'sub')=('C', 'z'); n=4",
        '  : name : n',
        '3 : eps  : 4',
    ]),
    'single_item_list_groupby': ('frame', {'groupby': ['grp'], 'stats': {'n': 'size'}, 'sortrows': {}, 'columns': ['n']},
     [
        '#',
        "# ('grp',)=('B',); n=2",
        '  : n  ',
        '0 : 3  ',
        '2 : 200',
        '#',
        "# ('grp',)=('A',); n=2",
        '  : n ',
        '1 : 10',
        '4 : 5 ',
        '#',
        "# ('grp',)=('C',); n=1",
        '  : n',
        '3 : 4',
    ]),
    'colsep_float_format_align': ('frame', {'groupby': 'grp', 'stats': {'n': 'max'}, 'sortrows': {}, 'columns': ['n'], 'colsep': ' ; ', 'float_format': '%.1f', 'align': 'right'},
     [
        '#',
        '# grp=B; n=200',
        '  ;   n',
        '0 ;   3',
        '2 ; 200',
        '#',
        '# grp=A; n=10',
        '  ;  n',
        '1 ; 10',
        '4 ;  5',
        '#',
        '# grp=C; n=4',
        '  ; n',
        '3 ; 4',
    ]),
    'literal_stat_value': ('frame', {'groupby': 'grp', 'stats': {'n': 7}, 'sortrows': {}, 'columns': ['n']},
     [
        '#',
        '# grp=B; n=7',
        '  : n  ',
        '0 : 3  ',
        '2 : 200',
        '#',
        '# grp=A; n=7',
        '  : n ',
        '1 : 10',
        '4 : 5 ',
        '#',
        '# grp=C; n=7',
        '  : n',
        '3 : 4',
    ]),
    'multilevel_column_labels': ('mi', {'groupby': ('m', 'grp'), 'stats': {('d', 'n'): 'sum'}, 'sortrows': {}, 'columns': [('d', 'name'), ('d', 'n')]},
     [
        '#',
        "# ('m', 'grp')=B; ('d', 'n')=203",
        "  : ('d', 'name') : ('d', 'n')",
        '0 : alpha|x       : 3         ',
        '2 : delta         : 200       ',
        '#',
        "# ('m', 'grp')=A; ('d', 'n')=15",
        "  : ('d', 'name') : ('d', 'n')",
        '1 : beta',
        'gamma    : 10        ',
        '4 :               : 5         ',
        '#',
        "# ('m', 'grp')=C; ('d', 'n')=4",
        "  : ('d', 'name') : ('d', 'n')",
        '3 : eps           : 4         ',
    ]),
}


@pytest.mark.parametrize('name', sorted(TO_BLOCKS))
def test_to_blocks_golden(name):
    which, kwargs, expected = TO_BLOCKS[name][0], TO_BLOCKS[name][1], TO_BLOCKS[name][2]
    df = frame() if which == 'frame' else multilevel_frame()
    assert rp.to_blocks(df, **kwargs) == text(expected)


# pins the genome-specific defaults; moves to NeighborhoodDF.to_blocks in the refactor
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
    assert rp.to_blocks(genome_frame()) == text(GENOME_DEFAULT_BLOCKS)


# ------------------------------------------------------------------- to_yaml
YAML_DOC = [
    '- parameters:',
    '    groupby: [grp]',
    '',
    '# -----',
    '- key: [C]',
    '  stats: |',
    '    | score | n |',
    '    | ----- | - |',
    '    | 3.0   | 4 |',  # changed: an int aggregate keeps its dtype (was 4.0)
    '  names: |',
    '    | name |',
    '    | ---- |',
    '    | eps  |',
    '  rows: |',
    '    | name | n |',
    '    | ---- | - |',
    '    | eps  | 4 |',
    '',
    '# -----',
    '- key: [A]',
    '  stats: |',
    '    | score  | n  |',
    '    | ------ | -- |',
    '    | 10.125 | 15 |',  # changed: was 15.0, int sum no longer upcast
    '  names: |',
    '    | name       |',
    '    | ---------- |',
    '    | beta gamma |',
    '    |            |',
    '  rows: |',
    '    | name       | n  |',
    '    | ---------- | -- |',
    '    | beta gamma | 10 |',
    '    |            | 5  |',
    '',
    '# -----',
    '- key: [B]',
    '  stats: |',
    '    | score | n   |',
    '    | ----- | --- |',
    '    | 2.25  | 203 |',  # changed: was 203.0, int sum no longer upcast
    '  names: |',
    '    | name     |',
    '    | -------- |',
    '    | delta    |',
    '    | alpha\\|x |',
    '  rows: |',
    '    | name     | n   |',
    '    | -------- | --- |',
    '    | delta    | 200 |',
    '    | alpha\\|x | 3   |',
]

YAML_KWARGS = dict(groupby='grp', columns=['name', 'n'], stats={'score': 'max', 'n': 'sum'},
                   sortrows={'n': False}, sortblocks={'n': True},
                   names=lambda block: block[['name']].drop_duplicates())


def test_to_yaml_golden():
    """A stats mapping (an int sum stays an int), a callable table, sortrows, sortblocks and columns."""
    assert yb.to_yaml(frame(), **YAML_KWARGS) == text(YAML_DOC) + '\n'


def test_to_yaml_buf_receives_the_same_text():
    buf = io.StringIO()
    assert yb.to_yaml(frame(), buf=buf, **YAML_KWARGS) is None
    assert buf.getvalue() == text(YAML_DOC) + '\n'


# Two grouping columns, no tables, no columns: everything but the groupby.
YAML_MINIMAL = [
    '- parameters:',
    '    groupby: [grp, sub]',
    '',
    '# -----',
    '- key: [B, x]',
    '  rows: |',
    '    | name     | score | n   |',
    '    | -------- | ----- | --- |',
    '    | alpha\\|x | 1.5   | 3   |',
    '    | delta    | 2.25  | 200 |',
    '',
    '# -----',
    '- key: [A, y]',
    '  rows: |',
    '    | name       | score  | n  |',
    '    | ---------- | ------ | -- |',
    '    | beta gamma |        | 10 |',
    '    |            | 10.125 | 5  |',
    '',
    '# -----',
    '- key: [C, z]',
    '  rows: |',
    '    | name | score | n |',
    '    | ---- | ----- | - |',
    '    | eps  | 3.0   | 4 |',
]


def test_to_yaml_minimal_golden():
    assert yb.to_yaml(frame(), groupby=['grp', 'sub']) == text(YAML_MINIMAL) + '\n'


# ------------------------------------------------------------------ read_yaml
# changed: the document read here (YAML_DOC) now writes the int sum as 4/15/203, so the
# 'n' / 'n_stats' columns it feeds are int64 rather than float64.
READ_YAML = {
    'merge_true': ({},
     {'columns': ['grp', 'name', 'n', 'name_names', 'score', 'n_stats'],
      'dtypes': ['str', 'str', 'int64', 'str', 'float64', 'int64'],
      'data': {'grp': ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C'],
               'name': ['beta gamma',
                        'beta gamma',
                        None,
                        None,
                        'delta',
                        'delta',
                        'alpha|x',
                        'alpha|x',
                        'eps'],
               'n': [10, 10, 5, 5, 200, 200, 3, 3, 4],
               'name_names': ['beta gamma',
                              None,
                              'beta gamma',
                              None,
                              'delta',
                              'alpha|x',
                              'delta',
                              'alpha|x',
                              'eps'],
               'score': [10.125, 10.125, 10.125, 10.125, 2.25, 2.25, 2.25, 2.25, 3.0],
               'n_stats': [15, 15, 15, 15, 203, 203, 203, 203, 4]}}),
    'merge_false': ({'merge': False},
     {'stats': {'columns': ['grp', 'score', 'n'],
                'dtypes': ['str', 'float64', 'int64'],
                'data': {'grp': ['C', 'A', 'B'],
                         'score': [3.0, 10.125, 2.25],
                         'n': [4, 15, 203]}},
      'names': {'columns': ['grp', 'name'],
                'dtypes': ['str', 'str'],
                'data': {'grp': ['C', 'A', 'A', 'B', 'B'],
                         'name': ['eps', 'beta gamma', None, 'delta', 'alpha|x']}},
      'rows': {'columns': ['grp', 'name', 'n'],
               'dtypes': ['str', 'str', 'int64'],
               'data': {'grp': ['C', 'A', 'A', 'B', 'B'],
                        'name': ['eps', 'beta gamma', None, 'delta', 'alpha|x'],
                        'n': [4, 10, 5, 200, 3]}}}),
    'merge_list': ({'merge': ['stats']},
     {'merged': {'columns': ['grp', 'score', 'n'],
                 'dtypes': ['str', 'float64', 'int64'],
                 'data': {'grp': ['C', 'A', 'B'],
                          'score': [3.0, 10.125, 2.25],
                          'n': [4, 15, 203]}},
      'names': {'columns': ['grp', 'name'],
                'dtypes': ['str', 'str'],
                'data': {'grp': ['C', 'A', 'A', 'B', 'B'],
                         'name': ['eps', 'beta gamma', None, 'delta', 'alpha|x']}},
      'rows': {'columns': ['grp', 'name', 'n'],
               'dtypes': ['str', 'str', 'int64'],
               'data': {'grp': ['C', 'A', 'A', 'B', 'B'],
                        'name': ['eps', 'beta gamma', None, 'delta', 'alpha|x'],
                        'n': [4, 10, 5, 200, 3]}}}),
    'merge_mapping': ({'merge': {'stats': 'n'}},
     {'merged': {'columns': ['grp', 'score', 'n'],
                 'dtypes': ['str', 'float64', 'int64'],
                 'data': {'grp': ['C', 'A', 'B'],
                          'score': [3.0, 10.125, 2.25],
                          'n': [4, 15, 203]}},
      'names': {'columns': ['grp', 'name'],
                'dtypes': ['str', 'str'],
                'data': {'grp': ['C', 'A', 'A', 'B', 'B'],
                         'name': ['eps', 'beta gamma', None, 'delta', 'alpha|x']}},
      'rows': {'columns': ['grp', 'name', 'n'],
               'dtypes': ['str', 'str', 'int64'],
               'data': {'grp': ['C', 'A', 'A', 'B', 'B'],
                        'name': ['eps', 'beta gamma', None, 'delta', 'alpha|x'],
                        'n': [4, 10, 5, 200, 3]}}}),
    'keep_with_drop': ({'keep': 'stats', 'drop': ['names']},
     {'merged': {'columns': ['grp', 'name', 'n', 'score', 'n_stats'],
                 'dtypes': ['str', 'str', 'int64', 'float64', 'int64'],
                 'data': {'grp': ['A', 'A', 'B', 'B', 'C'],
                          'name': ['beta gamma', None, 'delta', 'alpha|x', 'eps'],
                          'n': [10, 5, 200, 3, 4],
                          'score': [10.125, 10.125, 2.25, 2.25, 3.0],
                          'n_stats': [15, 15, 203, 203, 4]}},
      'stats': {'columns': ['grp', 'score', 'n'],
                'dtypes': ['str', 'float64', 'int64'],
                'data': {'grp': ['C', 'A', 'B'],
                         'score': [3.0, 10.125, 2.25],
                         'n': [4, 15, 203]}}}),
    'drop': ({'drop': ['stats', 'names']},
     {'columns': ['grp', 'name', 'n'],
      'dtypes': ['str', 'str', 'int64'],
      'data': {'grp': ['C', 'A', 'A', 'B', 'B'],
               'name': ['eps', 'beta gamma', None, 'delta', 'alpha|x'],
               'n': [4, 10, 5, 200, 3]}}),
    'keep_wins_over_drop': ({'keep': ['names'], 'drop': ['names', 'stats']},
     {'merged': {'columns': ['grp', 'name', 'n', 'name_names'],
                 'dtypes': ['str', 'str', 'int64', 'str'],
                 'data': {'grp': ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C'],
                          'name': ['beta gamma',
                                   'beta gamma',
                                   None,
                                   None,
                                   'delta',
                                   'delta',
                                   'alpha|x',
                                   'alpha|x',
                                   'eps'],
                          'n': [10, 10, 5, 5, 200, 200, 3, 3, 4],
                          'name_names': ['beta gamma',
                                         None,
                                         'beta gamma',
                                         None,
                                         'delta',
                                         'alpha|x',
                                         'delta',
                                         'alpha|x',
                                         'eps']}},
      'names': {'columns': ['grp', 'name'],
                'dtypes': ['str', 'str'],
                'data': {'grp': ['C', 'A', 'A', 'B', 'B'],
                         'name': ['eps', 'beta gamma', None, 'delta', 'alpha|x']}}}),
}


@pytest.mark.parametrize('name', sorted(READ_YAML))
def test_read_yaml_golden(name):
    kwargs, expected = READ_YAML[name]
    assert snap(yb.read_yaml(text(YAML_DOC) + '\n', **kwargs)) == expected


def test_read_yaml_reads_a_path_and_a_file_object(tmp_path):
    path = tmp_path / 'doc.yml'
    path.write_text(text(YAML_DOC) + '\n')
    expected = READ_YAML['drop'][1]
    assert snap(yb.read_yaml(str(path), drop=['stats', 'names'])) == expected
    with open(path) as fh:
        assert snap(yb.read_yaml(fh, drop=['stats', 'names'])) == expected


def test_read_yaml_unknown_table_raises():
    with pytest.raises(KeyError):
        yb.read_yaml(text(YAML_DOC) + '\n', keep='nope')


# ------------------------------------------------------------ md_table / read
MD_TABLE = '    | name       | score  | n   |\n    | ---------- | ------ | --- |\n    | alpha\\|x   | 1.5    | 3   |\n    | beta gamma |        | 10  |\n    | delta      | 2.25   | 200 |\n    | eps        | 3.0    | 4   |\n    |            | 10.125 | 5   |'
MD_TABLE_INDEXED = '| n   | name       |\n| --- | ---------- |\n| 3   | alpha\\|x   |\n| 10  | beta gamma |\n| 200 | delta      |\n| 4   | eps        |\n| 5   |            |'
MD_TABLE_EMPTY = '    | name | n |\n    | ---- | - |'


def test_md_table_golden():
    assert yb.md_table(frame()[['name', 'score', 'n']]) == MD_TABLE
    assert yb.md_table(frame()[['name', 'n']].set_index('n'), index=True, indent='') \
        == MD_TABLE_INDEXED
    assert yb.md_table(frame().iloc[:0][['name', 'n']]) == MD_TABLE_EMPTY


def test_md_table_escapes_pipes_and_flattens_newlines():
    out = yb.md_table(pd.DataFrame({'t': ['a|b', 'c\nd', 'e \n f']}), indent='')
    assert out.split('\n') == ['| t    |', '| ---- |', '| a\\|b |', '| c d  |', '| e f  |']


def test_read_md_table_round_trip():
    """The escaped pipe comes back as a pipe, the flattened newline as a space, the
    missing string as NaN (None here), numbers keep their numeric dtype."""
    got = yb.read_md_table(MD_TABLE)
    assert snap(got) == {'columns': ['name', 'score', 'n'],
     'dtypes': ['str', 'float64', 'int64'],
     'data': {'name': ['alpha|x', 'beta gamma', 'delta', 'eps', None],
              'score': [1.5, None, 2.25, 3.0, 10.125],
              'n': [3, 10, 200, 4, 5]}}
    got = yb.read_md_table(MD_TABLE_INDEXED)
    assert snap(got) == {'columns': ['n', 'name'],
     'dtypes': ['int64', 'str'],
     'data': {'n': [3, 10, 200, 4, 5], 'name': ['alpha|x', 'beta gamma', 'delta', 'eps', None]}}


def test_read_md_table_without_a_table_is_empty():
    assert yb.read_md_table('no table here').shape == (0, 0)
    assert yb.read_md_table('').shape == (0, 0)


# ------------------------------------------------------------------------ pad
PAD_SERIES = {'left': (['  a', 'bbb', ' cc'], ['..a', 'bbb', '.cc']),
     'right': (['a  ', 'bbb', 'cc '], ['a..', 'bbb', 'cc.']),
     'both': ([' a ', 'bbb', ' cc'], ['.a.', 'bbb', '.cc'])}


@pytest.mark.parametrize('side', ['left', 'right', 'both'])
def test_pad_series(side):
    s = pd.Series(['a', 'bbb', 'cc'])
    assert rpf.pad(s, side=side).tolist() == PAD_SERIES[side][0]
    assert rpf.pad(s, side=side, fillchar='.').tolist() == PAD_SERIES[side][1]


def test_pad_series_default_side_is_left():
    assert rpf.pad(pd.Series(['a', 'bbb', 'cc'])).tolist() == PAD_SERIES['left'][0]


def test_pad_converts_to_str_first():
    ints, floats = (['  1', ' 22', '333'], ['  1.5', '22.25'])
    assert rpf.pad(pd.Series([1, 22, 333])).tolist() == ints
    assert rpf.pad(pd.Series([1.5, 22.25])).tolist() == floats


def test_pad_dataframe_with_dict_side():
    """Columns missing from the dict default to 'left'. (A plain string `side` on a
    DataFrame is a known bug being fixed separately, so it is not pinned.)"""
    d = pd.DataFrame({'x': ['a', 'bbb', 'cc'], 'y': [1, 22, 333], 'z': ['q', 'rr', 's']})
    out = rpf.pad(d, side={'x': 'right', 'y': 'left'})
    assert out.to_dict('list') == {'x': ['a  ', 'bbb', 'cc '], 'y': ['  1', ' 22', '333'], 'z': [' q', 'rr', ' s']}
    assert [str(t) for t in out.dtypes] == ['str', 'str', 'str']
    assert rpf.pad(d, side={'x': 'both'}, fillchar='*').to_dict('list') == {'x': ['*a*', 'bbb', '*cc'], 'y': ['**1', '*22', '333'], 'z': ['*q', 'rr', '*s']}


# ------------------------------------------------------------------- to_color
TO_COLOR = {
    'default': ({},
     ['\x1b[38;5;033mA\x1b[m\x1b[38;5;168mC\x1b[m\x1b[38;5;127mD\x1b[m',
      '\x1b[38;5;033mW\x1b[m\x1b[38;5;124mK\x1b[m',
      '\x1b[39;5;000mX\x1b[m\x1b[39;5;000mB\x1b[m']),
    'bg': ({'style': 'bg'},
     ['\x1b[48;5;033mA\x1b[0m\x1b[48;5;168mC\x1b[0m\x1b[48;5;127mD\x1b[0m',
      '\x1b[48;5;033mW\x1b[0m\x1b[48;5;124mK\x1b[0m',
      '\x1b[49;5;000mX\x1b[0m\x1b[49;5;000mB\x1b[0m']),
    'dict_fg': ({'colors': {'A': '1', 'C': '22'}},
     ['\x1b[38;5;1mA\x1b[m\x1b[38;5;22mC\x1b[m\x1b[39;5;000mD\x1b[m',
      '\x1b[39;5;000mW\x1b[m\x1b[39;5;000mK\x1b[m',
      '\x1b[39;5;000mX\x1b[m\x1b[39;5;000mB\x1b[m']),
    'dict_background': ({'colors': {'A': '1', 'C': '22'}, 'style': 'background'},
     ['\x1b[48;5;1mA\x1b[0m\x1b[48;5;22mC\x1b[0m\x1b[49;5;000mD\x1b[0m',
      '\x1b[49;5;000mW\x1b[0m\x1b[49;5;000mK\x1b[0m',
      '\x1b[49;5;000mX\x1b[0m\x1b[49;5;000mB\x1b[0m']),
    'padding_left': ({'padding': 'left'},
     ['\x1b[38;5;033mA\x1b[m\x1b[38;5;168mC\x1b[m\x1b[38;5;127mD\x1b[m',
      '\x1b[39;5;000m \x1b[m\x1b[38;5;033mW\x1b[m\x1b[38;5;124mK\x1b[m',
      '\x1b[39;5;000m \x1b[m\x1b[39;5;000mX\x1b[m\x1b[39;5;000mB\x1b[m']),
    'padding_right_fillchar': ({'padding': 'right', 'fillchar': '-', 'colors': {'A': '1'}},
     ['\x1b[38;5;1mA\x1b[m\x1b[39;5;000mC\x1b[m\x1b[39;5;000mD\x1b[m',
      '\x1b[39;5;000mW\x1b[m\x1b[39;5;000mK\x1b[m\x1b[39;5;000m-\x1b[m',
      '\x1b[39;5;000mX\x1b[m\x1b[39;5;000mB\x1b[m\x1b[39;5;000m-\x1b[m']),
    'fallback_color': ({'colors': '7'},
     ['\x1b[39;5;7mA\x1b[m\x1b[39;5;7mC\x1b[m\x1b[39;5;7mD\x1b[m',
      '\x1b[39;5;7mW\x1b[m\x1b[39;5;7mK\x1b[m',
      '\x1b[39;5;7mX\x1b[m\x1b[39;5;7mB\x1b[m']),
}


@pytest.mark.parametrize('name', sorted(TO_COLOR))
def test_to_color_series(name):
    """Escape codes are pinned through repr(). The default colour table is the
    shipped aminoacids.yml; an unmapped character gets the 000 'default' code."""
    kwargs, expected = TO_COLOR[name]
    out = rpf.to_color(pd.Series(['ACD', 'WK', 'XB']), **kwargs)
    assert repr(out.tolist()) == repr(expected)


# ------------------------------------------------------- optimize_memory_usage
@pytest.mark.filterwarnings('ignore::Warning')   # select_dtypes('object') deprecation under pandas 3
def test_optimize_memory_usage_dtypes_and_values():
    """Ints and floats are downcast; text (object or str) is made a category when at
    most half of its values are distinct, and left alone otherwise."""
    df = pd.DataFrame({
        'i': [1, 2, 3, 4], 'big': [1, 2, 3, 70000], 'neg': [-1, 5, 0, 2],
        'f': [1.5, 2.5, np.nan, 4.0], 'lo': ['a', 'a', 'b', 'b'], 'hi': ['w', 'x', 'y', 'z'],
        'o_lo': pd.Series(['a', 'a', 'b', 'b'], dtype=object),
        'o_hi': pd.Series(['w', 'x', 'y', 'z'], dtype=object),
        'b': [True, False, True, True]})
    before = df.copy()
    out = rpf.optimize_memory_usage(df)
    assert {c: str(t) for c, t in out.dtypes.items()} == {
        'i': 'int8', 'big': 'int32', 'neg': 'int8', 'f': 'float32', 'lo': 'category', 'hi': 'str',
        'o_lo': 'category', 'o_hi': 'object', 'b': 'bool'}
    pd.testing.assert_frame_equal(out.astype(object), df.astype(object))
    pd.testing.assert_frame_equal(df, before)    # the input is not modified


# ------------------------------------------------------------ display options
@pytest.fixture
def display_options():
    names = ['max_rows', 'max_columns', 'max_colwidth', 'width']
    saved = {n: getattr(pd.options.display, n) for n in names}
    yield
    for n, v in saved.items():
        setattr(pd.options.display, n, v)


def test_show_display_options_shape():
    out = rpf.show_display_options()
    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == ['option', 'value']
    assert out['option'].is_monotonic_increasing        # dir() order
    shown = dict(zip(out['option'], out['value']))
    assert shown['max_rows'] == pd.options.display.max_rows
    assert shown['width'] == pd.options.display.width


def test_print_everything_defaults(display_options):
    rpf.print_everything()
    d = pd.options.display
    assert (d.max_rows, d.max_columns, d.max_colwidth, d.width) == (
        10000000000, 1000000000, 1000000000, 100000)


def test_print_everything_arguments_and_verbose(display_options, capsys):
    rpf.print_everything(max_rows=42, max_columns=7, max_colwidth=33, width=99, verbose=True)
    d = pd.options.display
    assert (d.max_rows, d.max_columns, d.max_colwidth, d.width) == (42, 7, 33, 99)
    captured = capsys.readouterr()
    assert captured.out == ''
    assert captured.err == 'Pandas options have been reset!\n'


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-q']))
