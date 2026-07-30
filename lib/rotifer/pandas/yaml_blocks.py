#!/usr/bin/env python3
"""Block output as hand-editable YAML carrying markdown tables.

A document is a YAML sequence. The first element holds the file-level parameters,
notably the grouping columns, which are therefore stored once rather than repeated
on every block. Every later element is one block: a `key` naming the group, plus
any number of further keys, each holding a markdown table as a literal scalar.

    - parameters:
        groupby: [c80e3]

    # ----------------------------------------------------------------------
    - key: [A]
      stats: |
        | lineage | pid |
        | ------- | --- |
        | 1       | 2   |
      rows: |
        | pid | plen | product       |
        | --- | ---- | ------------- |
        | p1  | 10   | prod one      |
        | p2  | 200  | has \\| a pipe |

Writing formats the text directly so the comments, blank lines and column padding
survive; reading is `yaml.safe_load` plus `pandas.read_csv`, so no structure is
parsed by hand. Padding is cosmetic: the separator regex strips surrounding
whitespace, so a hand-edited table with ragged pipes still parses.

Reading combines the tables it finds: those sharing a column signature have their
rows concatenated, and the resulting frames are merged on the group columns.
"""

import io
import os
import re

import yaml
import pandas as pd

from rotifer.pandas import to_string, to_blocks, _group_key, _py, _resolve_stats

INDENT = '    '
RESERVED = ('key', 'parameters')
SEPARATOR = '# ' + '-' * 5

# A pipe that is not backslash-escaped, plus any padding around it.
_CELL_SEP = r'\s*(?<!\\)\|\s*'
_RULE_RE = re.compile(r'^\s*\|[\s:|-]+\|\s*$')


# --------------------------------------------------------------------- helpers
def _flow(value):
    """Render a scalar or sequence as a compact YAML flow scalar, letting PyYAML
    handle quoting so values containing commas or brackets stay valid."""
    if isinstance(value, tuple):
        value = list(value)
    if isinstance(value, list):
        value = [_py(v) for v in value]
    else:
        value = _py(value)
    return yaml.safe_dump(value, default_flow_style=True, width=10 ** 6,
                          allow_unicode=True).strip()


def md_table(df, indent=INDENT, index=False, **kwargs):
    """Render a DataFrame as a markdown table, indented for a YAML block scalar.

    Pipes are escaped and newlines flattened before padding -- escaping afterwards
    would leave every affected column one character short, and an embedded newline
    would terminate the block scalar."""
    body = df.copy()
    for c in body.columns:
        s = body[c]
        if s.dtype == object or str(s.dtype).startswith('string'):
            body[c] = (s.astype('string')
                        .str.replace('|', r'\|', regex=False)
                        .str.replace(r'\s*\n\s*', ' ', regex=True))
    lines = to_string(body, sep=' | ', index=index, **kwargs).split('\n')
    if not len(body):
        lines = lines[:1]                   # to_string leaves a trailing blank row
    rule = ' | '.join('-' * len(c) for c in lines[0].split(' | '))
    return '\n'.join(f'{indent}| {ln} |' for ln in [lines[0], rule] + lines[1:])


def read_md_table(text):
    """Parse one markdown table. pandas does the tokenizing; the separator regex
    only tells it which pipes delimit cells."""
    rows = [ln for ln in text.split('\n') if ln.strip().startswith('|')]
    rows = [ln for ln in rows if not _RULE_RE.match(ln)]
    if not rows:
        return pd.DataFrame()
    frame = pd.read_csv(io.StringIO('\n'.join(rows)), sep=_CELL_SEP, engine='python',
                        skipinitialspace=True, dtype=None).iloc[:, 1:-1]
    return frame.apply(lambda s: s.str.replace(r'\|', '|', regex=False)
                       if s.dtype == object else s)


# --------------------------------------------------------------------- writing
def yaml_header(df, groupby='c80e3', stats=None, tables=None, name='stats', **kwargs):
    """to_blocks `header` callback: emits the block's metadata tables.

    `stats` becomes one markdown table; `tables` adds more as
    {name: DataFrame or callable(df) -> DataFrame}."""
    parts = []
    stats = _resolve_stats(df, stats or {})
    if stats:
        parts.append(f'  {name}: |\n' + md_table(pd.DataFrame([stats])))
    for label, spec in (tables or {}).items():
        frame = spec(df) if callable(spec) else spec
        parts.append(f'  {label}: |\n' + md_table(frame))
    return '\n'.join(parts)


def yaml_block_rows(df, groupby='c80e3', header=yaml_header, colsep=None,
                    columns=None, sample=None, stats=None, sortrows=None,
                    name='rows', **kwargs):
    """to_blocks `apply` callback: emits one block element of the sequence."""
    key = _group_key(df, groupby)   # _group_key normalizes the groupby itself
    key = list(key) if isinstance(key, tuple) else [key]

    statDict = _resolve_stats(df, stats or {})
    headerStr = header(df, groupby=groupby, stats=statDict, **kwargs)

    block = df.copy()
    if sample and sample < len(block):
        block = block.sample(sample)
    if sortrows:
        sortcols = [c for c in sortrows if c in block]
        if sortcols:
            block = block.sort_values(sortcols, ascending=[sortrows[c] for c in sortcols])
    if columns:
        block = block.filter(columns)

    parts = [SEPARATOR, f'- key: {_flow(key)}']
    if headerStr:
        parts.append(headerStr)
    parts.append(f'  {name}: |\n' + md_table(block))
    return pd.DataFrame({'blocks': ['\n'.join(parts)], **statDict})


def to_yaml(df, groupby='c80e3', buf=None, header=yaml_header,
            apply=yaml_block_rows, sep='\n\n', **kwargs):
    """Write `df` as the YAML block document described above.

    A blank line precedes each block separator; nothing else is padded out, so the
    envelope costs about 0.1% over the plain to_blocks rendering."""
    # Recorded for the reader only; `groupby` itself is passed on untouched, since
    # pandas reads a tuple as a single key rather than a list of columns.
    labels = list(groupby) if isinstance(groupby, (list, tuple)) else [groupby]
    head = f'- parameters:\n    groupby: {_flow(labels)}'

    blocks = to_blocks(df, groupby=groupby, header=header, apply=apply, sep=sep, **kwargs)
    text = head + sep + blocks + '\n'

    if buf is None:
        return text
    if isinstance(buf, (str, os.PathLike)):
        with open(buf, 'w') as fh:
            fh.write(text)
        return None
    buf.write(text)
    return None


# --------------------------------------------------------------------- reading
def read_yaml(source, merge=False, drop=None, how='outer', suffixes=('_stats', ''),
              groupby=None):
    """Read a document written by to_yaml.

    source   : YAML text, a path, or a file-like object.
    merge    : True folds every table into one DataFrame, merging on the group
               columns; False returns {name: DataFrame}.
    drop     : table names to discard before combining, e.g. drop=['stats'].
    how      : join used when merging tables of differing signatures.
    suffixes : disambiguates columns present in more than one table -- the default
               tags the block-level ones, since stats are conventionally named
               after the data columns they summarise.
    groupby  : overrides the grouping columns recorded in the parameters element.

    Group columns are injected into each table from the block `key`, so they are
    not repeated on every row of the document.
    """
    if hasattr(source, 'read'):
        text = source.read()
    elif isinstance(source, (str, os.PathLike)) and '\n' not in str(source) and os.path.exists(source):
        with open(source) as fh:
            text = fh.read()
    else:
        text = source

    doc = yaml.safe_load(text)
    if not isinstance(doc, list):
        raise ValueError('not a to_yaml document: top level is not a YAML sequence')

    params = {}
    for element in doc:
        if isinstance(element, dict) and 'parameters' in element:
            params = element['parameters'] or {}
            break
    groupcols = groupby if groupby is not None else params.get('groupby', [])
    groupcols = list(groupcols) if isinstance(groupcols, (list, tuple)) else [groupcols]

    drop = set(drop or ())
    collected = []
    for element in doc:
        if not isinstance(element, dict) or 'parameters' in element:
            continue
        key = element.get('key', [])
        key = list(key) if isinstance(key, (list, tuple)) else [key]
        for name, table in element.items():
            if name in RESERVED or name in drop or not isinstance(table, str):
                continue
            frame = read_md_table(table)
            for col, value in zip(groupcols, key):
                if col not in frame.columns:    # injected, not stored per row
                    frame.insert(0, col, value)
            collected.append((name, frame))

    if not collected:
        raise ValueError('no tables found in document')

    # same column signature -> concatenate rows; keep the first name as the label
    bysig = {}
    for name, frame in collected:
        bysig.setdefault(tuple(frame.columns), (name, []))[1].append(frame)
    frames = {name: pd.concat(parts, ignore_index=True) for name, parts in bysig.values()}
    if not merge:
        return frames

    # different signatures -> merge on the group columns, widest table first so the
    # per-row tables keep their granularity and block-level values broadcast onto them
    ordered = sorted(frames.values(), key=len, reverse=True)
    on = [c for c in groupcols if all(c in f.columns for f in ordered)]
    if not on:
        raise ValueError(f'no group columns shared by every table (looked for {groupcols})')
    merged = ordered[0]
    for right in ordered[1:]:
        merged = merged.merge(right, on=on, how=how, suffixes=(suffixes[1], suffixes[0]))
    return merged


read_blocks = read_yaml  # alias
