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
survive; reading is a safe YAML load plus a split on unescaped pipes, so no
structure is parsed by hand. Padding is cosmetic: surrounding whitespace is
stripped, so a hand-edited table with ragged pipes still parses.

Writing does not go through to_blocks. Escaping and string conversion happen once
for the whole frame, statistics are aggregated for every block in a single pass,
and each block is cut with an arrow `take`, so no per-block DataFrame is built and
no Python-level pass is made over the rows. That is worth roughly 3.5x on a frame
of 100k rows in 500 groups, rising to about 6x at 5000 groups, where the per-block
overhead had dominated. pyarrow is required here rather than optional.

Reading combines the tables by the key they are stored under: the same key across
blocks has its rows concatenated, and the resulting frames are merged on the group
columns.
"""

import io
import os
import re
from collections.abc import Mapping

import yaml
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc

from rotifer.pandas import _py
from rotifer.core import functions as rcf

CONFIG = rcf.loadConfig(__name__, defaults = {
    'indent': 4,
    'colsep': ' | ',
    'header': True,
    'ndf': {
        'annotation': lambda x: pd.DataFrame({'pfam':x.pfam.mode(), 'curated':x.pfam.mode() }),
        'stats': {
            'superkingdom': ('superkingdom', 'nunique'),
            'lineage': ('lineage', 'nunique'),
            'phylum': ('phylum', 'nunique'),
            'taxid': ('taxid', 'nunique'),
            'c80i70': ('c80i70', 'nunique'),
            'pid': ('pid', 'nunique'),
        },
        'dist': lambda x: x.groupby('pfam').c80i70.nunique().reset_index(),
        'columns': ['pid', 'c80i70', 'c100i100', 'plen', 'pfam', 'aravind', 'product', 'organism', 'lineage', 'classification'],
        'sortblocks': {
            'phylum': False,
            'c80i70': False
        },
        'sortrows': {
            'qc80e3': True,
            'superkingdom': True,
            'phylum': True,
            'lineage': True,
            'classification': True,
            'c80i70': True
        },
    },
})
INDENT = ' ' * CONFIG['indent']
MERGED = 'merged'                       # read_yaml returns the joined frame here
RESERVED = ('key', 'parameters', MERGED)
SEPARATOR = '# ' + '-' * 5

# libyaml where it was built, which reads these documents about 13x faster than
# the pure-Python scanner -- with block scalars this large, that is most of the
# cost of reading. Both loaders produce the same object.
LOADER = getattr(yaml, 'CSafeLoader', yaml.SafeLoader)

# A pipe that is not backslash-escaped: as a separator string for read_csv, and
# compiled for splitting rows directly.
_CELL_SEP = r'\s*(?<!\\)\|\s*'
_CELL_RE = re.compile(r'(?<!\\)\|')
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


def _escape_columns(df):
    """DataFrame -> one arrow string array per column, escaped and ready to pad.

    Pipes are escaped and newlines flattened *before* padding: escaping afterwards
    would leave every affected column one character short, and an embedded newline
    would terminate the YAML block scalar. Doing this once for a whole frame, then
    slicing blocks out of the result, avoids repeating it for every block."""
    out = []
    for c in df.columns:
        arr = pa.Array.from_pandas(df[c].astype('string').fillna(''), type=pa.string())
        arr = pc.replace_substring(arr, '|', r'\|')
        out.append(pc.replace_substring_regex(arr, r'\s*\n\s*', ' '))
    return out


def _md_arrow(cols, names, indent=INDENT):
    """Render pre-escaped arrow columns as a markdown table.

    Padding and the surrounding pipes are both applied by arrow kernels, so no
    Python-level pass is made over the rows."""
    padded, headers = [], []
    for arr, name in zip(cols, names):
        width = max(pc.max(pc.utf8_length(arr)).as_py() or 0, len(name))
        padded.append(pc.utf8_rpad(arr, width, padding=' '))
        headers.append(name.ljust(width))
    head = f'{indent}| ' + ' | '.join(headers) + ' |'
    rule = f'{indent}| ' + ' | '.join('-' * len(h) for h in headers) + ' |'
    if not padded or len(padded[0]) == 0:
        return f'{head}\n{rule}'
    body = pc.binary_join_element_wise(*padded, pa.scalar(' | '))
    body = pc.binary_join_element_wise(pa.scalar(f'{indent}| '), body, pa.scalar(' |'),
                                       pa.scalar(''))
    return '\n'.join([head, rule] + body.to_pylist())


def md_table(df, indent=INDENT, index=False):
    """Render a DataFrame as a markdown table, indented for a YAML block scalar."""
    if index:
        df = df.reset_index()
    return _md_arrow(_escape_columns(df), [str(c) for c in df.columns], indent=indent)


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


def _row_rank(df, sortrows):
    """Positions ranked by the sort columns. Ranking once and reordering each
    block's positions keeps rows sorted *within* a block without disturbing the
    order the blocks themselves appear in."""
    cols = [c for c in sortrows if c in df]
    if not cols:
        return None
    order = (df.reset_index(drop=True)
               .sort_values(cols, ascending=[sortrows[c] for c in cols], kind='stable')
               .index.to_numpy())
    rank = np.empty(len(df), dtype=np.int64)
    rank[order] = np.arange(len(df))
    return rank


def _group_agg(grouped, spec, table):
    """Aggregate one table's spec for every block in a single pass."""
    missing = [c for c in spec if c not in grouped.obj.columns]
    if missing:
        raise KeyError(f'{table} column {missing[0]!r} is not in the DataFrame')
    return grouped.agg(**{str(c): (c, f) for c, f in spec.items()})


def _split_tables(grouped, tables):
    """Sort each table spec into the two ways a table can be produced.

    A mapping is an aggregation, {column: aggregation}, evaluated once for every
    block; a callable receives each block as a DataFrame and returns its table.
    Aggregations are much the cheaper of the two, since they never materialise a
    per-block DataFrame."""
    aggregated, callbacks = {}, {}
    for table, spec in tables.items():
        if table in RESERVED:
            raise ValueError(f'table name {table!r} is reserved: {RESERVED[0]!r} and '
                             f'{RESERVED[1]!r} carry block metadata, and {MERGED!r} is '
                             f'where read_yaml returns the joined frame')
        if isinstance(spec, Mapping):
            aggregated[table] = _group_agg(grouped, spec, table)
        elif callable(spec):
            callbacks[table] = spec
        else:
            raise TypeError(f'table {table!r} must be a mapping of column to '
                            f'aggregation, or a callable, not {type(spec).__name__}')
    return aggregated, callbacks


def to_yaml(df, groupby='c80e3', buf=None, columns=None, sortrows=None,
            sortblocks=None, sample=None, name='rows', sep='\n\n',
            **tables):
    """Write `df` as the YAML block document described above.

    Column text is escaped and converted once for the whole frame, per-block
    aggregations are evaluated in a single pass, and each block is cut with an
    arrow `take`, so the default path builds no per-block DataFrame.

    columns    : columns of the main table. Defaults to everything but the
                 grouping columns, which the reader restores from the block key.
    name       : name of that main table, `rows` by default.
    sortrows   : {column: ascending} ordering rows within a block. Applied per
                 block, so it never changes the order of the blocks themselves.
    sortblocks : {column: ascending} ordering the blocks, addressing a grouping
                 column or any column produced by an aggregation table below.
    sample     : rows per block, sampled without replacement when a block is
                 larger.

    Every remaining keyword argument defines one further table, the keyword
    naming it in the document. Two forms are accepted:

      stats={'lineage': 'nunique', 'pid': 'nunique'}
          A mapping of column to aggregation. Evaluated for all blocks at once
          with a single groupby.agg, producing a one-row table per block. This is
          much the cheaper form, and `stats` is the conventional name for it.

      lineages=lambda block: block[['lineage']].drop_duplicates()
          A callable receiving each block as a DataFrame and returning the table
          to render. Flexible, but it materialises a DataFrame per block.

    So the whole shape of a block is caller-defined:

        to_yaml(df, groupby='c80e3',
                stats={'lineage': 'nunique'},
                span={'plen': 'max'},
                products=lambda b: b['product'].value_counts().reset_index())

    writes `stats`, `span` and `products` tables into every block, in that order,
    followed by the main table. The tables need have nothing in common: they are
    never joined to one another here. Columns produced by any aggregation table
    are visible to sortblocks, so blocks can be ordered by a statistic, and if two
    tables define the same column name the one named first decides. Names are
    yours to choose except `key` and `parameters`, which carry block metadata.

    Rows are brought back together by read_yaml, not here -- it keys them on the
    grouping columns, which it restores into every table from the block key.

    Columns named in sortrows or sortblocks that are not present are ignored, so a
    single ordering can be reused across frames that differ in shape.
    """
    labels = list(groupby) if isinstance(groupby, (list, tuple)) else [groupby]
    head = f'- parameters:\n    groupby: {_flow(labels)}'

    grouped = df.groupby(groupby, sort=False)
    aggregated, callbacks = _split_tables(grouped, tables)
    rank = _row_rank(df, sortrows) if sortrows else None

    body = df.filter(columns) if columns else df.drop(columns=labels, errors='ignore')
    cols = _escape_columns(body)                  # once, not per block
    names = [str(c) for c in body.columns]
    positions = {k: np.asarray(v) for k, v in grouped.indices.items()}

    # Group keys in first-appearance order. sortblocks may name a grouping column
    # or a column of any aggregation table; those are collected one column at a
    # time rather than by joining the tables, which are free to share nothing with
    # each other. Where two tables define the same column, the one named first
    # wins. All the frames come from the same grouping, so they align positionally.
    index = next(iter(aggregated.values())).index if aggregated else grouped.size().index
    keys = list(index)
    if sortblocks:
        frame = pd.DataFrame(index=index).reset_index()
        for col in sortblocks:
            if col in frame.columns:
                continue
            for table in aggregated.values():
                if col in table.columns:
                    frame[col] = table[col].to_numpy()
                    break
        by = [c for c in sortblocks if c in frame.columns]
        if by:
            frame = frame.sort_values(by, ascending=[sortblocks[c] for c in by], kind='stable')
            keys = [keys[i] for i in frame.index]

    rng = np.random.default_rng() if sample else None
    parts = [head]
    for key in keys:
        idx = positions[key]
        if rank is not None:
            idx = idx[np.argsort(rank[idx], kind='stable')]
        if sample and sample < len(idx):
            idx = np.sort(rng.choice(idx, sample, replace=False))
        taken = pa.array(idx)

        keylist = list(key) if isinstance(key, tuple) else [key]
        block = [SEPARATOR, f'- key: {_flow(keylist)}']
        for table in tables:                      # in the order the caller named them
            if table in aggregated:
                row = aggregated[table].loc[key]
                rendered = _md_arrow([pa.array([str(_py(v))]) for v in row],
                                     [str(c) for c in aggregated[table].columns])
            else:
                rendered = md_table(callbacks[table](df.take(idx)))
            block.append(f'  {table}: |\n' + rendered)
        block.append(f'  {name}: |\n' + _md_arrow([c.take(taken) for c in cols], names))
        parts.append('\n'.join(block))

    text = sep.join(parts) + '\n'

    if buf is None:
        return text
    if isinstance(buf, (str, os.PathLike)):
        with open(buf, 'w') as fh:
            fh.write(text)
        return None
    buf.write(text)
    return None


# --------------------------------------------------------------------- reading
def read_yaml(source, merge=True, drop=None, keep=None, how='outer'):
    """Read a document written by to_yaml.

    source : YAML text, a path, or a file-like object.
    how    : join used when merging tables held under different keys.
    merge  : what to fold together, and on what.

        True          the default: every table is merged on the grouping columns.
        None, False   nothing is merged.
        ['a', 'b']    only these tables are merged, on the grouping columns.
        {'a': 'pid'}  only these tables are merged, each joined on the columns
        {'a': [...]}  given for it rather than on the grouping columns.

    drop   : tables to leave out entirely.
    keep   : tables to return in their original, unmerged form.

    drop and keep each take one table name or an array-like of them, and **keep
    wins over drop**: a table named in both is returned. keep is also independent
    of merging, so a table can be folded into the merged frame and handed back
    untouched at the same time.

    Everything the document holds comes back unless it was dropped: tables left
    out of a partial merge are returned as they are, alongside the merged frame,
    which is keyed 'merged'. The return is that bare DataFrame when it is the only
    thing to return -- the usual case under the default merge=True -- and a
    {name: DataFrame} mapping otherwise. Naming a table that the document does not
    hold raises KeyError rather than passing unnoticed.

    The grouping columns come from the document's parameters element, so a
    document always describes itself. Their values are injected into each table
    from the block `key`, which is why they are not repeated on every row, and
    they are what the tables are keyed on when recombined.

    Tables are recombined by the key they are stored under, which is what makes
    two of them equivalent -- a `stats` table is never concatenated onto a `span`
    table just because their columns happen to match. Same key across blocks means
    the rows are concatenated; the resulting frames are then merged on the group
    columns, and a column reaching the result from more than one table is suffixed
    with that table's key, so aggregating `pid` under `stats` lands as `pid_stats`
    beside the row-level `pid`. The widest table keeps the unsuffixed names.
    """
    if hasattr(source, 'read'):
        text = source.read()
    elif isinstance(source, (str, os.PathLike)) and '\n' not in str(source) and os.path.exists(source):
        with open(source) as fh:
            text = fh.read()
    else:
        text = source

    doc = yaml.load(text, Loader=LOADER)
    if not isinstance(doc, list):
        raise ValueError('not a to_yaml document: top level is not a YAML sequence')

    params = {}
    for element in doc:
        if isinstance(element, dict) and 'parameters' in element:
            params = element['parameters'] or {}
            break
    groupcols = params.get('groupby', [])
    groupcols = list(groupcols) if isinstance(groupcols, (list, tuple)) else [groupcols]

    # keep wins over drop, so a table named in both survives to be returned.
    keep = [keep] if isinstance(keep, str) else list(keep or ())
    dropped = set([drop] if isinstance(drop, str) else list(drop or ())) - set(keep)

    # Gather each table's text as one string. The key a table is stored under is
    # what makes two of them equivalent, not their columns: a `stats` table and a
    # `span` table are different things whether or not they share column names.
    # Every block wrote the table with the same header on the first line and the
    # ---- rule on the second, so one header is kept and later blocks contribute
    # only their data rows.
    # Rows are cut into cells here and never rebuilt into text, so the frames are
    # assembled from lists at the end rather than re-parsed from a string. The
    # grouping columns ride along as extra cells on each row, so a row carries its
    # own key and nothing has to be realigned afterwards; a table already naming a
    # grouping column keeps the one it stores.
    rows, header, missing, escaped = {}, {}, {}, {}
    for element in doc:
        if not isinstance(element, dict) or 'parameters' in element:
            continue
        key = element.get('key', [])
        key = list(key) if isinstance(key, (list, tuple)) else [key]
        values = dict(zip(groupcols, key))
        for name, table in element.items():
            if name in RESERVED or name in dropped or not isinstance(table, str):
                continue
            head, _, rest = table.partition('\n')
            _, _, body = rest.partition('\n')           # the second line is the rule
            if name not in rows:
                stored = [c.strip() for c in head.split('|')][1:-1]
                missing[name] = [c for c in groupcols if c not in stored]
                header[name] = stored + missing[name]
                rows[name], escaped[name] = [], False
            extra = [_py(values.get(c)) for c in missing[name]]
            # An escaped pipe is rare, so only pay for the regex where one occurs.
            here = r'\|' in body
            escaped[name] |= here
            for line in body.split('\n'):
                if not line.strip():
                    continue
                cells = _CELL_RE.split(line) if here else line.split('|')
                rows[name].append([c.strip() or None for c in cells[1:-1]] + extra)

    if not rows:
        raise ValueError('no tables found in document')

    frames = {}
    for name, table in rows.items():
        frame = pd.DataFrame(table, columns=header[name])
        if escaped[name]:
            frame = frame.apply(lambda s: s.str.replace(r'\|', '|', regex=False)
                                if s.dtype == object else s)
        # read_csv used to type the columns; do it here instead, per column and
        # all-or-nothing, so an identifier that merely looks numeric is left alone
        # unless the whole column is numeric.
        for c in frame.columns:
            if frame[c].dtype == object:
                try:
                    frame[c] = pd.to_numeric(frame[c])
                except (ValueError, TypeError):
                    pass
        front = [c for c in groupcols if c in frame.columns]
        frames[name] = frame[front + [c for c in frame.columns if c not in front]]
    unknown = [name for name in keep if name not in frames]
    if unknown:
        raise KeyError(f'no table named {unknown[0]!r} in the document '
                       f'(have {sorted(frames)})')
    if merge is None or merge is False:
        return frames

    # Which tables to merge, and on what. True takes them all on the grouping
    # columns; a list picks some of them, still on the grouping columns; a mapping
    # picks some and says which columns each is joined on.
    if merge is True:
        # No order was given, so take the document in reverse: to_yaml appends the
        # main table after the aggregations, which puts it first here, where it
        # keeps its column names and sets the row granularity.
        joins = {name: groupcols for name in reversed(frames)}
    elif isinstance(merge, Mapping):
        joins = {name: [c] if isinstance(c, str) else list(c) for name, c in merge.items()}
    else:
        joins = {name: groupcols for name in ([merge] if isinstance(merge, str) else merge)}

    unknown = [name for name in joins if name not in frames]
    if unknown:
        raise KeyError(f'no table named {unknown[0]!r} in the document '
                       f'(have {sorted(frames)})')
    if not joins:
        raise ValueError('merge names no tables')

    # Tables are merged in the order merge gave them: the first is the frame the
    # rest are joined onto, so it is the one that keeps its column names unsuffixed
    # and sets the row granularity. A column arriving from a later table is
    # suffixed with the name that table has in the document, so the merged frame
    # stays self-describing.
    names = list(joins)
    merged = frames[names[0]]
    for name in names[1:]:
        right = frames[name]
        on = [c for c in joins[name] if c in merged.columns and c in right.columns]
        if not on:
            raise ValueError(f'table {name!r} shares none of {joins[name]} with the '
                             f'tables merged before it')
        merged = merged.merge(right, on=on, how=how, suffixes=('', f'_{name}'))

    # Everything the document held comes back: the merged frame, whatever a
    # partial merge left out, and any original asked for by keep. A lone merged
    # frame is handed back bare rather than wrapped in a mapping of one.
    result = {MERGED: merged}
    result.update({name: frame for name, frame in frames.items() if name not in joins})
    result.update({name: frames[name] for name in keep})
    return merged if len(result) == 1 else result


read_blocks = read_yaml  # alias
