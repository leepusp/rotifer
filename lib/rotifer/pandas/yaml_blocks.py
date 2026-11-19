#!/usr/bin/env python3
"""Block output as hand-editable YAML carrying markdown tables.

A document is a YAML sequence. The first element holds the file-level parameters,
notably the grouping columns, which are therefore stored once rather than repeated
on every block. Every later element is one block: a `key` naming the group, plus
any number of further keys, each holding a markdown table as a literal scalar.

    - parameters:
        groupby: [grp]

    # ----------------------------------------------------------------------
    - key: [A]
      stats: |
        | kinds | id |
        | ----- | -- |
        | 1     | 2  |
      rows: |
        | id | size | product       |
        | -- | ---- | ------------- |
        | p1 | 10   | prod one      |
        | p2 | 200  | has \\| a pipe |

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

import re
from collections.abc import Mapping

import yaml
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc

from ._labels import _py, _as_list
from .blocks import plan_blocks, MERGED, RESERVED
from .text import pad_columns
from ._io import _read_source, _write_sink
from rotifer.core import functions as rcf

config = rcf.loadConfig(__name__, defaults = {
    'indent': 4,
})
SEPARATOR = '# ' + '-' * 5

# libyaml where it was built, which reads these documents about 13x faster than
# the pure-Python scanner -- with block scalars this large, that is most of the
# cost of reading. Both loaders produce the same object.
LOADER = getattr(yaml, 'CSafeLoader', yaml.SafeLoader)

# A pipe that is not backslash-escaped, for splitting rows into cells.
_CELL_RE = re.compile(r'(?<!\\)\|')
_RULE_RE = re.compile(r'^\s*\|[\s:|-]+\|\s*$')
# Rows buffered as Python strings before being packed into an arrow array.
_CHUNK = 1 << 15


# --------------------------------------------------------------------- helpers
def _is_text(series):
    """True for a column still holding parsed text. Cells arrive as str, but the
    dtype carrying them depends on the pandas version: object before 3.0, str
    from 3.0 on, where text is no longer inferred as object."""
    return pd.api.types.is_string_dtype(series)


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


def _md_arrow(cols, names, indent=None):
    """Render pre-escaped arrow columns as a markdown table.

    Padding and the surrounding pipes are both applied by arrow kernels, so no
    Python-level pass is made over the rows. indent=None takes it from config."""
    if indent is None:
        indent = ' ' * config['indent']
    padded, widths = pad_columns(cols, names)
    headers = [name.ljust(w) for name, w in zip(names, widths)]
    head = f'{indent}| ' + ' | '.join(headers) + ' |'
    rule = f'{indent}| ' + ' | '.join('-' * len(h) for h in headers) + ' |'
    if not padded or len(padded[0]) == 0:
        return f'{head}\n{rule}'
    body = pc.binary_join_element_wise(*padded, pa.scalar(' | '))
    body = pc.binary_join_element_wise(pa.scalar(f'{indent}| '), body, pa.scalar(' |'),
                                       pa.scalar(''))
    return '\n'.join([head, rule] + body.to_pylist())


def md_table(df, indent=None, index=False):
    """Render a DataFrame as a markdown table, indented for a YAML block scalar."""
    if index:
        df = df.reset_index()
    return _md_arrow(_escape_columns(df), [str(c) for c in df.columns], indent=indent)


def read_md_table(text):
    """Parse one markdown table, cutting and typing cells exactly as read_yaml does."""
    rows = [ln for ln in text.split('\n') if ln.strip().startswith('|')]
    rows = [ln for ln in rows if not _RULE_RE.match(ln)]
    if not rows:
        return pd.DataFrame()
    table = _open_table(_header_cells(rows[0]), [])
    _add_rows(table, rows[1:], [])
    return _to_frames({'table': table}, [])['table']


# --------------------------------------------------------------------- writing


def to_yaml(df, groupby, buf=None, columns=None, sortrows=None,
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

      stats={'kind': 'nunique', 'id': 'nunique'}
          A mapping of column to aggregation. Evaluated for all blocks at once
          with a single groupby.agg, producing a one-row table per block. This is
          much the cheaper form, and `stats` is the conventional name for it.
          Named aggregation works too, {'n_kinds': ('kind', 'nunique')}, naming
          the output column apart from its source.

      kinds=lambda block: block[['kind']].drop_duplicates()
          A callable receiving each block as a DataFrame and returning the table
          to render. Flexible, but it materialises a DataFrame per block.

    So the whole shape of a block is caller-defined:

        to_yaml(df, groupby='grp',
                stats={'kind': 'nunique'},
                span={'size': 'max'},
                products=lambda b: b['product'].value_counts().reset_index())

    writes `stats`, `span` and `products` tables into every block, in that order,
    followed by the main table. The tables need have nothing in common: they are
    never joined to one another here. Columns produced by any aggregation table
    are visible to sortblocks, so blocks can be ordered by a statistic, and if two
    tables define the same column name the one named first decides. Names are
    yours to choose except `key` and `parameters`, which carry block metadata,
    and `merged`, where read_yaml returns the joined frame.

    Rows are brought back together by read_yaml, not here -- it keys them on the
    grouping columns, which it restores into every table from the block key.

    Columns named in sortrows or sortblocks that are not present are ignored, so a
    single ordering can be reused across frames that differ in shape.
    """
    rng = np.random.default_rng() if sample else None
    plan = plan_blocks(df, groupby, sortrows=sortrows, sortblocks=sortblocks,
                       sample=sample, tables=tables, rng=rng)
    head = f'- parameters:\n    groupby: {_flow(plan.groupby)}'

    body = df.filter(columns) if columns else df.drop(columns=plan.groupby, errors='ignore')
    cols = _escape_columns(body)                  # once, not per block
    names = [str(c) for c in body.columns]

    # Aggregate values are read column by column so each keeps its own dtype: a
    # row taken across columns would upcast an int next to a float.
    agg = {table: ([str(c) for c in frame.columns],
                   [frame.iloc[:, j].to_numpy() for j in range(frame.shape[1])])
           for table, frame in plan.tables.items()}

    parts = [head]
    for i, key in enumerate(plan.keys):
        taken = pa.array(plan.positions[key])

        keylist = list(key) if isinstance(key, tuple) else [key]
        block = [SEPARATOR, f'- key: {_flow(keylist)}']
        for table in tables:                      # in the order the caller named them
            if table in agg:
                heads, values = agg[table]
                rendered = _md_arrow([pa.array([str(_py(v[i]))]) for v in values], heads)
            else:
                rendered = md_table(plan.callbacks[table](df.take(plan.positions[key])))
            block.append(f'  {table}: |\n' + rendered)
        block.append(f'  {name}: |\n' + _md_arrow([c.take(taken) for c in cols], names))
        parts.append('\n'.join(block))

    text = sep.join(parts) + '\n'

    return _write_sink(text, buf)


# --------------------------------------------------------------------- reading
class _Table:
    """One table being collected: cells accumulate column by column, packed into
    arrow string arrays every _CHUNK rows so no list of row lists ever exists."""
    __slots__ = ('header', 'nstored', 'buffer', 'chunks', 'extras', 'escaped', 'nrows')

    def __init__(self, stored, extra_names):
        self.header = stored + extra_names
        self.nstored = len(stored)
        self.buffer = [[] for _ in stored]           # pending cells per stored column
        self.chunks = [[] for _ in stored]           # packed arrow arrays per column
        self.extras = [[] for _ in extra_names]      # group columns, as Python values
        self.escaped = False
        self.nrows = 0


def _header_cells(line):
    """Column names of a markdown header line, outer pipes dropped."""
    return [c.strip() for c in line.split('|')][1:-1]


def _flush(table):
    """Pack the buffered cells of every stored column into arrow arrays."""
    for buf, chunks in zip(table.buffer, table.chunks):
        if buf:
            chunks.append(pa.array(buf, type=pa.string()))
            buf.clear()


def _open_table(stored, groupcols):
    """A collector for a table stored with these columns; group columns missing from
    it are carried as extra columns."""
    return _Table(stored, [c for c in groupcols if c not in stored])


def _add_rows(table, lines, extra):
    """Cut `lines` into cells and append them column-wise; `extra` holds one value
    per extra column, repeated on every row."""
    width = table.nstored
    escaped = any(r'\|' in ln for ln in lines)       # rare, so only then pay for the regex
    table.escaped |= escaped
    cut = _CELL_RE.split if escaped else (lambda ln: ln.split('|'))
    rows = []
    for line in lines:
        if not line.strip():
            continue
        cells = cut(line)[1:-1]
        if len(cells) != width:
            if len(cells) > width:
                raise ValueError(f'row has {len(cells)} cells, header has {width}: {line!r}')
            cells += [''] * (width - len(cells))
        rows.append(cells)
    if not rows:
        return
    for buf, col in zip(table.buffer, zip(*rows)):
        buf.extend([c.strip() or None for c in col])
    for values, v in zip(table.extras, extra):
        values.extend([v] * len(rows))
    table.nrows += len(rows)
    if len(table.buffer[0] if table.buffer else ()) >= _CHUNK:
        _flush(table)


def _load_document(source):
    """Source (text, path or file-like) -> the YAML sequence it holds."""
    doc = yaml.load(_read_source(source), Loader=LOADER)
    if not isinstance(doc, list):
        raise ValueError('not a to_yaml document: top level is not a YAML sequence')
    return doc


def _parameters(doc):
    """The grouping columns the document's parameters element names."""
    params = {}
    for element in doc:
        if isinstance(element, dict) and 'parameters' in element:
            params = element['parameters'] or {}
            break
    return _as_list(params.get('groupby', []))


def _collect_rows(doc, groupcols, dropped):
    """Document -> {name: _Table} with every block's rows gathered under the key
    they are stored under; one header is kept and later blocks add data rows only."""
    tables = {}
    for element in doc:
        if not isinstance(element, dict) or 'parameters' in element:
            continue
        values = dict(zip(groupcols, _as_list(element.get('key', []))))
        for name, text in element.items():
            if name in RESERVED or name in dropped or not isinstance(text, str):
                continue
            head, _, rest = text.partition('\n')
            _, _, body = rest.partition('\n')          # the second line is the rule
            if name not in tables:
                tables[name] = _open_table(_header_cells(head), groupcols)
            table = tables[name]
            names = table.header[table.nstored:]
            _add_rows(table, body.split('\n'), [_py(values.get(c)) for c in names])
    return tables


def _numeric(column):
    """Arrow string column -> Series, numeric only if every present value parses."""
    try:
        return pd.Series(pc.cast(column, pa.int64()).to_numpy())
    except (pa.ArrowInvalid, pa.ArrowNotImplementedError):
        pass
    series = pd.Series(column, dtype='str')
    try:
        return pd.to_numeric(series)
    except (ValueError, TypeError):
        return series


def _text_column(table, j):
    """Stored column j of a table -> typed Series (unescaped, numeric if it can be)."""
    _flush(table)
    if not table.nrows:
        return pd.Series([], dtype='int64')     # what to_numeric made of an empty column
    column = pa.chunked_array(table.chunks[j], type=pa.string())
    if column.null_count == len(column):
        return pd.Series([None] * len(column), dtype=object)
    if table.escaped and pc.any(pc.match_substring(column, '\\|')).as_py():
        column = pc.replace_substring(column, '\\|', '|')
    return _numeric(column)


def _to_frames(collected, groupcols):
    """{name: _Table} -> {name: DataFrame}, group columns first, typed per column."""
    frames = {}
    for name, table in collected.items():
        cols = [_text_column(table, j) for j in range(table.nstored)]
        cols += [pd.Series(values) if values else pd.Series([], dtype='int64')
                 for values in table.extras]
        frame = pd.DataFrame(dict(enumerate(cols)))
        frame.columns = table.header
        front = [c for c in groupcols if c in frame.columns]
        frames[name] = frame[front + [c for c in frame.columns if c not in front]]
    return frames


def _check_tables(names, frames):
    """KeyError naming the first of `names` the document does not hold."""
    unknown = [name for name in names if name not in frames]
    if unknown:
        raise KeyError(f'no table named {unknown[0]!r} in the document '
                       f'(have {sorted(frames)})')


def _merge_plan(merge, frames, groupcols):
    """merge argument -> ordered {table: [join columns]}; {} means no merge."""
    if merge is None or merge is False:
        return {}
    if merge is True:
        # No order was given, so take the document in reverse: to_yaml appends the
        # main table after the aggregations, which puts it first here, where it
        # keeps its column names and sets the row granularity.
        plan = {name: groupcols for name in reversed(frames)}
    elif isinstance(merge, Mapping):
        plan = {name: _as_list(c) for name, c in merge.items()}
    else:
        plan = {name: groupcols for name in _as_list(merge)}
    _check_tables(plan, frames)
    if not plan:
        raise ValueError('merge names no tables')
    return plan


def _merge(frames, plan, how):
    """Join the planned tables in order onto the first, suffixing later columns with
    the name of the table they came from."""
    names = list(plan)
    merged = frames[names[0]]
    for name in names[1:]:
        right = frames[name]
        on = [c for c in plan[name] if c in merged.columns and c in right.columns]
        if not on:
            raise ValueError(f'table {name!r} shares none of {plan[name]} with the '
                             f'tables merged before it')
        merged = merged.merge(right, on=on, how=how, suffixes=('', f'_{name}'))
    return merged


def _assemble(merged, frames, plan, keep):
    """The merged frame, whatever the plan left out and any original asked for by
    keep; a lone merged frame comes back bare. No merge returns every table."""
    if merged is None:
        return frames
    result = {MERGED: merged}
    result.update({name: frame for name, frame in frames.items() if name not in plan})
    result.update({name: frames[name] for name in keep})
    return merged if len(result) == 1 else result


def read_yaml(source, merge=True, drop=None, keep=None, how='outer'):
    """Read a document written by to_yaml.

    source : YAML text, a path, or a file-like object.
    how    : join used when merging tables held under different keys.
    merge  : what to fold together, and on what.

        True          the default: every table is merged on the grouping columns.
        None, False   nothing is merged.
        ['a', 'b']    only these tables are merged, on the grouping columns.
        {'a': 'id'}  only these tables are merged, each joined on the columns
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
    with that table's key, so aggregating `id` under `stats` lands as `id_stats`
    beside the row-level `id`. The widest table keeps the unsuffixed names.
    """
    doc = _load_document(source)
    groupcols = _parameters(doc)
    keep = _as_list(keep)                         # keep wins over drop
    dropped = set(_as_list(drop)) - set(keep)

    collected = _collect_rows(doc, groupcols, dropped)
    if not collected:
        raise ValueError('no tables found in document')
    frames = _to_frames(collected, groupcols)
    _check_tables(keep, frames)
    plan = _merge_plan(merge, frames, groupcols)
    merged = _merge(frames, plan, how) if plan else None
    return _assemble(merged, frames, plan, keep)


read_blocks = read_yaml  # alias
