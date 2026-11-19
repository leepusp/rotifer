#!/usr/bin/env python3
"""Block reports as hand-editable YAML documents holding markdown tables.

:func:`to_yaml` writes the same blocks as :func:`rotifer.pandas.to_blocks`, which
are planned by :func:`rotifer.pandas.blocks.plan_blocks`, but in a format that can
be edited by hand and read back with :func:`read_yaml`. :func:`md_table` and
:func:`read_md_table` convert a single DataFrame to and from one markdown table.

A document is a YAML sequence. The first element holds the file-level parameters,
notably the grouping columns, which are stored once and not repeated on every
block. Each later element is one block: a ``key`` listing the group's values, one
per grouping column, plus any number of further entries, each a markdown table
stored as a literal block scalar. The names ``key``, ``parameters`` and ``merged``
are reserved and cannot name a table.

    - parameters:
        groupby: [grp]

    # -----
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

A literal pipe inside a cell is written ``\\|``. Padding is cosmetic: whitespace
around cells is stripped when reading, so a hand-edited table with ragged pipes
still parses. Writing formats the text directly, which keeps the comments, blank
lines and column padding. Reading is a safe YAML load followed by a split on
unescaped pipes.

Writing does not go through :func:`rotifer.pandas.to_blocks`. Text is escaped and
converted once for the whole frame, aggregations are computed for all blocks in
one pass and each block is cut out with an arrow ``take``, so no DataFrame is
built per block. This needs pyarrow.

Reading combines the tables by the name they are stored under: the rows of a table
of the same name are concatenated across blocks, and the resulting tables are
merged on the grouping columns.
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
def _flow(value):
    """Render a scalar or sequence as a one-line YAML flow value, quoted by PyYAML when needed."""
    if isinstance(value, tuple):
        value = list(value)
    if isinstance(value, list):
        value = [_py(v) for v in value]
    else:
        value = _py(value)
    return yaml.safe_dump(value, default_flow_style=True, width=10 ** 6,
                          allow_unicode=True).strip()


def _escape_columns(df):
    """Convert each column of `df` to an arrow string array, escaping pipes and flattening newlines.

    This happens before padding, so column widths include the escape characters,
    and once for the whole frame, so blocks can be cut from the result.
    """
    out = []
    for c in df.columns:
        arr = pa.Array.from_pandas(df[c].astype('string').fillna(''), type=pa.string())
        arr = pc.replace_substring(arr, '|', r'\|')
        out.append(pc.replace_substring_regex(arr, r'\s*\n\s*', ' '))
    return out


def _md_arrow(cols, names, indent=None):
    """Render pre-escaped arrow columns as an indented markdown table, using arrow kernels.

    `indent` is a string; None takes the width from the ``indent`` setting.
    """
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
    """Render a DataFrame as a padded markdown table.

    The table has a header line, a rule line and one line per row. Missing values
    become empty cells, a pipe inside a cell is written ``\\|`` and a newline is
    replaced by a space.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to render.
    indent : str, optional
        Prefix of every line, for embedding the table in a YAML block scalar. Four
        spaces (the ``indent`` setting) when None.
    index : bool, default False
        Write the index as the first column.

    Returns
    -------
    str
        The table, without a trailing newline.

    See Also
    --------
    read_md_table : Parse such a table back into a DataFrame.
    rotifer.pandas.to_string : Aligned plain text instead of markdown.

    Examples
    --------
    >>> import pandas as pd
    >>> print(md_table(pd.DataFrame({'id': ['p1', 'p2'], 'kb': [10, 200]}), indent=''))
    | id | kb  |
    | -- | --- |
    | p1 | 10  |
    | p2 | 200 |
    """
    if index:
        df = df.reset_index()
    return _md_arrow(_escape_columns(df), [str(c) for c in df.columns], indent=indent)


def read_md_table(text):
    """Parse one markdown table into a DataFrame.

    Lines that do not start with a pipe and the rule line under the header are
    ignored. Cells are stripped of surrounding whitespace and ``\\|`` is read as a
    pipe. Empty cells become missing values. A column is converted to numbers only
    when every non-empty cell in it parses as one. A row with fewer cells than the
    header is padded with missing values.

    Parameters
    ----------
    text : str
        The table, indented or not.

    Returns
    -------
    pandas.DataFrame
        One column per header cell. An empty frame when `text` holds no table.

    Raises
    ------
    ValueError
        A row has more cells than the header.

    See Also
    --------
    md_table : Write a DataFrame as such a table.
    read_yaml : Read a whole YAML document of tables.

    Examples
    --------
    >>> read_md_table('''
    ... | id | kb  |
    ... | -- | --- |
    ... | p1 | 10  |
    ... | p2 | 200 |
    ... ''')
       id   kb
    0  p1   10
    1  p2  200
    """
    rows = [ln for ln in text.split('\n') if ln.strip().startswith('|')]
    rows = [ln for ln in rows if not _RULE_RE.match(ln)]
    if not rows:
        return pd.DataFrame()
    table = _open_table(_header_cells(rows[0]), [])
    _add_rows(table, rows[1:], [])
    return _to_frames({'table': table}, [])['table']


# --------------------------------------------------------------------- writing


_PLAIN = re.compile(r'[A-Za-z_][A-Za-z0-9_]*\Z')
_RESOLVED = frozenset('yes no true false on off null'.split())


def _flow_keys(keylist):
    """`_flow` for a key list, skipping PyYAML when every part is a plain word."""
    for v in keylist:
        if not (type(v) is str and _PLAIN.match(v) and v.lower() not in _RESOLVED):
            return _flow(keylist)
    return '[' + ', '.join(keylist) + ']'


def _one_row_table(indent, names, cells):
    """A header, rule and one padded row, from names and ready-made cell strings."""
    widths = [max(len(n), len(c)) for n, c in zip(names, cells)]
    head = [n.ljust(w) for n, w in zip(names, widths)]
    return (f'{indent}| ' + ' | '.join(head) + ' |\n'
            + f'{indent}| ' + ' | '.join('-' * w for w in widths) + ' |\n'
            + f'{indent}| ' + ' | '.join(c.ljust(w) for c, w in zip(cells, widths)) + ' |')


def _rows_table(indent, names, name_widths, rows, lengths, pos):
    """The main table of one block: the rows at `pos`, padded to the block's widths."""
    pos = pos.tolist()
    widths = [max(nw, int(ln[pos].max())) if pos else nw
              for nw, ln in zip(name_widths, lengths)]
    lead = f'{indent}| '
    lines = [lead + ' | '.join(n.ljust(w) for n, w in zip(names, widths)) + ' |',
             lead + ' | '.join('-' * w for w in widths) + ' |']
    if names and pos:
        line = lead + ' | '.join(f'%-{w}s' for w in widths) + ' |'
        cols = [[r[p] for p in pos] for r in rows]
        lines += [line % cells for cells in zip(*cols)]
    return '\n'.join(lines)


def to_yaml(df, groupby, buf=None, columns=None, sortrows=None,
            sortblocks=None, sample=None, random_state=None, name='rows', sep='\n\n',
            **tables):
    """Write `df` as a YAML document of blocks holding markdown tables.

    The document format is described in the module docstring. Each block holds
    one table per extra keyword (aggregations or callables, in the order given),
    followed by the main table of the block's rows.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to write.
    groupby : label or list of labels
        Column(s) defining the blocks. They are stored once, in the document
        parameters, and as the key of each block.
    buf : str, path or file-like, optional
        Where to write. The text is returned when None.
    columns : list, optional
        Columns of the main table. All but the grouping columns when None.
    sortrows : dict, optional
        ``{column: ascending}``. Orders the rows inside each block.
    sortblocks : dict, optional
        ``{column: ascending}``. Orders the blocks. It may name a grouping column
        or a column produced by an aggregation table.
    sample : int, optional
        Keep at most this many rows per block, chosen at random.
    random_state : int or numpy.random.Generator, optional
        Seed for `sample`. A fresh generator is used when None.
    name : str, default 'rows'
        Name of the main table in each block.
    sep : str, default '\\n\\n'
        Text between the parameters and the first block and between blocks.
    **tables
        One further table per keyword, named by the keyword. Two forms are
        accepted. A mapping, ``stats={'kind': 'nunique'}`` or the named form
        ``stats={'n_kinds': ('kind', 'nunique')}``, is aggregated for all blocks
        in one pass and gives a one-row table per block; a 2-tuple value is always
        read as ``(column, aggregation)``. A callable,
        ``kinds=lambda block: block[['kind']].drop_duplicates()``, receives each
        block as a DataFrame and returns the table to render.

    Returns
    -------
    str or None
        The document, or None when `buf` is given.

    Raises
    ------
    KeyError
        A grouping column or a column named by an aggregation is missing.
    ValueError
        A table uses a reserved name, or an aggregation fails.
    TypeError
        A table is neither a mapping nor a callable.

    See Also
    --------
    read_yaml : Read the document back.
    rotifer.pandas.to_blocks : The same blocks as aligned plain text.
    rotifer.pandas.blocks.plan_blocks : Ordering, sampling and aggregation rules.

    Notes
    -----
    Pipes inside cells are escaped as ``\\|`` and newlines are replaced by a space,
    so a cell cannot break the table or the YAML block. Missing values become empty
    cells. Aggregated values are written with ``str``.

    Blocks appear in order of first appearance of their key unless `sortblocks` is
    given. `sortrows` never moves blocks. `sample` is drawn per block before
    `sortrows` is applied. Columns named in either that do not exist are ignored.
    Rows whose group key is missing (NaN or None) are not written.

    The names ``key`` and ``parameters`` carry block metadata and ``merged`` is where
    :func:`read_yaml` returns the joined frame, so none of them can name a table.
    Tables are not joined to each other here; :func:`read_yaml` recombines them on
    the grouping columns. A callable table builds a DataFrame for every block and is
    much slower than a mapping.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'g': ['A', 'A', 'B'], 'id': ['p1', 'p2', 'p3'],
    ...                    'kb': [10, 20, 30]})
    >>> print(to_yaml(df, 'g', stats={'n': ('id', 'nunique')}))
    - parameters:
        groupby: [g]
    <BLANKLINE>
    # -----
    - key: [A]
      stats: |
        | n |
        | - |
        | 2 |
      rows: |
        | id | kb |
        | -- | -- |
        | p1 | 10 |
        | p2 | 20 |
    <BLANKLINE>
    # -----
    - key: [B]
      stats: |
        | n |
        | - |
        | 1 |
      rows: |
        | id | kb |
        | -- | -- |
        | p3 | 30 |
    <BLANKLINE>
    """
    plan = plan_blocks(df, groupby, sortrows=sortrows, sortblocks=sortblocks,
                       sample=sample, tables=tables, rng=np.random.default_rng(random_state))
    head = f'- parameters:\n    groupby: {_flow(plan.groupby)}'

    body = df.filter(columns) if columns is not None else df.drop(columns=plan.groupby, errors='ignore')
    cols = _escape_columns(body)                  # once, not per block
    names = [str(c) for c in body.columns]

    # Aggregate values are read column by column so each keeps its own dtype: a
    # row taken across columns would upcast an int next to a float. Each is turned
    # into its text once, for all blocks.
    indent = ' ' * config['indent']
    agg = {table: ([str(c) for c in frame.columns],
                   [[str(_py(v)) for v in frame.iloc[:, j].tolist()]
                    for j in range(frame.shape[1])])
           for table, frame in plan.tables.items()}

    # The main table is cut from Python strings, converted once.
    rows = [c.to_pylist() for c in cols]
    heads = [len(n) for n in names]
    lengths = [np.fromiter(map(len, r), dtype=np.int64, count=len(r)) for r in rows]

    parts = [head]
    for i, key in enumerate(plan.keys):
        pos = plan.positions[key]

        keylist = list(key) if isinstance(key, tuple) else [key]
        block = [SEPARATOR, f'- key: {_flow_keys(keylist)}']
        for table in tables:                      # in the order the caller named them
            if table in agg:
                heads_t, values = agg[table]
                rendered = _one_row_table(indent, heads_t, [v[i] for v in values])
            else:
                rendered = md_table(plan.callbacks[table](df.take(pos)))
            block.append(f'  {table}: |\n' + rendered)
        block.append(f'  {name}: |\n' + _rows_table(indent, names, heads, rows, lengths, pos))
        parts.append('\n'.join(block))

    text = sep.join(parts) + '\n'

    return _write_sink(text, buf)


# --------------------------------------------------------------------- reading
class _Table:
    """Cells of one table being collected, kept column by column and packed into arrow arrays every _CHUNK rows."""
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
    """Column names of a markdown header line."""
    return [c.strip() for c in line.split('|')][1:-1]


def _flush(table):
    """Pack the cells buffered for each stored column into arrow arrays."""
    for buf, chunks in zip(table.buffer, table.chunks):
        if buf:
            chunks.append(pa.array(buf, type=pa.string()))
            buf.clear()


def _open_table(stored, groupcols):
    """Start collecting a table with these stored columns; grouping columns it lacks become extra columns."""
    return _Table(stored, [c for c in groupcols if c not in stored])


def _add_rows(table, lines, extra):
    """Split `lines` into cells and append them column-wise, repeating the `extra` values on every row.

    Short rows are padded with empty cells; a row with too many cells raises ValueError.
    """
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
    """Load the YAML sequence from text, a path or a file-like object; ValueError if it is not a sequence."""
    doc = yaml.load(_read_source(source), Loader=LOADER)
    if not isinstance(doc, list):
        raise ValueError('not a to_yaml document: top level is not a YAML sequence')
    return doc


def _parameters(doc):
    """Grouping columns named by the parameters element of the document, as a list."""
    params = {}
    for element in doc:
        if isinstance(element, dict) and 'parameters' in element:
            params = element['parameters'] or {}
            break
    return _as_list(params.get('groupby', []))


def _collect_rows(doc, groupcols, dropped):
    """Gather the rows of every block into one _Table per table name, skipping the `dropped` names.

    The header comes from the first block that holds the table; later blocks add rows only.
    """
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
    """Convert an arrow string column to a Series, numeric only if every present value parses."""
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
    """Stored column `j` of a table as a Series, with pipes unescaped and numbers converted where possible."""
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
    """Build a DataFrame per collected table, with the grouping columns first."""
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
    """Raise KeyError naming the first of `names` that is not a table in `frames`."""
    unknown = [name for name in names if name not in frames]
    if unknown:
        raise KeyError(f'no table named {unknown[0]!r} in the document '
                       f'(have {sorted(frames)})')


def _merge_plan(merge, frames, groupcols):
    """Turn the `merge` argument of read_yaml into an ordered {table: [join columns]}; {} means no merge."""
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
    """Merge the planned tables in order onto the first.

    Columns that collide get the suffix ``_<table>``. Raises ValueError if a table
    shares no join column with those before it.
    """
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
    """Choose what read_yaml returns: the merged frame, tables left out of the plan and those in `keep`.

    A lone merged frame is returned bare; with no merge every table is returned.
    """
    if merged is None:
        return frames
    result = {MERGED: merged}
    result.update({name: frame for name, frame in frames.items() if name not in plan})
    result.update({name: frames[name] for name in keep})
    return merged if len(result) == 1 else result


def read_yaml(source, merge=True, drop=None, keep=None, how='outer'):
    """Read a YAML document written by :func:`to_yaml` into DataFrames.

    Parameters
    ----------
    source : str, path or file-like
        YAML text, a path to a file, or an open file.
    merge : bool, list, str or dict, default True
        Which tables to merge, and on what:

        * ``True``: merge every table on the grouping columns.
        * ``None`` or ``False``: merge nothing.
        * ``['a', 'b']`` or ``'a'``: merge only these tables, on the grouping columns.
        * ``{'a': 'id'}`` or ``{'a': ['id', ...]}``: merge only these tables, each
          on the columns given for it.
    drop : str or list, optional
        Tables to leave out entirely.
    keep : str or list, optional
        Tables to return as stored, without merging, in addition to any merged
        frame. If a table is in both `drop` and `keep`, `keep` wins.
    how : str, default 'outer'
        Join type used by :meth:`pandas.DataFrame.merge`.

    Returns
    -------
    pandas.DataFrame or dict
        A single DataFrame when it is the only result, which is the usual case with
        ``merge=True``. Otherwise a ``{name: DataFrame}`` mapping: the merged frame
        under the key ``'merged'``, the tables left out of a partial merge, and
        those requested with `keep`. With no merge, every table not dropped.
        Each frame has the grouping columns first.

    Raises
    ------
    ValueError
        The document is not a YAML sequence, holds no tables, a row has more cells
        than its header, `merge` names no tables, or a merged table shares no
        join column with the ones before it.
    KeyError
        `merge` or `keep` names a table the document does not hold.

    See Also
    --------
    to_yaml : Write such a document.
    read_md_table : Read one markdown table.

    Notes
    -----
    The grouping columns are listed in the document parameters. Their values are
    taken from each block's key and added to every table, so they need not appear
    in the stored tables.

    Tables are recombined by the name they are stored under: tables of the same
    name are concatenated across blocks, so a ``stats`` table is never joined to
    a ``rows`` table by accident of matching columns. The resulting tables are then
    merged on the grouping columns, or on the columns given in `merge`. When a
    column other than a join column exists in more than one table, those from later
    tables get the suffix ``_<table>``, so ``id`` aggregated in ``stats`` becomes
    ``id_stats``. With ``merge=True`` the tables are merged in reverse document
    order, so the main table comes first and keeps its column names.

    Each column is converted to numbers only if every non-empty cell parses as a
    number, otherwise it stays text. Empty cells become missing values.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'g': ['A', 'A', 'B'], 'id': ['p1', 'p2', 'p3'],
    ...                    'kb': [10, 20, 30]})
    >>> text = to_yaml(df, 'g', stats={'n': ('id', 'nunique')})
    >>> read_yaml(text)
       g  id  kb  n
    0  A  p1  10  2
    1  A  p2  20  2
    2  B  p3  30  1
    >>> tables = read_yaml(text, merge=False)
    >>> sorted(tables)
    ['rows', 'stats']
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
