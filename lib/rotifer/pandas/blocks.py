"""Block reports: group a table by one or more columns and render each group as a block.

A block report splits the rows of a DataFrame by the values of one or more
grouping columns. Each group becomes a block: a header line naming the group
(optionally followed by per-block statistics) and then the rows of the group.

Two writers produce block reports. :func:`to_blocks`, defined here, writes
aligned plain text. :func:`rotifer.pandas.yaml_blocks.to_yaml` writes a YAML
document whose blocks hold markdown tables, which :func:`rotifer.pandas.yaml_blocks.read_yaml`
reads back. Both call :func:`plan_blocks`, so the order of the blocks, the order
and sampling of the rows and the per-block tables are decided in one place and
are the same in both formats.

The text format is a header of two lines, ``#`` and ``# <groupby>=<key>`` with any
statistics appended as ``; <name>=<value>``, followed by the block's rows as
columns separated by ``colsep``. Blocks are separated by a newline. For two groups
``A`` (2 rows) and ``B`` (1 row) with a statistic ``n`` the report reads::

    #
    # g=A; n=2
      : id : kb
    0 : p1 : 10
    1 : p2 : 20
    #
    # g=B; n=1
      : id : kb
    2 : p3 : 30

This module also holds the pieces the two writers share: :class:`BlockPlan`,
:func:`plan_blocks`, :func:`default_header` and the names reserved for block
metadata in YAML documents.
"""

import types
import warnings
import numpy as np
import pandas as pd
from collections.abc import Mapping
from dataclasses import dataclass

from ._labels import _as_label, _fmt_label, _group_key, _group_labels, _py
from ._io import _write_sink
from .text import to_string, _check_align, _text_columns, _join_rows

def _resolve_stats(df, stats: Mapping):
    """Evaluate {column: aggregation name} on `df`, returning {column: value}.

    Used only by the deprecated per-block path. Array-like column labels become
    tuples. Raises KeyError for a missing column. A string the column has no
    attribute for is passed through unchanged, so resolving twice is harmless and
    a misspelled aggregation on an existing column is not detected.
    """
    out = {}
    for c, f in stats.items():
        c = _as_label(c)
        if c not in df:
            raise KeyError(f'stats column {c!r} is not in the DataFrame')
        val = f
        if isinstance(f, str) and hasattr(df[c], f):
            val = getattr(df[c], f)
            val = val() if isinstance(val, types.MethodType) else val
        out[c] = val
    return out

MERGED = 'merged'                       # read_yaml returns the joined frame here
RESERVED = ('key', 'parameters', MERGED)


@dataclass(frozen=True)
class BlockPlan:
    """Blocks of a grouped table, decided once for the whole frame.

    Returned by :func:`plan_blocks` and consumed by the text and YAML writers. The
    fields below are also available as attributes; the dataclass is frozen.

    Parameters
    ----------
    groupby : list
        Grouping column labels, in the order given.
    keys : list
        Group keys in block order. A key is a tuple with one value per grouping
        column, or a bare value when `groupby` named a single column as a scalar
        label.
    positions : dict
        ``{key: numpy.ndarray}`` of integer row positions in the original frame,
        already ordered by ``sortrows`` and reduced by ``sample``.
    tables : dict
        ``{name: pandas.DataFrame}`` of per-block aggregates. Row ``i`` of every
        frame belongs to ``keys[i]``. Tables come in the order the caller named
        them.
    callbacks : dict
        ``{name: callable}`` for table specifications that are functions. Each
        callable takes one block as a DataFrame and returns that block's table.
    """
    groupby: list
    keys: list
    positions: dict
    tables: dict
    callbacks: dict


def _row_rank(df, sortrows):
    """Rank of every row under the sortrows columns, or None if none is present.

    Ranking once and reordering each block's positions by rank sorts the rows
    inside a block without changing the order of the blocks.
    """
    cols = [c for c in sortrows if c in df]
    if not cols:
        return None
    order = (df[cols].reset_index(drop=True)
               .sort_values(cols, ascending=[sortrows[c] for c in cols], kind='stable')
               .index.to_numpy())
    rank = np.empty(len(df), dtype=np.int64)
    rank[order] = np.arange(len(df))
    return rank


def _group_agg(grouped, spec, table):
    """Evaluate one aggregation table for every block in a single groupby pass.

    `spec` maps an output name to an aggregation (name or callable) of the column
    of the same name, to a named aggregation ``(column, aggregation)``, or to a
    plain value repeated for every block. A 2-tuple is always read as
    ``(column, aggregation)``. Raises KeyError for a missing column and ValueError,
    naming the table and column, for an aggregation that fails. Returns a frame
    indexed by group key with columns in spec order.
    """
    obj = grouped.obj
    entries = []
    for out, val in spec.items():
        col, agg = val if isinstance(val, tuple) and len(val) == 2 else (out, val)
        col = _as_label(col)
        if col not in obj.columns:
            raise KeyError(f'{table} column {col!r} is not in the DataFrame')
        entries.append((_as_label(out), col, agg))

    named = {f'_{i}': (col, agg) for i, (_, col, agg) in enumerate(entries)
             if isinstance(agg, str) or callable(agg)}
    if named:
        try:
            frame = grouped.agg(**named)
        except Exception as error:
            for i, (_, col, agg) in enumerate(entries):
                if f'_{i}' in named:
                    try:
                        grouped.agg(**{f'_{i}': named[f'_{i}']})
                    except Exception as one:
                        raise ValueError(f'table {table!r}: cannot aggregate column '
                                         f'{col!r} with {agg!r}') from one
            raise
    else:
        frame = pd.DataFrame(index=grouped.size().index)
    for i, (_, _, agg) in enumerate(entries):
        if f'_{i}' not in named:
            frame[f'_{i}'] = agg
    frame = frame[[f'_{i}' for i in range(len(entries))]]
    frame.columns = pd.Index([out for out, _, _ in entries], tupleize_cols=False)
    return frame


def _split_tables(grouped, tables):
    """Split table specifications into aggregations and callables.

    A mapping is evaluated once for all blocks; a callable receives each block as
    a DataFrame, which is much slower. Raises ValueError for a reserved table name
    and TypeError for any other kind of specification.
    """
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


def plan_blocks(df, groupby, *, sortrows=None, sortblocks=None, sample=None,
                tables=None, rng=None):
    """Group `df` once and decide everything about its blocks except their text.

    This is the planning step shared by :func:`to_blocks` and
    :func:`rotifer.pandas.yaml_blocks.to_yaml`. It decides which rows belong to
    which block, the order of blocks and of rows within them, which rows are
    sampled and the values of every aggregation table.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to split into blocks.
    groupby : label or list of labels
        Column(s) defining the blocks. A list of several columns gives tuple keys.
    sortrows : dict, optional
        ``{column: ascending}``. Orders the rows inside each block.
    sortblocks : dict, optional
        ``{column: ascending}``. Orders the blocks. It may name a grouping column
        or a column produced by any aggregation table.
    sample : int, optional
        Keep at most this many rows in each block.
    tables : dict, optional
        ``{name: specification}``, one per-block table each. A specification is
        either a mapping (an aggregation table) or a callable (see Notes).
    rng : numpy.random.Generator, optional
        Source of randomness for `sample`. A fresh generator is used when None.

    Returns
    -------
    BlockPlan
        Keys, row positions, aggregation tables and callables, ready for a writer.

    Raises
    ------
    KeyError
        A grouping column, or a column named by an aggregation table, is not in
        `df`.
    ValueError
        A table uses a reserved name, or one of its aggregations fails.
    TypeError
        A table specification is neither a mapping nor a callable.

    See Also
    --------
    to_blocks : Text report built from this plan.
    rotifer.pandas.yaml_blocks.to_yaml : YAML report built from this plan.

    Notes
    -----
    Blocks come out in order of first appearance of their key, unless
    `sortblocks` reorders them. If two aggregation tables define the same column,
    the table named first decides the order.

    `sortrows` orders rows inside each block and never moves a block. `sample` is
    drawn per block, before `sortrows` is applied, so the sampled rows are then
    sorted like any other. Columns named in `sortrows` or `sortblocks` that do not
    exist are ignored.

    Rows whose group key is missing (NaN or None) belong to no block, as with
    :meth:`pandas.DataFrame.groupby` by default.

    An aggregation table accepts two forms: ``{column: aggregation}``, which
    aggregates the column of that name, and ``{name: (column, aggregation)}``
    which names the output separately. A 2-tuple value is always read as
    ``(column, aggregation)``. An aggregation is a method name such as
    ``'nunique'`` or a callable. A value that is neither is repeated unchanged
    in every block. All aggregation tables are computed with one groupby pass.

    A callable table receives each block as a DataFrame and returns that block's
    table. It is more flexible, but it builds a DataFrame for every block, so its
    cost grows with the number of blocks.

    The names ``key`` and ``parameters`` (block metadata in YAML documents) and
    ``merged`` (where :func:`rotifer.pandas.yaml_blocks.read_yaml` returns the
    joined frame) are reserved and cannot name a table.
    """
    tables = tables or {}
    labels, scalar = _group_labels(df, groupby)
    grouped = df.groupby(labels[0] if scalar else labels, sort=False)
    aggregated, callbacks = _split_tables(grouped, tables)

    def norm(key):
        return key if scalar or isinstance(key, tuple) else (key,)

    index = next(iter(aggregated.values())).index if aggregated else grouped.size().index
    keys = [norm(k) for k in index]
    positions = {norm(k): np.asarray(v) for k, v in grouped.indices.items()}

    # sortblocks columns are collected one at a time rather than by joining the
    # tables, which are free to share nothing with each other. All frames come
    # from the same grouping, so they align positionally.
    order = np.arange(len(keys))
    if sortblocks:
        by, values = [], {}
        for col in sortblocks:
            if col in labels:
                i = labels.index(col)
                values[len(by)] = [k[i] if not scalar else k for k in keys]
            else:
                found = next((t for t in aggregated.values() if col in t.columns), None)
                if found is None:
                    continue
                values[len(by)] = found[col].to_numpy()
            by.append(col)
        if by:
            frame = pd.DataFrame(values).sort_values(
                list(range(len(by))), ascending=[sortblocks[c] for c in by], kind='stable')
            order = frame.index.to_numpy()
    keys = [keys[i] for i in order]
    aggregated = {name: t.iloc[order] for name, t in aggregated.items()}
    aggregated = {name: aggregated[name] for name in tables if name in aggregated}

    rank = _row_rank(df, sortrows) if sortrows else None
    if sample and rng is None:
        rng = np.random.default_rng()
    for key in keys:
        idx = positions[key]
        if sample and sample < len(idx):
            idx = np.sort(rng.choice(idx, sample, replace=False))
        if rank is not None:
            idx = idx[np.argsort(rank[idx], kind='stable')]
        positions[key] = idx

    return BlockPlan(groupby=labels, keys=keys, positions=positions,
                     tables=aggregated, callbacks=callbacks)


REMOVAL = '2026-11-30'

# Options of text.to_string that to_blocks forwards; anything else is a mistake.
_TO_STRING_OPTIONS = ('index', 'na_rep', 'float_format', 'align', 'insert')

def _deprecated(what, instead, stacklevel=2):
    """Warn with a FutureWarning that `what` goes away on REMOVAL; `instead` says what to use."""
    # stacklevel counts from the deprecated function: 2 is its caller.
    warnings.warn(f'{what} is deprecated and will be removed on {REMOVAL} (issue #TBD); {instead}.',
                  FutureWarning, stacklevel=stacklevel + 1)


def default_header(key, groupby, stats=None, **kwargs):
    """Return the header of one block.

    .. deprecated::
       Will be removed on 2026-11-30. Passing the block as a DataFrame,
       ``default_header(df, groupby=..., stats=...)``, is deprecated: use
       :func:`default_header` with the group key and a ``{name: value}`` dict.

    The header has two lines: ``#`` and ``# <groupby>=<key>`` followed by
    ``; <name>=<value>`` for each statistic. :func:`to_blocks` calls
    ``header(key, groupby, stats)``, and a replacement header takes the same three
    arguments and returns the text to put above the rows of the block.

    Parameters
    ----------
    key : object or tuple
        Group key of the block: a bare value, or a tuple with one value per
        grouping column. A DataFrame selects the deprecated call form.
    groupby : label or tuple of labels
        Grouping column label(s), shown to the left of the ``=``.
    stats : dict, optional
        ``{name: value}`` computed for this block. In the deprecated call form it
        is ``{column: aggregation name}``, evaluated on the DataFrame.
    **kwargs
        Ignored. They are accepted so that the deprecated call form keeps working.

    Returns
    -------
    str
        The header text, without a trailing newline.

    Warns
    -----
    FutureWarning
        When `key` is a DataFrame.

    See Also
    --------
    to_blocks : Calls the header once per block.

    Examples
    --------
    >>> print(default_header('A', 'g', {'n': 2}))
    #
    # g=A; n=2
    """
    if isinstance(key, pd.DataFrame):
        _deprecated('default_header(df, ...) with a DataFrame',
                    'pass the group key and a {name: value} dict: '
                    'default_header(key, groupby, stats)', stacklevel=2)
        return _legacy_header(key, groupby, stats)
    tail = ''.join(f'; {_fmt_label(c)}={v}' for c, v in (stats or {}).items())
    return f'#\n# {_fmt_label(groupby)}={_fmt_label(key)}{tail}'


def _legacy_header(df, groupby, stats):
    """default_header for the deprecated call form: the block as a DataFrame."""
    groupby = _as_label(groupby)
    key = _group_key(df, groupby)
    labels, scalar = _group_labels(df, groupby)
    if not scalar:
        # Key arity follows label arity: pandas unwraps a one-label grouping, and
        # an empty or unknown group has no key at all.
        key = key if isinstance(key, tuple) else (key,) * len(labels) if key is None else (key,)
    tail = ''.join(f'; {_fmt_label(c)}={v}' for c, v in _resolve_stats(df, stats or {}).items())
    return f'#\n# {_fmt_label(groupby)}={_fmt_label(key)}{tail}'


def default_block_rows(df, groupby, header=default_header, colsep=" : ",
                          columns=None, sample=None, stats=None, sortrows=None, **kwargs):
    """Render one block as a one-row DataFrame whose ``blocks`` column holds its text.

    .. deprecated::
       Will be removed on 2026-11-30. :func:`to_blocks` builds the blocks itself;
       customise them with its `header`, `columns` and `stats` arguments.

    Parameters
    ----------
    df : pandas.DataFrame
        Rows of one block.
    groupby : label or list of labels
        Grouping column label(s), passed to `header`.
    header : callable, default default_header
        Called as ``header(df, groupby=..., stats=...)``.
    colsep : str, default ' : '
        Separator between columns.
    columns : list, optional
        Columns to show; all when None.
    sample : int, optional
        Keep at most this many rows, chosen at random.
    stats : dict, optional
        ``{column: aggregation name}`` evaluated on the block and passed to `header`.
    sortrows : dict, optional
        ``{column: ascending}`` ordering the rows.
    **kwargs
        Passed to `header`; those that are options of :func:`rotifer.pandas.to_string`
        are passed to it as well.

    Returns
    -------
    pandas.DataFrame
        One row: the ``blocks`` column with the text, then one column per statistic.

    Warns
    -----
    FutureWarning
        Always.
    """
    _deprecated('default_block_rows', 'to_blocks builds the blocks itself now', stacklevel=2)
    statDict = _resolve_stats(df, stats or {})
    if header is default_header:
        # Build the header through the private helper: calling default_header(df, ...)
        # would emit a second FutureWarning from inside the library.
        headerStr = _legacy_header(df, groupby, statDict)
    else:
        headerStr = header(df, groupby=groupby, stats=statDict, **kwargs)

    block = df.copy()
    if sample and sample < len(block):
        block = block.sample(sample)
    if sortrows:
        sortcols = [c for c in sortrows if c in block]
        if sortcols:
            block = block.sort_values(sortcols, ascending=[sortrows[c] for c in sortcols])

    if columns is not None:
        block = block.filter(columns)
    own = {k: v for k, v in kwargs.items() if k in _TO_STRING_OPTIONS}
    block = [headerStr + "\n" + to_string(block, sep=colsep, **own)]
    return pd.DataFrame({'blocks': block, **statDict})


def _call_header(header, key, label, values):
    """Call ``header(key, label, values)``, turning a failure of a hook written for the
    old ``(df, groupby, stats)`` signature into an explanatory TypeError."""
    if header is default_header:
        return header(key, label, values)
    try:
        return header(key, label, values)
    except (TypeError, AttributeError) as error:
        raise TypeError('header hooks are now called as header(key, groupby, stats); '
                        'key is the group key and stats a {name: value} dict, no longer '
                        'the block DataFrame') from error


def to_blocks(df, groupby, *, columns=None, stats=None, sortrows=None,
              sortblocks=None, sample=None, header=default_header, colsep=' : ',
              sep='\n', buf=None, apply=None, random_state=None, **to_string_kwargs):
    """Write `df` as text blocks, one per group: a header, then the rows of the group.

    .. deprecated::
       Will be removed on 2026-11-30. The `apply` argument is deprecated: build
       the blocks with `header`, `columns` and `stats` instead.

    Each block starts with the header produced by `header` (by default
    :func:`default_header`: ``#`` and ``# <groupby>=<key>; <stat>=<value>``) and
    continues with the block's rows, formatted by :func:`rotifer.pandas.to_string`.
    Blocks are joined with `sep`.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to write.
    groupby : label or list of labels
        Column(s) to group by. Blocks come out in order of first appearance of
        their key unless `sortblocks` is given.
    columns : list, optional
        Columns shown in the rows of each block. All columns except the grouping
        ones when None.
    stats : dict, optional
        Per-block statistics shown in the header: ``{column: aggregation}`` or
        ``{name: (column, aggregation)}``. A 2-tuple value is always read as
        ``(column, aggregation)``. Computed once for all blocks.
    sortrows : dict, optional
        ``{column: ascending}``. Orders the rows inside each block.
    sortblocks : dict, optional
        ``{column: ascending}``. Orders the blocks. It may name a grouping column
        or a column of `stats`.
    sample : int, optional
        Rows per block. A block with more rows is sampled at random.
    header : callable, default default_header
        Called as ``header(key, groupby, stats) -> str`` for each block. `key` is
        the group key (a bare value or a tuple), `groupby` the grouping label(s)
        and `stats` the ``{name: value}`` dict of that block (empty without
        `stats`). It returns the text placed above the rows.
    colsep : str, default ' : '
        Separator between the columns of a row.
    sep : str, default '\\n'
        Separator between blocks.
    buf : str, path or file-like, optional
        Where to write. The text is returned when None.
    apply : callable, optional
        Deprecated hook returning a DataFrame with a ``blocks`` column. It runs
        the original group-by-group algorithm, and all other keywords are passed
        to it.
    random_state : int or numpy.random.Generator, optional
        Seed for `sample`. A fresh generator is used when None.
    **to_string_kwargs
        Options of :func:`rotifer.pandas.to_string` applied to the rows: ``index``,
        ``na_rep``, ``float_format``, ``align`` and ``insert``.

    Returns
    -------
    str or None
        The report as text, or None when `buf` is given.

    Raises
    ------
    TypeError
        An unknown keyword was passed (unless `apply` is used), or `header` is
        written for the old ``(df, groupby, stats)`` signature.
    KeyError
        A grouping column or a column used by `stats` is missing.

    Warns
    -----
    FutureWarning
        When `apply` is given.

    See Also
    --------
    rotifer.pandas.yaml_blocks.to_yaml : The same blocks as YAML with markdown tables.
    default_header : The default `header`.
    plan_blocks : The grouping, ordering and sampling shared by both writers.

    Notes
    -----
    Rows whose group key is missing (NaN or None) belong to no block. `colsep`
    separates the columns inside a row and `sep` separates whole blocks. The
    ordering and sampling rules are those of :func:`plan_blocks`: `sample` is drawn
    per block before `sortrows` is applied.

    Examples
    --------
    >>> import pandas as pd
    >>> from rotifer.pandas.blocks import to_blocks
    >>> df = pd.DataFrame({'g': ['A', 'A', 'B'], 'id': ['p1', 'p2', 'p3'],
    ...                    'kb': [10, 20, 30]})
    >>> print(to_blocks(df, 'g', stats={'n': ('id', 'nunique')}))
    #
    # g=A; n=2
      : id : kb
    0 : p1 : 10
    1 : p2 : 20
    #
    # g=B; n=1
      : id : kb
    2 : p3 : 30
    """
    import pyarrow as pa
    if apply is None:
        unknown = [k for k in to_string_kwargs if k not in _TO_STRING_OPTIONS]
        if unknown:
            raise TypeError(f'to_blocks() got an unexpected keyword argument {unknown[0]!r}')
    else:
        _deprecated('to_blocks(apply=...)', 'to_blocks now builds the blocks itself; '
                    'customise them with header, columns and stats', stacklevel=2)
        return _to_blocks_apply(df, groupby, buf, header, apply, colsep, sortblocks, sep,
                                {k: v for k, v in dict(columns=columns, stats=stats,
                                                       sortrows=sortrows, sample=sample).items()
                                 if v is not None} | to_string_kwargs)

    plan = plan_blocks(df, groupby, sortrows=sortrows, sortblocks=sortblocks, sample=sample,
                       tables={'stats': stats} if stats else None,
                       rng=np.random.default_rng(random_state))
    label = _as_label(groupby)

    body = df if columns is None else df.filter(columns)
    body = body.drop(columns=plan.groupby, errors='ignore')
    opts = {k: to_string_kwargs[k] for k in _TO_STRING_OPTIONS if k in to_string_kwargs}
    align = opts.get('align', 'left')
    _check_align(align)
    fast = 'insert' not in opts or not opts['insert']
    if fast:
        names, arrays = _text_columns(body, opts.get('na_rep', ' '), opts.get('float_format'),
                                      opts.get('index', True))

    stat_table = plan.tables.get('stats')
    if stat_table is not None:
        stat_names = list(stat_table.columns)
        stat_values = [stat_table.iloc[:, j].tolist() for j in range(stat_table.shape[1])]

    blocks = []
    for i, key in enumerate(plan.keys):
        idx = plan.positions[key]
        values = {} if stat_table is None else {
            c: _py(v[i]) for c, v in zip(stat_names, stat_values)}
        if fast:
            rows = _join_rows(names, arrays, colsep, align, taken=pa.array(idx))
        else:
            rows = to_string(body.take(idx), sep=colsep, **opts)
        blocks.append(_call_header(header, key, label, values) + "\n" + rows)

    return _write_sink(sep.join(blocks), buf)


def _to_blocks_apply(df, groupby, buf, header, apply, colsep, sortblocks, sep, kwargs):
    """Deprecated path of to_blocks: call the user-supplied `apply` once per group."""
    labels, scalar = _group_labels(df, groupby)
    results = []
    # The caller was warned once; nested deprecated hooks must not warn per block.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FutureWarning)
        for key, block in df.groupby(labels[0] if scalar else labels, sort=False):
            block = block.drop(columns=labels)
            object.__setattr__(block, 'name', key)
            results.append(apply(block, groupby=groupby, header=header, colsep=colsep, **kwargs))
    blocks = pd.concat(results, ignore_index=True) if results else pd.DataFrame({'blocks': []})

    if sortblocks:
        sortcols = [c for c in sortblocks if c in blocks]
        blocks = blocks.sort_values(sortcols, ascending=[sortblocks[c] for c in sortcols])

    return _write_sink(sep.join(blocks.blocks), buf)
