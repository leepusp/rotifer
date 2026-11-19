"""Block output: the shared block engine (plan_blocks), default_header and to_blocks."""

import types
import warnings
import numpy as np
import pandas as pd
from collections.abc import Mapping
from dataclasses import dataclass

from ._labels import _as_label, _fmt_label, _group_key, _group_labels, _py
from ._io import _write_sink
from .text import to_string, pad_columns

def _resolve_stats(df, stats: Mapping):
    """{col: aggfunc_name} -> {col: value}. `stats` is any dictionary-like object
    (anything implementing the Mapping protocol, i.e. `.items()`), mapping a
    column label to an aggregation name. Array-like keys are normalized to tuples,
    so multi-level column labels stay hashable and round-trippable.

    Raises KeyError when a column is missing, rather than silently emitting the
    aggregation name in place of a value. Idempotent: a string the column has no
    attribute for is passed through unchanged, so re-resolving an already-resolved
    dict is a no-op -- which is also why a misspelled aggregation on an existing
    column cannot be flagged here, being indistinguishable from an aggregation
    that returned a string."""
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
    """Everything a block writer needs, computed once for the whole frame.

    groupby   : list of grouping column labels.
    keys      : group keys in block order. A key is a bare value when one column
                was named, otherwise a tuple with one value per label.
    positions : {key: row positions}, already ordered by sortrows and sampled.
    tables    : {name: DataFrame} of per-block aggregates, row i belonging to
                keys[i], in the order the caller named the tables.
    callbacks : {name: callable(block DataFrame) -> DataFrame} for table specs
                that need the block itself.
    """
    groupby: list
    keys: list
    positions: dict
    tables: dict
    callbacks: dict


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
    """Aggregate one table's spec for every block in a single pass.

    spec maps an output name to one of
      * an aggregation (name or callable) of the column of the same name,
        {'species': 'nunique'};
      * a named aggregation, {'n_species': ('species', 'nunique')};
      * a plain value, repeated for every block, {'n': 7}.
    Returns a DataFrame indexed by group key, columns in spec order."""
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
    frame = grouped.agg(**named) if named else pd.DataFrame(index=grouped.size().index)
    for i, (_, _, agg) in enumerate(entries):
        if f'_{i}' not in named:
            frame[f'_{i}'] = agg
    frame = frame[[f'_{i}' for i in range(len(entries))]]
    frame.columns = pd.Index([out for out, _, _ in entries], tupleize_cols=False)
    return frame


def _split_tables(grouped, tables):
    """Sort each table spec into the two ways a table can be produced.

    A mapping is an aggregation, evaluated once for every block; a callable
    receives each block as a DataFrame and returns its table. Aggregations are
    much the cheaper of the two, since they never materialise a per-block
    DataFrame."""
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
    """Group `df` once and decide everything about the blocks except their text.

    Blocks are ordered by first appearance of their key, unless sortblocks
    ({column: ascending}) reorders them; it may name a grouping column or a column
    of any aggregation table (the first table defining a column wins). sortrows
    orders rows inside each block and never moves a block. sample keeps at most
    that many rows per block, drawn with `rng` (a numpy Generator, fresh by
    default) and then sorted like any other rows. Columns named in sortrows or
    sortblocks that do not exist are ignored.

    tables is {name: spec}; see _group_agg for the mapping forms, a callable
    spec lands in BlockPlan.callbacks."""
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
    # stacklevel counts from the deprecated function: 2 is its caller.
    warnings.warn(f'{what} is deprecated and will be removed on {REMOVAL}; {instead}.',
                  FutureWarning, stacklevel=stacklevel + 1)


def default_header(key, groupby, stats=None, **kwargs):
    """Header of one block: ``#`` then ``# <groupby>=<key>; <stat>=<value>; ...``.

    key     : the block's group key, a bare value or a tuple matching `groupby`.
    groupby : the grouping label, or the tuple of labels.
    stats   : {name: value} for this block, or None for none.

    to_blocks calls ``header(key, groupby, stats)``; a replacement header takes
    the same three arguments and returns the text above the block's rows.

    Calling it the old way, ``default_header(df, groupby=..., stats=...)`` with the
    block's DataFrame, still works but is deprecated."""
    if isinstance(key, pd.DataFrame):
        _deprecated('default_header(df, ...) with a DataFrame',
                    'pass the group key and a {name: value} dict: '
                    'default_header(key, groupby, stats)', stacklevel=2)
        return _legacy_header(key, groupby, stats)
    tail = ''.join(f'; {_fmt_label(c)}={v}' for c, v in (stats or {}).items())
    return f'#\n# {_fmt_label(groupby)}={_fmt_label(key)}{tail}'


def _legacy_header(df, groupby, stats):
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
    """Deprecated: the old per-block hook of to_blocks, returning a one-row
    DataFrame whose `blocks` column holds the block text."""
    _deprecated('default_block_rows', 'to_blocks builds the blocks itself now', stacklevel=2)
    statDict = _resolve_stats(df, stats or {})
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
    block = [headerStr + "\n" + to_string(block, sep=colsep, **kwargs)]
    return pd.DataFrame({'blocks': block, **statDict})


def _string_columns(body, na_rep, float_format, index):
    """Text of every cell of `body`, converted once for the whole frame so a block
    is later cut out with an arrow `take`. Mirrors text.to_string cell for cell."""
    import pyarrow as pa
    columns = []
    if index:
        columns.append((str(body.index.name or ""), body.index.to_series().astype("string")))
    for c in body.columns:
        s = body[c]
        if float_format is not None and pd.api.types.is_float_dtype(s.dtype):
            text = s.map(lambda v: na_rep if pd.isna(v) else float_format % v).astype("string")
        else:
            text = s.astype("string")
        columns.append((str(c), text.fillna(na_rep)))
    names = [name for name, _ in columns]
    return names, [pa.Array.from_pandas(text, type=pa.string()) for _, text in columns]


def _render_rows(names, arrays, taken, sep, align):
    """The rows `taken` of pre-converted columns, padded and joined like to_string."""
    import pyarrow as pa
    import pyarrow.compute as pc
    just = str.rjust if align == "right" else str.ljust
    padded, widths = pad_columns([a.take(taken) for a in arrays], names, align=align)
    header = sep.join(just(name, w) for name, w in zip(names, widths))
    if not padded:
        return header
    lines = pc.binary_join_element_wise(*padded, pa.scalar(sep))
    return header + "\n" + "\n".join(lines.to_pylist())


def to_blocks(df, groupby, *, columns=None, stats=None, sortrows=None,
              sortblocks=None, sample=None, header=default_header, colsep=' : ',
              sep='\n', buf=None, apply=None, **to_string_kwargs):
    """Write `df` as text blocks, one per group: a header, then the block's rows.

    groupby    : column label, or list of labels, to group by. Blocks come out in
                 order of first appearance.
    columns    : columns shown in each block's rows; every column but the grouping
                 ones when None.
    stats      : {column: aggregation} or {name: (column, aggregation)} computed
                 once per block and handed to `header`.
    sortrows   : {column: ascending} ordering the rows inside each block.
    sortblocks : {column: ascending} ordering the blocks; it may name a grouping
                 column or a stats column.
    sample     : rows per block, drawn at random when a block is larger.
    header     : callable ``header(key, groupby, stats) -> str``, where stats is a
                 {name: value} dict for the block. This replaces the old
                 ``header(df, groupby=..., stats=...)`` signature, and a header
                 written for it must be changed.
    colsep, sep: between columns in a row; between blocks.
    buf        : path or file-like to write to; the text is returned when None.
    apply      : deprecated, see below.

    Any other keyword must be one of the to_string options: index, na_rep,
    float_format, align, insert. Anything else raises TypeError.

    apply, the old hook returning a DataFrame with a `blocks` column, is
    deprecated and runs the previous group-by-group algorithm."""
    import pyarrow as pa
    unknown = [k for k in to_string_kwargs if k not in _TO_STRING_OPTIONS]
    if unknown:
        raise TypeError(f'to_blocks() got an unexpected keyword argument {unknown[0]!r}')

    if apply is not None:
        _deprecated('to_blocks(apply=...)', 'to_blocks now builds the blocks itself; '
                    'customise them with header, columns and stats', stacklevel=2)
        return _to_blocks_apply(df, groupby, buf, header, apply, colsep, sortblocks, sep,
                                {k: v for k, v in dict(columns=columns, stats=stats,
                                                       sortrows=sortrows, sample=sample).items()
                                 if v is not None} | to_string_kwargs)

    plan = plan_blocks(df, groupby, sortrows=sortrows, sortblocks=sortblocks, sample=sample,
                       tables={'stats': stats} if stats else None)
    label = _as_label(groupby)

    body = df if columns is None else df.filter(columns)
    body = body.drop(columns=plan.groupby, errors='ignore')
    opts = {k: to_string_kwargs[k] for k in _TO_STRING_OPTIONS if k in to_string_kwargs}
    fast = 'insert' not in opts or not opts['insert']
    if fast:
        names, arrays = _string_columns(body, opts.get('na_rep', ' '), opts.get('float_format'),
                                        opts.get('index', True))

    stat_table = plan.tables.get('stats')
    if stat_table is not None:
        stat_names = list(stat_table.columns)
        stat_values = [stat_table.iloc[:, j].to_numpy() for j in range(stat_table.shape[1])]

    blocks = []
    for i, key in enumerate(plan.keys):
        idx = plan.positions[key]
        values = {} if stat_table is None else {
            c: _py(v[i]) for c, v in zip(stat_names, stat_values)}
        if fast:
            rows = _render_rows(names, arrays, pa.array(idx), colsep, opts.get('align', 'left'))
        else:
            rows = to_string(body.take(idx), sep=colsep, **opts)
        blocks.append(header(key, label, values) + "\n" + rows)

    return _write_sink(sep.join(blocks), buf)


def _to_blocks_apply(df, groupby, buf, header, apply, colsep, sortblocks, sep, kwargs):
    """The old algorithm for a user-supplied `apply`: one call per group."""
    labels, scalar = _group_labels(df, groupby)
    results = []
    for key, block in df.groupby(labels[0] if scalar else labels, sort=False):
        block = block.drop(columns=labels)
        object.__setattr__(block, 'name', key)
        results.append(apply(block, groupby=groupby, header=header, colsep=colsep, **kwargs))
    blocks = pd.concat(results, ignore_index=True) if results else pd.DataFrame({'blocks': []})

    if sortblocks:
        sortcols = [c for c in sortblocks if c in blocks]
        blocks = blocks.sort_values(sortcols, ascending=[sortblocks[c] for c in sortcols])

    return _write_sink(sep.join(blocks.blocks), buf)
