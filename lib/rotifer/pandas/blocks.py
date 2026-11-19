"""Block output: default_header, default_block_rows and to_blocks."""

import types
import numpy as np
import pandas as pd
from collections.abc import Mapping
from dataclasses import dataclass

from ._labels import _as_label, _fmt_label, _group_key, _group_labels
from ._io import _write_sink
from .text import to_string

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
        {'lineage': 'nunique'};
      * a named aggregation, {'n_lineages': ('lineage', 'nunique')};
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


def _to_string_arrow(columns, sep, align, just):
    """Pad and join with pyarrow kernels. ~1.7x faster than the pandas path."""
    import pyarrow as pa
    import pyarrow.compute as pc

    pad = pc.utf8_lpad if align == "right" else pc.utf8_rpad

    padded, headers = [], []
    for name, body in columns:
        arr = pa.Array.from_pandas(body, type=pa.string())
        width = max(pc.max(pc.utf8_length(arr)).as_py() or 0, len(name))
        padded.append(pad(arr, width, padding=" "))
        headers.append(just(name, width))

    header = sep.join(headers)
    if not padded:
        return header
    lines = pc.binary_join_element_wise(*padded, pa.scalar(sep))
    return header + "\n" + "\n".join(lines.to_pylist())

def _to_string_pandas(columns, sep, align, just):
    """Fallback for environments without pyarrow."""
    acc = "rjust" if align == "right" else "ljust"

    padded, headers = [], []
    for name, body in columns:
        w = body.str.len().max()
        width = max(int(w) if pd.notna(w) else 0, len(name))
        padded.append(getattr(body.str, acc)(width).to_numpy().tolist())
        headers.append(just(name, width))

    header = sep.join(headers)
    if not padded:
        return header
    return header + "\n" + "\n".join([sep.join(row) for row in zip(*padded)])

def to_string(df, sep=" : ", index=True, na_rep=" ", insert=None,
                      float_format=None, align="left", **kwargs):
    just = str.rjust if align == "right" else str.ljust

    if insert:
        df = df.copy()
        for data in insert:
            if len(data) != 3:
                raise ValueError(f"insert entries must be (loc, column, value) tuples, got {data!r}")
            df.insert(loc=data[0], column=data[1], value=data[2])

    columns = []
    if index:
        columns.append((str(df.index.name or ""), df.index.to_series().astype("string")))
    for c in df.columns:
        s = df[c]
        if float_format is not None and pd.api.types.is_float_dtype(s.dtype):
            body = s.map(lambda v: na_rep if pd.isna(v) else float_format % v).astype("string")
        else:
            body = s.astype("string")
        columns.append((str(c), body.fillna(na_rep)))

    try:
        return _to_string_arrow(columns, sep, align, just)
    except ImportError:
        return _to_string_pandas(columns, sep, align, just)

def default_header(df, groupby='c80e3', stats=None, **kwargs):
    if stats is None:
        stats = {'lineage': 'nunique', 'block_id': 'nunique'}
    groupby = _as_label(groupby)
    key = _group_key(df, groupby)
    # A tuple groupby that is not itself a column label is a *sequence* of labels,
    # so pad the key back to matching arity (_group_key unwraps 1-tuples to match
    # pandas) and a parser can zip the two tuples together.
    if isinstance(groupby, tuple) and groupby not in df.columns and not isinstance(key, tuple):
        key = (key,)
    statStr = _resolve_stats(df, stats) # idempotent
    statStr = "; ".join(f'{_fmt_label(c)}={v}' for c, v in statStr.items())
    return f'#\n# {_fmt_label(groupby)}={_fmt_label(key)}; {statStr}'

def default_block_rows(df, groupby='c80e3', header=default_header, colsep=" : ",
                          columns=None, sample=None, stats=None, sortrows=None, **kwargs):
    if columns is None:
        columns = ['pid', 'c80i70', 'c100i100', 'plen', 'pfam', 'aravind',
                    'product', 'organism', 'lineage', 'classification']
    if sortrows is None:
        sortrows = {'lineage': True, 'classification': True, 'organism': True}
    if stats is None:
        stats = {'lineage':'nunique', 'block_id':'nunique', 'c80i70':'nunique', 'pid':'nunique'}

    statDict = _resolve_stats(df, stats)
    headerStr = header(df, groupby=groupby, stats=statDict, **kwargs)

    block = df.copy()
    if sample and sample < len(block):
        block = block.sample(sample)
    if sortrows:
        sortcols = [c for c in sortrows if c in block]
        if sortcols:
            block = block.sort_values(sortcols, ascending=[sortrows[c] for c in sortcols])

    block = block.filter(columns)
    block = [headerStr + "\n" + to_string(block, sep=colsep, **kwargs)]
    return pd.DataFrame({'blocks': block, **statDict})

def to_blocks(df, groupby='c80e3', buf=None, header=default_header,
              apply=default_block_rows, colsep=" : ",
              sortblocks=None, sep="\n", **kwargs):
    blocks = df.groupby(groupby, sort=False)
    blocks = blocks.apply(apply, groupby=groupby, header=header, colsep=colsep,
                          include_groups=False, **kwargs)

    if sortblocks:
        sortcols = [c for c in sortblocks if c in blocks]
        blocks = blocks.sort_values(sortcols, ascending=[sortblocks[c] for c in sortcols])

    blocks = sep.join(blocks.blocks)

    return _write_sink(blocks, buf)
