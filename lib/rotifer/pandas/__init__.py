import types
import pandas as pd
from collections.abc import Mapping
from . import *

def _group_key(df, groupby):
    """Group key inside a GroupBy.apply callback (df.name), or derived from the
    grouping columns when called standalone."""
    single = not isinstance(groupby, (list, tuple))
    groupby = (groupby,) if single else tuple(groupby) # tuple: a list is unhashable below
    if len(df):
        if not single and groupby in df.columns: # one column labelled by a tuple (multi-level columns)
            return df[groupby].iloc[0]
        elif all(c in df for c in groupby):
            key = tuple(df[c].iloc[0] for c in groupby)
            return key[0] if len(key) == 1 else key # pandas unwraps single-column groupings
    return getattr(df, 'name', None) # set by GroupBy.apply on each group

def _as_label(x):
    """Normalize an array-like column label (or list of labels) to a tuple, so it
    is hashable and survives `repr` -> `ast.literal_eval` round-tripping. Scalars,
    including strings, pass through untouched."""
    if isinstance(x, str) or not hasattr(x, '__iter__'):
        return x
    return tuple(_as_label(e) for e in x)

def _py(v):
    """numpy/pandas scalar -> builtin, so repr() yields a Python literal."""
    return v.item() if hasattr(v, 'item') else v

def _fmt_label(x):
    """Render a label or group key for the block header. Tuples are emitted as
    Python literals so a parser can recover them with ast.literal_eval; numpy
    scalars are demoted to builtins first, since repr(np.int64(1)) is not a
    literal. Scalars stay bare, keeping the common `# c80e3=A` form readable."""
    if isinstance(x, tuple):
        return repr(tuple(_fmt_label(e) if isinstance(e, tuple) else _py(e) for e in x))
    return str(_py(x))

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
            if len(data) == 3:
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
    from rotifer.pandas import functions as rpf

    blocks = df.groupby(groupby, sort=False, group_keys=False)
    blocks = blocks.apply(apply, groupby=groupby, header=header, colsep=colsep, **kwargs)

    if sortblocks:
        sortcols = [c for c in sortblocks if c in blocks]
        blocks = blocks.sort_values(sortcols, ascending=[sortblocks[c] for c in sortcols])

    blocks = sep.join(blocks.blocks)

    if buf is None:
        return blocks
    if isinstance(buf, (str, os.PathLike)):
        with open(buf, "w") as f:
            f.write(blocks)
        return None
    buf.write(blocks) # file-like object; caller closes it
    return None

def to_network(df, target=['pfam'], ftype=['CDS'], interaction=True, ignore = [], strand = True, separator="+", separator_regex=False, stats=None, replace=None, replace_regex=False):
    import numpy as np
    if isinstance(ftype, str):
        ftype = [ftype]
    if stats is None:
        stats = dict(weight=('block_id', 'count'), blocks=('block_id', 'nunique'))

    if strand:
        w = df.neighbors(before=len(df), after=len(df), strand='same')
        w['rid'] = list(range(1,len(w)+1))
        w.rid = w.rid * w.strand
        w.sort_values(['rid'], inplace=True)
    else:
        w = df.copy()

    # Building the source column
    w = w.query(f'type.isin({ftype})').reset_index(drop=True)
    w['source'] = w[target[0]]
    for col in target[1:]:
        w['source'] = np.where(w['source'].isna(), w[col], w['source'])

    # Parse target columns
    if separator is not None:
        w['source'] = w['source'].fillna("?").str.split(separator, regex=separator_regex)
        w = w.explode(column='source')

    # Rename components and cleanup
    w.source = w.source.replace(replace, regex=replace_regex)
    if ignore:
        w = w[~w.source.isin(ignore)].copy()

    # Building target data
    w['tpid'] = w['pid'].shift(-1)
    w['tblock_id'] = w['block_id'].shift(-1)
    w['target'] = w['source'].shift(-1)

    # Selecting same block rows
    w = w[w.block_id == w.tblock_id].copy()

    # Fix source target order when not restricted to the same strand
    if not strand and w.source > w.target:
        w.loc[reverse,['source','target']] = w.loc[reverse, ['target', 'source']].values

    if interaction:
        w['interaction'] = np.where(w.pid == w.tpid, 'fusion', 'neighbor')
        w = w.groupby(['source', 'target', 'interaction'])
    else:
        w = w.groupby(['source', 'target'])

    w = w.agg(**stats).reset_index()
    return w
