"""Helpers that normalise column labels and group keys."""

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
