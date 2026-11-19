"""Private helpers that normalise column labels and group keys for block output."""

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
    literal. Scalars stay bare, keeping the common `# grp=A` form readable."""
    if isinstance(x, tuple):
        return repr(tuple(_fmt_label(e) if isinstance(e, tuple) else _py(e) for e in x))
    return str(_py(x))

def _as_list(x):
    """None -> [], a string or non-iterable scalar -> [x], any other iterable -> list(x)."""
    if x is None:
        return []
    if isinstance(x, str) or not hasattr(x, '__iter__'):
        return [x]
    return list(x)


def _group_labels(df, groupby):
    """Normalise a groupby argument to ``(labels, scalar)``.

    ``labels`` is the list of grouping column labels. ``scalar`` says the caller
    named ONE column, so block keys are bare values rather than tuples: a plain
    label, or a tuple that is itself a column of ``df`` (multi-level columns). A
    list, or a tuple that is not a column, is a sequence of labels and always
    yields tuple keys of matching arity, even for one label."""
    groupby = _as_label(groupby) if isinstance(groupby, tuple) else groupby
    if isinstance(groupby, tuple):
        if groupby in df.columns:
            return [groupby], True
        return list(groupby), False
    if isinstance(groupby, list):
        return [_as_label(c) for c in groupby], False
    return [_as_label(groupby)], True
