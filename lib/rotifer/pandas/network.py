"""Neighbourhood network built from a DataFrame (will move to rotifer.graph)."""

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
