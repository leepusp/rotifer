__doc__ = """
Older helpers for regions in one-dimensional coordinates.

- ``complement`` lists the gaps between the regions of each sequence.
- ``filter_nonoverlapping_regions`` and ``architecture`` are deprecated
  wrappers around ``rotifer.interval.architecture.select_regions``; they
  keep their old arguments and output for existing callers.
- ``config`` holds their per-format presets (``config['hmmer']``, ...),
  which can be overridden by the user's rotifer configuration.
"""

import math
import warnings
from collections import OrderedDict

import numpy as np
import pandas as pd

import rotifer
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)
config = loadConfig(__name__, defaults = {
    'default': {
        'reference': ['query'],
        'collapse': 'hit',
        'start': 'qstart',
        'end': 'qend',
        'maximum_overlap': 0.4,
        'criteria': OrderedDict({'sum_probs':False, 'probability':False, 'score':False, 'region_length':False}),
    },
    'archtsv': {
        'reference': ['ID'],
        'collapse': 'domain',
        'start': 'start',
        'end': 'end',
        'maximum_overlap': 0.4,
        'criteria': OrderedDict({'evalue':True, 'region_length':False}),
    },
    'hhsuite': {
        'reference': ['query'],
        'collapse': 'hit',
        'start': 'qstart',
        'end': 'qend',
        'maximum_overlap': 0.4,
        'criteria': OrderedDict({'sum_probs':False, 'probability':False, 'score':False, 'region_length':False}),
    },
    'hmmer': {
        'reference': ['sequence'],
        'collapse': 'model',
        'start': 'estart',
        'end': 'eend',
        'maximum_overlap': 0.4,
        'criteria': OrderedDict({'evalue':True, 'region_length':False}),
    },
})

_REMOVAL = "2026-11-30 (issue #173)"
"""When the deprecated wrappers go away, and the issue tracking it."""

def _select(df, reference, start, end, maximum_overlap, criteria):
    """
    Run ``select_regions`` with the old arguments and rules.

    Parameters
    ----------
    df : pandas.DataFrame
        Regions, one per row. It is not modified.
    reference : str or list of str
        Column(s) naming the sequence of each region.
    start, end : str
        Columns with the first and last position of each region.
    maximum_overlap : int or float
        Below 1, a fraction of the shorter region; from 1 up, a number of
        residues (a non-integer number is rounded down, since overlaps are
        whole residues).
    criteria : dict or None
        Ranking columns mapped to True (smallest first) or False (largest
        first). None means ``config['default']['criteria']``.

    Returns
    -------
    pandas.DataFrame
        The kept rows, as ``select_regions`` returns them with
        ``nested=False`` (a region inside a kept region is always
        dropped, as these functions always did).
    """
    from rotifer.interval.architecture import select_regions
    if isinstance(reference, str):
        reference = [reference]
    if criteria is None:
        criteria = config['default']['criteria']
    if maximum_overlap < 1:
        overlap = float(maximum_overlap)
    else:
        overlap = int(math.floor(maximum_overlap))
    return select_regions(df, reference=list(reference), start=start, end=end,
                          criteria=dict(criteria), overlap=overlap, nested=False)

def filter_nonoverlapping_regions(df,
        reference       = config['default']['reference'],
        start           = config['default']['start'],
        end             = config['default']['end'],
        maximum_overlap = config['default']['maximum_overlap'],
        criteria        = None,
        *args, **kwargs
    ):
    """
    Select the best set of non-overlapping regions of each sequence.

    Deprecated: use ``rotifer.interval.architecture.select_regions`` with
    ``nested=False``, which this function calls.

    Regions are ranked by ``criteria`` within each sequence, then by the
    smallest start, then by input order, and taken best first; a region is
    dropped when it overlaps a kept region by more than
    ``maximum_overlap`` or lies inside it.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per region. It is not modified.
    reference : str or list of str, default ['query']
        Column(s) naming the sequence of each region.
    start : str, default 'qstart'
        Column with the first position of each region.
    end : str, default 'qend'
        Column with the last position of each region.
    maximum_overlap : int or float, default 0.4
        Below 1, the largest fraction of the shorter region two regions
        may share; from 1 up, the largest number of shared residues.
    criteria : dict, optional
        Ranking columns, in order of precedence, mapped to True to prefer
        the smallest values or False to prefer the largest.
        ``'region_length'`` is ``end - start + 1`` unless ``df`` has such
        a column. None means ``config['default']['criteria']`` (largest
        sum_probs, probability, score and region length).
    *args, **kwargs
        Ignored. They let callers unpack a ``config`` preset, which also
        has a ``collapse`` key.

    Returns
    -------
    pandas.DataFrame
        The kept rows of ``df`` with all its columns and index labels,
        sorted by ``reference``, ``start`` and ``end``.

    Warns
    -----
    FutureWarning
        Always: this function will be removed.

    See Also
    --------
    rotifer.interval.architecture.select_regions : the function to use
        instead.

    Examples
    --------
    >>> import warnings
    >>> import pandas as pd
    >>> df = pd.DataFrame({'sequence': ['p1', 'p1'], 'model': ['A', 'B'],
    ...                    'estart': [1, 50], 'eend': [60, 100], 'evalue': [1e-10, 1e-5]})
    >>> with warnings.catch_warnings():
    ...     warnings.simplefilter('ignore', FutureWarning)
    ...     filter_nonoverlapping_regions(df, **config['hmmer'])['model'].tolist()
    ['A', 'B']
    """
    warnings.warn(
        "rotifer.interval.utils.filter_nonoverlapping_regions is deprecated and will be removed on "
        + _REMOVAL + "; use rotifer.interval.architecture.select_regions(..., nested=False) instead",
        FutureWarning, stacklevel=2)
    return _select(df, reference, start, end, maximum_overlap, criteria)

def architecture(df,
        reference       = config['default']['reference'],
        collapse        = config['default']['collapse'],
        start           = config['default']['start'],
        end             = config['default']['end'],
        maximum_overlap = config['default']['maximum_overlap'],
        criteria        = None,
    ):
    """
    Join the best non-overlapping regions of each sequence into one string.

    Deprecated: use ``rotifer.interval.architecture.domain2architecture``
    for architecture tables, or ``select_regions`` followed by a groupby.

    Regions are selected as ``filter_nonoverlapping_regions`` does, then
    the values of ``collapse`` of each sequence are joined with ``+`` in
    coordinate order. Nested regions are never kept, so there are no
    brackets.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per region. It is not modified.
    reference : str or list of str, default ['query']
        Column(s) naming the sequence of each region.
    collapse : str, default 'hit'
        Column with the region names to join.
    start, end, maximum_overlap, criteria
        As in ``filter_nonoverlapping_regions``.

    Returns
    -------
    pandas.DataFrame
        The ``reference`` column(s) and ``architecture`` (str), one row per
        sequence, sorted by ``reference``, with a 0..n-1 index. Empty, with
        these columns, when ``df`` is empty.

    Warns
    -----
    FutureWarning
        Always: this function will be removed.

    See Also
    --------
    rotifer.interval.architecture.domain2architecture : architecture table
        with regions and E-values.

    Examples
    --------
    >>> import warnings
    >>> import pandas as pd
    >>> df = pd.DataFrame({'sequence': ['p1', 'p1', 'p2'], 'model': ['A', 'B', 'C'],
    ...                    'estart': [1, 50, 5], 'eend': [60, 100, 80],
    ...                    'evalue': [1e-10, 1e-5, 1e-4]})
    >>> with warnings.catch_warnings():
    ...     warnings.simplefilter('ignore', FutureWarning)
    ...     architecture(df, **config['hmmer'])
      sequence architecture
    0       p1          A+B
    1       p2            C
    """
    warnings.warn(
        "rotifer.interval.utils.architecture is deprecated and will be removed on "
        + _REMOVAL + "; use rotifer.interval.architecture.domain2architecture, "
        "or select_regions and a groupby, instead",
        FutureWarning, stacklevel=2)
    if isinstance(reference, str):
        reference = [reference]
    reference = list(reference)
    kept = _select(df, reference, start, end, maximum_overlap, criteria)
    if kept.empty:
        return kept[reference].assign(architecture=pd.Series(dtype=object)).reset_index(drop=True)
    return (kept.groupby(reference)[collapse]
                .agg(lambda x: "+".join(str(y) for y in x))
                .rename('architecture')
                .reset_index())

def complement(df, include=True, exclude=False, add=False, reference='ID', start='start', end='end',
               length='length'):
    """
    List the gaps between the regions of each sequence.

    Each row of ``df`` is a region of the sequence named in ``reference``,
    from ``start`` to ``end`` (inclusive, 1-based). Overlapping and nested
    regions are merged first, so each position is either covered or in
    exactly one gap. Gaps are the stretches before the first region (from
    position 1), between regions, and after the last region up to the
    sequence length, when the length is known.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per region. It is not modified. Rows with a missing start
        or end, or with ``end < start``, cover nothing.
    include : bool or list of str, default True
        Columns of the result. True: the columns of ``df``, in its order
        (columns not computed here are missing values). A list: these
        columns, in this order. False: only the computed columns
        (``reference``, ``start``, ``end`` and the gap length).
    exclude : list of str or False, default False
        Columns to drop from the result. Names not in the result are
        ignored.
    add : dict or False, default False
        Columns to set on every gap, as ``{column: value}``, e.g.
        ``{'domain': 'unk'}``. Applied after ``include`` and ``exclude``.
    reference : str or list of str, default 'ID'
        Column(s) naming the sequence of each region.
    start, end : str, default 'start', 'end'
        Columns with the first and last position of each region.
    length : str, default 'length'
        Column with the length of each region's sequence (the first
        non-missing value of each sequence is used). When ``df`` has it,
        the gap after the last region is listed, and the same column in
        the result holds each gap's length. When it does not, there is no
        trailing gap and no length column.

    Returns
    -------
    pandas.DataFrame
        One row per gap, sorted by ``reference`` and ``start``, with a
        0..n-1 index; ``start`` and ``end`` are int64. Sequences with no
        gap have no row. Empty when there are no gaps.

    See Also
    --------
    rotifer.devel.alpha.sequence.complete.complete : adds these gaps to
        the regions.

    Examples
    --------
    >>> import pandas as pd
    >>> regions = pd.DataFrame({'ID': ['s1', 's1', 's1'], 'start': [10, 20, 50],
    ...                         'end': [30, 25, 80], 'length': [100, 100, 100]})
    >>> complement(regions)
       ID  start  end  length
    0  s1      1    9       9
    1  s1     31   49      19
    2  s1     81  100      20
    """
    if isinstance(reference, str):
        reference = [reference]
    reference = list(reference)
    n = len(df)
    has_length = length in df.columns

    # One code per sequence, in sorted order of the reference columns
    if n:
        seq = df.groupby(reference, sort=True, dropna=False).ngroup().to_numpy()
    else:
        seq = np.zeros(0, dtype=np.int64)
    nseq = int(seq.max()) + 1 if n else 0
    first_row = np.zeros(nseq, dtype=np.int64)
    first_row[seq[::-1]] = np.arange(n)[::-1]

    first = pd.to_numeric(df[start], errors='coerce').to_numpy(dtype='float64')
    last = pd.to_numeric(df[end], errors='coerce').to_numpy(dtype='float64')
    valid = ~np.isnan(first) & ~np.isnan(last) & (last >= first)
    regions = pd.DataFrame({'seq': seq[valid], 'start': first[valid].astype(np.int64),
                            'end': last[valid].astype(np.int64)})
    regions = regions.sort_values(['seq', 'start'], kind='stable')

    # Gaps before and between regions: a region starts after everything
    # before it, in the same sequence, has ended.
    reach = regions.groupby('seq')['end'].cummax()
    before = reach.groupby(regions['seq']).shift(1).fillna(0).astype(np.int64)
    inner = regions['start'] > before + 1
    pieces = [pd.DataFrame({'seq': regions['seq'][inner], 'start': (before + 1)[inner],
                            'end': regions['start'][inner] - 1})]

    # Gap after the last region, when the sequence length is known
    if has_length:
        size = pd.Series(pd.to_numeric(df[length], errors='coerce').to_numpy()).groupby(seq).first()
        covered = reach.groupby(regions['seq']).max().reindex(range(nseq), fill_value=0)
        size = size.reindex(range(nseq))
        tail = size.notna() & (covered + 1 <= size)
        pieces.append(pd.DataFrame({'seq': np.flatnonzero(tail.to_numpy()),
                                    'start': (covered + 1)[tail].to_numpy(),
                                    'end': size[tail].to_numpy().astype(np.int64)}))

    gaps = pd.concat(pieces, ignore_index=True).astype({'seq': np.int64, 'start': np.int64, 'end': np.int64})
    gaps = gaps.sort_values(['seq', 'start'], kind='stable').reset_index(drop=True)

    result = df[reference].iloc[first_row[gaps['seq'].to_numpy()]].reset_index(drop=True)
    result[start] = gaps['start'].to_numpy()
    result[end] = gaps['end'].to_numpy()
    if has_length:
        result[length] = (gaps['end'] - gaps['start'] + 1).to_numpy()

    if isinstance(include, (list, tuple)):
        result = result.reindex(columns=list(include))
    elif include:
        result = result.reindex(columns=df.columns)

    if isinstance(exclude, (list, tuple, set)):
        result = result.drop(columns=[c for c in exclude if c in result.columns])

    if isinstance(add, dict):
        for column, value in add.items():
            result[column] = value

    return result
