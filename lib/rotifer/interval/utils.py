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

def complement(df, include=True, exclude=False, add=False, reference='ID', start='start', end='end'):
    
    '''
    Given a DataFrame describing regions in one-dimensional
    reference systems, return a similar DataFrame describing
    the set of regions required to fully cover all the length
    of the references.

    The input table is expected to include one or more columns
    describing the reference object and a pair of integer columns
    corresponding to the start and end positions of each region.

    If a column with the length of each interval is supplied, the
    last interval and the length of each fragment will be included
    in the returned dataframe.

    Parameters
    ----------
    df : Pandas dataframe
       The dataframe must have at least a reference column and coordinates
    include: list
    exclude: list
    add: dict
    reference : (list of) strings, default ['sequence']
        Name of the column(s) storing sequence identifiers
    start : string, default 'start'
        Name of the column with the first coordinate of each interval
    end : string, default 'end'
        Name of the column with the last coordinate of each interval
    '''
    if isinstance(reference, str):
        reference = [ reference ]

    regions = df.sort_values([*reference, start, end])
    if 'length' not in regions.columns:
        regions['length'] = regions.ID.map(regions.set_index(reference).end.to_dict())
        regions.loc[regions.length.isna(), 'length'] = regions.loc[regions.length.isna(), 'end']

    first = pd.DataFrame({
        'ID':     regions.ID.tolist(),
        'start':  np.where(regions.ID != regions.ID.shift(1),1,regions.end.shift(1) + 1).tolist(),
        'end':    (regions.start - 1).tolist(),
    }).astype({'start':int}).eval('length = end - start + 1').query('end > 0')

    uncovered = (regions.groupby(['ID','length']).end.max()+1)
    uncovered = uncovered.reset_index()
    uncovered = uncovered.query('end < length')
    uncovered = uncovered.rename({'end':'start'}, axis=1)
    uncovered = uncovered.eval('end = length').eval('length = end - start + 1')
    uncovered = pd.concat([ first, uncovered ]).sort_values(['ID','start','end'])

    if isinstance(include, list):
        uncovered = uncovered.filter(include).reindex(df.columns, axis=1)
    elif include:
        uncovered = uncovered.reindex(df.columns, axis=1)

    if isinstance(exclude, list):
        uncovered = uncovered.drop(exclude, axis=1)

    if isinstance(add, dict):
        for keys in values.keys():
            uncovered[keys] = values[keys]

    return uncovered

