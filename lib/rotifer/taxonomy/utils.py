"""Reduce NCBI classifications to short lineages of preferred clades.

An NCBI classification lists every clade from the top of the tree down to
the organism, joined by ``"; "``::

    Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; ...

Most of those clades are too fine or too obscure to group proteins by.
A *lineage* keeps only the clades on a list of preferred taxa, lowercased
and joined by ``">"``::

    bacteria>pseudomonadota>gammaproteobacteria>enterobacterales

The list ships with rotifer as ``share/rotifer/taxonomy/taxonomy.txt``.
:class:`rotifer.genome.data.NeighborhoodDF` (``update_lineage=``) and the
NCBI, UniProt and ete3 taxonomy cursors in :mod:`rotifer.db` fill their
``lineage`` column with :func:`lineage`.

Functions
---------
lineage
    Reduce classifications to lineages.
load_preferred_taxa
    Read the list of preferred taxa.
loadPreferredTaxa
    Deprecated name of :func:`load_preferred_taxa`.
"""

import warnings

import numpy as np
import pandas as pd

from rotifer.core.functions import findDataFiles

DEFAULT_TAXA = ':taxonomy.taxonomy.txt'
"""Where :func:`load_preferred_taxa` looks by default, in
:func:`rotifer.core.functions.findDataFiles` notation."""

REMOVAL = '2026-11-30'
_TRACKING_ISSUE = '160'


def _deprecated(what, instead, stacklevel=2):
    """Warn that `what` will be removed on REMOVAL.

    Parameters
    ----------
    what : str
        The deprecated name, as the user wrote it.
    instead : str
        What to use instead, as a phrase that completes the sentence.
    stacklevel : int, default 2
        Counted from the deprecated function: 2 points the warning at the
        line that called it.

    Warns
    -----
    FutureWarning
        Always.
    """
    warnings.warn(f'{what} is deprecated and will be removed on {REMOVAL} '
                  f'(issue #{_TRACKING_ISSUE}); {instead}.',
                  FutureWarning, stacklevel=stacklevel + 1)


def _same(x):
    """Return `x` unchanged: the `case` used when it is None."""
    return x


def _missing(value):
    """Tell whether a scalar is None, NaN, NaT or ``pd.NA``."""
    if isinstance(value, str):
        return False
    try:
        flag = pd.isna(value)
    except (TypeError, ValueError):
        return False
    return flag is True or isinstance(flag, np.bool_) and bool(flag)


def load_preferred_taxa(path=DEFAULT_TAXA):
    """Read the list of preferred taxa.

    Each line of the file names one clade. Lines that start with ``#``
    and blank lines are skipped. Names are stripped of surrounding
    whitespace and lowercased.

    Parameters
    ----------
    path : str, default ``':taxonomy.taxonomy.txt'``
        A file path, or a name in :func:`rotifer.core.functions.findDataFiles`
        notation. The default finds ``taxonomy/taxonomy.txt`` under
        ``~/.rotifer/share``, under ``$ROTIFER_DATA`` and under rotifer's
        ``share/rotifer``, and reads **every** copy it finds: a user file
        adds names to the shipped list, it does not replace it.

    Returns
    -------
    set of str
        The lowercased clade names. It is empty when no file is found.

    See Also
    --------
    lineage : Uses this set when it is not given one.

    Examples
    --------
    >>> taxa = load_preferred_taxa()
    >>> {'bacteria', 'pseudomonadota', 'enterobacterales'} <= taxa
    True
    >>> 'escherichia' in taxa
    False
    """
    taxa = set()
    for filename in findDataFiles(path, all=True):
        with open(filename, 'rt') as handle:
            for line in handle:
                name = line.strip()
                if name and not name.startswith('#'):
                    taxa.add(name.lower())
    return taxa


def loadPreferredTaxa():
    """Read the default list of preferred taxa (deprecated).

    Deprecated: use :func:`load_preferred_taxa`. This name will be
    removed on 2026-11-30.

    Returns
    -------
    set of str
        Same as ``load_preferred_taxa()``.

    Warns
    -----
    FutureWarning
        Always.
    """
    _deprecated('loadPreferredTaxa', 'use load_preferred_taxa', stacklevel=2)
    return load_preferred_taxa()


def lineage(data, preferred_taxa=None, column='classification', insep='; ', outsep='>', case=str.lower):
    """Reduce classifications to lineages of preferred clades.

    Each classification is split on `insep`. Each name is stripped of
    surrounding whitespace and passed through `case`, and the names found
    in `preferred_taxa` are joined with `outsep`, in their original order.
    A classification with no preferred clade gives the empty string.

    Parameters
    ----------
    data : pandas.Series, pandas.DataFrame or iterable of str
        The classifications, one per row or item. For a DataFrame, the
        `column` column is used. Any other iterable (list, tuple, generator)
        is read item by item. A single ``str`` is rejected rather than read
        as a sequence of characters. Values that are not strings and not
        missing are converted with ``str``. The input is not modified.
    preferred_taxa : iterable of str, optional
        The clades to keep. Each name is stripped and passed through `case`
        before matching, so ``['Bacteria']`` matches ``'bacteria'``.
        None (the default) reads :func:`load_preferred_taxa`. An empty
        iterable keeps nothing, so every lineage is ``''``.
    column : str, default 'classification'
        The column that holds the classifications when `data` is a DataFrame.
    insep : str, default '; '
        Separator between the names of a classification.
    outsep : str, default '>'
        Separator between the names of a lineage.
    case : callable or None, default str.lower
        Applied to every name, on both sides, before matching; the names in
        the result are its output. None compares and returns names as they
        are, so case matters.

    Returns
    -------
    pandas.Series or list of str
        A :class:`pandas.Series` when `data` is a Series or a
        :class:`pandas.DataFrame`, otherwise a list.

        - The Series has one value per input row, in input order. Its index
          has the input's labels and name (a RangeIndex comes back as an
          int64 Index), and its name is the input Series' name or `column`.
          Its dtype is pandas' default string dtype (``str``); an empty
          input gives an empty ``object`` Series.
        - The list has one item per input item, in input order.
        - A missing classification (None, NaN or ``pd.NA``) gives a missing
          lineage: NaN in a Series, None in a list.

    Raises
    ------
    TypeError
        If `data` is a ``str``.
    KeyError
        If `data` is a DataFrame without `column`.

    See Also
    --------
    load_preferred_taxa : The default `preferred_taxa`.
    rotifer.genome.data.NeighborhoodDF : ``update_lineage=`` adds a
        ``lineage`` column with this function.

    Notes
    -----
    Each distinct classification is reduced once, and the results are
    mapped back onto the rows. Tables with many rows per organism, such as
    gene neighbourhoods, cost about as much as their number of organisms.

    Examples
    --------
    >>> import pandas as pd
    >>> lineage(['Bacteria; Pseudomonadota; Gammaproteobacteria; '
    ...          'Enterobacterales; Enterobacteriaceae; Escherichia'])
    ['bacteria>pseudomonadota>gammaproteobacteria>enterobacterales']

    A Series keeps its index, and a missing classification stays missing:

    >>> s = pd.Series(['Archaea; Methanobacteriota; Methanobacteria', None,
    ...                'Eukaryota; Metazoa; Chordata; Mammalia'], index=[10, 5, 7])
    >>> lineage(s)
    10                archaea>methanobacteria
    5                                     NaN
    7     eukaryota>metazoa>chordata>mammalia
    dtype: str

    Your own list of clades, in any case:

    >>> lineage(['Bacteria; Bacillota; Bacilli'], preferred_taxa=['Bacteria', 'Bacilli'])
    ['bacteria>bacilli']
    """
    if isinstance(data, str):
        raise TypeError('data must be a Series, a DataFrame or an iterable of '
                        'str, not a single str')
    if case is None:
        case = _same
    if preferred_taxa is None:
        preferred_taxa = load_preferred_taxa()
    keep = {case(str(name).strip()) for name in preferred_taxa}

    def kept(value):
        names = (case(name.strip()) for name in value.split(insep))
        return [name for name in names if name in keep]

    def reduce(value):
        return outsep.join(kept(value if isinstance(value, str) else str(value)))

    if isinstance(data, pd.DataFrame):
        data = data[column]
    if isinstance(data, pd.Series):
        index = data.index
        if len(data) == 0:
            return pd.Series([], dtype=object, index=index, name=data.name)
        if isinstance(index, pd.RangeIndex):
            index = pd.Index(index.to_numpy(), name=index.name)
        # str values first, as the old code did, so 1 and 1.0 stay apart;
        # missing values stay missing and factorize gives them code -1
        codes, uniques = pd.factorize(data.astype(str).to_numpy(dtype=object))
        table = np.empty(len(uniques) + 1, dtype=object)
        table[:-1] = [reduce(value) for value in uniques]
        table[-1] = None
        return pd.Series(table[codes], index=index, name=data.name, dtype='str')
    return [None if _missing(item) else reduce(item) for item in data]
