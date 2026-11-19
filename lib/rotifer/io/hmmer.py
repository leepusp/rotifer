__doc__ = """
Parser for HMMER 3 reports (hmmscan, hmmsearch, phmmer, jackhmmer).

Each domain hit becomes one row of a domain table with the columns
listed in ``COLUMNS``:

========= ==============================================================
sequence  sequence identifier
model     model (profile) identifier
start     first aligned residue in the sequence
end       last aligned residue in the sequence
evalue    domain's independent E-value (i-Evalue)
cov       fraction of the query covered by the alignment (the sequence
          for hmmscan and jackhmmer, the model for hmmsearch)
qstart    first aligned model state
qend      last aligned model state
iteration jackhmmer round
score     domain score
talilen   length of the aligned region in the target
qlen      query length
estart    envelope start in the sequence
eend      envelope end in the sequence
program   HMMER program name
========= ==============================================================

``parse_hmmer`` streams rows as strings, exactly as written by the
``hmmer2table`` command; ``read_hmmer`` returns a pandas DataFrame.
"""

import re
import rotifer
from rotifer.io import strip_suffix

logger = rotifer.logging.getLogger(__name__)

COLUMNS = ['sequence', 'model', 'start', 'end', 'evalue', 'cov', 'qstart', 'qend',
           'iteration', 'score', 'talilen', 'qlen', 'estart', 'eend', 'program']

DTYPES = {'start': int, 'end': int, 'evalue': float, 'cov': float, 'qstart': int,
          'qend': int, 'iteration': int, 'score': float, 'talilen': int, 'qlen': int,
          'estart': int, 'eend': int}

CLEANERS = {
    'version':   lambda name: re.sub(r'\.\d+$', '', name),
    'profiledb': lambda name: re.sub(r'(_[^_]+)?\.\d+$', '', name),
}
"""Pre-defined identifier cleaners, selected by name in ``clean``."""

_TARGETS = {'model': 'model', 'query': 'model', 'sequence': 'sequence', 'hit': 'sequence'}

_PROGRAM = re.compile(r'^#\s(\w+)\s+::')
_ROUND = re.compile(r'@@\s+Round:\s+(\d+)')
_QUERY = re.compile(r'^Query:\s+(\S+)\s+.*\[\S=(\d+)')
_HIT = re.compile(r'^>> (\S+)')
_DOMAIN = re.compile(r'^\s+\d+\s+[?!]\s+')

def _cleaners(clean):
    """
    Resolve ``clean`` into one function per target.

    Parameters
    ----------
    clean : dict or None
        Maps ``model`` (alias ``query``) or ``sequence`` (alias ``hit``) to
        the name of a pre-defined cleaner (see :data:`CLEANERS`) or to a
        function taking and returning a string. None means no cleaning.

    Returns
    -------
    dict
        ``{'model' or 'sequence': function}``.

    Raises
    ------
    ValueError
        If a target or a cleaner name is unknown, or if one target is given
        twice (``model`` and ``query``, or ``sequence`` and ``hit``).
    """
    functions = {}
    for target, cleaner in (clean or {}).items():
        if target not in _TARGETS:
            raise ValueError(f"unknown clean target {target!r}: use model or sequence")
        if not callable(cleaner):
            if cleaner not in CLEANERS:
                raise ValueError(f"unknown cleaner {cleaner!r}: use one of {', '.join(CLEANERS)} or a function")
            cleaner = CLEANERS[cleaner]
        if _TARGETS[target] in functions:
            raise ValueError(f"two cleaners for {_TARGETS[target]!r}: use either model or query, "
                             "and either sequence or hit")
        functions[_TARGETS[target]] = cleaner
    return functions

def parse_hmmer(handle, evalue=None, clean=None, force_basename=None, filename=None):
    """
    Parse a HMMER 3 report and yield one row per domain.

    Rows are lists of strings in :data:`COLUMNS` order, written as
    ``hmmer2table`` writes them. Domains are read from the "Domain
    annotation for each ..." sections. For hmmscan the query is the
    sequence and the hits are models; for hmmsearch, phmmer and jackhmmer
    it is the reverse, and ``sequence`` and ``model`` are swapped so that
    ``sequence`` is always the protein. ``cov`` is the fraction of the query
    covered (the sequence for hmmscan, the model otherwise), with two
    decimals. ``talilen`` is ``end - start``, one less than the aligned
    length, as the Perl tool wrote it.

    Parameters
    ----------
    handle : iterable of str
        Lines of a HMMER report (an open file, for example). Several
        concatenated reports are read in turn.
    evalue : float, optional
        Skip domains whose independent E-value is above this value (the
        cutoff is inclusive). None keeps every domain.
    clean : dict, optional
        Identifier cleaners: ``model`` (alias ``query``) or ``sequence``
        (alias ``hit``) mapped to the name of a pre-defined cleaner (see
        :data:`CLEANERS`) or to a function taking and returning a string.
    force_basename : list of str, optional
        Replace the model name with the report's file name without the
        first of these suffixes (regular expressions) that matches.
    filename : str, optional
        Name of the report; required by ``force_basename``.

    Yields
    ------
    list of str
        One domain, in :data:`COLUMNS` order.

    Raises
    ------
    ValueError
        If ``force_basename`` is given without ``filename``, or ``clean``
        is invalid (see :func:`_cleaners`).

    Notes
    -----
    A domain line without 16 fields is still read, and a warning is logged.

    Examples
    --------
    >>> from rotifer.io.hmmer import parse_hmmer
    >>> with open("test/data/golden/inputs/hmmscan.out") as handle:
    ...     next(parse_hmmer(handle))
    ['NP_414542.1', 'Leader_Thr', '1', '21', '1.2e-10', '1.00', '1', '23', '1', '41.6', '20', '21', '1', '21', 'hmmscan']
    """
    if force_basename and filename is None:
        raise ValueError("force_basename needs the report's filename")
    cleaners = _cleaners(clean)
    query = hit = qlen = program = None
    iteration = 1
    in_table = False
    for line in handle:
        line = line.rstrip('\n')
        if line == '//':
            iteration = 1
            continue
        match = _PROGRAM.match(line)
        if match:
            program = match.group(1)
            continue
        match = _ROUND.search(line)
        if match:
            iteration = int(match.group(1))
            continue
        match = _QUERY.match(line)
        if match:
            query, qlen = match.group(1), int(match.group(2))
            in_table = False
            continue
        match = _HIT.match(line)
        if match:
            hit = match.group(1)
            continue
        if 'Domain annotation for each' in line:
            in_table = True
            continue
        if not in_table or not _DOMAIN.match(line):
            continue

        row = line.split()
        if evalue is not None and float(row[5]) > float(evalue):
            continue
        if len(row) != 16:
            logger.warning(f"Unexpected number of columns at {line}")

        # hmmscan searches a sequence against models; the others, the reverse
        if program == 'hmmscan':
            coverage = (int(row[10]) - int(row[9]) + 1) / qlen
            ids = {'sequence': query, 'model': hit}
        else:
            coverage = (int(row[7]) - int(row[6]) + 1) / qlen
            ids = {'sequence': hit, 'model': query}
        for target, cleaner in cleaners.items():
            ids[target] = cleaner(ids[target])
        if force_basename:
            ids['model'] = strip_suffix(filename, force_basename)

        yield [ids['sequence'], ids['model'], row[9], row[10], row[5], '%.2f' % coverage,
               row[6], row[7], str(iteration), row[2], str(int(row[10]) - int(row[9])),
               str(qlen), row[12], row[13], program]

def read_hmmer(source, evalue=None, clean=None, force_basename=None):
    """
    Load HMMER 3 reports into a domain table.

    Parameters
    ----------
    source : str, os.PathLike, list of them, or file object
        Report file name(s) (``-`` is standard input; ``.gz`` and ``.bz2``
        are decompressed) or an open report. A single string is one file.
    evalue, clean, force_basename
        See :func:`parse_hmmer`. With ``force_basename`` the model name comes
        from each file's name.

    Returns
    -------
    pandas.DataFrame
        One row per domain, in report order, with :data:`COLUMNS`:
        ``sequence``, ``model`` and ``program`` are str; ``evalue``, ``cov``
        and ``score`` float64; the others int64. The frame is empty, with
        these columns and dtypes, when no domain passes.

    Raises
    ------
    FileNotFoundError
        If a report does not exist.
    ValueError
        See :func:`parse_hmmer`.

    Examples
    --------
    >>> from rotifer.io.hmmer import read_hmmer
    >>> df = read_hmmer("test/data/golden/inputs/hmmscan.out", clean={"model": "version"})
    >>> df.iloc[:2, :6]
          sequence          model  start  end        evalue   cov
    0  NP_414542.1     Leader_Thr      1   21  1.200000e-10  1.00
    1  NP_414543.1  Homoserine_dh    296  353  7.300000e+03  0.07
    """
    import pandas as pd
    from rotifer.io import each_source
    rows = []
    for name, handle in each_source(source):
        rows.extend(parse_hmmer(handle, evalue=evalue, clean=clean,
                                force_basename=force_basename, filename=name))
    return pd.DataFrame(rows, columns=COLUMNS).astype(DTYPES)
