__doc__ = """
Parser for BLAST+ pairwise text reports (blastp, psiblast, rpsblast...).

Each high-scoring pair (HSP) becomes one row of a table with the
columns in ``COLUMNS``:

============== ==========================================================
hit_name       hit (subject) identifier
query_name     query identifier
hit_start      first aligned residue in the hit
hit_end        last aligned residue in the hit
evalue         HSP E-value
query_coverage fraction of the query covered by the HSP
query_start    first aligned residue in the query
query_end      last aligned residue in the query
iteration      search round (always 1, see below)
bits           HSP bit score
hit_length     hit length
============== ==========================================================

By default overlapping HSPs of the same query and hit are merged:
when two HSPs overlap by more than ``overlap`` (a fraction of both, or a
number of residues) in both sequences, only the one with the lowest
E-value is kept; if the HSP kept first was already below the inclusion
``threshold``, its coordinates are widened instead.

With ``swap=True`` query and hit columns are exchanged, which turns
RPS-BLAST reports (protein queries, domain hits) into domain tables for
``rotifer.interval.architecture``.

This reproduces the former Perl ``blast2table`` (BioPerl based),
including its quirks: each PSI-BLAST round is a separate result and is
reported as iteration 1; rows are written after the whole report is read.

``iter_blast`` streams rows as strings, exactly as written by the
``blast2table`` command; ``read_blast`` returns a pandas DataFrame.
"""

import re
import rotifer
from rotifer.io import strip_suffix, written_as_fraction

logger = rotifer.logging.getLogger(__name__)

COLUMNS = ['hit_name', 'query_name', 'hit_start', 'hit_end', 'evalue', 'query_coverage',
           'query_start', 'query_end', 'iteration', 'bits', 'hit_length']

AVAILABLE_COLUMNS = COLUMNS + ['query_length', 'hit_coverage']
"""Columns that may be selected with ``columns``."""

DTYPES = {'hit_start': int, 'hit_end': int, 'evalue': float, 'query_coverage': float,
          'query_start': int, 'query_end': int, 'iteration': int, 'bits': float,
          'hit_length': int, 'query_length': int, 'hit_coverage': float}

_SWAPPED = ['name', 'start', 'end', 'coverage', 'length']

# Report parsing

_QUERY = re.compile(r'^Query=\s*(.*)$')
_LENGTH = re.compile(r'^Length\s*=\s*(\d+)')
_HIT = re.compile(r'^>\s*(\S+)\s*(.*)$')
_SCORE = re.compile(r'^\s*Score\s*=\s*([\d.e+-]+)\s+bits.*Expect(?:\(\d+\))?\s*=\s*([^,\s]+)')
_ALIGNMENT = re.compile(r'^(Query|Sbjct)\s+(\d+)\s+\S+\s+(\d+)\s*$')

def parse_blast_text(handle):
    """
    Read a BLAST+ pairwise text report into nested dicts.

    Parameters
    ----------
    handle : iterable of str
        Lines of a report from blastp, psiblast, rpsblast or another BLAST+
        program, in the default pairwise format (``-outfmt 0``).

    Yields
    ------
    dict
        One per result: one per query, and one per round for PSI-BLAST.
        Keys: ``query_name``, ``query_description``, ``query_length`` (int
        or None) and ``hits``, a list of dicts with ``name``,
        ``description``, ``length`` and ``hsps``. HSPs are dicts with
        ``bits`` and ``evalue`` (strings, as printed) and ``query_start``,
        ``query_end``, ``hit_start``, ``hit_end`` (ints, 1-based, start
        at most end on both strands).

    Examples
    --------
    >>> from rotifer.io.blast import parse_blast_text
    >>> with open("test/data/golden/inputs/psiblast.out") as handle:
    ...     result = next(parse_blast_text(handle))
    >>> result["query_name"], result["query_length"], result["hits"][0]["hsps"][0]["hit_end"]
    ('NP_414555.1', 638, 638)
    """
    result = hit = hsp = None
    reading = None  # 'query' or 'hit' while reading a multi-line description
    for line in handle:
        line = line.rstrip('\n')
        match = _QUERY.match(line)
        if match:
            if result is not None:
                yield result
            words = match.group(1).split(None, 1)
            result = {'query_name': words[0] if words else '', 'query_length': None,
                      'query_description': words[1] if len(words) > 1 else '', 'hits': []}
            hit = hsp = None
            reading = 'query'
            continue
        if result is None:
            continue
        match = _LENGTH.match(line)
        if match:
            if reading == 'query' and result['query_length'] is None:
                result['query_length'] = int(match.group(1))
            elif reading == 'hit':
                hit['length'] = int(match.group(1))
            reading = None
            continue
        if reading and line.strip():
            target = result if reading == 'query' else hit
            key = 'query_description' if reading == 'query' else 'description'
            target[key] = (target[key] + ' ' + line.strip()).strip()
            continue
        match = _HIT.match(line)
        if match:
            hit = {'name': match.group(1), 'description': match.group(2).strip(), 'length': None, 'hsps': []}
            result['hits'].append(hit)
            hsp = None
            reading = 'hit'
            continue
        if hit is None:
            continue
        match = _SCORE.match(line)
        if match:
            hsp = {'bits': match.group(1), 'evalue': match.group(2), 'query_start': None,
                   'query_end': None, 'hit_start': None, 'hit_end': None}
            hit['hsps'].append(hsp)
            continue
        match = _ALIGNMENT.match(line)
        if match and hsp is not None:
            prefix = 'query' if match.group(1) == 'Query' else 'hit'
            first, last = int(match.group(2)), int(match.group(3))
            start, end = hsp[f'{prefix}_start'], hsp[f'{prefix}_end']
            hsp[f'{prefix}_start'] = min(first, last) if start is None else min(start, first, last)
            hsp[f'{prefix}_end'] = max(first, last) if end is None else max(end, first, last)
    if result is not None:
        yield result

# Identifier cleaners: (name, 'query' or 'hit', hit dict) -> new name

def _genbank2gi(name, kind=None, hit=None):
    """Return the GI number of a ``gi|...`` identifier; other names as they are."""
    if not re.match(r'^gi\|[^|]+\|\w+\|+[^|]+', name):
        return name
    return re.split(r'\|+', name)[1]

def _genbank2acc(name, kind=None, hit=None):
    """Return the unversioned accession of a ``gi|...|ref|ACC.1|`` identifier; other names as they are."""
    if not re.match(r'^gi\|[^|]+\|\w+\|+[^|]+', name):
        return name
    return re.sub(r'\.\d+$', '', re.split(r'\|+', name)[3])

def _profiledb(name, kind=None, hit=None):
    """Drop a trailing version (``.1``) from hit names; query names as they are."""
    return re.sub(r'\.\d+$', '', name) if kind == 'hit' else name

def _cdd(name, kind=None, hit=None):
    """Take hit names from the second field of CDD descriptions, and GI numbers for queries."""
    if kind == 'query':
        return _genbank2gi(name)
    fields = hit['description'].split(', ')
    return fields[1] if len(fields) > 1 else name

CLEANERS = {'cdd': _cdd, 'genbank2gi': _genbank2gi, 'genbank2acc': _genbank2acc, 'profiledb': _profiledb}
"""Pre-defined identifier cleaners, selected by name in ``clean``."""

def _cleaners(clean):
    """
    Resolve ``clean`` into a list of cleaner functions.

    Parameters
    ----------
    clean : list or None
        Names from :data:`CLEANERS`, functions ``f(name, kind, hit)``, or
        ``none``/``0`` to turn cleaning off (an empty list is returned).

    Returns
    -------
    list of callable

    Raises
    ------
    ValueError
        If a name is not in :data:`CLEANERS`.
    """
    functions = []
    for cleaner in clean or []:
        if callable(cleaner):
            functions.append(cleaner)
        elif str(cleaner).lower() in ('none', '0'):
            return []
        elif cleaner in CLEANERS:
            functions.append(CLEANERS[cleaner])
        else:
            raise ValueError(f"unknown cleaner {cleaner!r}: use one of {', '.join(CLEANERS)}, none or a function")
    return functions

# Table building

def _format_evalue(text):
    """Write an E-value as blast2table does: ``%.1f`` above 1, ``0`` for 0, ``%.3e`` otherwise."""
    text = text.rstrip(',')
    if text.startswith('e'):
        text = '1' + text
    value = float(text)
    return '%.1f' % value if value > 1 else '0' if value == 0 else '%.3e' % value

def _format_bits(text):
    """Write a bit score as blast2table does: no decimals above 99, one otherwise."""
    value = float(text)
    return '%.0f' % value if value > 99 else '%.1f' % value

def _swap(row):
    """Return a copy of a row with its query and hit columns exchanged."""
    row = dict(row)
    for column in _SWAPPED:
        row[f'query_{column}'], row[f'hit_{column}'] = row.get(f'hit_{column}'), row.get(f'query_{column}')
    return row

def iter_blast_records(handle, evalue=10, threshold=0.01, overlap='0.4', nonredundant=True,
                       clean=None, rename_query=None, force_basename=None, filename='-'):
    """
    Parse a BLAST+ text report and yield one dict per HSP.

    This is the table :func:`iter_blast` writes, before columns are chosen
    or swapped. Names are cleaned and renamed, HSPs above ``evalue`` are
    dropped, and, with ``nonredundant``, overlapping HSPs of the same query
    and hit are merged by :func:`rotifer.io.search._merge_hsps` once the
    whole report is read. PSI-BLAST rounds are separate results of the
    report, and HSPs of different rounds are merged like any others; every
    row has ``iteration`` 1, as in the former Perl tool.

    Parameters
    ----------
    handle : iterable of str
        Lines of a BLAST+ pairwise text report.
    evalue, threshold, overlap, nonredundant, clean, rename_query, force_basename, filename
        See :func:`iter_blast`.

    Yields
    ------
    dict
        Keys ``query_name``, ``query_length``, ``iteration``, ``hit_name``,
        ``hit_length``, ``evalue`` and ``bits`` (strings formatted as
        blast2table prints them), ``query_start``, ``query_end``,
        ``hit_start``, ``hit_end`` (ints) and ``query_coverage``,
        ``hit_coverage`` (strings with three decimals, or ``NA``). Merged
        regions come in the order they were first seen, and keep the
        coverages of the HSP their E-value comes from.

    Raises
    ------
    ValueError
        If ``clean`` names an unknown cleaner.

    Examples
    --------
    >>> from rotifer.io.blast import iter_blast_records
    >>> with open("test/data/golden/inputs/psiblast.out") as handle:
    ...     row = next(iter_blast_records(handle))
    >>> row["hit_name"], row["evalue"], row["bits"], row["query_coverage"]
    ('NP_414555.1', '0', '1290', '1.000')
    """
    cleaners = _cleaners(clean)
    rows, keys = [], []
    for result in parse_blast_text(handle):
        query_name, query_length = result['query_name'], result['query_length']
        if re.match(r'^Query_\d+', query_name) and result['query_description']:
            query_name = result['query_description'].split()[0]
        if force_basename:
            query_name = strip_suffix(filename, force_basename)
        elif rename_query is not None:
            query_name = rename_query(query_name, filename) if callable(rename_query) else rename_query
        elif not query_name:
            logger.warning(f"Missing query name in report {filename}! Using filename...")
            query_name = filename
        match = re.match(r'^[^/]+/(\d+)-(\d+)$', query_name)
        region_start = int(match.group(1)) if match else None
        for cleaner in cleaners:
            query_name = cleaner(query_name, 'query', None)
        for hit in result['hits']:
            if not hit['hsps']:
                continue
            hit_name = hit['name']
            for cleaner in cleaners:
                hit_name = cleaner(hit_name, 'hit', hit)
            if re.match(r'^gnl\|BL_ORD_ID\|\d+$', hit_name) and hit['description']:
                hit_name = hit['description'].split()[0]
            for hsp in hit['hsps']:
                row = {'query_name': query_name, 'query_length': query_length, 'iteration': 1,
                       'hit_name': hit_name, 'hit_length': hit['length']}
                row['evalue'] = _format_evalue(hsp['evalue'])
                if evalue is not None and float(row['evalue']) > float(evalue):
                    continue
                row['bits'] = _format_bits(hsp['bits'])
                row['query_start'], row['query_end'] = hsp['query_start'], hsp['query_end']
                if region_start is not None:
                    row['query_start'] += region_start - 1
                    row['query_end'] += region_start - 1
                row['hit_start'], row['hit_end'] = hsp['hit_start'], hsp['hit_end']
                row['query_coverage'] = '%.3f' % ((row['query_end'] - row['query_start'] + 1) / query_length) \
                    if query_length else 'NA'
                row['hit_coverage'] = '%.3f' % ((row['hit_end'] - row['hit_start'] + 1) / hit['length']) \
                    if hit['length'] else 'NA'

                if not nonredundant:
                    yield row
                    continue
                rows.append(row)
                keys.append((query_name, hit_name))
    if not nonredundant or not rows:
        return
    from rotifer.io.search import _merge_hsps
    groups = {}
    for position, key in enumerate(keys):
        groups.setdefault(key, []).append(position)
    merged = _merge_hsps(
        groups.values(),
        [row['hit_start'] for row in rows], [row['hit_end'] for row in rows],
        [row['query_start'] for row in rows], [row['query_end'] for row in rows],
        [float(row['evalue']) for row in rows],
        float(threshold), float(overlap), written_as_fraction(overlap))
    for _, source, hit_start, hit_end, query_start, query_end in merged:
        yield {**rows[source], 'hit_start': hit_start, 'hit_end': hit_end,
               'query_start': query_start, 'query_end': query_end}

def iter_blast(handle, columns=COLUMNS, swap=False, **kwargs):
    """
    Parse a BLAST+ text report and yield one row of strings per HSP.

    The rows are what the ``blast2table`` command writes.

    Parameters
    ----------
    handle : iterable of str
        Lines of a BLAST+ pairwise text report.
    columns : list of str, default ``COLUMNS``
        Output columns, from :data:`AVAILABLE_COLUMNS`.
    swap : bool, default False
        Exchange the query and hit columns (names, starts, ends, coverages
        and lengths). RPS-BLAST reports (proteins searched against domain
        models) then give a domain table with the protein in ``hit_name``.
    evalue : float, default 10
        Drop HSPs whose E-value, rounded as printed, is above this value.
        None keeps them all.
    threshold : float, default 0.01
        Inclusion threshold of the merge: a kept region with a lower
        E-value is widened by overlapping HSPs instead of replaced.
    overlap : str or number, default '0.4'
        How much two HSPs must share, on both the query and the hit, to be
        merged: a fraction of each HSP's length if written with a decimal
        point (``'0.4'``, ``0.4``, ``'1.0'``), else a number of residues
        (``'1'``, ``10``); see :func:`rotifer.io.written_as_fraction`.
    nonredundant : bool, default True
        Merge redundant HSPs. False yields every HSP.
    clean : list, optional
        Identifier cleaners, applied in order: names from :data:`CLEANERS`
        (``cdd`` takes domain names from RPS-BLAST hit descriptions;
        ``none`` or ``0`` turns cleaning off) or functions
        ``f(name, kind, hit)`` with ``kind`` ``'query'`` or ``'hit'``.
    rename_query : str or callable, optional
        New query name, or a function ``f(name, filename)`` returning it.
    force_basename : list of str, optional
        Use the report's file name, without the first of these suffixes
        (regular expressions) that matches, as the query name.
    filename : str, default '-'
        Name of the report, for ``force_basename``, ``rename_query`` and
        warnings.

    Yields
    ------
    list of str
        One HSP, in ``columns`` order; missing values are empty strings.

    Raises
    ------
    ValueError
        If ``clean`` names an unknown cleaner.

    Examples
    --------
    >>> from rotifer.io.blast import iter_blast
    >>> with open("test/data/golden/inputs/psiblast.out") as handle:
    ...     next(iter_blast(handle))
    ['NP_414555.1', 'NP_414555.1', '1', '638', '0', '1.000', '1', '638', '1', '1290', '638']
    """
    for row in iter_blast_records(handle, **kwargs):
        if swap:
            row = _swap(row)
        yield ['' if row.get(column) is None else str(row[column]) for column in columns]

def read_blast(source, columns=COLUMNS, swap=False, **kwargs):
    """
    Load BLAST+ text reports into a table, one row per (merged) HSP.

    Parameters
    ----------
    source : str, os.PathLike, list of them, or file object
        Report file name(s) (``-`` is standard input; ``.gz`` and ``.bz2``
        are decompressed) or an open report. A single string is one file.
    columns : list of str, default ``COLUMNS``
        Output columns, from :data:`AVAILABLE_COLUMNS`.
    swap : bool, default False
        Exchange query and hit columns (see :func:`iter_blast`).
    **kwargs
        ``evalue``, ``threshold``, ``overlap``, ``nonredundant``, ``clean``,
        ``rename_query`` and ``force_basename``, as in :func:`iter_blast`.

    Returns
    -------
    pandas.DataFrame
        The rows :func:`iter_blast` yields for every report, in order, with
        ``columns``. Columns in :data:`DTYPES` are numeric (``NA`` becomes
        NaN): coordinates, lengths and ``iteration`` int64, ``evalue``,
        ``bits`` and the coverages float64, also when the table is empty.

    Raises
    ------
    FileNotFoundError
        If a report does not exist.
    ValueError
        If ``clean`` names an unknown cleaner.

    See Also
    --------
    rotifer.io.search.read_blast_tabular : the same merge on BLAST's
        tabular output, into the standard search table.

    Examples
    --------
    >>> from rotifer.io.blast import read_blast
    >>> domains = read_blast("test/data/golden/inputs/rpsblast.out", swap=True, clean=["cdd"])
    >>> domains.loc[:1, ["hit_name", "query_name", "hit_start", "hit_end", "evalue"]]
          hit_name                      query_name  hit_start  hit_end  evalue
    0  NP_414543.1                            thrA          1      820     0.0
    1  NP_414543.1  Aspartokinase_Aspartate_kinase          3      461     0.0
    """
    import pandas as pd
    from rotifer.io import each_source
    rows = []
    for name, handle in each_source(source):
        rows.extend(iter_blast(handle, columns=columns, swap=swap, filename=name, **kwargs))
    df = pd.DataFrame(rows, columns=columns)
    for column in columns:
        if column in DTYPES:
            df[column] = pd.to_numeric(df[column], errors='coerce')
            if DTYPES[column] is float:
                # an empty column would otherwise come back as int64
                df[column] = df[column].astype('float64')
    return df
