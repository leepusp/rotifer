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
    Parse a BLAST+ pairwise text report.

    Yields one dict per result (one per query, and one per round for
    PSI-BLAST) with keys ``query_name``, ``query_description``,
    ``query_length`` and ``hits``. Hits are dicts with ``name``,
    ``description``, ``length`` and ``hsps``; HSPs are dicts with
    ``bits``, ``evalue`` (strings, as in the report), ``query_start``,
    ``query_end``, ``hit_start`` and ``hit_end`` (integers).
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
    if not re.match(r'^gi\|[^|]+\|\w+\|+[^|]+', name):
        return name
    return re.split(r'\|+', name)[1]

def _genbank2acc(name, kind=None, hit=None):
    if not re.match(r'^gi\|[^|]+\|\w+\|+[^|]+', name):
        return name
    return re.sub(r'\.\d+$', '', re.split(r'\|+', name)[3])

def _profiledb(name, kind=None, hit=None):
    return re.sub(r'\.\d+$', '', name) if kind == 'hit' else name

def _cdd(name, kind=None, hit=None):
    if kind == 'query':
        return _genbank2gi(name)
    fields = hit['description'].split(', ')
    return fields[1] if len(fields) > 1 else name

CLEANERS = {'cdd': _cdd, 'genbank2gi': _genbank2gi, 'genbank2acc': _genbank2acc, 'profiledb': _profiledb}
"""Pre-defined identifier cleaners, selected by name in ``clean``."""

def _cleaners(clean):
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
    text = text.rstrip(',')
    if text.startswith('e'):
        text = '1' + text
    value = float(text)
    return '%.1f' % value if value > 1 else '0' if value == 0 else '%.3e' % value

def _format_bits(text):
    value = float(text)
    return '%.0f' % value if value > 99 else '%.1f' % value

def _same_hsp(a, b, overlap):
    """
    True if two HSPs overlap in both the query and the hit.
    """
    fraction = written_as_fraction(overlap)
    overlaps = 0
    for kind in ('hit', 'query'):
        length = min(a[f'{kind}_end'], b[f'{kind}_end']) - max(a[f'{kind}_start'], b[f'{kind}_start']) + 1
        if fraction:
            limit = float(overlap)
            if length > limit * (a[f'{kind}_end'] - a[f'{kind}_start'] + 1) and \
               length > limit * (b[f'{kind}_end'] - b[f'{kind}_start'] + 1):
                overlaps += 1
        elif length > float(overlap):
            overlaps += 1
    return overlaps == 2

def _swap(row):
    row = dict(row)
    for column in _SWAPPED:
        row[f'query_{column}'], row[f'hit_{column}'] = row.get(f'hit_{column}'), row.get(f'query_{column}')
    return row

def iter_blast_records(handle, evalue=10, threshold=0.01, overlap='0.4', nonredundant=True,
                       clean=None, rename_query=None, force_basename=None, filename='-'):
    """
    Parse a BLAST+ text report and yield one dict per HSP (see
    ``iter_blast`` for the parameters), before any column swap.
    """
    cleaners = _cleaners(clean)
    parsed, index = [], {}
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
                key = (query_name, hit_name)
                for i in index.get(key, []):
                    if not _same_hsp(parsed[i], row, overlap):
                        continue
                    if float(parsed[i]['evalue']) < float(threshold):
                        # Already significant: keep the E-value, widen the coordinates
                        for kind in ('hit', 'query'):
                            if row[f'{kind}_end'] - row[f'{kind}_start'] > \
                               parsed[i][f'{kind}_end'] - parsed[i][f'{kind}_start']:
                                parsed[i][f'{kind}_start'] = row[f'{kind}_start']
                                parsed[i][f'{kind}_end'] = row[f'{kind}_end']
                    elif float(row['evalue']) < float(parsed[i]['evalue']):
                        parsed[i] = row
                    break
                else:
                    parsed.append(row)
                    index.setdefault(key, []).append(len(parsed) - 1)
    yield from parsed

def iter_blast(handle, columns=COLUMNS, swap=False, **kwargs):
    """
    Parse a BLAST+ text report and yield one row (a list of strings) per
    HSP, with the given ``columns``.

    Parameters
    ----------
    handle : iterable of str
        Lines of a BLAST+ pairwise text report.
    columns : list of str
        Output columns (see ``AVAILABLE_COLUMNS``).
    swap : bool, default False
        Exchange query and hit columns.
    evalue : float, default 10
        Maximum HSP E-value.
    threshold : float, default 0.01
        Inclusion threshold used by PSI-BLAST (see the module docstring).
    overlap : str or number, default '0.4'
        Overlap above which HSPs are redundant: a fraction if written with
        a decimal point, else a number of residues.
    nonredundant : bool, default True
        Merge redundant HSPs.
    clean : list
        Identifier cleaners, applied in order: names of pre-defined
        cleaners (see ``CLEANERS``; ``cdd`` takes domain names from
        RPS-BLAST hit descriptions) or functions ``f(name, kind, hit)``
        where ``kind`` is ``'query'`` or ``'hit'``.
    rename_query : str or function, optional
        New query name, or function ``f(name, filename)``.
    force_basename : list of str, optional
        Use the report's file name, without these suffixes, as query name.
    filename : str
        Name of the report.
    """
    for row in iter_blast_records(handle, **kwargs):
        if swap:
            row = _swap(row)
        yield ['' if row.get(column) is None else str(row[column]) for column in columns]

def read_blast(source, columns=COLUMNS, swap=False, **kwargs):
    """
    Load BLAST+ text reports into a table.

    Parameters
    ----------
    source : str, list of str or file object
        Report file name(s) or an open report.
    columns, swap, **kwargs
        See ``iter_blast``.

    Returns
    -------
    pandas.DataFrame

    Examples
    --------
    >>> from rotifer.io.blast import read_blast
    >>> domains = read_blast("seqs.rpsblast.out", swap=True, clean=['cdd'])
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
    return df
