__doc__ = """
Parser for Phobius predictions (default, long output format).

Signal peptides and transmembrane segments become rows of a domain
table with the columns in ``COLUMNS``: sequence accession, region type
(``SP`` or ``TM``), start, end and an arbitrary E-value, so they can be
merged with domain tables from HMMER or BLAST (see
``rotifer.interval.architecture``).

``parse_phobius`` streams rows as strings, exactly as written by the
``phobius2table`` command; ``read_phobius`` returns a pandas DataFrame.
"""

import re
from rotifer.io import split_fields

COLUMNS = ['sequence', 'domain', 'start', 'end', 'evalue']
DTYPES = {'start': int, 'end': int, 'evalue': float}

_ID = re.compile(r'^ID\s+(\S+)')
_FEATURE = re.compile(r'^FT\s+(TRANSMEM|SIGNAL)')
_REGION = re.compile(r'[:/|](\d+)[-.]{1,2}(\d+)$')

def accession(identifier):
    """
    Extract the accession from a sequence identifier.

    Parameters
    ----------
    identifier : str
        A sequence identifier as Phobius prints it.

    Returns
    -------
    str
        The identifier without a region suffix (``/10-200``, ``:10..200``,
        ``|10-200``). A plain identifier is returned as it is; for a
        pipe-separated one (``gi|123|ref|NP_1.1|``, ``sp|P12345|NAME``) the
        second field is returned. An empty string when the first or second
        field is empty.

    Examples
    --------
    >>> from rotifer.io.phobius import accession
    >>> accession("sp|P12345|NAME"), accession("NP_414542.1/10-200")
    ('P12345', 'NP_414542.1')
    """
    identifier = _REGION.sub('', identifier)
    fields = [field.strip() for field in split_fields(identifier, r'\|')]
    if len(fields) == 1:
        return fields[0]
    if len(fields) > 1 and fields[0] and fields[1]:
        return fields[1]
    return ''

def parse_phobius(handle, evalue='0'):
    """
    Parse a Phobius report and yield one row per signal peptide or TM segment.

    Parameters
    ----------
    handle : iterable of str
        Lines of a Phobius report in the default (long) format. Only ``ID``
        lines and ``FT`` lines for ``SIGNAL`` and ``TRANSMEM`` are read.
    evalue : str or float, default '0'
        Text written in the ``evalue`` column of every row. Phobius gives no
        E-value; this one lets the rows pass E-value filters when they are
        merged with HMMER or BLAST domains.

    Yields
    ------
    list of str
        ``[sequence, domain, start, end, evalue]``, with ``domain`` ``SP``
        or ``TM`` and the sequence's accession from :func:`accession`.

    Examples
    --------
    >>> from rotifer.io.phobius import parse_phobius
    >>> with open("test/data/golden/inputs/phobius.out") as handle:
    ...     next(parse_phobius(handle))
    ['NP_414546.1', 'SP', '1', '23', '0']
    """
    evalue = str(evalue)
    sequence = ''
    for line in handle:
        match = _ID.match(line)
        if match:
            sequence = accession(match.group(1))
            continue
        if _FEATURE.match(line):
            fields = line.split()
            kind = 'SP' if 'SIGNAL' in fields[1] else 'TM'
            yield [sequence, kind, fields[2], fields[3], evalue]

def read_phobius(source, evalue=0):
    """
    Load Phobius reports into a domain table.

    Parameters
    ----------
    source : str, os.PathLike, list of them, or file object
        Report file name(s) (``-`` is standard input; ``.gz`` and ``.bz2``
        are decompressed) or an open report. A single string is one file.
    evalue : float, default 0
        Value of the ``evalue`` column (see :func:`parse_phobius`).

    Returns
    -------
    pandas.DataFrame
        One row per signal peptide or TM segment, in report order, with
        :data:`COLUMNS`: ``sequence`` and ``domain`` (str), ``start`` and
        ``end`` (int64), ``evalue`` (float64). Empty, with these columns,
        when the reports have no such segment.

    Raises
    ------
    FileNotFoundError
        If a report does not exist.

    Examples
    --------
    >>> from rotifer.io.phobius import read_phobius
    >>> read_phobius("test/data/golden/inputs/phobius.out", evalue=0.0101).head(2)
          sequence domain  start  end  evalue
    0  NP_414546.1     SP      1   23  0.0101
    1  NP_414548.1     TM     12   30  0.0101
    """
    import pandas as pd
    from rotifer.io import each_source
    rows = []
    for name, handle in each_source(source):
        rows.extend(parse_phobius(handle, evalue=evalue))
    return pd.DataFrame(rows, columns=COLUMNS).astype(DTYPES)
