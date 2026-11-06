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

    Region suffixes (``/10-200``, ``:10..200``) are removed. Plain
    identifiers are returned as they are; for pipe-separated identifiers
    (``gi|123|ref|NP_1.1|``, ``sp|P12345|NAME``) the second field is
    returned. Returns an empty string when nothing can be extracted.
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
    Parse a Phobius report and yield one row (a list of strings) per
    signal peptide or transmembrane segment, in ``COLUMNS`` order.

    Parameters
    ----------
    handle : iterable of str
        Lines of a Phobius report in the default (long) format.
    evalue : str or float, default '0'
        Value written in the E-value column of every row.
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
    source : str, list of str or file object
        Report file name(s) or an open report.
    evalue : float, default 0
        Value of the E-value column.

    Returns
    -------
    pandas.DataFrame with ``COLUMNS``

    Examples
    --------
    >>> from rotifer.io.phobius import read_phobius
    >>> tm = read_phobius("seqs.phobius.out", evalue=0.0101)
    """
    import pandas as pd
    from rotifer.io import each_source
    rows = []
    for name, handle in each_source(source):
        rows.extend(parse_phobius(handle, evalue=evalue))
    return pd.DataFrame(rows, columns=COLUMNS).astype(DTYPES)
