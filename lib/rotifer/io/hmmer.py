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
    Turn ``{target: cleaner}`` into ``{'model'|'sequence': function}``.
    """
    functions = {}
    for target, cleaner in (clean or {}).items():
        if target not in _TARGETS:
            raise ValueError(f"unknown clean target {target!r}: use model or sequence")
        if not callable(cleaner):
            if cleaner not in CLEANERS:
                raise ValueError(f"unknown cleaner {cleaner!r}: use one of {', '.join(CLEANERS)} or a function")
            cleaner = CLEANERS[cleaner]
        functions.setdefault(_TARGETS[target], cleaner)
    return functions

def parse_hmmer(handle, evalue=None, clean=None, force_basename=None, filename=None):
    """
    Parse one HMMER 3 report and yield one row (a list of strings) per
    domain, in ``COLUMNS`` order.

    Parameters
    ----------
    handle : iterable of str
        Lines of a HMMER report (an open file, for example).
    evalue : float, optional
        Skip domains whose independent E-value is above this cutoff.
    clean : dict, optional
        Identifier cleaners, mapping ``model`` (alias ``query``) or
        ``sequence`` (alias ``hit``) to the name of a pre-defined cleaner
        (see ``CLEANERS``) or to a function taking and returning a string.
    force_basename : list of str, optional
        Replace the model name by the report's file name without these
        suffixes (regular expressions).
    filename : str, optional
        Name of the report, required by ``force_basename``.
    """
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
    source : str, list of str or file object
        Report file name(s) (``.gz`` and ``.bz2`` are decompressed) or an
        open report.
    evalue, clean, force_basename
        See ``parse_hmmer``.

    Returns
    -------
    pandas.DataFrame with ``COLUMNS``

    Examples
    --------
    >>> from rotifer.io.hmmer import read_hmmer
    >>> df = read_hmmer("seqs.hmmscan.out", clean={'model': 'version'})
    """
    import pandas as pd
    from rotifer.io import each_source
    rows = []
    for name, handle in each_source(source):
        rows.extend(parse_hmmer(handle, evalue=evalue, clean=clean,
                                force_basename=force_basename, filename=name))
    return pd.DataFrame(rows, columns=COLUMNS).astype(DTYPES)
