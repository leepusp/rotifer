# HH-suite parser

__author__ = "Robson Francisco de Souza, Rodolfo Alvarenga Ribeiro, Gianlucca G. Nicastro"
__version__ = 0.002

# Libraries
import os
import sys
import pandas as pd
from glob import glob
from collections import namedtuple

# Internal data
default_columns = ['query','hit','qstart','qend','evalue','qcov','hstart','hend','aligned_cols','score',
                'qalilen','qlen','hlen','probability','identity','similarity','sum_probs','hcov','pvalue']

# Column name to hhr_alignment attribute or function of the record
_COL2ATTR = {
    'query':   'query_id',
    'hit':     'template_id',
    'qstart':  lambda x: x.start[0],
    'qend':    lambda x: x.end[0],
    'qcov':    lambda x: (x.end[0] - x.start[0] + 1) / x.query_length,
    'hstart':  lambda x: x.start[1],
    'hend':    lambda x: x.end[1],
    'qlen':    'query_length',
    'qalilen': lambda x: x.end[0] - x.start[0] + 1,
    'hlen':    'template_length',
    'hcov':    lambda x: (x.end[1] - x.start[1] + 1) / x.template_length if x.template_length > 0 else 0,
}
_RECORD_FIELDS = ['query_id', 'query_length', 'query_neff', 'template_id', 'template_length',
                  'template_info', 'template_neff', 'query_ali', 'template_ali', 'start', 'end']
_ALLOWED = list(default_columns) + _RECORD_FIELDS
# Documented dtypes; the other fields keep what parse_result gives
_COLUMN_TYPES = {}
_COLUMN_TYPES.update({c: 'str' for c in ('query', 'hit')})
_COLUMN_TYPES.update({c: 'int64' for c in ('qstart', 'qend', 'hstart', 'hend', 'qalilen', 'qlen', 'hlen', 'aligned_cols')})
_COLUMN_TYPES.update({c: 'float64' for c in ('qcov', 'hcov', 'probability', 'evalue', 'pvalue', 'score',
                                             'similarity', 'sum_probs', 'identity')})

# HH-suite output parser
def hhsuite2pandas(infile, columns=default_columns):
    """
    Load the hits of HH-suite reports into a table, one row per hit.

    Each ``.hhr`` report (written by hhblits, hhsearch or hhalign with
    ``-o``) lists the hits of one query; every hit becomes a row, in the
    order of the files and, within a file, of the report.

    Parameters
    ----------
    infile : str, os.PathLike or list of them
        One report, or a list of reports. A single string is one file name,
        not a list of characters.
    columns : list of str, default ``default_columns``
        Columns of the table, in this order. Allowed: every name in
        ``default_columns`` (``query``, ``hit``, ``qstart``, ``qend``,
        ``evalue``, ``qcov``, ``hstart``, ``hend``, ``aligned_cols``,
        ``score``, ``qalilen``, ``qlen``, ``hlen``, ``probability``,
        ``identity``, ``similarity``, ``sum_probs``, ``hcov``, ``pvalue``)
        and the fields of the records :func:`parse_result` returns
        (``query_id``, ``query_length``, ``query_neff``, ``template_id``,
        ``template_length``, ``template_info``, ``template_neff``,
        ``query_ali``, ``template_ali``, ``start``, ``end``). A name given
        twice gives one column. The list is not modified.

    Returns
    -------
    pandas.DataFrame
        One row per hit, with a default RangeIndex and these columns when
        requested:

        * ``query``, ``hit`` (str): query and hit (template) identifiers.
        * ``qstart``, ``qend``, ``hstart``, ``hend`` (int64): first and last
          aligned columns of the query and of the hit, 1-based, inclusive.
        * ``qalilen`` (int64): ``qend - qstart + 1``.
        * ``qlen``, ``hlen`` (int64): query and hit lengths (match columns).
        * ``qcov`` (float64): ``qalilen / qlen``.
        * ``hcov`` (float64): ``(hend - hstart + 1) / hlen``, 0 when ``hlen``
          is 0.
        * ``probability`` (float64): HH-suite's probability, in percent
          (99.9, not 0.999).
        * ``evalue``, ``pvalue``, ``score``, ``similarity``, ``sum_probs``
          (float64): as printed in the report.
        * ``identity`` (float64): fraction of identical residues, 0 to 1
          (the report prints a percentage).
        * ``aligned_cols`` (int64): number of aligned columns.

        Other fields keep the type :func:`parse_result` gives them
        (``query_ali`` and ``template_ali`` are str). The frame is empty,
        with the requested columns and these dtypes, when the reports have
        no hits.

    Raises
    ------
    ValueError
        If ``columns`` names a column that is not listed above. Nothing is
        read in that case.
    FileNotFoundError
        If a report does not exist.
    HHRFormatError
        If a report's alignment coordinates cannot be read.

    See Also
    --------
    read_hhr : the same table from a directory of reports.

    Notes
    -----
    Earlier versions returned ``score`` as text (``'120.5'``) and
    ``probability`` truncated to an integer. ``bin/hhsuite2table`` still
    prints them that way.

    Examples
    --------
    >>> from rotifer.io.hhsuite import hhsuite2pandas
    >>> df = hhsuite2pandas("test/io/data/sample.hhr", columns=["hit", "probability", "score", "pvalue"])
    >>> df
             hit  probability  score        pvalue
    0  PF00001.1         99.9  120.5  2.200000e-24
    1  PF00002.2         50.0   20.1  1.000000e-04
    """
    unknown = [c for c in columns if c not in _ALLOWED]
    if unknown:
        raise ValueError("Unknown column(s) %s; allowed: %s"
                         % (", ".join(map(str, unknown)), ", ".join(_ALLOWED)))
    columns = list(columns)
    data = {c: [] for c in columns}
    for f in _expand(infile):
        for x in read_result(f):
            for c in columns:
                attr = _COL2ATTR.get(c, c)
                data[c].append(attr(x) if callable(attr) else getattr(x, attr))
    dtypes = {c: _COLUMN_TYPES[c] for c in columns if c in _COLUMN_TYPES}
    if any(data.values()):
        df = pd.DataFrame(data=data)
        return df.astype(dtypes)
    df = pd.DataFrame({c: pd.Series([], dtype=dtypes.get(c, object)) for c in columns})
    return df

def _expand(source, suffix=".hhr"):
    """
    Turn a report, a list of reports or a directory into a list of files.

    Parameters
    ----------
    source : str, os.PathLike, list or tuple
        A directory, a single report, or a list of reports. A single string
        is one path, not a list of characters.
    suffix : str, default ".hhr"
        Only used for a directory: its files whose names end with
        ``suffix`` are listed.

    Returns
    -------
    list of str
        A list or tuple is returned as a list, in its order; a directory
        gives its matching files sorted by name; any other path gives a
        one-element list (it is not checked for existence).
    """
    if isinstance(source, (list, tuple)):
        return [os.fspath(f) for f in source]
    source = os.fspath(source)
    if os.path.isdir(source):
        return [os.path.join(source, f) for f in sorted(os.listdir(source)) if f.endswith(suffix)]
    return [source]

"""
Parser for hhr result files created with hhblits|hhsearch|hhalign -o <hhr_file>
"""
# __author__ = 'Markus Meier (markus.meier@mpibpc.mpg.de)'
# __version__ = '1.0'
# __license__ = "GPL-3"

class HHRFormatError(Exception):
    """
    Raised when the alignment coordinates of an HH-suite report are not numbers.

    Parameters
    ----------
    value : str
        What could not be converted; the message is prefixed with ``ERROR: ``.
    """
    def __init__(self, value):
        self.value = "ERROR: "+value

    def __str__(self):
        return repr(self.value)

def parse_result(lines):
    """
    Parse the lines of one HH-suite report into one record per hit.

    Parameters
    ----------
    lines : iterable of str
        Lines of an ``.hhr`` report, with their line endings.

    Returns
    -------
    list of collections.namedtuple
        One ``hhr_alignment`` per hit, in report order, with the fields
        ``query_id``, ``query_length``, ``query_neff``, ``template_id``,
        ``template_length``, ``template_info``, ``template_neff``,
        ``query_ali`` and ``template_ali`` (str), ``start`` and ``end``
        (tuples of query and template coordinates), ``probability``,
        ``evalue``, ``score``, ``aligned_cols``, ``identity``,
        ``similarity``, ``sum_probs`` and ``pvalue``.

    Raises
    ------
    HHRFormatError
        If an alignment coordinate cannot be converted to an integer.
    """
    results = []
    hhr_alignment = namedtuple('hhr_alignment', ['query_id', 'query_length', 'query_neff',
                                                 'template_id', 'template_length', 'template_info',
                                                 'template_neff', 'query_ali', 'template_ali',
                                                 'start', 'end', 'probability', 'evalue', 'score',
                                                 'aligned_cols', 'identity', 'similarity', 'sum_probs','pvalue'])
    query_id = None
    query_length = None
    query_neff = None
    query_seq = []
    template_id = None
    template_length = None
    template_seq = []
    template_info = None
    query_start = None
    query_end = None
    template_start = None
    template_end = None
    probability = None
    evalue = None
    pvalue = None
    pvalueDict = dict()
    score = None
    identity = None
    similarity = None
    template_neff = None
    sum_probs = None
    aligned_cols = None

    skipped_ali_tags = ["ss_dssp", "ss_pred", "Consensus"]

    is_alignment_section = False
    is_hit_table = False

    for line in lines:
        if "P-value" in line:
            is_hit_table = True
            continue
        if is_hit_table:
            if line == "\n" or len(line) == 0:
                is_hit_table = False
                continue
            tokens = line.split()
            hitID = int(tokens[0])
            pvalue = float(line[49:56].strip())
            pvalueDict[hitID] = pvalue
        if line.startswith("No") and ("P-value" not in line):
            tokens = line.split()
            hitID = int(tokens[1])

        if(line.startswith("Query")):
            query_id = line.split()[1]
        elif(line.startswith("Match_columns")):
            query_length = int(line.split()[1])
        elif(line.startswith("Neff")):
            query_neff = float(line.split()[1])
        elif(is_alignment_section and (line.startswith("No") or line.startswith("Done!"))):
            if query_start is not None:
                # At 'No N', hitID already points to the next hit, N; at
                # 'Done!' it is still the hit being closed
                pvalue = pvalueDict[hitID if line.startswith("Done!") else hitID - 1]
                result = hhr_alignment(query_id, query_length, query_neff,
                                       template_id, template_length, template_info, template_neff,
                                       "".join(query_seq), "".join(template_seq), (query_start, template_start),
                                       (query_end, template_end), probability, evalue, score,
                                       aligned_cols, identity, similarity, sum_probs, pvalue)
                results.append(result)
            template_id = None
            template_info = None
            query_seq = []
            template_seq = []

            query_start = None
            query_end = None
            template_start = None
            template_end = None
        elif(line.startswith("Probab")):
            tokens = line.split()
            probability = float(tokens[0].split("=")[1])
            evalue = float(tokens[1].split("=")[1])
            score = float(tokens[2].split("=")[1])
            aligned_cols = int(tokens[3].split("=")[1])
            identity = float(tokens[4].split("=")[1].replace("%", "")) / 100.0
            similarity = float(tokens[5].split("=")[1])
            sum_probs = float(tokens[6].split("=")[1])
            if(len(tokens) > 7):
                template_neff = float(tokens[7].split("=")[1])
            continue
        elif(line.startswith(">")):
            is_alignment_section = True
            template_id = line[1:].split()[0]
            template_info = line
            is_hit_table = False
        elif(line.startswith("Q")):
            tokens = line.split()
            if (tokens[1] != 'Consensus'):
                continue
            #if(tokens[1] in skipped_ali_tags):
            #   continue

            try:
                token_2 = tokens[2].replace("(", "").replace(")", "")
                token_2 = int(token_2)
            except ValueError as e:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of query alignment").format(tokens[2])) from e

            if query_start is None:
                query_start = token_2
            query_start = min(query_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except ValueError as e:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of query alignment").format(tokens[4])) from e

            if query_end is None:
                query_end = token_4
            query_end = max(query_end, token_4)
            query_seq.append(tokens[3])
        elif(line.startswith("T")):
            tokens = line.split()
            if (tokens[1] != 'Consensus'):
                continue
            #if(tokens[1] in skipped_ali_tags):
            #   continue
            template_seq.append(tokens[3])

            try:
                token_2 = tokens[2].replace("(", "").replace(")", "")
                token_2 = int(token_2)
            except ValueError as e:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of template alignment").format(tokens[2])) from e

            if template_start is None:
                template_start = token_2
            template_start = min(template_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except ValueError as e:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of template alignment").format(tokens[4])) from e

            if template_end is None:
                template_end = token_4
            template_end = max(template_end, token_4)

            try:
                token_5 = tokens[5].replace("(", "").replace(")", "")
                token_5 = int(token_5)
            except ValueError as e:
                raise HHRFormatError(("Converting failure of template length ({}) "
                                      "in template alignment").format(tokens[5])) from e
            template_length = token_5

    if(template_id is not None and query_start is not None):
        pvalue = pvalueDict[hitID]
        result = hhr_alignment(query_id, query_length, query_neff,
                               template_id, template_length, template_info, template_neff,
                               "".join(query_seq), "".join(template_seq), (query_start, template_start),
                               (query_end, template_end), probability, evalue, score,
                               aligned_cols, identity, similarity, sum_probs, pvalue)
        results.append(result)

    return results

def read_result(input_file):
    """
    Read one HH-suite report file and parse it.

    Parameters
    ----------
    input_file : str or os.PathLike
        Path of an ``.hhr`` report. Undecodable bytes are ignored.

    Returns
    -------
    list of collections.namedtuple
        The records :func:`parse_result` returns.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    HHRFormatError
        If an alignment coordinate cannot be read.
    """
    with open(input_file, errors='ignore') as fh:
        lines = fh.readlines()
        return parse_result(lines)

# Parse output tables
def read_hhr(indir, suffix=".hhr"):
    """
    Load HH-suite reports into one table: a report, a list or a directory.

    Parameters
    ----------
    indir : str, os.PathLike, list or tuple
        A report, a list of reports, or a directory. For a directory, every
        file whose name ends with ``suffix`` is read, in name order.
    suffix : str, default ".hhr"
        File name ending of the reports in a directory.

    Returns
    -------
    pandas.DataFrame
        The table :func:`hhsuite2pandas` returns with ``default_columns``,
        for all the reports together.

    Raises
    ------
    FileNotFoundError
        If a report does not exist.

    See Also
    --------
    parse_hhr : one table per report.

    Examples
    --------
    >>> from rotifer.io import hhsuite
    >>> hhsuite.read_hhr("test/io/data/sample.hhr")[["hit", "qstart", "qend", "probability"]]
             hit  qstart  qend  probability
    0  PF00001.1       3    47         99.9
    1  PF00002.2      10    29         50.0
    """
    return hhsuite2pandas(_expand(indir, suffix))

def parse_hhr(indir, suffix=".hhr"):
    """
    Yield one table per HH-suite report: a report, a list or a directory.

    Parameters
    ----------
    indir : str, os.PathLike, list or tuple
        A report, a list of reports, or a directory. For a directory, every
        file whose name ends with ``suffix`` is read, in name order. A single
        report gives a single table.
    suffix : str, default ".hhr"
        File name ending of the reports in a directory.

    Yields
    ------
    pandas.DataFrame
        For each report, the table :func:`hhsuite2pandas` returns with
        ``default_columns``.

    Raises
    ------
    FileNotFoundError
        If a report does not exist (raised when that report is reached).

    Examples
    --------
    >>> from rotifer.io import hhsuite
    >>> [len(table) for table in hhsuite.parse_hhr("test/io/data/sample.hhr")]
    [2]
    """
    for infile in _expand(indir, suffix):
        yield hhsuite2pandas(infile)
