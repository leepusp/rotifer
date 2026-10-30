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

# HH-suite output parser
def hhsuite2pandas(infile, columns=default_columns):
    '''
    Parse a list of HH-suite output files.
    '''
    import types
    if not isinstance(infile,list):
        infile = [infile]

    # Initialize hash for dataframe
    data = {}
    for c in columns: data[c] = []

    # Column name to hhr_alignment attribute or lambda function
    col2attr = { 
            'query':  'query_id',
            'hit':    'template_id',
            'qstart': lambda x: x.start[0],
            'qend':   lambda x: x.end[0],
            'qcov':   lambda x: (x.end[0] - x.start[0] + 1) / x.query_length,
            'hstart': lambda x: x.start[1],
            'hend':   lambda x: x.end[1],
            'qlen':   'query_length',
            'qalilen': lambda x: x.end[0] - x.start[0] + 1,
            'hlen':   'template_length',
            'hcov':   lambda x: (x.end[1] - x.start[1] + 1) / x.template_length if x.template_length > 0 else 0
            }

    # Transfer data from hhr_alignment's list of tuples to a dictionary of lists
    for f in infile:
        for x in read_result(f):
            for c in columns:
                if c in col2attr:
                    if type(col2attr[c]) is types.LambdaType:
                        data[c].append(col2attr[c](x))
                    elif hasattr(x,col2attr[c]):
                        data[c].append(getattr(x,col2attr[c]))
                    else:
                        sys.stderr.write("Column "+c+" is not a function or hhr_alignment attribute!\n")
                elif hasattr(x,c):
                    data[c].append(getattr(x,c))
                else:
                    sys.stderr.write("Column "+c+" is unknown\n")

    df = pd.DataFrame(data=data)
    df.probability = df.probability.astype(int)
    df.score = df['score'].map(lambda x: '%-.1f' % x)
    return df

"""
Parser for hhr result files created with hhblits|hhsearch|hhalign -o <hhr_file>
"""
# __author__ = 'Markus Meier (markus.meier@mpibpc.mpg.de)'
# __version__ = '1.0'
# __license__ = "GPL-3"

class HHRFormatError(Exception):
    def __init__(self, value):
        self.value = "ERROR: "+value

    def __str__(self):
        return repr(self.value)

def parse_result(lines):
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
                pvalue = pvalueDict[hitID-1]
                result = hhr_alignment(query_id, query_length, query_neff,
                                       template_id, template_length, template_info, template_neff,
                                       query_seq, template_seq, (query_start, template_start),
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
            except:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of query alignment").format(tokens[2]))

            if query_start is None:
                query_start = token_2
            query_start = min(query_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of query alignment").format(tokens[4]))

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
            except:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of template alignment").format(tokens[2]))

            if template_start is None:
                template_start = token_2
            template_start = min(template_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of template alignment").format(tokens[4]))

            if template_end is None:
                template_end = token_4
            template_end = max(template_end, token_4)

            try:
                token_5 = tokens[5].replace("(", "").replace(")", "")
                token_5 = int(token_5)
            except:
                raise HHRFormatError(("Converting failure of template length ({}) "
                                      "in template alignment").format(tokens[5]))
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
    with open(input_file, errors='ignore') as fh:
        lines = fh.readlines()
        return parse_result(lines)

# Parse output tables
def read_hhr(indir, suffix=".hhr"):
    '''
    Parse **all** HH-suite (``*.hhr``) output files located in a
    directory or, alternatively, a list of .hhr files

    Parameters
    ----------
    indir : string
      A HH-suite results file, a list of files or a directory with results
    suffix: string, default .hhr
      The suffix, a.k.a extension, for the HH-suite results in indir

    Returns
    -------
      A Pandas dataframe

    Examples
    --------
    Parsing all ``*.hhr`` files at once

    >>> from rotifer.io import hhsuite
    >>> df = hhsuite.read_hhr(hhrdir)
    '''
    if not isinstance(indir,list):
        if os.path.exists(indir) and os.path.isdir(indir):
            indir = os.path.realpath(indir)
            indir = glob(f'{indir}/*{suffix}')
    return hhsuite2pandas(indir)

def parse_hhr(indir, suffix=".hhr"):
    '''
    Parse HH-suite (``*.hhr``) output files.

    Parameters
    ----------
    indir : string
      A HH-suite results file, a list of files or a directory with results
    suffix: string, default .hhr
      The suffix, a.k.a extension, for the HH-suite results in indir

    Returns
    -------
      A generator of Pandas dataframes

    Examples
    --------
    Iterating over each file

    >>> from rotifer.io import hhsuite
    >>> for hhr in hhsuite.parse_hhr(hhrdir):
    ...     hhr.to_csv("hhr.tsv", sep="\\t", index=False)
    '''
    if not isinstance(indir,list):
        if os.path.exists(indir) and os.path.isdir(indir):
            indir = os.path.realpath(indir)
            indir = glob(f'{indir}/*{suffix}')
    for infile in indir:
        yield hhsuite2pandas(infile)
