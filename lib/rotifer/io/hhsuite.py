# HH-suite parser

__author__ = "Robson Francisco de Souza, Rodolfo Alvarenga Ribeiro"
__version__ = 0.002

# Libraries
import os
import sys
import pandas as pd
from glob import glob
from rotifer.core.functions import import_path

# Add path to rotifer libraries to Python's library search list
hhr = import_path(os.path.realpath(__file__).split("rotifer")[0] + "rotifer/bin/hhsuite2table")

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
    return hhr.hhsuite2pandas(indir)

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
        yield hhr.hhsuite2pandas(infile)
