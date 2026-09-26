__doc__ = """
Access biological data stored on the local machine.

This package groups the cursors that never touch the network:
:mod:`rotifer.db.local.easel` fetches sequences from FASTA files
indexed by Easel's ``esl-sfetch`` and :mod:`rotifer.db.local.ete3`
queries the local copy of the NCBI Taxonomy database managed by the
ETE toolkit.

Configuration
-------------
The default sequence database path is read from the user
configuration (key ``local_database_path``) and falls back to
``fadb/nr/nr`` under the directory named by the ``ROTIFER_DATA``
environment variable.
"""

import re
import os
import types
import typing
import sqlite3
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from io import StringIO

import rotifer
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
from rotifer.db.core import BaseCursor
from rotifer.genome.data import NeighborhoodDF
from rotifer.genome.utils import seqrecords_to_dataframe
import rotifer.devel.beta.sequence as rdbs
logger = rotifer.logging.getLogger(__name__)

# Defaults
config = loadConfig(__name__.replace('rotifer.',':'), defaults = {
    'local_database_path': os.path.join(GlobalConfig['data'],"fadb","nr","nr"),
})
