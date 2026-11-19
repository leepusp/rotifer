# Frozen copy of lib/rotifer/taxonomy/utils.py as of master 0d57415, kept
# byte for byte below this header. test_lineage_equivalence.py runs it next
# to the current module to prove the rewrite returns the same bytes on every
# input the old code handled. Delete it together with the loadPreferredTaxa
# shim, on the shim's removal date.
#!/usr/bin/env python3

# Data here is a child of Pandas DataFrame

import pandas as pd
from rotifer.core.functions import findDataFiles

def loadPreferredTaxa():
    preferred_taxa = []
    for p in findDataFiles(':taxonomy.taxonomy.txt', all=True):
        p = open(p,'rt')
        preferred_taxa.extend([ s.strip("\n").lower() for s in p.readlines() if s[0] != '#' ])
        p.close()
    preferred_taxa = set(preferred_taxa)
    return preferred_taxa

# Extract list of important lineages from taxonomy
def lineage(data=None, preferred_taxa=None, column="classification", insep="; ", outsep=">", case=str.lower):
    '''
    Transform detailed linearized taxonomies into reduced lineages using a set of preferred taxa.

    Available arguments:
      data           : may be Pandas Series or DataFrames
      preferred_taxa : list of clades that should be kept
      column         : for dataframes, name of the column containing the taxonomies
      insep          : clade name separator in input data
      outsep         : clade name separator in the returned stream
      case           : subroutine to handle case independent matching, set to None to disable
    '''
    def same(x):
        return(x)
    if case == None:
        case = same
    if not preferred_taxa:
        preferred_taxa = loadPreferredTaxa()
    if isinstance(data,pd.DataFrame):
        data = data[column]
    if isinstance(data,pd.Series):
        lineage = data.astype(str).str.split(insep).explode()
        lineage = lineage.groupby(level=0).agg(lambda x: outsep.join([ case(y) for y in x if case(y) in preferred_taxa ]))
    else:
        lineage = [ outsep.join([ case(y) for y in x.split(insep) if case(y) in preferred_taxa ]) for x in data ]
    return lineage

