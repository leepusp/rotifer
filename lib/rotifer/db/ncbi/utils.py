"""
Helper functions for identical protein group (IPG) reports.

These utilities operate on the dataframes produced by
:class:`rotifer.db.ncbi.entrez.IPGCursor` and are used by the gene
neighborhood cursors to choose, for each queried protein, the best
genomic sequence to download.
"""

import types

def best_ipgs(ipgs):
    """
    Select the first ranked row of each identical protein group.

    Rows are ranked by their original position in the report (the
    ``order`` column), so the selected row corresponds to the entry
    NCBI listed first for each IPG.

    Parameters
    ----------
    ipgs : pandas.DataFrame
        Identical protein group reports.

    Returns
    -------
    pandas.DataFrame
        One row per IPG identifier.
    """
    best = ipgs.sort_values(['id','order'], ascending=[True,True])
    best = best.drop_duplicates(['id'], keep='first')
    return best

def ipgs_to_dicts(ipgs):
    """
    Split IPG rows into assembly based and nucleotide based maps.

    IPGs that name at least one genome assembly are grouped by
    assembly accession; the remaining IPGs are grouped by nucleotide
    accession.

    Parameters
    ----------
    ipgs : pandas.DataFrame
        Identical protein group reports.

    Returns
    -------
    tuple of dict
        A pair ``(assemblies, nucleotides)``. In both dictionaries
        each key maps to a dictionary from protein accession
        (``pid``) to IPG representative.
    """
    if len(ipgs) == 0:
        return dict(), dict()

    # Split
    assemblies = ipgs[ipgs.assembly.notna()].id.unique().tolist()
    nucleotides = ipgs[~ipgs.id.isin(assemblies)]
    assemblies = ipgs[ipgs.id.isin(assemblies)]

    # By assembly
    if len(assemblies) > 0:
        assemblies = assemblies.assembly.unique().tolist()
        assemblies = ipgs[ipgs.assembly.isin(assemblies)]
        assemblies = assemblies.filter(['assembly','pid','representative'])
        assemblies = assemblies.drop_duplicates(ignore_index=True)
        assemblies = assemblies.groupby('assembly').apply(lambda x: x.set_index('pid').representative.to_dict())
        assemblies = assemblies.to_dict()
    else:
        assemblies = dict()

    # By nucleotide
    if len(nucleotides) > 0:
        nucleotides = nucleotides.nucleotide.unique().tolist()
        nucleotides = ipgs[ipgs.nucleotide.isin(nucleotides)]
        nucleotides = nucleotides.filter(['nucleotide','pid','representative'])
        nucleotides = nucleotides.drop_duplicates(ignore_index=True)
        nucleotides = nucleotides.groupby('nucleotide').apply(lambda x: x.set_index('pid').representative.to_dict())
        nucleotides = nucleotides.to_dict()
    else:
        nucleotides = dict()

    return assemblies, nucleotides


