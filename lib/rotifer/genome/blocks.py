"""Presets for writing a NeighborhoodDF as block reports.

:meth:`rotifer.genome.data.NeighborhoodDF.to_blocks` and
:meth:`rotifer.genome.data.NeighborhoodDF.to_yaml` apply the presets defined
here. Each preset is a plain dictionary of keyword arguments for
:func:`rotifer.pandas.to_blocks` or :func:`rotifer.pandas.to_yaml`. A keyword given
to either method replaces the preset's value for that call. The dictionaries can
also be copied and edited in Python.

:mod:`rotifer.pandas` knows nothing about gene neighbourhoods, so the choices that
make a neighbourhood table readable live here: group by the sequence cluster
``c80e3``, show the protein, cluster, domain, product and taxonomy columns, and
summarise each block by how many distinct lineages, proteins and clusters it
holds. The columns named are those of a NeighborhoodDF that has been annotated
with sequence clusters (``c80e3``, ``c80i70``, ``c100i100``) and with taxonomy
(``organism``, ``lineage``, ``classification`` and, for the YAML preset,
``superkingdom``, ``phylum`` and ``taxid``). A missing column raises KeyError when it is named by ``groupby`` or ``stats``,
and is ignored when it is named by ``columns``, ``sortrows`` or ``sortblocks``.

The presets are not read through :func:`rotifer.core.functions.loadConfig`, because
the YAML overrides it offers cannot express the lambdas the YAML preset needs.
"""

import pandas as pd

# Columns shown in each block: protein id, clusters at 80%/70% and 100%/100%, length,
# domains, product, and the taxonomy of the source organism.
_COLUMNS = ['pid', 'c80i70', 'c100i100', 'plen', 'pfam', 'aravind', 'product',
            'organism', 'lineage', 'classification']

# Preset for NeighborhoodDF.to_blocks.
NEIGHBORHOOD_BLOCKS = dict(
    groupby='c80e3',            # one block per sequence cluster
    columns=_COLUMNS,
    # Header statistics: number of distinct lineages, neighborhoods, clusters, proteins.
    stats={'lineage': 'nunique', 'block_id': 'nunique', 'c80i70': 'nunique',
           'pid': 'nunique'},
    # Rows in a block ordered by taxonomy, so related organisms sit together.
    sortrows={'lineage': True, 'classification': True, 'organism': True},
)

# Preset for NeighborhoodDF.to_yaml.
NEIGHBORHOOD_YAML = dict(
    groupby='c80e3',            # one block per sequence cluster
    columns=_COLUMNS,
    # Rows ordered by taxonomy, then by 70% cluster. 'qc80e3' matches no column of a
    # NeighborhoodDF, so that entry is ignored.
    sortrows={'qc80e3': True, 'superkingdom': True, 'phylum': True, 'lineage': True,
              'classification': True, 'c80i70': True},
    # Blocks with the most phyla first, then those with the most 70% clusters.
    sortblocks={'phylum': False, 'c80i70': False},
    # Callable table: most common Pfam domain(s) of the block.
    annotation=lambda x: pd.DataFrame({'pfam': x.pfam.mode(), 'curated': x.pfam.mode()}),
    # Per-block counts of distinct taxa, clusters and proteins; sortblocks reads them.
    stats={
        'superkingdom': ('superkingdom', 'nunique'),
        'lineage': ('lineage', 'nunique'),
        'phylum': ('phylum', 'nunique'),
        'taxid': ('taxid', 'nunique'),
        'c80i70': ('c80i70', 'nunique'),
        'pid': ('pid', 'nunique'),
    },
    # Callable table: number of 70% clusters carrying each Pfam domain in the block.
    dist=lambda x: x.groupby('pfam').c80i70.nunique().reset_index(),
)
