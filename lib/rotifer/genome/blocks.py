"""Genome presets for rotifer.pandas.to_blocks and rotifer.pandas.yaml_blocks.to_yaml.

rotifer.pandas knows nothing about neighborhoods, so the choices that make a
neighborhood table readable (what to group by, which columns to show, how to sort
and summarise a block) live here. NeighborhoodDF.to_blocks and NeighborhoodDF.to_yaml
apply them; any keyword passed to those methods replaces the preset's value.

The presets are plain module-level dictionaries, so they can be copied and edited
in Python. They are deliberately not read through rotifer.core.functions.loadConfig:
the YAML overrides it offers cannot express the lambdas the YAML preset needs.
"""

import pandas as pd

_COLUMNS = ['pid', 'c80i70', 'c100i100', 'plen', 'pfam', 'aravind', 'product',
            'organism', 'lineage', 'classification']

NEIGHBORHOOD_BLOCKS = dict(
    groupby='c80e3',
    columns=_COLUMNS,
    stats={'lineage': 'nunique', 'block_id': 'nunique', 'c80i70': 'nunique',
           'pid': 'nunique'},
    sortrows={'lineage': True, 'classification': True, 'organism': True},
)

NEIGHBORHOOD_YAML = dict(
    groupby='c80e3',
    columns=_COLUMNS,
    sortrows={'qc80e3': True, 'superkingdom': True, 'phylum': True, 'lineage': True,
              'classification': True, 'c80i70': True},
    sortblocks={'phylum': False, 'c80i70': False},
    annotation=lambda x: pd.DataFrame({'pfam': x.pfam.mode(), 'curated': x.pfam.mode()}),
    stats={
        'superkingdom': ('superkingdom', 'nunique'),
        'lineage': ('lineage', 'nunique'),
        'phylum': ('phylum', 'nunique'),
        'taxid': ('taxid', 'nunique'),
        'c80i70': ('c80i70', 'nunique'),
        'pid': ('pid', 'nunique'),
    },
    dist=lambda x: x.groupby('pfam').c80i70.nunique().reset_index(),
)
