"""
Genome files, feature tables and gene neighbourhoods.

The modules work as a chain:

:mod:`rotifer.genome.io`
    Read GenBank, GFF and other Biopython formats into sequence records
    with assembly accessions (:func:`~rotifer.genome.io.parse`).
:mod:`rotifer.genome.utils`
    Turn those records into a table with one row per feature
    (:func:`~rotifer.genome.utils.seqrecords_to_dataframe`).
:mod:`rotifer.genome.data`
    :class:`~rotifer.genome.data.NeighborhoodDF`, the table class, with
    methods to find neighbourhoods, trim blocks, draw compact views,
    write reports and build networks.
:mod:`rotifer.genome.blocks`
    The presets behind ``NeighborhoodDF.to_blocks`` and ``to_yaml``.

:func:`~rotifer.genome.io.parse` and
:func:`~rotifer.genome.io.GenbankAssemblyFromFilename` are re-exported
here for older code that imports them from the package.
"""

from rotifer.genome.io import parse
from rotifer.genome.io import GenbankAssemblyFromFilename
