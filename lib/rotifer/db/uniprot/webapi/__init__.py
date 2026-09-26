"""
Query UniProt's REST web services.

:mod:`rotifer.db.uniprot.webapi.idmapping` implements the identifier
mapping workflow: submit a batch of accessions, poll until the job
finishes and download the results.

Note
----
The ``local_database_path`` configuration default below is not used
anywhere in this package; see ``docs/OPEN_QUESTIONS.md``.
"""

import os
import rotifer
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

Extra API parameters
--------------------
UniProt's API accepts many options these classes know nothing about.
Any keyword a cursor does not recognise is passed straight through as a
query parameter, so they all remain available:

>>> fc = webapi.FastaCursor(includeIsoform=True)                # doctest: +SKIP
>>> sc = webapi.SearchCursor(database='uniprotkb',              # doctest: +SKIP
...                          fields='accession,id,length')

Cursors
-------
:class:`SequenceCursor`
    Entries with their annotation, as ``Bio.SeqRecord`` objects.
:class:`FastaCursor`
    The same entries as plain sequences.
:class:`SearchCursor`
    Any resource, as a dataframe, by identifier or by query.
:class:`ProteomeCursor`
    Proteome descriptions.
:class:`TaxonomyCursor`
    Taxonomy records.
:class:`MappingCursor`
    Identifier translation, in either direction, read off the
    cross-references an entry carries.

See Also
--------
rotifer.db.uniprot : delegators combining this backend with the others
rotifer.db.ncbi : the equivalent cursors for NCBI
"""

import rotifer
from rotifer.db.uniprot.webapi.core import (
    API_URL,
    RESOURCES,
    Resource,
    build_session,
    config,
    detect_resource,
    group_by_resource,
)
from rotifer.db.uniprot.webapi.entries import FastaCursor, SequenceCursor
from rotifer.db.uniprot.webapi.search import (
    ProteomeCursor,
    SearchCursor,
    TaxonomyCursor,
)
from rotifer.db.uniprot.webapi.idmapping import MappingCursor

logger = rotifer.logging.getLogger(__name__)

__all__ = [
    'API_URL',
    'FastaCursor',
    'MappingCursor',
    'ProteomeCursor',
    'RESOURCES',
    'Resource',
    'SearchCursor',
    'SequenceCursor',
    'TaxonomyCursor',
    'build_session',
    'config',
    'detect_resource',
    'group_by_resource',
]
