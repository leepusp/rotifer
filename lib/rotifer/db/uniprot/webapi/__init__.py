"""
Read UniProt through its REST API.

This subpackage is the online backend of :mod:`rotifer.db.uniprot`. It
covers the resources UniProt publishes -- UniProtKB, UniParc, UniRef,
proteomes and taxonomy -- behind cursors that follow the same
interface as the NCBI ones, so that switching source is mostly a
matter of switching import.

Masking the resources
---------------------
Every cursor here takes ``database='auto'`` by default, which decides
from each identifier which resource should answer for it. One call may
therefore mix identifiers of different kinds:

>>> from rotifer.db.uniprot import webapi                       # doctest: +SKIP
>>> sc = webapi.SequenceCursor()                                # doctest: +SKIP
>>> records = sc.fetchall(['P00750','UPI0000000001','UP000000625'])  # doctest: +SKIP

Detection is by identifier syntax first and by asking the API second,
so the common cases cost no extra request. Naming a resource, e.g.
``database='uniref'``, skips detection entirely.

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
