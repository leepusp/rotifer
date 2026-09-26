"""
Sequence cursors for UniProt's REST API.

The cursors here answer the question the NCBI cursors answer for
GenBank: given identifiers, return sequence records. What they add is
that the identifiers need not all name the same thing. A single call
may mix UniProtKB accessions, UniParc UPIs, UniRef cluster names and
whole proteomes; each group is fetched from the endpoint that knows
about it and the records come back together.

>>> from rotifer.db.uniprot import webapi                  # doctest: +SKIP
>>> sc = webapi.SequenceCursor()                           # doctest: +SKIP
>>> records = sc.fetchall(['P00750','UPI0000000001'])      # doctest: +SKIP

See Also
--------
rotifer.db.ncbi.entrez.SequenceCursor : the same idea, for NCBI
rotifer.db.uniprot.webapi.core : resource detection and HTTP plumbing
"""

# Import external modules
import io
import types
from Bio import SeqIO

# Import rotifer modules
import rotifer
import rotifer.db.methods
from rotifer.db.uniprot.webapi import core
logger = rotifer.logging.getLogger(__name__)


def _record_ids(record):
    """
    Every identifier by which a record might have been requested.

    A caller may name a UniProtKB entry by accession, by entry name or
    by a secondary accession, and the identifier that comes back
    depends on the format: FASTA carries ``sp|P00750|TPA_HUMAN`` while
    the flat file carries ``P00750``. Reporting all of them is what
    keeps an entry from being counted as missing merely because it was
    asked for by another of its names.

    Parameters
    ----------
    record : Bio.SeqRecord.SeqRecord

    Returns
    -------
    set of str
    """
    ids = set()
    for value in (getattr(record, 'id', None), getattr(record, 'name', None)):
        if not value:
            continue
        ids.add(value)
        # FASTA from UniProtKB: db|accession|entry_name
        if '|' in value:
            ids.update([ x for x in value.split('|') if x ])
    annotations = getattr(record, 'annotations', None) or {}
    for key in ('accessions', 'proteome'):
        value = annotations.get(key)
        if isinstance(value, str):
            ids.add(value)
        elif isinstance(value, (list, tuple, set)):
            ids.update([ str(x) for x in value ])
    return ids


class BaseSequenceCursor(rotifer.db.methods.SequenceCursor, core.BaseUniProtWebCursor):
    """
    Shared behaviour of the UniProt sequence cursors.

    Subclasses differ only in the formats they ask for, so everything
    that decides *which* request to send lives here.

    Parameters
    ----------
    database : str, default 'auto'
        UniProt resource to query, or ``auto`` to decide from each
        identifier. See
        :class:`rotifer.db.uniprot.webapi.core.BaseUniProtWebCursor`.
    **kwargs
        Extra UniProt query parameters, e.g. ``fields`` or
        ``includeIsoform``.
    """

    #: Resource name mapped to (API format, Bio.SeqIO format).
    _formats = {
        'uniprotkb': ('fasta', 'fasta'),
        'uniparc':   ('fasta', 'fasta'),
        'uniref':    ('fasta', 'fasta'),
        'proteomes': ('fasta', 'fasta'),
    }

    def getids(self, obj, *args, **kwargs):
        """
        Extract identifiers from parsed records.

        Parameters
        ----------
        obj : Bio.SeqRecord.SeqRecord, list or None
            Records produced by :meth:`parser`.

        Returns
        -------
        set of str
        """
        if isinstance(obj, types.NoneType):
            return set()
        if not isinstance(obj, (list, tuple)):
            obj = [obj]
        ids = set()
        for record in obj:
            ids.update(_record_ids(record))
        return ids

    def _chunks(self, targets):
        """
        Turn identifiers into the requests that will answer for them.

        Parameters
        ----------
        targets : iterable of str
            Identifiers to fetch.

        Returns
        -------
        list of dict
            One entry per request, carrying the resource name, the
            path, the query parameters and the Bio.SeqIO format the
            reply should be parsed with.
        """
        groups, unknown = core.group_by_resource(
            targets, database=self.database, probe=self.probe,
            session=self.session, timeout=self.timeout,
        )
        if unknown:
            self.update_missing(unknown, 'Unrecognised UniProt identifier', retry=False)

        chunks = []
        for name, accessions in groups.items():
            resource = self.resource(name)
            if name not in self._formats:
                # No sequences to be had here, by any route. A
                # proteome has none of its own either, but it expands
                # into entries that do, so the test is what this
                # cursor can fetch rather than what the resource
                # itself stores.
                self.update_missing(accessions, f'Resource {name} has no sequences',
                                    retry=False, final=True)
                continue
            api_format, seqio_format = self._formats[name]

            if name == 'proteomes':
                # A proteome is not an entry but a query: its
                # sequences are UniProtKB entries. One request each,
                # so every record can be attributed to its proteome.
                for upid in accessions:
                    chunks.append({
                        'resource': name,
                        'path': 'uniprotkb/stream',
                        'params': self.query_parameters(query=f'proteome:{upid}', format=api_format),
                        'format': seqio_format,
                        'proteome': upid,
                    })
            elif resource.batch_path:
                chunks.append({
                    'resource': name,
                    'path': resource.batch_path,
                    'params': self.query_parameters(**{resource.batch_param: ",".join(accessions),
                                                       'format': api_format}),
                    'format': seqio_format,
                    'proteome': None,
                })
            else:
                # No batch endpoint, so ask the search index for all
                # of them at once instead of one request each.
                query = " OR ".join([ f'{resource.query_field}:{x}' for x in accessions ])
                chunks.append({
                    'resource': name,
                    'path': f'{name}/stream',
                    'params': self.query_parameters(query=query, format=api_format),
                    'format': seqio_format,
                    'proteome': None,
                })
        return chunks

    def fetcher(self, accession, *args, **kwargs):
        """
        Send one request per resource and collect the replies.

        Parameters
        ----------
        accession : iterable of str
            Identifiers to fetch.

        Returns
        -------
        list of dict
            The chunks built by :meth:`_chunks`, each with the reply
            text added under ``text``.
        """
        targets = self.parse_ids(accession)
        stream = []
        for chunk in self._chunks(targets):
            reply = self.get(chunk['path'], **chunk['params'])
            chunk = dict(chunk)
            chunk['text'] = reply.text
            stream.append(chunk)
        return stream

    def parser(self, stream, accession, *args, **kwargs):
        """
        Parse the replies collected by :meth:`fetcher`.

        Parameters
        ----------
        stream : list of dict
            Chunks carrying the reply text.
        accession : iterable of str
            Identifiers that were requested.

        Returns
        -------
        list of Bio.SeqRecord.SeqRecord
        """
        if isinstance(stream, types.NoneType):
            return []
        records = []
        for chunk in stream:
            text = chunk.get('text') or ''
            if not text.strip():
                continue
            for record in SeqIO.parse(io.StringIO(text), chunk['format']):
                if chunk.get('proteome'):
                    record.annotations['proteome'] = chunk['proteome']
                records.append(record)
        return records


class FastaCursor(BaseSequenceCursor):
    """
    Fetch UniProt sequences in FASTA.

    Every resource that has sequences can produce FASTA, so this
    cursor covers UniProtKB, UniParc, UniRef and proteomes with one
    format.

    Parameters
    ----------
    database : str, default 'auto'
        UniProt resource, or ``auto``.
    **kwargs
        Extra UniProt query parameters.

    Examples
    --------
    >>> from rotifer.db.uniprot import webapi                   # doctest: +SKIP
    >>> fc = webapi.FastaCursor()                               # doctest: +SKIP
    >>> records = fc.fetchall(['P00750','UniRef50_P00750'])     # doctest: +SKIP

    Ask for isoforms as well, using a parameter of UniProt's API:

    >>> fc = webapi.FastaCursor(includeIsoform=True)            # doctest: +SKIP
    """
    pass


class SequenceCursor(BaseSequenceCursor):
    """
    Fetch UniProt entries with their annotation.

    UniProtKB entries are taken as flat files, which carry the
    features and cross-references FASTA leaves out, and are parsed
    with Biopython's ``swiss`` parser. The other resources publish no
    such format, so they fall back to FASTA; a mixed query therefore
    returns annotated records for UniProtKB and plain ones for the
    rest.

    Parameters
    ----------
    database : str, default 'auto'
        UniProt resource, or ``auto``.
    **kwargs
        Extra UniProt query parameters.

    Examples
    --------
    >>> from rotifer.db.uniprot import webapi         # doctest: +SKIP
    >>> sc = webapi.SequenceCursor()                  # doctest: +SKIP
    >>> record = sc.fetchall('P00750')                # doctest: +SKIP

    See Also
    --------
    FastaCursor : the same entries without annotation
    """

    _formats = {
        'uniprotkb': ('txt', 'swiss'),
        'uniparc':   ('fasta', 'fasta'),
        'uniref':    ('fasta', 'fasta'),
        'proteomes': ('txt', 'swiss'),
    }
