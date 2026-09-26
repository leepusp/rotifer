"""
Tabular cursors for UniProt's REST API.

Where :mod:`rotifer.db.uniprot.webapi.entries` returns sequences, the
cursors here return rows: one per entry, as a
:class:`pandas.DataFrame`. That covers the resources that have no
sequences of their own -- proteomes and taxonomy -- and gives every
resource a way to be queried by anything other than an identifier.

:class:`SearchCursor` is the general case. It accepts identifiers like
any other cursor, but also takes a raw UniProt query through
:meth:`SearchCursor.search`, and passes any keyword it does not
recognise on to the API, so options this package never heard of remain
reachable.

Result sets are walked with the cursor token UniProt returns in the
``Link`` header, which is what UniProt recommends: offsets are not
supported beyond the first pages, and ``/stream`` is offered for
result sets too large to page through at all.

See Also
--------
rotifer.db.uniprot.webapi.core : resource detection and HTTP plumbing
"""

# Import external modules
import types
import pandas as pd

# Import rotifer modules
import rotifer
from rotifer.db.uniprot.webapi import core
logger = rotifer.logging.getLogger(__name__)


class BaseTableCursor(core.BaseUniProtWebCursor):
    """
    Base class of the cursors that return dataframes.

    Parameters
    ----------
    database : str, default 'auto'
        UniProt resource to query, or ``auto`` to decide from each
        identifier.
    page_size : int, optional
        Number of entries per page. Defaults to the ``page_size``
        configuration entry.
    stream : bool, default False
        Fetch whole result sets from UniProt's ``/stream`` endpoint
        instead of paging. Faster for large results, but UniProt
        returns them in one response, so memory grows with the result.
    **kwargs
        Extra UniProt query parameters, e.g. ``fields``.
    """

    _format = 'json'

    def __init__(self, database='auto', page_size=None, stream=False, *args, **kwargs):
        super().__init__(database=database, *args, **kwargs)
        self.page_size = page_size or core.config['page_size'] or core._defaults['page_size']
        self.use_stream = stream

    _reserved = core.BaseUniProtWebCursor._reserved | frozenset({'page_size', 'stream'})

    def empty(self):
        """
        Build an empty result.

        Returns
        -------
        pandas.DataFrame
        """
        return pd.DataFrame()

    def getids(self, obj, *args, **kwargs):
        """
        Extract identifiers from a result dataframe.

        Every column that can hold an entry's own identifier is
        scanned, because a mixed query returns rows from resources
        that name that column differently.

        Parameters
        ----------
        obj : pandas.DataFrame or None

        Returns
        -------
        set of str
        """
        if isinstance(obj, types.NoneType) or not isinstance(obj, pd.DataFrame) or obj.empty:
            return set()
        columns = { r.id_key for r in core.RESOURCES.values() if r.id_key }
        columns.add('query_id')
        ids = set()
        for column in columns:
            if column in obj.columns:
                ids.update(obj[column].dropna().astype(str))
        return ids

    def _collect(self, path, params):
        """
        Run one search and return every row it produced.

        Parameters
        ----------
        path : str
            Resource name, e.g. ``uniprotkb``.
        params : dict
            Query parameters.

        Returns
        -------
        list of dict
            The ``results`` arrays of every page, concatenated.
        """
        rows = []
        if self.use_stream:
            replies = [self.get(f'{path}/stream', **params)]
        else:
            replies = self.pages(f'{path}/search', **params)
        for reply in replies:
            payload = reply.json()
            rows.extend(payload.get('results', []) or [])
        return rows

    def fetcher(self, accession, *args, **kwargs):
        """
        Query every resource the identifiers belong to.

        Parameters
        ----------
        accession : iterable of str
            Identifiers to look up.

        Returns
        -------
        list of dict
            One entry per resource, carrying its rows.
        """
        targets = self.parse_ids(accession)
        groups, unknown = core.group_by_resource(
            targets, database=self.database, probe=self.probe,
            session=self.session, timeout=self.timeout,
        )
        if unknown:
            self.update_missing(unknown, 'Unrecognised UniProt identifier', retry=False)

        stream = []
        for name, accessions in groups.items():
            resource = self.resource(name)
            query = " OR ".join([ f'{resource.query_field}:{x}' for x in accessions ])
            params = self.query_parameters(query=query, format='json', size=self.page_size)
            if self.use_stream:
                params.pop('size', None)
            stream.append({'resource': name, 'rows': self._collect(name, params)})
        return stream

    def parser(self, stream, accession, *args, **kwargs):
        """
        Turn the collected rows into one dataframe.

        The JSON UniProt returns is deeply nested, so it is flattened
        one level: nested objects become dotted column names and lists
        are left as they are.

        Parameters
        ----------
        stream : list of dict
            Output of :meth:`fetcher`.
        accession : iterable of str
            Identifiers that were requested.

        Returns
        -------
        pandas.DataFrame
        """
        if isinstance(stream, types.NoneType):
            return self.empty()
        frames = []
        for chunk in stream:
            rows = chunk.get('rows') or []
            if not rows:
                continue
            frame = pd.json_normalize(rows)
            frame['resource'] = chunk['resource']
            frames.append(frame)
        if not frames:
            return self.empty()
        return pd.concat(frames, ignore_index=True)

    def fetchall(self, accessions, *args, **kwargs):
        """
        Fetch every identifier as a single dataframe.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to look up.

        Returns
        -------
        pandas.DataFrame
        """
        frames = [ x for x in self.fetchone(accessions, *args, **kwargs)
                   if isinstance(x, pd.DataFrame) and not x.empty ]
        if not frames:
            return self.empty()
        return pd.concat(frames, ignore_index=True)


class SearchCursor(BaseTableCursor):
    """
    Search any UniProt resource.

    Used like the other cursors it looks identifiers up, but it also
    exposes UniProt's query language directly, which is the only way
    to ask questions that are not about a known identifier.

    Parameters
    ----------
    database : str, default 'auto'
        UniProt resource. ``auto`` works when looking identifiers up;
        :meth:`search` needs a resource named explicitly, since a raw
        query carries no identifier to detect one from.
    page_size : int, optional
        Entries per page.
    stream : bool, default False
        Use ``/stream`` instead of paging.
    **kwargs
        Extra UniProt query parameters.

    Examples
    --------
    Look up identifiers of several kinds at once:

    >>> from rotifer.db.uniprot import webapi              # doctest: +SKIP
    >>> sc = webapi.SearchCursor()                         # doctest: +SKIP
    >>> df = sc.fetchall(['P00750','UP000005640'])         # doctest: +SKIP

    Run a query of your own, restricted to the fields you want:

    >>> sc = webapi.SearchCursor(database='uniprotkb',     # doctest: +SKIP
    ...                          fields='accession,id,length')
    >>> df = sc.search('taxonomy_id:9606 AND reviewed:true')   # doctest: +SKIP
    """

    def search(self, query, database=None, **kwargs):
        """
        Run a UniProt query and return every match.

        Parameters
        ----------
        query : str
            A query in UniProt's query language, e.g.
            ``taxonomy_id:9606 AND reviewed:true``.
        database : str, optional
            Resource to search. Defaults to this cursor's
            ``database``, which must then not be ``auto``.
        **kwargs
            Extra query parameters for this call only.

        Returns
        -------
        pandas.DataFrame

        Raises
        ------
        ValueError
            If no resource was named, here or on the cursor.
        """
        database = database or self.database
        if database == 'auto':
            raise ValueError(
                'search() needs a resource: a raw query carries no identifier '
                'to detect one from. Pass database=, or set it on the cursor.'
            )
        if database not in core.RESOURCES:
            raise ValueError(f'Unknown UniProt resource {database!r}: '
                             f'expected one of {sorted(core.RESOURCES)}')
        params = self.query_parameters(query=query, format='json', size=self.page_size, **kwargs)
        if self.use_stream:
            params.pop('size', None)
        rows = self._collect(database, params)
        return self.parser([{'resource': database, 'rows': rows}], [])


class ProteomeCursor(BaseTableCursor):
    """
    Fetch UniProt proteome descriptions.

    This is the UniProt counterpart of
    :class:`rotifer.db.ncbi.GenomeCursor`: it describes whole
    proteomes rather than the sequences in them. Use
    :class:`rotifer.db.uniprot.webapi.SequenceCursor` with the same
    identifiers to get the sequences.

    Parameters
    ----------
    **kwargs
        Extra UniProt query parameters.

    Examples
    --------
    >>> from rotifer.db.uniprot import webapi        # doctest: +SKIP
    >>> pc = webapi.ProteomeCursor()                 # doctest: +SKIP
    >>> df = pc.fetchall('UP000005640')              # doctest: +SKIP
    """

    def __init__(self, database='proteomes', *args, **kwargs):
        super().__init__(database=database, *args, **kwargs)


class TaxonomyCursor(BaseTableCursor):
    """
    Fetch UniProt taxonomy records.

    Parameters
    ----------
    **kwargs
        Extra UniProt query parameters.

    Examples
    --------
    >>> from rotifer.db.uniprot import webapi        # doctest: +SKIP
    >>> tc = webapi.TaxonomyCursor()                 # doctest: +SKIP
    >>> df = tc.fetchall([9606, 562])                # doctest: +SKIP

    See Also
    --------
    rotifer.db.ncbi.TaxonomyCursor : the same idea, for NCBI
    """

    def __init__(self, database='taxonomy', *args, **kwargs):
        super().__init__(database=database, *args, **kwargs)
