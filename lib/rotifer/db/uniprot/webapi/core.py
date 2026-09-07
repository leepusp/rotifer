"""
Shared plumbing for UniProt's REST API.

This module holds what every web cursor in
:mod:`rotifer.db.uniprot.webapi` needs and none of them should repeat:
the description of each UniProt resource, the rules that decide which
resource an identifier belongs to, an HTTP session that retries the way
UniProt asks clients to, and a cursor base class that turns all of it
into the ``fetcher``/``parser`` pair the rotifer cursor machinery calls.

Resources and masking
---------------------
UniProt publishes its data as several resources -- UniProtKB, UniParc,
UniRef, proteomes, taxonomy -- each with its own endpoint, its own
identifier syntax and its own idea of what an entry is. The cursors
here are meant to hide that: a caller passes identifiers and gets
records back, whether the identifiers name UniProtKB accessions,
UniParc UPIs or a whole proteome. :func:`detect_resource` is what makes
that possible, and :attr:`RESOURCES` is what it consults.

Retrieval follows UniProt's own recommendations. Batches go to the
endpoints built for them (``/uniprotkb/accessions``, ``/uniparc/upis``)
rather than to one request per identifier; result sets are walked with
the opaque cursor UniProt returns in the ``Link`` header rather than by
computing offsets; and large result sets can be taken from ``/stream``,
which is what UniProt offers instead of deep pagination.

See Also
--------
rotifer.db.uniprot.webapi.entries : sequence and entry cursors
rotifer.db.uniprot.webapi.search : generic paged search
rotifer.db.uniprot.webapi.idmapping : the asynchronous mapping service
"""

# Import external modules
import re
import types
import typing
import requests
from requests.adapters import HTTPAdapter, Retry

# Import rotifer modules
import rotifer
import rotifer.db.core
import rotifer.db.parallel
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

#: Root of UniProt's REST service.
API_URL = "https://rest.uniprot.org"

# Configuration
#
# The name is spelled out rather than derived from ``__name__`` so that
# every module of this subpackage shares one configuration file,
# ``<user_path>/db/uniprot/webapi.yml``, and so that importing this
# module from the package ``__init__`` cannot become circular.
_defaults = {
    'batch_size': 100,
    'maxgetitem': 100,
    'threads': 5,
    'timeout': 60,
    'tries': 3,
    'page_size': 500,
}
config = loadConfig(':db.uniprot.webapi', defaults=_defaults)


class Resource:
    """
    One UniProt resource and the few facts the cursors need about it.

    Parameters
    ----------
    name : str
        Resource name, which is also the first path segment of its
        endpoints, e.g. ``uniprotkb``.
    query_field : str
        Name of the search field that matches an entry's own
        identifier, used to build queries for resources that have no
        batch endpoint.
    batch_path : str, optional
        Path of the endpoint that accepts several identifiers at once.
        ``None`` when the resource has none.
    batch_param : str, optional
        Name of the query parameter carrying those identifiers.
    sequences : bool, default False
        Whether entries of this resource carry sequences and can
        therefore be asked for in FASTA.
    id_key : str, optional
        Key holding an entry's own identifier in the JSON this
        resource returns. Each resource spells it differently, and
        tabular cursors need it to tell which entries came back.
    """

    def __init__(self, name, query_field, batch_path=None, batch_param=None,
                 sequences=False, id_key=None):
        self.name = name
        self.query_field = query_field
        self.batch_path = batch_path
        self.batch_param = batch_param
        self.sequences = sequences
        self.id_key = id_key

    def __repr__(self):
        return f'Resource({self.name!r})'


#: Every resource these cursors can address, keyed by name.
RESOURCES = {
    'uniprotkb': Resource('uniprotkb', 'accession', 'uniprotkb/accessions', 'accessions',
                          sequences=True, id_key='primaryAccession'),
    'uniparc':   Resource('uniparc',   'upi',       'uniparc/upis',         'upis',
                          sequences=True, id_key='uniParcId'),
    'uniref':    Resource('uniref',    'id',        sequences=True, id_key='id'),
    'proteomes': Resource('proteomes', 'upid',      id_key='id'),
    'taxonomy':  Resource('taxonomy',  'tax_id',    id_key='taxonId'),
}

# Identifier syntax, in the order it must be tried.
#
# UniRef comes before UniProtKB because a UniRef identifier ends with
# the accession of its representative, and a proteome identifier is
# tested before UniParc only for readability: UP is followed by digits
# in one and by an I in the other, so they cannot collide.
_PATTERNS = (
    ('uniref',    re.compile(r'^UniRef(100|90|50)_\S+$')),
    ('proteomes', re.compile(r'^UP[0-9]{9}$')),
    ('uniparc',   re.compile(r'^UPI[0-9A-F]{10}$')),
    # UniProt's own accession syntax, allowing an isoform suffix.
    ('uniprotkb', re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|'
                             r'[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$')),
    # UniProtKB entry names, e.g. TPA_HUMAN.
    ('uniprotkb', re.compile(r'^[A-Z0-9]{1,10}_[A-Z0-9]{1,5}$')),
    ('taxonomy',  re.compile(r'^[0-9]+$')),
)

#: Order in which the API is asked about an identifier no pattern matched.
_PROBE_ORDER = ('uniprotkb', 'uniparc', 'uniref', 'proteomes', 'taxonomy')


def build_session(tries=None, backoff=0.25):
    """
    Build an HTTP session that retries the way UniProt asks clients to.

    UniProt returns 429 when a client is going too fast and the usual
    5xx family when a node is briefly unavailable; both are worth
    retrying after a pause, and nothing else is.

    Parameters
    ----------
    tries : int, optional
        Total number of attempts per request. Defaults to the
        ``tries`` configuration entry.
    backoff : float, default 0.25
        Backoff factor between attempts.

    Returns
    -------
    requests.Session
    """
    if isinstance(tries, types.NoneType):
        tries = config['tries'] or _defaults['tries']
    retries = Retry(
        total=tries,
        backoff_factor=backoff,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(['GET', 'POST']),
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session


def detect_resource(accession, probe=True, session=None, timeout=None):
    """
    Decide which UniProt resource an identifier belongs to.

    Identifier syntax is checked first, because it costs nothing and
    answers almost every case. Only when no pattern matches is UniProt
    itself asked, by requesting the identifier from each resource in
    turn until one acknowledges it.

    Parameters
    ----------
    accession : str
        A UniProt identifier of any kind.
    probe : bool, default True
        Whether to ask the API when no pattern matches. Set to False
        to keep detection entirely offline.
    session : requests.Session, optional
        Session used for probing. One is built when not given.
    timeout : int, optional
        Timeout of each probe request.

    Returns
    -------
    str or None
        The resource name, or ``None`` when nothing recognised the
        identifier.

    Examples
    --------
    >>> from rotifer.db.uniprot.webapi import core
    >>> core.detect_resource('P00750', probe=False)
    'uniprotkb'
    >>> core.detect_resource('UPI0000000001', probe=False)
    'uniparc'
    >>> core.detect_resource('UniRef50_P00750', probe=False)
    'uniref'
    >>> core.detect_resource('UP000005640', probe=False)
    'proteomes'
    """
    accession = str(accession).strip()
    if not accession:
        return None
    for name, pattern in _PATTERNS:
        if pattern.match(accession):
            return name
    if not probe:
        return None

    if isinstance(session, types.NoneType):
        session = build_session()
    if isinstance(timeout, types.NoneType):
        timeout = config['timeout'] or _defaults['timeout']
    for name in _PROBE_ORDER:
        url = f'{API_URL}/{name}/{requests.utils.quote(accession)}'
        try:
            reply = session.get(url, params={'format': 'json'}, timeout=timeout)
        except Exception:
            logger.debug(f'Probing {url} failed', exc_info=1)
            continue
        if reply.status_code == 200:
            return name
    return None


def group_by_resource(accessions, database='auto', probe=True, session=None, timeout=None):
    """
    Split identifiers into the resources that can answer for them.

    This is what lets one cursor accept a mixture of UniProtKB
    accessions, UniParc UPIs and proteome identifiers: each group is
    then fetched from its own endpoint and the results are merged.

    Parameters
    ----------
    accessions : iterable of str
        Identifiers to sort.
    database : str, default 'auto'
        Name of a resource to assign every identifier to, bypassing
        detection. ``auto`` detects each one.
    probe : bool, default True
        Passed to :func:`detect_resource`.
    session : requests.Session, optional
        Passed to :func:`detect_resource`.
    timeout : int, optional
        Passed to :func:`detect_resource`.

    Returns
    -------
    tuple of (dict, set)
        A dictionary mapping resource names to sorted lists of
        identifiers, and the set of identifiers no resource claimed.

    Raises
    ------
    ValueError
        If `database` names an unknown resource.
    """
    if database != 'auto' and database not in RESOURCES:
        raise ValueError(f'Unknown UniProt resource {database!r}: '
                         f'expected "auto" or one of {sorted(RESOURCES)}')
    groups = dict()
    unknown = set()
    for accession in accessions:
        if database == 'auto':
            name = detect_resource(accession, probe=probe, session=session, timeout=timeout)
        else:
            name = database
        if isinstance(name, types.NoneType):
            unknown.add(accession)
            continue
        groups.setdefault(name, []).append(accession)
    return { k: sorted(v) for k, v in groups.items() }, unknown


class BaseUniProtWebCursor(rotifer.db.parallel.SimpleParallelProcessCursor):
    """
    Base class of every cursor that reads UniProt's REST API.

    It adds three things to the generic parallel cursor: a session
    that survives being sent to a worker process, the extra query
    parameters a caller passed as keyword arguments, and the two ways
    UniProt offers to walk a result set.

    Parameters
    ----------
    database : str, default 'auto'
        Which UniProt resource to query. ``auto`` decides from each
        identifier, which is what lets one cursor serve UniProtKB,
        UniParc, UniRef and proteomes at once. Naming a resource
        explicitly skips detection.
    probe : bool, default True
        Whether detection may ask the API about identifiers no
        pattern matched. False keeps detection offline.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, optional
        Number of attempts per request.
    batch_size : int, optional
        Number of identifiers per parallel batch.
    threads : int, optional
        Number of worker processes. Kept low by default: UniProt asks
        that clients not open many connections at once.
    timeout : int, optional
        Timeout of each request, in seconds.
    **kwargs
        Any other keyword argument is passed to UniProt as a query
        parameter, so that options this class knows nothing about --
        ``fields``, ``includeIsoform``, ``compressed`` and the rest --
        remain available. See
        https://www.uniprot.org/help/api_queries.

    Attributes
    ----------
    parameters : dict
        The extra query parameters collected from ``**kwargs``.
    """

    #: Format asked of the API. Overridden by subclasses.
    _format = 'json'

    def __init__(
            self,
            database='auto',
            probe=True,
            progress=True,
            tries=config['tries'] or _defaults['tries'],
            batch_size=config['batch_size'] or _defaults['batch_size'],
            threads=config['threads'] or _defaults['threads'],
            timeout=config['timeout'] or _defaults['timeout'],
            *args, **kwargs):
        # Every remaining keyword is an API query parameter, so they
        # must be taken out before the cursor machinery sees them.
        self.parameters = { k: v for k, v in kwargs.items() if k not in self._reserved }
        for key in self.parameters:
            kwargs.pop(key, None)
        super().__init__(progress=progress, tries=tries, batch_size=batch_size,
                         threads=threads, *args, **kwargs)
        self.database = database
        self.probe = probe
        self.timeout = timeout
        self.maxgetitem = config['maxgetitem'] or _defaults['maxgetitem']
        self._session = None

        # A request UniProt refuses is not worth repeating, and an
        # identifier it does not know is not worth asking twice.
        self.giveup.update(["400 Client Error", "404 Client Error"])

    #: Keywords consumed by the cursor machinery rather than by UniProt.
    _reserved = frozenset({
        'progress', 'tries', 'batch_size', 'threads', 'timeout',
        'database', 'probe', 'readers', 'writers',
    })

    # A requests.Session holds sockets and cannot be pickled, but
    # fetchone() sends the cursor to worker processes. Dropping it here
    # lets each process build its own on first use.
    def __getstate__(self):
        state = self.__dict__.copy()
        state['_session'] = None
        return state

    @property
    def session(self):
        """
        HTTP session of this process, built on first use.

        Returns
        -------
        requests.Session
        """
        if isinstance(getattr(self, '_session', None), types.NoneType):
            self._session = build_session(tries=self.tries)
        return self._session

    def resource(self, name):
        """
        Look up one entry of :attr:`RESOURCES`.

        Parameters
        ----------
        name : str
            Resource name.

        Returns
        -------
        Resource
        """
        return RESOURCES[name]

    def query_parameters(self, **overrides):
        """
        Merge the caller's extra parameters with this call's own.

        Parameters given here win, so that a cursor can set the format
        it knows how to parse while still letting the caller add
        ``fields`` or any other option.

        Parameters
        ----------
        **overrides
            Parameters set by the calling method.

        Returns
        -------
        dict
        """
        params = dict(self.parameters)
        params.update({ k: v for k, v in overrides.items() if not isinstance(v, types.NoneType) })
        return params

    def get(self, path, **params):
        """
        Run one GET request against the API.

        Parameters
        ----------
        path : str
            Path below :data:`API_URL`, without a leading slash.
        **params
            Query parameters.

        Returns
        -------
        requests.Response

        Raises
        ------
        requests.HTTPError
            If the reply carries an error status.
        """
        url = f'{API_URL}/{path}'
        reply = self.session.get(url, params=params, timeout=self.timeout)
        reply.raise_for_status()
        return reply

    def pages(self, path, **params):
        """
        Walk a result set by following UniProt's ``Link`` header.

        UniProt paginates with an opaque cursor rather than with
        offsets, and returns the next page's whole URL in a ``Link``
        header. Following it is what UniProt recommends over computing
        ``offset`` values, which deep result sets do not support.

        Parameters
        ----------
        path : str
            Path of the search endpoint, e.g. ``uniprotkb/search``.
        **params
            Query parameters of the first request.

        Yields
        ------
        requests.Response
            One response per page.
        """
        reply = self.get(path, **params)
        yield reply
        url = _next_link(reply.headers)
        while url:
            reply = self.session.get(url, timeout=self.timeout)
            reply.raise_for_status()
            yield reply
            url = _next_link(reply.headers)

    def stream(self, path, **params):
        """
        Fetch a whole result set from UniProt's ``/stream`` endpoint.

        ``/stream`` returns every match in one response and is what
        UniProt offers instead of paging through very large result
        sets.

        Parameters
        ----------
        path : str
            Path of the stream endpoint, e.g. ``uniprotkb/stream``.
        **params
            Query parameters.

        Returns
        -------
        requests.Response
        """
        return self.get(path, **params)


_NEXT_LINK = re.compile(r'<(.+)>; *rel="next"')


def _next_link(headers):
    """
    Extract the next page's URL from a ``Link`` header.

    Parameters
    ----------
    headers : mapping
        Response headers.

    Returns
    -------
    str or None
        The URL, or ``None`` when this was the last page.
    """
    link = headers.get("Link")
    if not link:
        return None
    match = _NEXT_LINK.match(link)
    return match.group(1) if match else None
