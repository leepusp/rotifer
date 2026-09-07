"""
UniProt's identifier mapping service.

UniProt does not answer mapping queries from a URL: a job is submitted,
polled until it finishes and only then read, possibly over several
pages. :class:`MappingCursor` wraps that exchange so it looks like
every other cursor, and returns the same three column table as
:class:`rotifer.db.uniprot.mirror.MappingCursor` and
:class:`rotifer.db.uniprot.clickhouse.MappingCursor`, so the three
can be used interchangeably or stacked behind one delegator.

The module level functions below are the thin wrappers around each
step of the exchange. They predate the cursor and are kept because
they are convenient on their own.

See Also
--------
rotifer.db.uniprot.mirror.MappingCursor : same table, from a local file
rotifer.db.uniprot.clickhouse.MappingCursor : same table, indexed
"""

# Import external modules
import re
import csv
import time
import json
import zlib
import types
import pandas as pd
import requests
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
from requests.adapters import HTTPAdapter, Retry

# Import rotifer modules
import rotifer
import rotifer.db.methods
from rotifer.db.uniprot.webapi import core
logger = rotifer.logging.getLogger(__name__)

POLLING_INTERVAL = 3
API_URL = core.API_URL

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


class MappingCursor(rotifer.db.methods.MappingCursor, core.BaseUniProtWebCursor):
    """
    Map identifiers between UniProt and the databases it references.

    Results are returned as the same three column table the other
    identifier mapping cursors produce: the accession that was
    matched, the name of the database the mapping points at, and the
    identifier in that database.

    Parameters
    ----------
    polling_interval : int, default 3
        Seconds between polls of a running job.
    **kwargs
        Extra UniProt query parameters.

    Attributes
    ----------
    columns : list of str
        ``['accession', 'id_type', 'id']``.

    Examples
    --------
    Map GenBank protein identifiers onto UniProtKB:

    >>> from rotifer.db.uniprot import webapi            # doctest: +SKIP
    >>> mc = webapi.MappingCursor()                        # doctest: +SKIP
    >>> mc.fetchall(['BAE76179.1'], source=['EMBL-CDS'],   # doctest: +SKIP
    ...             target=mc.UNIPROTKB)

    Go the other way, from UniProtKB to RefSeq:

    >>> mc.fetchall(['P00750'], source=mc.UNIPROTKB,       # doctest: +SKIP
    ...             target=['RefSeq'])
    """

    # Names a delegator shares with its backends, which describe where
    # to look rather than what to ask, and must not reach UniProt as
    # query parameters.
    _reserved = core.BaseUniProtWebCursor._reserved | frozenset({
        'polling_interval', 'path', 'engine', 'host', 'port', 'dbname',
        'table', 'release', 'id_type', 'source', 'target',
    })

    def __init__(self, polling_interval=POLLING_INTERVAL, *args, **kwargs):
        # Identifiers are named by the caller, not detected from their
        # syntax, so resource detection has nothing to add here.
        kwargs.setdefault('database', 'uniprotkb')
        kwargs.setdefault('probe', False)
        super().__init__(*args, **kwargs)
        self.polling_interval = polling_interval
        self._databases = None
        # The service takes large batches, and every batch costs a
        # submission and a poll, so few and large is the cheap shape.
        self.maxgetitem = 5000

    def fetcher(self, accession, source=None, target=None, *args, **kwargs):
        """
        Submit a mapping job per pair of databases and read the results.

        UniProt's service maps between one pair at a time, so a query
        naming several databases becomes several jobs. Each costs a
        submission and at least one poll, which is why the ends must be
        named: there are 94 possible targets, and enumerating them
        would be minutes of polling for a single query.

        Parameters
        ----------
        accession : iterable of str
            Identifiers to translate.
        source, target : str or list of str
            The two ends of the mapping. Neither may be omitted.

        Returns
        -------
        list of dict
            One entry per result row, carrying the pair it came from.

        Raises
        ------
        ValueError
            If either end is left open.
        """
        targets = sorted(self.parse_ids(accession))
        if not targets:
            return []
        source = self.parse_databases(source)
        target = self.parse_databases(target)
        if not source or not target:
            raise ValueError(
                'UniProt maps between one pair of databases at a time, so both '
                'source and target must be named: this backend cannot answer '
                '"every database".'
            )

        rows = []
        for from_db in source:
            for to_db in target:
                reply = self.session.post(
                    f'{API_URL}/idmapping/run',
                    data={'from': self._service_name(from_db),
                          'to': self._service_name(to_db),
                          'ids': ",".join(targets)},
                    timeout=self.timeout,
                )
                reply.raise_for_status()
                job = reply.json()['jobId']
                self._wait(job)
                for row in self._results(job):
                    row = dict(row)
                    row['_source_type'] = from_db
                    row['_target_type'] = to_db
                    rows.append(row)
        return rows

    def databases(self):
        """
        Name the databases UniProt's mapping service can translate.

        The service publishes its own vocabulary, so it is asked
        rather than guessed, and the answer is remembered for the
        life of the cursor. Names are reported in the spelling used by
        ``idmapping.dat``, so that they can be compared with what the
        local backends hold.

        Only databases usable at *both* ends are reported: the service
        marks each as readable, writable or both, and a mapping needs
        one of each. The list is smaller than it looks -- Pfam, for
        one, is neither, so no mapping job can reach it even though
        UniProt entries carry Pfam cross-references.

        Returns
        -------
        set of str or None
            None when the service cannot be reached, so that a network
            failure narrows nothing.
        """
        if not isinstance(getattr(self, '_databases', None), types.NoneType):
            return self._databases
        try:
            reply = self.session.get(f'{API_URL}/configure/idmapping/fields', timeout=self.timeout)
            reply.raise_for_status()
            payload = reply.json()
        except Exception:
            logger.debug('Could not list the databases of the mapping service', exc_info=1)
            return None
        names = set()
        for group in payload.get('groups', []) or []:
            for item in group.get('items', []) or []:
                name = item.get('name')
                if name and item.get('from') and item.get('to'):
                    names.add(self._local_name(name))
        self._databases = names
        return self._databases

    def _local_name(self, database):
        """
        Spell a service database name the way ``idmapping.dat`` does.

        Parameters
        ----------
        database : str
            Name as the service reports it.

        Returns
        -------
        str
        """
        return self._LOCAL_NAMES.get(database, database)

    def _service_name(self, database):
        """
        Spell a database name the way the mapping service expects.

        The service and ``idmapping.dat`` disagree on some names, and
        on the accession itself, which the service calls a field of
        UniProtKB rather than a database.

        Parameters
        ----------
        database : str
            Database name.

        Returns
        -------
        str
        """
        if database == self.UNIPROTKB:
            return 'UniProtKB_AC-ID'
        return self._SERVICE_NAMES.get(database, database)

    #: Names that differ between ``idmapping.dat`` and the service.
    _SERVICE_NAMES = {
        'RefSeq': 'RefSeq_Protein',
        'EMBL-CDS': 'EMBL-GenBank-DDBJ_CDS',
        'EMBL': 'EMBL-GenBank-DDBJ',
    }

    #: The same correspondence, read the other way.
    _LOCAL_NAMES = { v: k for k, v in _SERVICE_NAMES.items() }

    def _wait(self, job):
        """
        Poll a job until its results are ready.

        Parameters
        ----------
        job : str
            Job identifier returned by the service.

        Raises
        ------
        RuntimeError
            If the job reports an error status.
        """
        while True:
            reply = self.session.get(f'{API_URL}/idmapping/status/{job}', timeout=self.timeout)
            reply.raise_for_status()
            payload = reply.json()
            status = payload.get('jobStatus')
            if isinstance(status, types.NoneType) or status in ('FINISHED',):
                return
            if status in ('RUNNING', 'NEW', 'QUEUED'):
                time.sleep(self.polling_interval)
                continue
            raise RuntimeError(f'UniProt mapping job {job} failed: {status}')

    def _results(self, job):
        """
        Read every page of a finished job.

        Parameters
        ----------
        job : str
            Job identifier.

        Returns
        -------
        list of dict
        """
        rows = []
        params = {'format': 'json', 'size': self.page_size}
        for reply in self.pages(f'idmapping/results/{job}', **params):
            payload = reply.json()
            rows.extend(payload.get('results', []) or [])
        return rows

    #: Entries per page when reading results.
    page_size = 500

    def parser(self, stream, accession, *args, **kwargs):
        """
        Turn mapping results into this cursor's table.

        Parameters
        ----------
        stream : list of dict
            Output of :meth:`fetcher`.
        accession : iterable of str
            Identifiers that were requested.

        Returns
        -------
        pandas.DataFrame
            The columns listed in
            :attr:`~rotifer.db.methods.MappingCursor.columns`.
        """
        if isinstance(stream, types.NoneType) or not stream:
            return self.empty()
        rows = []
        for entry in stream:
            queried = entry.get('from')
            found = entry.get('to')
            # UniProtKB targets arrive as whole entries, everything
            # else as a bare identifier.
            if isinstance(found, dict):
                found = found.get('primaryAccession') or found.get('id')
            if isinstance(queried, types.NoneType) or isinstance(found, types.NoneType):
                continue
            source_type = entry.get('_source_type')
            target_type = entry.get('_target_type')
            # The accession is whichever end is UniProtKB; when neither
            # is, the service does not report the one it joined through.
            if source_type == self.UNIPROTKB:
                accession_value = str(queried)
            elif target_type == self.UNIPROTKB:
                accession_value = str(found)
            else:
                accession_value = None
            rows.append({
                'source': str(queried), 'source_type': source_type,
                'accession': accession_value,
                'target': str(found), 'target_type': target_type,
            })
        if not rows:
            return self.empty()
        return pd.DataFrame(rows, columns=self.columns)

# ----------------------------------------------------------------------
# The original function level interface, kept for callers that use it.
# ----------------------------------------------------------------------

def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        logger.error(response.text)
        raise


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                logger.debug(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    logger.debug(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)


def get_data_frame_from_tsv_results(tsv_results):
    reader = csv.DictReader(tsv_results, delimiter="\t", quotechar='"')
    t = pd.DataFrame(list(reader))
    return t


def genbank_to_uniprot(from_db="EMBL-GenBank-DDBJ_CDS", to_db="UniProtKB", ids=["BAE76179.1"]):
    job_id = submit_id_mapping(
        from_db=from_db,
        to_db=to_db,
        ids=ids)
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_stream(link+"?compressed=true&fields=accession%2Cxref_pdb%2Cxref_alphafolddb%2C&format=tsv")
    r = get_data_frame_from_tsv_results(results)
    return(r)


def AF_link(id_list=None):
    r = genbank_to_uniprot(ids=id_list)
    r['urlAF'] = "https://alphafold.ebi.ac.uk/files/AF-" + r['AlphaFoldDB'].str.split(';', expand=True)[0] + "-F1-model_v4.pdb"
    r.loc[r.urlAF == "https://alphafold.ebi.ac.uk/files/AF--F1-model_v4.pdb", 'urlAF'] = None
    r['urlPDB'] = "https://files.rcsb.org/download/" + r['PDB'].str.split(';', expand=True)[0 ]+ ".pdb"
    r.loc[r.urlPDB == "https://files.rcsb.org/download/.pdb", 'urlPDB'] = None
    r = r[r.urlAF.notnull()].reset_index()
    return r


def af_to_seq(seqobj):
    r = seqobj.copy()
    u = AF_link(r.df.id.to_list())
    for x in range(0, len(u.From)):
        r = r.add_pdb(pdb_id=u.From[x], pdb_file=u.urlAF[x])
    return r
