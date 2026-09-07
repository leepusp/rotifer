"""
UniProt's identifier mapping service.

UniProt does not answer mapping queries from a URL: a job is submitted,
polled until it finishes and only then read, possibly over several
pages. :class:`IdMappingCursor` wraps that exchange so it looks like
every other cursor, and returns the same three column table as
:class:`rotifer.db.uniprot.mirror.IdMappingCursor` and
:class:`rotifer.db.uniprot.clickhouse.IdMappingCursor`, so the three
can be used interchangeably or stacked behind one delegator.

The module level functions below are the thin wrappers around each
step of the exchange. They predate the cursor and are kept because
they are convenient on their own.

See Also
--------
rotifer.db.uniprot.mirror.IdMappingCursor : same table, from a local file
rotifer.db.uniprot.clickhouse.IdMappingCursor : same table, indexed
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


class IdMappingCursor(rotifer.db.methods.IdMappingCursor, core.BaseUniProtWebCursor):
    """
    Map identifiers between UniProt and the databases it references.

    Results are returned as the same three column table the other
    identifier mapping cursors produce: the accession that was
    matched, the name of the database the mapping points at, and the
    identifier in that database.

    Parameters
    ----------
    from_db : str, default 'UniProtKB_AC-ID'
        Database the queried identifiers belong to, in the spelling
        UniProt's service expects, e.g. ``RefSeq_Protein`` or
        ``EMBL-GenBank-DDBJ_CDS``.
    to_db : str, default 'UniProtKB'
        Database to map to. Reported in the ``id_type`` column.
    column : str, default 'accession'
        Which column the queried identifiers are matched against, kept
        for compatibility with the other identifier mapping cursors.
        ``accession`` means the query holds UniProtKB accessions,
        ``id`` that it holds identifiers of ``from_db``.
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
    >>> ic = webapi.IdMappingCursor(                     # doctest: +SKIP
    ...     from_db='EMBL-GenBank-DDBJ_CDS', to_db='UniProtKB', column='id')
    >>> df = ic.fetchall(['BAE76179.1'])                 # doctest: +SKIP

    Go the other way, from UniProtKB to RefSeq:

    >>> ic = webapi.IdMappingCursor(to_db='RefSeq_Protein')   # doctest: +SKIP
    >>> df = ic.fetchall(['P00750'])                          # doctest: +SKIP
    """

    _reserved = core.BaseUniProtWebCursor._reserved | frozenset({
        'from_db', 'to_db', 'column', 'polling_interval',
    })

    def __init__(
            self,
            from_db='UniProtKB_AC-ID',
            to_db='UniProtKB',
            column='accession',
            polling_interval=POLLING_INTERVAL,
            *args, **kwargs):
        # Identifiers are named by the service, not detected from
        # their syntax, so resource detection has nothing to add here.
        kwargs.setdefault('database', 'uniprotkb')
        kwargs.setdefault('probe', False)
        super().__init__(*args, **kwargs)
        self.from_db = from_db
        self.to_db = to_db
        self.column = column
        self.polling_interval = polling_interval
        # The service takes large batches, and every batch costs a
        # submission and a poll, so few and large is the cheap shape.
        self.maxgetitem = 5000

    def fetcher(self, accession, *args, **kwargs):
        """
        Submit a mapping job and read its results.

        Parameters
        ----------
        accession : iterable of str
            Identifiers to map.

        Returns
        -------
        list of dict
            The ``results`` array returned by the service.
        """
        targets = sorted(self.parse_ids(accession))
        if not targets:
            return []
        reply = self.session.post(
            f'{API_URL}/idmapping/run',
            data={'from': self.from_db, 'to': self.to_db, 'ids': ",".join(targets)},
            timeout=self.timeout,
        )
        reply.raise_for_status()
        job = reply.json()['jobId']
        self._wait(job)
        return self._results(job)

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
        Turn mapping results into the three column table.

        Parameters
        ----------
        stream : list of dict
            Output of :meth:`fetcher`.
        accession : iterable of str
            Identifiers that were requested.

        Returns
        -------
        pandas.DataFrame
            Columns ``accession``, ``id_type`` and ``id``.
        """
        if isinstance(stream, types.NoneType) or not stream:
            return self.empty()
        rows = []
        for entry in stream:
            source = entry.get('from')
            target = entry.get('to')
            # UniProtKB targets arrive as whole entries, everything
            # else as a bare identifier.
            if isinstance(target, dict):
                target = target.get('primaryAccession') or target.get('id')
            if isinstance(source, types.NoneType) or isinstance(target, types.NoneType):
                continue
            rows.append({'accession': str(source), 'id_type': self.to_db, 'id': str(target)})
        if not rows:
            return self.empty()
        return pd.DataFrame(rows, columns=self.columns)

    def getids(self, obj, *args, **kwargs):
        """
        Report which queried identifiers were mapped.

        The queried identifiers are always in the ``accession``
        column, whatever database they belong to, because that is the
        column the service echoes back.

        Parameters
        ----------
        obj : pandas.DataFrame or None

        Returns
        -------
        set of str
        """
        if isinstance(obj, types.NoneType):
            return set()
        if isinstance(obj, pd.DataFrame):
            if obj.empty or 'accession' not in obj.columns:
                return set()
            return set(obj['accession'].dropna().astype(str))
        return super().getids(obj)


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
