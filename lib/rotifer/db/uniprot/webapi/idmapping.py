"""
Translate identifiers through UniProt's entry cross-references.

A UniProtKB entry lists the databases it is cross-referenced to, so
fetching entries answers a mapping in either direction and needs no
database named at either end. :class:`MappingCursor` does that, and
returns the same table as
:class:`rotifer.db.uniprot.mirror.MappingCursor` and
:class:`rotifer.db.uniprot.clickhouse.MappingCursor`, so the three can
be stacked behind one delegator.

UniProt also runs an asynchronous mapping service, which translates
between one pair of databases per job. The module level functions
below drive it and are kept because they are useful on their own, but
the cursor no longer uses them: a job must be submitted and polled
before it returns anything, which costs seconds where reading an entry
costs one, and it cannot be asked for every database at once.

What entries carry is not what ``idmapping.dat`` carries. A curator
records cross-references to Pfam, GO, InterPro and the like, which the
flat file does not hold; the file holds identifiers derived from the
sequence, such as UniRef clusters, UniParc and CRC64 checksums, which
are not cross-references and are not here. Neither source contains the
other.

See Also
--------
rotifer.db.uniprot.mirror.MappingCursor : the flat file
rotifer.db.uniprot.clickhouse.MappingCursor : the same, indexed
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
    Translate identifiers through UniProt's entry cross-references.

    A UniProtKB entry lists every database it is cross-referenced to,
    so one request answers a mapping in either direction: the entries
    are fetched, and the cross-references they carry are read off. That
    is what lets this backend answer an open ended query, which
    UniProt's asynchronous mapping service cannot: the service
    translates between one pair of databases per job, and there are
    ninety-odd of them.

    It is also the faster route. Fetching entries is a plain paged GET,
    while a mapping job must be submitted, polled and only then read,
    which costs seconds before any work is done.

    What it carries is not what ``idmapping.dat`` carries. Entries
    cross-reference the databases a curator recorded -- Pfam, GO,
    InterPro, AlphaFoldDB among them -- while the flat file also holds
    identifiers derived from the sequence, such as UniRef clusters,
    UniParc and CRC64 checksums, which are not cross-references and do
    not appear here. Neither source contains the other, which is why a
    delegator asks both.

    Parameters
    ----------
    **kwargs
        Extra UniProt query parameters.

    Attributes
    ----------
    columns : list of str
        ``['source', 'source_type', 'accession', 'target',
        'target_type']``.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.MappingCursor : the same query, from an indexed copy
    rotifer.db.uniprot.mirror.MappingCursor : the same query, from the flat file

    Examples
    --------
    Every cross-reference of an accession, which needs no database named:

    >>> from rotifer.db.uniprot import webapi                  # doctest: +SKIP
    >>> mc = webapi.MappingCursor()                            # doctest: +SKIP
    >>> mc.fetchall(['Q6GZX4'], source=mc.UNIPROTKB)           # doctest: +SKIP

    Which accession an external identifier belongs to, without saying
    which database it comes from:

    >>> mc.fetchall(['YP_031579.1'], target=mc.UNIPROTKB)      # doctest: +SKIP

    Between two databases:

    >>> mc.fetchall(['AAT09660.1'], source=['EMBL'],           # doctest: +SKIP
    ...             target=['RefSeq','Pfam'])
    """

    # Names a delegator shares with its backends, which describe where
    # to look rather than what to ask, and must not reach UniProt as
    # query parameters.
    _reserved = core.BaseUniProtWebCursor._reserved | frozenset({
        'path', 'engine', 'host', 'port', 'dbname', 'table', 'release',
        'id_type', 'source', 'target', 'polling_interval',
    })

    def __init__(self, *args, **kwargs):
        # Identifiers are named by the caller, not detected from their
        # syntax, so resource detection has nothing to add here.
        kwargs.setdefault('database', 'uniprotkb')
        kwargs.setdefault('probe', False)
        super().__init__(*args, **kwargs)
        self._databases = None
        self._fields = None
        # The batch endpoint takes many accessions at once, and a
        # reverse query is one OR per identifier, so keep both within
        # what a URL and the server will accept.
        self.maxgetitem = 100

    def databases(self):
        """
        Name the databases UniProt cross-references entries to.

        UniProt publishes the list, so it is asked rather than
        guessed, and remembered for the life of the cursor. The names
        are the ones entries themselves report, which is what makes
        them comparable with the databases the local backends hold.

        Returns
        -------
        set of str or None
            None when the list cannot be fetched, so that a network
            failure narrows nothing.
        """
        if not isinstance(self._databases, types.NoneType):
            return self._databases
        try:
            reply = self.get('configure/uniprotkb/allDatabases')
            self._databases = { x['name'] for x in reply.json() if x.get('name') }
        except Exception:
            logger.debug('Could not list UniProt cross-referenced databases', exc_info=1)
            return None
        return self._databases

    def fields(self):
        """
        Map each database to the result field that returns it.

        Asking for named fields instead of whole entries is worth the
        lookup: a batch of five entries is 3 KiB with the fields named
        and 147 KiB without. UniProt labels each field with the
        database it carries, so the correspondence is read from the
        API rather than guessed from the spelling, which would be
        wrong for a sixth of them.

        Returns
        -------
        dict
            Database name to field name. Empty when the list cannot be
            fetched, which simply means whole entries are requested.
        """
        if not isinstance(self._fields, types.NoneType):
            return self._fields
        self._fields = {}
        try:
            reply = self.get('configure/uniprotkb/result-fields')
            for group in reply.json():
                for field in group.get('fields', []) or []:
                    name, label = field.get('name'), field.get('label')
                    if name and label and name.startswith('xref_'):
                        self._fields[label] = name
        except Exception:
            logger.debug('Could not list UniProt result fields', exc_info=1)
        return self._fields

    def _requested_fields(self, target):
        """
        Choose the fields one request should ask for.

        Parameters
        ----------
        target : list of str or None
            Databases wanted, or None for every one of them.

        Returns
        -------
        str or None
            A comma separated field list, or None to ask for whole
            entries, which is necessary when every database is wanted
            or when one of them has no field of its own.
        """
        if isinstance(target, types.NoneType):
            return None
        known = self.fields()
        wanted = [ x for x in target if x != self.UNIPROTKB ]
        if not wanted or not known:
            return None
        if any(x not in known for x in wanted):
            return None
        return ",".join(['accession'] + [ known[x] for x in wanted ])

    def fetcher(self, accession, source=None, target=None, *args, **kwargs):
        """
        Fetch the entries that can answer the query.

        Identifiers that are UniProtKB accessions are taken from the
        batch endpoint; everything else is looked up through the
        cross-reference index, which finds an entry from an identifier
        without being told which database it belongs to.

        Parameters
        ----------
        accession : iterable of str
            Identifiers to translate.
        source, target : str or list of str, optional
            The two ends of the mapping. Either may be omitted.

        Returns
        -------
        list of dict
            The entries found.
        """
        targets = sorted(self.parse_ids(accession))
        if not targets:
            return []
        source = self.parse_databases(source)
        target = self.parse_databases(target)
        fields = self._requested_fields(target)

        entries = []
        if isinstance(source, types.NoneType) or self.UNIPROTKB in source:
            params = self.query_parameters(accessions=",".join(targets), format='json', fields=fields)
            try:
                entries.extend(self.get('uniprotkb/accessions', **params).json().get('results', []) or [])
            except Exception as error:
                # An identifier that is not an accession makes the whole
                # batch a bad request, which says nothing about the
                # others: the reverse lookup below still covers them.
                logger.debug(f'Batch accession lookup failed: {error}')

        others = None if isinstance(source, types.NoneType) else [ x for x in source if x != self.UNIPROTKB ]
        if isinstance(source, types.NoneType) or others:
            query = " OR ".join([ f'xref:"{x}"' for x in targets ])
            params = self.query_parameters(query=query, format='json',
                                           fields=fields, size=self.page_size)
            for reply in self.pages('uniprotkb/search', **params):
                entries.extend(reply.json().get('results', []) or [])
        return entries

    #: Entries per page of a reverse lookup.
    page_size = 100

    def parser(self, stream, accession, source=None, target=None, *args, **kwargs):
        """
        Read the mappings off the entries fetched.

        Each entry carries the accession and the cross-references it
        has. The identifiers asked about are matched against both, so
        an entry found through one of its cross-references reports
        that cross-reference as the source rather than the accession.

        Parameters
        ----------
        stream : list of dict
            Entries returned by :meth:`fetcher`.
        accession : iterable of str
            Identifiers that were requested.
        source, target : str or list of str, optional
            The two ends of the mapping.

        Returns
        -------
        pandas.DataFrame
            The columns listed in
            :attr:`~rotifer.db.methods.MappingCursor.columns`.
        """
        if isinstance(stream, types.NoneType) or not stream:
            return self.empty()
        queried = self.parse_ids(accession)
        source = self.parse_databases(source)
        target = self.parse_databases(target)

        rows = []
        for entry in stream:
            accession_value = entry.get('primaryAccession')
            if not accession_value:
                continue
            xrefs = [ (x.get('database'), x.get('id'))
                      for x in entry.get('uniProtKBCrossReferences', []) or []
                      if x.get('database') and x.get('id') ]

            # Which of the identifiers asked about this entry answers,
            # and by which name.
            sources = []
            if isinstance(source, types.NoneType) or self.UNIPROTKB in source:
                names = {accession_value} | set(entry.get('secondaryAccessions', []) or [])
                for name in sorted(names & queried):
                    sources.append((name, self.UNIPROTKB))
            for database, value in xrefs:
                if value in queried and (isinstance(source, types.NoneType) or database in source):
                    sources.append((value, database))
            if not sources:
                continue

            # What this entry offers at the other end.
            wanted = []
            if not isinstance(target, types.NoneType) and self.UNIPROTKB in target:
                wanted.append((accession_value, self.UNIPROTKB))
            for database, value in xrefs:
                if isinstance(target, types.NoneType) or database in target:
                    wanted.append((value, database))

            for source_value, source_type in sources:
                for target_value, target_type in wanted:
                    rows.append({
                        'source': str(source_value), 'source_type': str(source_type),
                        'accession': str(accession_value),
                        'target': str(target_value), 'target_type': str(target_type),
                    })
        if not rows:
            return self.empty()
        return pd.DataFrame(rows, columns=self.columns).drop_duplicates().reset_index(drop=True)


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
