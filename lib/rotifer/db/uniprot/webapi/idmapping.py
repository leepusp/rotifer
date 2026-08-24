"""
Map database identifiers using UniProt's ID mapping REST service.

The lower-level functions in this module (``submit_id_mapping``
through ``get_id_mapping_results_stream``) are adapted from the
client script UniProt publishes alongside its ID mapping REST API
documentation and implement the submit/poll/fetch job workflow used
by that service. :func:`genbank_to_uniprot`, :func:`AF_link` and
:func:`af_to_seq` are rotifer-specific wrappers built on top of them
to cross-reference GenBank/RefSeq proteins with AlphaFold and PDB
structures through UniProt.

Note
----
This module makes live network requests; none of its examples are
executed as doctests.
"""

import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry


POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    """
    Raise for HTTP errors, printing the JSON error body first.

    Parameters
    ----------
    response : requests.Response
        Response to check.

    Raises
    ------
    requests.HTTPError
    """
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids):
    """
    Submit an identifier mapping job to UniProt.

    Parameters
    ----------
    from_db : str
        Source database name, as accepted by UniProt's ID mapping
        API (for example ``EMBL-GenBank-DDBJ_CDS``).
    to_db : str
        Target database name (for example ``UniProtKB``).
    ids : list of str
        Identifiers to map.

    Returns
    -------
    str
        Identifier of the submitted job.
    """
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    """
    Extract the URL of the next results page from a Link header.

    Parameters
    ----------
    headers : mapping
        Response headers.

    Returns
    -------
    str or None
        The next page's URL, or None if there is no ``Link`` header
        with a ``rel="next"`` entry.
    """
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    """
    Block until a submitted ID mapping job finishes.

    Polls the job status every ``POLLING_INTERVAL`` seconds.

    Parameters
    ----------
    job_id : str
        Identifier returned by :func:`submit_id_mapping`.

    Returns
    -------
    bool
        True once the job is finished and has results or failed
        identifiers to report.

    Raises
    ------
    Exception
        If the job reaches a status other than ``RUNNING`` that is
        not a normal completion.
    """
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    """
    Iterate over the remaining pages of a paginated results response.

    Parameters
    ----------
    batch_response : requests.Response
        Response for the first page already fetched.
    file_format : str
        One of ``json``, ``tsv``, ``xlsx`` or ``xml``.
    compressed : bool
        Whether the response bodies are gzip-compressed.

    Yields
    ------
    object
        Decoded contents of each subsequent page, as returned by
        :func:`decode_results`.
    """
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    """
    Append one results page to the accumulated results.

    Parameters
    ----------
    all_results : dict or list
        Results accumulated so far.
    batch_results : dict or list
        Results of the page being merged in.
    file_format : str
        One of ``json``, ``tsv``, ``xlsx`` or ``xml``. For ``json``
        the ``results`` and ``failedIds`` keys are concatenated in
        place; for other formats a new combined list is returned.

    Returns
    -------
    dict or list
        The updated accumulated results.
    """
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
    """
    Fetch the results URL of a finished ID mapping job.

    Parameters
    ----------
    job_id : str
        Identifier returned by :func:`submit_id_mapping`.

    Returns
    -------
    str
        URL to pass to :func:`get_id_mapping_results_search` or
        :func:`get_id_mapping_results_stream`.
    """
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    """
    Decode a results response body according to its file format.

    Parameters
    ----------
    response : requests.Response
        Response carrying the results.
    file_format : str
        One of ``json``, ``tsv``, ``xlsx`` or ``xml``.
    compressed : bool
        Whether the response body is gzip-compressed.

    Returns
    -------
    object
        A dict for ``json``, a list of lines for ``tsv``, a
        single-element list of bytes for ``xlsx``, a single-element
        list of str for ``xml``, or the raw response text otherwise.
    """
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
    """
    Extract the XML namespace URI from an element's tag.

    Parameters
    ----------
    element : xml.etree.ElementTree.Element
        Element whose tag carries the namespace, in ``{uri}tag``
        form.

    Returns
    -------
    str
        The namespace URI, or an empty string if the tag has none.
    """
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    """
    Merge the ``entry`` elements of several UniProt XML pages.

    Parameters
    ----------
    xml_results : list of str
        XML document text for each results page.

    Returns
    -------
    bytes
        Serialized XML document combining every page's entries under
        the root element of the first page.
    """
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    """
    Print a one-line progress message for a fetched results page.

    Parameters
    ----------
    batch_index : int
        Zero-based index of the page just fetched.
    size : int
        Page size, in results.
    total : int
        Total number of results expected.
    """
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    """
    Fetch all pages of a finished ID mapping job through the search endpoint.

    Parameters
    ----------
    url : str
        Results URL returned by :func:`get_id_mapping_results_link`.

    Returns
    -------
    object
        The combined results, decoded as described in
        :func:`decode_results`; XML results are merged into a single
        document by :func:`merge_xml_results`.
    """
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
    """
    Fetch a finished ID mapping job through the streaming endpoint.

    Unlike :func:`get_id_mapping_results_search`, streaming returns
    the whole result set in one request and is preferred for large
    jobs.

    Parameters
    ----------
    url : str
        Results URL returned by :func:`get_id_mapping_results_link`.

    Returns
    -------
    object
        The decoded results, as described in :func:`decode_results`.
    """
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


import csv
import pandas as pd

def get_data_frame_from_tsv_results(tsv_results):
    """
    Parse TSV-formatted ID mapping results into a dataframe.

    Parameters
    ----------
    tsv_results : iterable of str
        TSV lines, such as the ``tsv``-format output of
        :func:`decode_results` (the header row is used as column
        names).

    Returns
    -------
    pandas.DataFrame
    """
    reader = csv.DictReader(tsv_results, delimiter="\t", quotechar='"')
    t = pd.DataFrame(list(reader))
    return t

def genbank_to_uniprot(from_db="EMBL-GenBank-DDBJ_CDS", to_db="UniProtKB", ids=["BAE76179.1"]):
    """
    Map protein accessions to UniProt entries and their PDB/AlphaFold cross-references.

    Submits an ID mapping job, waits for it to finish and returns
    the ``accession``, ``xref_pdb`` and ``xref_alphafolddb`` fields
    for each mapped identifier.

    Parameters
    ----------
    from_db : str, default 'EMBL-GenBank-DDBJ_CDS'
        Source database name, as accepted by UniProt's ID mapping API.
    to_db : str, default 'UniProtKB'
        Target database name.
    ids : list of str
        Identifiers to map.

    Returns
    -------
    pandas.DataFrame
        One row per mapped identifier.
    """
    job_id = submit_id_mapping(
        from_db=from_db,
        to_db=to_db,
        ids=ids)
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_stream(link+"?compressed=true&fields=accession%2Cxref_pdb%2Cxref_alphafolddb%2C&format=tsv")
    r = get_data_frame_from_tsv_results(results)
    #if r.AlphaFold:
    #    r['urlAF'] = [ "https://alphafold.ebi.ac.uk/files/AF-" + x + "-F1-model_v4.pdb" for x in r.loc[r.AlphaFoldDB.str.split(';', expand=True)[0] == r.Entry].Entry ]
    #if r.PDB:
    #    r['urlpdb'] = [ "https://files.rcsb.org/download/" + x + ".pdb" for x in r.PDB.str.split(';', expand=True)[0] ]
    return(r)

def AF_link(id_list=None):
    """
    Build AlphaFold and PDB download URLs for a list of protein accessions.

    Parameters
    ----------
    id_list : list of str, optional
        Protein accessions, passed to :func:`genbank_to_uniprot`.

    Returns
    -------
    pandas.DataFrame
        Rows with an AlphaFold structure, with two added columns:
        ``urlAF`` (AlphaFold model URL) and ``urlPDB`` (RCSB PDB
        download URL, when a PDB cross-reference is present).
    """
    r = genbank_to_uniprot(ids=id_list)
    r['urlAF'] = "https://alphafold.ebi.ac.uk/files/AF-" + r['AlphaFoldDB'].str.split(';', expand=True)[0] + "-F1-model_v4.pdb"
    r.loc[r.urlAF == "https://alphafold.ebi.ac.uk/files/AF--F1-model_v4.pdb", 'urlAF'] = None
    r['urlPDB'] = "https://files.rcsb.org/download/" + r['PDB'].str.split(';', expand=True)[0 ]+ ".pdb"
    r.loc[r.urlPDB == "https://files.rcsb.org/download/.pdb", 'urlPDB'] = None
    r = r[r.urlAF.notnull()].reset_index()
    return r

def af_to_seq(seqobj):
    """
    Attach AlphaFold structures to a rotifer sequence object.

    Parameters
    ----------
    seqobj : rotifer.devel.beta.sequence.sequence
        Sequence object whose ``df.id`` column lists protein
        accessions to map through :func:`AF_link`.

    Returns
    -------
    rotifer.devel.beta.sequence.sequence
        A copy of ``seqobj`` with one PDB structure added per mapped
        accession, via ``add_pdb``.
    """
    r = seqobj.copy()
    u = AF_link(r.df.id.to_list())
    for x in range(0,len(u.From)):
        r = r.add_pdb(pdb_id=u.From[x], pdb_file=u.urlAF[x])
    return r


'''
job_id = submit_id_mapping(
    from_db="EMBL-GenBank-DDBJ_CDS",
    to_db="UniProtKB",
    ids=["BAE76179.1","AAC73502.1"])
if check_id_mapping_results_ready(job_id):
    link = get_id_mapping_results_link(job_id)
    results = get_id_mapping_results_search(link)
print(results)

def get_data_frame_from_tsv_results(tsv_results):
    reader = csv.DictReader(tsv_results, delimiter="\t", quotechar='"')
    return pd.DataFrame(list(reader))
job_id = submit_id_mapping(
    from_db="EMBL-GenBank-DDBJ_CDS",
    to_db="UniProtKB",
    ids=["BAE76179.1","AAC73502.1"])
if check_id_mapping_results_ready(job_id):
    link = get_id_mapping_results_link(job_id)
    results = get_id_mapping_results_stream(link+"?compressed=true&fields=accession%2Cxref_alphafolddb%2C&format=tsv")
r = get_data_frame_from_tsv_results(results)
r['urlAF'] = [ "https://alphafold.ebi.ac.uk/files/AF-" + x + "-F1-model_v4.pdb" for x in r.loc[r.AlphaFoldDB.str.split(';', expand=True)[0] == r.Entry].Entry]

'''
