#!/usr/bin/env python3
"""Tests for rotifer.db.uniprot.webapi.

Every test here runs offline. The cursors exist to talk to UniProt, but the
parts that are easy to get wrong are the parts that decide *what* to ask and
*how to read the answer*, and both can be checked without a network:

  * detection, which is pure string matching once `probe=False` rules out the
    API fallback
  * request planning, i.e. which endpoint and parameters each kind of
    identifier turns into
  * parsing, fed with captured payloads rather than live ones

The routing tests carry two regressions worth stating outright:

  * a proteome has no sequence of its own but expands into UniProtKB entries,
    so it must still produce a request. Gating that on `Resource.sequences`
    (False for proteomes) made the branch dead code and returned zero records
    *silently* -- which a test that only checks for exceptions will pass.
  * a taxonomy id in a sequence query has no sequences by any route, so it
    must be registered as missing rather than quietly dropped. Those two
    outcomes are opposite and a single wrong predicate produced both.

Run under pytest, or standalone: python test/db/uniprot/test_webapi.py
"""

import io
import pickle

import pandas as pd
import pytest
from Bio import SeqIO

from rotifer.db.uniprot import webapi
from rotifer.db.uniprot.webapi import core, entries, search
from rotifer.db.uniprot.webapi import idmapping as rdui


# --------------------------------------------------------------- detection

#: identifier -> resource it must be recognised as.
KNOWN = {
    'P00750': 'uniprotkb',            # classic accession
    'Q6GZX4': 'uniprotkb',
    'A0A0A0MS99': 'uniprotkb',        # the 10 character form
    'P00750-2': 'uniprotkb',          # isoform suffix
    'TPA_HUMAN': 'uniprotkb',         # entry name rather than accession
    'UPI0000000001': 'uniparc',
    'UniRef50_P00750': 'uniref',
    'UniRef90_P00750': 'uniref',
    'UniRef100_P00750': 'uniref',
    'UP000005640': 'proteomes',
    '9606': 'taxonomy',
}


@pytest.mark.parametrize('accession,resource', sorted(KNOWN.items()))
def test_detect_resource(accession, resource):
    assert core.detect_resource(accession, probe=False) == resource


def test_uniref_is_tested_before_uniprotkb():
    """A UniRef name ends with its representative's accession, so pattern order
    is what keeps UniRef50_P00750 from being read as the accession P00750."""
    assert core.detect_resource('UniRef50_P00750', probe=False) == 'uniref'
    assert core.detect_resource('P00750', probe=False) == 'uniprotkb'


def test_proteome_and_uniparc_do_not_collide():
    """Both start with UP; digits and an I are what separate them."""
    assert core.detect_resource('UP000005640', probe=False) == 'proteomes'
    assert core.detect_resource('UPI0000000001', probe=False) == 'uniparc'


@pytest.mark.parametrize('accession', ['', '   ', 'not an id', '???', 'P0'])
def test_unrecognised_identifiers_are_not_guessed(accession):
    """Offline, an identifier no pattern matches has no answer, and inventing
    one would send a request that cannot succeed."""
    assert core.detect_resource(accession, probe=False) is None


def test_group_by_resource_splits_a_mixed_query():
    groups, unknown = core.group_by_resource(
        ['P00750', 'Q6GZX4', 'UPI0000000001', 'UP000005640', '???'], probe=False)
    assert groups == {
        'uniprotkb': ['P00750', 'Q6GZX4'],
        'uniparc': ['UPI0000000001'],
        'proteomes': ['UP000005640'],
    }
    assert unknown == {'???'}


def test_naming_a_resource_skips_detection():
    """An explicit database is an assertion by the caller, so identifiers that
    no pattern would recognise still go to that resource."""
    groups, unknown = core.group_by_resource(['whatever'], database='uniref', probe=False)
    assert groups == {'uniref': ['whatever']} and not unknown


def test_unknown_resource_is_rejected():
    with pytest.raises(ValueError, match='Unknown UniProt resource'):
        core.group_by_resource(['P00750'], database='swissprot', probe=False)


# ------------------------------------------------------------ request plan

def plan(cursor, accessions):
    """Request plan as {resource: [params, ...]}, ignoring paths."""
    out = {}
    for chunk in cursor._chunks(cursor.parse_ids(accessions)):
        out.setdefault(chunk['resource'], []).append(chunk)
    return out


def test_batch_endpoints_are_used_where_they_exist():
    """UniProt publishes endpoints that take many identifiers at once; using
    them is the difference between one request and one per accession."""
    chunks = plan(webapi.FastaCursor(progress=False, probe=False), ['P00750', 'P12345'])
    assert list(chunks) == ['uniprotkb']
    chunk, = chunks['uniprotkb']
    assert chunk['path'] == 'uniprotkb/accessions'
    assert chunk['params']['accessions'] == 'P00750,P12345'


def test_resources_without_a_batch_endpoint_use_one_search():
    """UniRef has no batch endpoint, so the fallback must still be a single
    request rather than one per identifier."""
    chunks = plan(webapi.FastaCursor(progress=False, probe=False),
                  ['UniRef50_P00750', 'UniRef90_P00750'])
    chunk, = chunks['uniref']
    assert chunk['path'] == 'uniref/stream'
    assert chunk['params']['query'] == 'id:UniRef50_P00750 OR id:UniRef90_P00750'


def test_a_proteome_expands_into_a_uniprotkb_query():
    """Regression: a proteome has no sequences of its own, but it is not a
    dead end -- it expands into the UniProtKB entries it contains."""
    chunks = plan(webapi.SequenceCursor(progress=False, probe=False), ['UP000000625'])
    chunk, = chunks['proteomes']
    assert chunk['path'] == 'uniprotkb/stream'
    assert chunk['params']['query'] == 'proteome:UP000000625'
    assert chunk['proteome'] == 'UP000000625'


def test_one_request_per_proteome():
    """Records carry the proteome they came from, which only works if each
    proteome was asked for on its own."""
    chunks = plan(webapi.SequenceCursor(progress=False, probe=False),
                  ['UP000000625', 'UP000005640'])
    assert sorted(c['proteome'] for c in chunks['proteomes']) == ['UP000000625', 'UP000005640']


def test_taxonomy_in_a_sequence_query_is_reported_missing():
    """Regression, and the mirror image of the proteome case: taxonomy has no
    sequences by any route, so it must be registered rather than dropped."""
    cursor = webapi.FastaCursor(progress=False, probe=False)
    chunks = plan(cursor, ['P00750', '9606'])
    assert 'taxonomy' not in chunks
    assert '9606' in cursor.missing_ids()


def test_unrecognised_identifiers_are_registered():
    cursor = webapi.FastaCursor(progress=False, probe=False)
    plan(cursor, ['P00750', 'NOT_AN_ID_AT_ALL'])
    assert 'NOT_AN_ID_AT_ALL' in cursor.missing_ids()


def test_sequence_and_fasta_cursors_differ_only_in_format():
    """SequenceCursor asks for the flat file, which carries the annotation
    FASTA leaves out; both ask the same endpoint."""
    fasta, = plan(webapi.FastaCursor(progress=False, probe=False), ['P00750'])['uniprotkb']
    full, = plan(webapi.SequenceCursor(progress=False, probe=False), ['P00750'])['uniprotkb']
    assert fasta['params']['format'] == 'fasta' and fasta['format'] == 'fasta'
    assert full['params']['format'] == 'txt' and full['format'] == 'swiss'
    assert fasta['path'] == full['path']


# -------------------------------------------------------- extra parameters

def test_unknown_keywords_become_query_parameters():
    """UniProt accepts far more options than these classes model, so anything
    unrecognised has to reach the API rather than raise."""
    cursor = webapi.FastaCursor(progress=False, probe=False, includeIsoform=True,
                                fields='accession,id')
    assert cursor.parameters == {'includeIsoform': True, 'fields': 'accession,id'}
    chunk, = plan(cursor, ['P00750'])['uniprotkb']
    assert chunk['params']['includeIsoform'] is True
    assert chunk['params']['fields'] == 'accession,id'


def test_cursor_keywords_are_not_sent_to_uniprot():
    cursor = webapi.FastaCursor(progress=False, probe=False, threads=2, batch_size=7)
    assert cursor.parameters == {}


def test_the_calls_own_parameters_win():
    """A cursor sets the format it can parse; a caller must not be able to
    replace it by accident, but may still add parameters of their own."""
    cursor = webapi.FastaCursor(progress=False, probe=False, format='xml', fields='accession')
    merged = cursor.query_parameters(format='fasta')
    assert merged == {'format': 'fasta', 'fields': 'accession'}


# ------------------------------------------------------------------ paging

@pytest.mark.parametrize('headers,expected', [
    ({}, None),
    ({'Link': '<https://rest.uniprot.org/x?cursor=abc>; rel="next"'},
     'https://rest.uniprot.org/x?cursor=abc'),
    ({'Link': '<https://rest.uniprot.org/x>; rel="prev"'}, None),
])
def test_next_link(headers, expected):
    """UniProt pages with an opaque cursor handed back in a Link header, which
    is what it recommends following instead of computing offsets."""
    assert core._next_link(headers) == expected


# ----------------------------------------------------------------- parsing

FASTA = (">sp|P00750|TPA_HUMAN Tissue-type plasminogen activator\n"
         "MDAMKRGLCCVLLLCGAVFVSPSQEIHARFRRGAR\n"
         ">UPI0000000001 status=active\n"
         "MAAAAAAAAAAK\n")

SWISS = """ID   TEST_HUMAN              Reviewed;          12 AA.
AC   P00001; P00002;
DE   RecName: Full=Test protein;
OS   Homo sapiens (Human).
OC   Eukaryota; Metazoa.
OX   NCBI_TaxID=9606;
FT   CHAIN           1..12
FT                   /note="Test chain"
SQ   SEQUENCE   12 AA;  1234 MW;  0000000000000000 CRC64;
     MAAAAAAAAAAK
//
"""


def test_parser_reads_each_chunk_in_its_own_format():
    """A mixed query returns several formats at once, so the format travels
    with the reply rather than being a property of the cursor."""
    cursor = webapi.SequenceCursor(progress=False, probe=False)
    records = cursor.parser([
        {'resource': 'uniprotkb', 'format': 'fasta', 'text': FASTA, 'proteome': None},
        {'resource': 'uniprotkb', 'format': 'swiss', 'text': SWISS, 'proteome': None},
    ], [])
    assert [r.id for r in records] == ['sp|P00750|TPA_HUMAN', 'UPI0000000001', 'P00001']


def test_parser_tags_records_with_their_proteome():
    """The identifier that was asked for is nowhere in the reply, so it has to
    be attached, or the proteome counts as missing."""
    cursor = webapi.SequenceCursor(progress=False, probe=False)
    records = cursor.parser(
        [{'resource': 'proteomes', 'format': 'swiss', 'text': SWISS, 'proteome': 'UP000000625'}], [])
    assert [r.annotations['proteome'] for r in records] == ['UP000000625']
    assert 'UP000000625' in cursor.getids(records)


@pytest.mark.parametrize('stream', [None, [], [{'format': 'fasta', 'text': ''}]])
def test_parser_survives_an_empty_reply(stream):
    assert webapi.FastaCursor(progress=False, probe=False).parser(stream, []) == []


def test_record_ids_cover_every_name_an_entry_answers_to():
    """FASTA names an entry db|accession|entry_name while the flat file names
    it by accession alone. Reporting only one of them would count an entry as
    missing whenever it was asked for by the other."""
    records = list(SeqIO.parse(io.StringIO(SWISS), 'swiss'))
    assert entries._record_ids(records[0]) >= {'P00001', 'P00002', 'TEST_HUMAN'}
    fasta = list(SeqIO.parse(io.StringIO(FASTA), 'fasta'))
    assert entries._record_ids(fasta[0]) >= {'sp|P00750|TPA_HUMAN', 'P00750', 'TPA_HUMAN', 'sp'}


# ----------------------------------------------------------- id mapping

def test_id_mapping_returns_the_shared_three_column_table():
    """The mirror and clickhouse cursors return this table, so this one must
    too, or the three cannot be stacked behind one delegator."""
    cursor = webapi.IdMappingCursor(to_db='RefSeq_Protein', progress=False)
    frame = cursor.parser([{'from': 'P00750', 'to': 'NP_000921.1'}], [])
    assert list(frame.columns) == ['accession', 'id_type', 'id']
    assert frame.to_dict('records') == [
        {'accession': 'P00750', 'id_type': 'RefSeq_Protein', 'id': 'NP_000921.1'}]


def test_id_mapping_unwraps_uniprotkb_entries():
    """Mapping *to* UniProtKB returns whole entries where other targets return
    a bare identifier."""
    cursor = webapi.IdMappingCursor(to_db='UniProtKB', progress=False)
    frame = cursor.parser([{'from': 'BAE76179.1', 'to': {'primaryAccession': 'P00750'}}], [])
    assert frame['id'].tolist() == ['P00750']


def test_id_mapping_reports_what_was_asked_for():
    """getids answers 'which of my queries came back', so it reads the column
    the service echoes the query into, not the one holding the answers."""
    cursor = webapi.IdMappingCursor(progress=False)
    frame = pd.DataFrame([{'accession': 'P00750', 'id_type': 'RefSeq_Protein', 'id': 'NP_000921.1'}])
    assert cursor.getids(frame) == {'P00750'}


@pytest.mark.parametrize('rows', [[], [{'from': 'P00750'}], [{'to': 'X'}]])
def test_id_mapping_empty_results_keep_the_columns(rows):
    cursor = webapi.IdMappingCursor(progress=False)
    frame = cursor.parser(rows, [])
    assert frame.empty and list(frame.columns) == ['accession', 'id_type', 'id']


# --------------------------------------------------------------- tabular

def test_table_cursor_flattens_and_labels_rows():
    cursor = webapi.SearchCursor(progress=False, probe=False)
    frame = cursor.parser([
        {'resource': 'uniprotkb', 'rows': [{'primaryAccession': 'P00750',
                                            'organism': {'taxonId': 9606}}]},
        {'resource': 'proteomes', 'rows': [{'id': 'UP000005640'}]},
    ], [])
    assert frame['organism.taxonId'].tolist()[0] == 9606
    assert sorted(frame['resource']) == ['proteomes', 'uniprotkb']
    assert cursor.getids(frame) == {'P00750', 'UP000005640'}


def test_search_needs_a_named_resource():
    """A raw query carries no identifier, so there is nothing to detect a
    resource from and guessing one would query the wrong index."""
    with pytest.raises(ValueError, match='needs a resource'):
        webapi.SearchCursor(progress=False).search('reviewed:true')


def test_proteome_and_taxonomy_cursors_pin_their_resource():
    assert webapi.ProteomeCursor(progress=False).database == 'proteomes'
    assert webapi.TaxonomyCursor(progress=False).database == 'taxonomy'


# ------------------------------------------------------------- processes

def test_cursors_survive_being_sent_to_a_worker():
    """fetchone() ships the cursor to worker processes, and a requests.Session
    cannot be pickled, so it must not travel with it."""
    cursor = webapi.FastaCursor(progress=False, probe=False, fields='accession')
    assert cursor.session is not None          # force one into existence
    revived = pickle.loads(pickle.dumps(cursor))
    assert revived._session is None
    assert revived.parameters == {'fields': 'accession'}
    assert revived.session is not None         # and each process builds its own


if __name__ == '__main__':
    import sys
    failures = 0
    tests = [(name, obj) for name, obj in sorted(globals().items())
             if name.startswith('test_') and callable(obj)
             and not hasattr(obj, 'pytestmark')]
    for name, test in tests:
        try:
            test()
            ok = True
        except Exception as exc:
            ok, name = False, f'{name} ({type(exc).__name__}: {exc})'
        failures += not ok
        print('  %-58s %s' % (name.split(' (')[0], 'ok' if ok else 'FAILED'))
        if not ok:
            print('      %s' % name.split(' (', 1)[1][:-1])
    print('ALL PASS' if not failures else '%d FAILURES' % failures)
    sys.exit(1 if failures else 0)
