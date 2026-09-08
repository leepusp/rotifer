#!/usr/bin/env python3
"""Tests for the UniProt sequence and taxonomy delegators.

These two are thin over their backend today -- only the web service has
sequences or taxonomy -- so what is worth testing is not the delegation but
the two things the delegator itself decides: which identifiers a result
counts as answering for, and what shape the answer comes back in.

Both matter more than they look. An entry asked for by accession comes back
from FASTA as `sp|P00750|TPA_HUMAN`, and a delegator that does not recognise
that as P00750 reports a found entry as missing and asks the next backend for
it again. And a taxonomy frame is only interchangeable with the NCBI one if
the six shared columns are built the same way, down to which ancestors are
left out; the expected values here were taken from ete3's own output for the
same taxa, so a drift in either direction fails.

Nothing here touches the network: backends are injected, and the reference
values are literals.

Run under pytest, or standalone: python test/db/uniprot/test_delegators.py
"""

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import rotifer.db.core
from rotifer.db import uniprot


class Backend(rotifer.db.core.BaseCursor):
    """A backend holding fixed results, which records what it was asked."""

    def __init__(self, results=(), **kwargs):
        super().__init__(progress=False)
        self._results = list(results)
        self.asked = []
        self.built_with = kwargs

    def getids(self, obj, *args, **kwargs):
        return set()

    def fetchone(self, accessions, *args, **kwargs):
        wanted = self.parse_ids(accessions)
        self.asked.append(wanted)
        for result, ids in self._results:
            if wanted.intersection(ids):
                yield result
        self.update_missing(wanted - self.answers(), 'Not found.', retry=False)

    def answers(self):
        found = set()
        for _, ids in self._results:
            found.update(ids)
        return found


def sequences(*backends):
    """A SequenceCursor over the given backends, in order."""

    class Cursor(uniprot.SequenceCursor):
        def reset_cursors(self):
            self.cursors = dict(getattr(self, '_backends', {}))

    cursor = Cursor.__new__(Cursor)
    cursor._backends = { f'b{i}': b for i, b in enumerate(backends) }
    uniprot.SequenceCursor.__init__(
        cursor, readers=[f'b{i}' for i in range(len(backends))], progress=False)
    return cursor


def taxonomy(*backends):
    """A TaxonomyCursor over the given backends, in order."""

    class Cursor(uniprot.TaxonomyCursor):
        def reset_cursors(self):
            self.cursors = dict(getattr(self, '_backends', {}))

    cursor = Cursor.__new__(Cursor)
    cursor._backends = { f'b{i}': b for i, b in enumerate(backends) }
    uniprot.TaxonomyCursor.__init__(
        cursor, readers=[f'b{i}' for i in range(len(backends))], progress=False)
    return cursor


def fasta(header, residues='MKV'):
    return SeqRecord(Seq(residues), id=header, name=header, description='')


# ================================================================ sequences

def test_a_fasta_header_answers_for_the_accession_inside_it():
    """The one that matters. UniProt FASTA names an entry
    `sp|P00750|TPA_HUMAN`, so a delegator reading only `.id` would report an
    entry it just received as missing, and ask the next backend for it."""
    cursor = sequences()
    assert 'P00750' in cursor.getids(fasta('sp|P00750|TPA_HUMAN'))


def test_the_entry_name_answers_too():
    assert 'TPA_HUMAN' in sequences().getids(fasta('sp|P00750|TPA_HUMAN'))


def test_a_secondary_accession_answers_for_the_entry():
    """A flat file lists the accessions an entry has had. Asking by an old
    one must not look like an entry that does not exist."""
    record = fasta('P00750')
    record.annotations = {'accessions': ['P00750', 'Q15837']}
    assert sequences().getids(record) >= {'P00750', 'Q15837'}


def test_records_can_arrive_in_a_list():
    found = sequences().getids([fasta('P00750'), fasta('UPI0000000001')])
    assert found >= {'P00750', 'UPI0000000001'}


def test_nothing_answers_for_nothing():
    assert sequences().getids(None) == set()


def test_a_backend_is_only_asked_for_what_is_still_missing():
    first = Backend([(fasta('P00750'), {'P00750'})])
    second = Backend([(fasta('UPI0000000001'), {'UPI0000000001'})])
    cursor = sequences(first, second)
    cursor.fetchall(['P00750', 'UPI0000000001'])
    assert second.asked == [{'UPI0000000001'}]


def test_an_identifier_no_backend_has_is_reported_missing():
    cursor = sequences(Backend([(fasta('P00750'), {'P00750'})]))
    cursor.fetchall(['P00750', 'P99999999'])
    assert cursor.missing_ids() == {'P99999999'}


def test_the_resource_to_ask_reaches_the_backend():
    """database='auto' is what lets one call mix UniProtKB, UniParc and
    UniRef identifiers, so it has to arrive where detection happens."""
    cursor = uniprot.SequenceCursor(readers=[], progress=False)
    assert cursor.database == 'auto'
    named = uniprot.SequenceCursor(readers=[], database='uniref', progress=False)
    assert named.database == 'uniref'


def test_the_query_is_reported_on_but_a_lookup_is_not(capsys):
    """fetchall is a job; mc[x] is a lookup. Same rule as the mapping
    cursors, and progress is shared with the backends, so clearing it here
    silences theirs too."""
    cursor = sequences(Backend([(fasta('P00750'), {'P00750'})]))
    cursor.progress = True
    cursor['P00750']
    assert capsys.readouterr().err == ''
    assert cursor.progress is True
    cursor.fetchall(['P00750'])
    assert 'uniprot' in capsys.readouterr().err


# ================================================================ taxonomy

#: One taxon as UniProt returns it, cut down to the fields used here.
HUMAN = {
    'taxonId': 9606,
    'scientificName': 'Homo sapiens',
    'mnemonic': 'HUMAN',
    'lineage': [
        {'scientificName': 'cellular organisms', 'taxonId': 131567, 'hidden': True},
        {'scientificName': 'Eukaryota', 'taxonId': 2759, 'rank': 'domain'},
        {'scientificName': 'Opisthokonta', 'taxonId': 33154, 'hidden': True},
        {'scientificName': 'Metazoa', 'taxonId': 33208, 'rank': 'kingdom'},
    ],
}

#: What ete3 produces for the same taxon, which is what the shared columns
#: must reproduce. Trimmed to the ancestors above.
ETE3_CLASSIFICATION = 'Eukaryota; Opisthokonta; Metazoa; Homo sapiens'


def uniprot_frame(*rows):
    return pd.DataFrame(list(rows))


def test_the_shared_columns_are_all_there():
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert set(taxonomy().taxcols) <= set(frame.columns)


def test_the_classification_matches_the_one_ete3_builds():
    """This is what makes the two taxonomy cursors interchangeable, so it is
    checked against ete3's own output rather than against itself."""
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert frame.classification.iloc[0] == ETE3_CLASSIFICATION


def test_the_nodes_ncbi_leaves_out_are_left_out():
    """`cellular organisms` and `root` are in UniProt's lineage and not in
    NCBI's, so keeping them would shift every column that follows."""
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert 'cellular organisms' not in frame.classification.iloc[0]


def test_a_hidden_clade_is_kept():
    """UniProt marks Opisthokonta hidden, but ete3 lists it, so filtering on
    that flag would be the obvious wrong way to do the line above."""
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert 'Opisthokonta' in frame.classification.iloc[0]


def test_the_taxon_itself_closes_its_own_classification():
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert frame.classification.iloc[0].endswith('Homo sapiens')


def test_the_identifier_and_name_are_renamed_not_moved():
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert frame.taxid.iloc[0] == '9606'
    assert frame.organism.iloc[0] == 'Homo sapiens'
    assert frame.taxonId.iloc[0] == 9606         # and UniProt's own kept


def test_the_superkingdom_is_the_first_clade():
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert frame.superkingdom.iloc[0] == 'Eukaryota'


def test_the_ancestors_are_kept_whole():
    """The flattened string carries names only. UniProt gives a taxid and a
    rank per ancestor, which would be lost if `lineage` were overwritten in
    place."""
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert frame.lineage_taxa.iloc[0][1]['taxonId'] == 2759


def test_what_uniprot_adds_survives():
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert frame.mnemonic.iloc[0] == 'HUMAN'


def test_the_shared_columns_come_first():
    frame = taxonomy().standardize(uniprot_frame(HUMAN))
    assert frame.columns[:6].tolist() == taxonomy().taxcols


def test_a_backend_that_already_speaks_this_way_is_left_alone():
    """ete3 produces these columns itself. Rebuilding them would replace its
    lineage string with one derived from a lineage column it does not have."""
    ete3 = pd.DataFrame([{
        'taxid': '9606', 'organism': 'Homo sapiens', 'superkingdom': 'Eukaryota',
        'lineage': 'eukaryota>metazoa', 'classification': ETE3_CLASSIFICATION,
        'alternative_taxids': '9606',
    }])
    frame = taxonomy().standardize(ete3)
    assert frame.lineage.iloc[0] == 'eukaryota>metazoa'
    assert 'lineage_taxa' not in frame.columns


def test_a_taxon_with_no_ancestors_is_still_classified():
    frame = taxonomy().standardize(uniprot_frame(
        {'taxonId': 1, 'scientificName': 'root'}))
    assert frame.classification.iloc[0] == 'root'


def test_an_empty_answer_still_has_the_columns():
    """So that reading a column need not be guarded by a check that anything
    matched at all."""
    empty = taxonomy().fetchall([999999999])
    assert empty.empty
    assert empty.columns.tolist() == taxonomy().taxcols


def test_rows_answer_for_the_taxon_they_describe():
    cursor = taxonomy()
    assert cursor.getids(uniprot_frame(HUMAN)) == {'9606'}


def test_a_retired_taxid_answers_under_its_replacement():
    """A backend that followed a merge reports both, and without the second
    the query would look unanswered and be passed on."""
    frame = pd.DataFrame([{'taxid': '9606', 'alternative_taxids': '63221,9606'}])
    assert taxonomy().getids(frame) >= {'63221', '9606'}


def test_results_are_gathered_into_one_frame():
    backend = Backend([
        (uniprot_frame(HUMAN), {'9606'}),
        (uniprot_frame({'taxonId': 562, 'scientificName': 'Escherichia coli',
                        'lineage': [{'scientificName': 'Bacteria'}]}), {'562'}),
    ])
    frame = taxonomy(backend).fetchall([9606, 562])
    assert sorted(frame.taxid) == ['562', '9606']
    assert frame.columns[:6].tolist() == taxonomy().taxcols


def test_a_number_is_accepted_where_a_taxid_is_expected():
    backend = Backend([(uniprot_frame(HUMAN), {'9606'})])
    assert len(taxonomy(backend).fetchall(9606)) == 1


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-q']))
