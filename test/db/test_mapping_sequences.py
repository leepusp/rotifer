#!/usr/bin/env python3
"""Tests for looking a sequence up through the mapping cursors.

A sequence is not an identifier, but it has one: UniProt publishes a CRC64 of
the residues in idmapping.dat like any other cross-reference, so a sequence can
be asked about by checksum. That is worth having because it finds an entry
whatever the sequence has been called, and finds every entry that carries the
same residues -- one query here returns two accessions for one sequence.

Nothing in this file needs a network or a server: the conversion, the routing
and the bookkeeping are all local, and only the lookup itself is not.

The checksum is Biopython's, which produces the same digits as UniProt with a
`CRC-` prefix. That the prefix is dropped is checked against values taken from
the real table, since a checksum that is almost right finds nothing at all and
looks exactly like a sequence UniProt does not have.

Run under pytest, or standalone: python test/db/test_mapping_sequences.py
"""

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import rotifer.db.core
import rotifer.db.methods


class Cursor(rotifer.db.methods.MappingCursor, rotifer.db.core.BaseCursor):
    """A cursor that only parses: enough to see what a query becomes."""

    def fetchone(self, accessions, source=None, target=None, *args, **kwargs):
        self.asked = (sorted(self.parse_ids(accessions)), source, target)
        return iter(())


#: Q6GZX4 and its checksum, as idmapping.dat has them.
Q6GZX4 = (
    'MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS'
    'EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD'
    'AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL'
    'EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD'
    'SFRKIYTDLGWKFTPL'
)
Q6GZX4_CRC64 = 'B4840739BF7D4121'


def record(name, residues):
    return SeqRecord(Seq(residues), id=name, description='')


# ----------------------------------------------------------- the checksum

def test_the_checksum_matches_what_uniprot_publishes():
    """Taken from the real table. A checksum that is almost right finds
    nothing and looks like a sequence UniProt does not have."""
    assert rotifer.db.methods.MappingCursor.checksum(Q6GZX4) == Q6GZX4_CRC64


def test_the_biopython_prefix_is_not_part_of_it():
    value = rotifer.db.methods.MappingCursor.checksum(Q6GZX4)
    assert not value.startswith('CRC-') and len(value) == 16


def test_case_and_whitespace_do_not_change_it():
    """FASTA that has been through an editor should still find its entry."""
    plain = rotifer.db.methods.MappingCursor.checksum(Q6GZX4)
    assert rotifer.db.methods.MappingCursor.checksum('  ' + Q6GZX4.lower() + ' ') == plain


# ------------------------------------------------------ recognising them

@pytest.mark.parametrize('value', ['P00750', 'UPI0000000001', 42, None])
def test_an_identifier_is_not_mistaken_for_a_sequence(value):
    assert rotifer.db.methods.MappingCursor.is_sequence(value) is False


def test_a_biopython_record_is_a_sequence():
    assert rotifer.db.methods.MappingCursor.is_sequence(record('x', Q6GZX4)) is True


def test_an_object_holding_many_sequences_is_one_too():
    """rotifer's own sequence object keeps them in a dataframe, so a single
    object stands for as many sequences as it holds."""
    class Many:
        df = pd.DataFrame({'id': ['a', 'b'], 'sequence': [Q6GZX4, 'MKV']})

    assert rotifer.db.methods.MappingCursor.is_sequence(Many()) is True
    found = rotifer.db.methods.MappingCursor.sequence_checksums(Many())
    assert len(found) == 2
    assert Q6GZX4_CRC64 in found


# --------------------------------------------------- what a query becomes

def test_a_sequence_is_asked_about_by_its_checksum():
    cursor = Cursor(progress=False)
    cursor.fetchall(record('mine', Q6GZX4), target=['RefSeq'])
    asked, source, target = cursor.asked
    assert asked == [Q6GZX4_CRC64]


def test_naming_sources_still_looks_where_sequences_live():
    """Asking for RefSeq sources and then passing a sequence would look past
    it: a sequence is only ever found under its checksum."""
    cursor = Cursor(progress=False)
    cursor.fetchall(record('mine', Q6GZX4), source=['RefSeq'], target=['RefSeq'])
    assert 'CRC64' in cursor.asked[1]


def test_an_open_source_is_left_alone():
    """None already looks everywhere, so there is nothing to add."""
    cursor = Cursor(progress=False)
    cursor.fetchall(record('mine', Q6GZX4), target=['RefSeq'])
    assert cursor.asked[1] is None


def test_a_query_of_identifiers_alone_is_untouched():
    cursor = Cursor(progress=False)
    cursor.fetchall(['P00750'], source=['RefSeq'], target=['RefSeq'])
    assert cursor.asked == (['P00750'], ['RefSeq'], ['RefSeq'])


def test_sequences_and_identifiers_can_be_mixed():
    cursor = Cursor(progress=False)
    cursor.fetchall([record('mine', Q6GZX4), 'P00750'], target=['RefSeq'])
    assert sorted(cursor.asked[0]) == sorted([Q6GZX4_CRC64, 'P00750'])


def test_two_names_for_one_sequence_are_one_query():
    """The same residues under two names are one sequence, and asking twice
    would only make the table answer twice."""
    cursor = Cursor(progress=False)
    cursor.fetchall([record('one', Q6GZX4), record('two', Q6GZX4)], target=['RefSeq'])
    assert cursor.asked[0] == [Q6GZX4_CRC64]


# ------------------------------------------------------ finding the way back

def test_the_names_the_sequences_came_under_are_remembered():
    """A result names a sequence by its checksum, because that is what the
    data holds. Without this a caller could not say which of their sequences
    a row is about."""
    cursor = Cursor(progress=False)
    cursor.fetchall(record('mine', Q6GZX4), target=['RefSeq'])
    assert cursor.checksums.to_dict('records') == [
        {'source': Q6GZX4_CRC64, 'sequence': 'mine'}]


def test_both_names_of_one_sequence_are_kept():
    cursor = Cursor(progress=False)
    cursor.fetchall([record('one', Q6GZX4), record('two', Q6GZX4)], target=['RefSeq'])
    assert sorted(cursor.checksums.sequence) == ['one', 'two']


def test_the_table_merges_onto_a_result():
    """It is shaped to be joined on the source column, which is what a caller
    has to do to get from a row back to their sequence."""
    cursor = Cursor(progress=False)
    cursor.fetchall(record('mine', Q6GZX4), target=['RefSeq'])
    result = pd.DataFrame([{'source': Q6GZX4_CRC64, 'source_type': 'CRC64',
                            'accession': 'Q6GZX4', 'target': 'YP_031579.1',
                            'target_type': 'RefSeq'}])
    joined = result.merge(cursor.checksums, on='source')
    assert joined.sequence.tolist() == ['mine']


def test_nothing_is_remembered_when_no_sequence_was_given():
    cursor = Cursor(progress=False)
    cursor.fetchall(['P00750'], target=['RefSeq'])
    assert cursor.checksums.empty
    assert cursor.checksums.columns.tolist() == ['source', 'sequence']


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-q']))
