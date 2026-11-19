"""
Tests for rotifer.interval.architecture.select_regions and read_regions.

What has to hold, and why:

- A fraction overlap (float, or text with a decimal point or exponent) is
  a share of the shorter region, a whole number (int, or text such as
  '10') is a number of residues, and both limits are exclusive: this is
  how the former Perl tools read ``-o``.
- Regions are ranked by E-value, then length, then start; by_length puts
  length first and the start before the E-value. Explicit criteria take
  ascending and descending order. Remaining ties keep input order;
  kept regions with the same coordinates are listed best first.
- ``nested`` and ``ignore_overlap`` decide whether a region inside a kept
  region survives, and ``'all'`` keeps everything.
- Output keeps the input rows as they were (text, dtypes, index labels),
  sorted by sequence and numeric start and end, and the input is never
  modified.
- The vectorised rounds give the same result as checking the regions one
  at a time, which a straightforward reference implementation confirms on
  seeded random tables.
- read_regions filters and renames like the command, and keeps the text.

Fixtures: small DataFrames built by ``_regions`` (sequence, domain, start,
end, evalue) and text tables written to ``tmp_path``.

Run under pytest: PYTHONPATH=lib pytest test/interval/test_select_regions.py
"""

import gzip
import io

import numpy as np
import pandas as pd
import pytest

from rotifer.interval.architecture import INPUT_COLUMNS, read_regions, select_regions

def _regions(*rows, sequence='p1'):
    """
    DataFrame of regions from ``(domain, start, end, evalue)`` tuples, all
    in one sequence unless a fifth item names another.

    Parameters
    ----------
    *rows : tuple
        Domain, start, end, E-value and optionally the sequence.
    sequence : str, default 'p1'
        Sequence of the rows that do not name one.

    Returns
    -------
    pandas.DataFrame
        Columns sequence, domain, start, end, evalue.
    """
    return pd.DataFrame([[row[4] if len(row) > 4 else sequence, *row[:4]] for row in rows],
                        columns=INPUT_COLUMNS)

def _names(df):
    """
    Domain names of a result, in order.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of ``select_regions``.

    Returns
    -------
    list of str
    """
    return df['domain'].tolist()

@pytest.mark.parametrize("overlap,expected", [
    (0.4, ['A']),        # 10 shared > 0.4 * 20
    (0.5, ['A', 'B']),   # 10 shared is not more than 0.5 * 20
    (9, ['A']),
    (10, ['A', 'B']),    # exclusive limit in residues
    ('0.4', ['A']),
    ('0.5', ['A', 'B']),
    ('10', ['A', 'B']),
    ('9', ['A']),
    ('1.0', ['A', 'B']),  # a fraction: everything shorter than the whole region
    ('1e-1', ['A']),      # exponent means a fraction
    (np.int64(10), ['A', 'B']),
    (np.float64(0.4), ['A']),
])
def test_overlap_kinds_and_exclusive_limits(overlap, expected):
    """
    Fractions, residues and text overlaps, with exclusive limits.

    Parameters
    ----------
    overlap : int, float or str
        Overlap argument.
    expected : list of str
        Domains kept. A (1..100) shares 10 residues with B (91..110).
    """
    df = _regions(('A', 1, 100, 1e-10), ('B', 91, 110, 1e-5))
    assert _names(select_regions(df, overlap=overlap)) == expected

def test_string_one_point_zero_is_fraction_and_ten_is_residues():
    """
    '1.0' is a fraction (1.0 of the shorter region) while '10' is ten residues.
    """
    df = _regions(('A', 1, 100, 1e-10), ('B', 51, 70, 1e-5))   # B has 20 residues inside A
    assert _names(select_regions(df, overlap='1.0')) == ['A', 'B']   # 20 is not more than 1.0 * 20
    assert _names(select_regions(df, overlap='10')) == ['A']    # 20 > 10 residues

def test_bool_overlap_is_a_number():
    """
    True counts as one residue; booleans are not special.
    """
    df = _regions(('A', 1, 10, 1e-10), ('B', 10, 20, 1e-5))
    assert _names(select_regions(df, overlap=True)) == ['A', 'B']   # 1 shared is not > 1

def test_default_ranking_evalue_length_start():
    """
    Lowest E-value first, then longest, then smallest start.
    """
    df = _regions(('late', 50, 60, 1e-5), ('long', 20, 70, 1e-5), ('best', 1, 5, 1e-9))
    kept = select_regions(df, overlap=0)
    assert _names(kept) == ['best', 'long']   # 'late' lies in 'long' and has same E-value, shorter

def test_by_length_order_and_start_before_evalue():
    """
    by_length ranks length, then start, then E-value (not E-value before start).
    """
    df = _regions(('X', 10, 29, 1e-1), ('Y', 1, 20, 1e-9))   # same length, Y starts first, X worse E-value
    assert _names(select_regions(df, by_length=True)) == ['Y']
    df = _regions(('X', 1, 20, 1e-1), ('Y', 1, 20, 1e-9))     # same length and start: E-value decides
    assert _names(select_regions(df, by_length=True)) == ['Y']
    df = _regions(('S', 1, 15, 1e-20), ('L', 5, 40, 1e-1))    # longest wins whatever the E-value
    assert _names(select_regions(df, by_length=True)) == ['L']
    assert _names(select_regions(df)) == ['S']

def test_input_order_breaks_complete_ties():
    """
    Regions equal on every criterion keep input order.
    """
    df = _regions(('B', 1, 20, 1e-5), ('A', 1, 20, 1e-5))
    assert _names(select_regions(df)) == ['B']

def test_same_coordinates_keep_rank_order():
    """
    Kept regions with the same start and end come out best first, not in
    input order, as the former Perl tool sorted its kept regions.
    """
    df = _regions(('W', 1, 20, 1e-2), ('V', 1, 20, 1e-9))
    assert _names(select_regions(df, ignore_overlap={'all'})) == ['V', 'W']

def test_ignore_all_needs_no_domain_column():
    """
    ignore_overlap={'all'} does not look at domain names, so the table
    needs no domain column.
    """
    df = _regions(('A', 1, 100, 1e-10), ('B', 20, 40, 1e-5)).drop(columns='domain')
    assert len(select_regions(df, ignore_overlap={'all'})) == 2

def test_explicit_criteria_ascending_and_descending():
    """
    True prefers small values, False large ones; the order of the dict is the precedence.
    """
    df = _regions(('A', 1, 20, 1e-5), ('B', 5, 30, 1e-9))
    df['score'] = [10, 3]
    assert _names(select_regions(df, criteria={'score': False})) == ['A']
    assert _names(select_regions(df, criteria={'score': True})) == ['B']
    assert _names(select_regions(df, criteria={'evalue': True})) == ['B']
    assert _names(select_regions(df, criteria={'evalue': False})) == ['A']

def test_region_length_computed_and_from_a_column():
    """
    'region_length' is end - start + 1 unless the table has a column of that name.
    """
    df = _regions(('A', 1, 20, 1e-5), ('B', 5, 40, 1e-5))
    assert _names(select_regions(df, criteria={'region_length': False})) == ['B']
    df['region_length'] = [100, 1]   # the column wins for ranking
    assert _names(select_regions(df, criteria={'region_length': False})) == ['A']
    df['region_length'] = [100, 1]
    assert _names(select_regions(df, overlap=0.5, criteria={'region_length': True})) == ['B']

def test_criteria_and_by_length_are_exclusive():
    """
    Giving both criteria and by_length=True is an error.
    """
    with pytest.raises(ValueError):
        select_regions(_regions(('A', 1, 2, 1e-1)), criteria={'evalue': True}, by_length=True)

@pytest.mark.parametrize("nested,expected", [(True, ['A', 'B']), (False, ['A'])])
def test_nested(nested, expected):
    """
    nested=False drops a region contained in a kept one whatever the overlap.

    Parameters
    ----------
    nested : bool
        Argument of select_regions.
    expected : list of str
        Domains kept with a 1000-residue overlap allowance.
    """
    df = _regions(('A', 1, 100, 1e-10), ('B', 20, 40, 1e-5))
    assert _names(select_regions(df, overlap=1000, nested=nested)) == expected

def test_ignore_overlap_by_name_and_all():
    """
    Ignored domains are always kept and never remove others; 'all' keeps everything.
    """
    df = _regions(('A', 1, 100, 1e-10), ('B', 20, 40, 1e-5), ('C', 30, 50, 1e-3))
    assert _names(select_regions(df)) == ['A']
    assert _names(select_regions(df, ignore_overlap=['B'])) == ['A', 'B']        # B kept, C removed by A
    assert _names(select_regions(df, ignore_overlap={'A'})) == ['A', 'B']   # A removes nobody; B removes C
    assert _names(select_regions(df, ignore_overlap={'all'})) == ['A', 'B', 'C']
    assert _names(select_regions(df, ignore_overlap='B')) == ['A', 'B']
    assert _names(select_regions(df, ignore_overlap=['B'], nested=False)) == ['A', 'B']

def test_domain_auto_detection():
    """
    ignore_overlap looks in 'domain', else in 'model'; an explicit name overrides.
    """
    df = _regions(('A', 1, 100, 1e-10), ('B', 20, 40, 1e-5))
    models = df.rename(columns={'domain': 'model'})
    assert select_regions(models, ignore_overlap=['B'])['model'].tolist() == ['A', 'B']
    both = df.assign(model=['B', 'A'])   # 'domain' wins
    assert _names(select_regions(both, ignore_overlap=['B'])) == ['A', 'B']
    assert _names(select_regions(both, ignore_overlap=['A'], domain='model')) == ['A', 'B']

def test_several_reference_columns_and_missing_references():
    """
    Regions only conflict within the same combination of reference values;
    missing values form their own group.
    """
    df = pd.DataFrame({'q': ['a', 'a', 'a', None, None], 'h': ['x', 'x', 'y', 'x', 'x'],
                       'start': [1, 5, 5, 1, 5], 'end': [50, 60, 60, 50, 60],
                       'evalue': [1e-9, 1e-5, 1e-5, 1e-9, 1e-5], 'domain': list('ABCDE')})
    kept = select_regions(df, reference=['q', 'h'])
    assert _names(kept) == ['A', 'C', 'D']
    one = select_regions(df, reference='q')
    assert _names(one) == ['A', 'D']

def test_output_sorted_by_reference_then_numeric_coordinates():
    """
    Rows come out by sequence, then numeric start and end (text numbers too).
    """
    df = _regions(('late', '100', '120', '1e-9', 'p2'), ('b', '9', '30', '1e-9', 'p1'),
                  ('a', '10', '25', '1e-9', 'p1'), ('top', '2', '5', '1e-9', 'p1'))
    kept = select_regions(df, overlap=100)
    assert _names(kept) == ['top', 'b', 'a', 'late']

def test_string_tables_keep_their_text():
    """
    Text columns are ranked as Perl numbers and returned unchanged.
    """
    df = _regions(('A', '1', '60', '1e-10'), ('B', '50', '69', '1E-5'), ('C', '70', '120', '0.001'))
    kept = select_regions(df)
    assert kept['evalue'].tolist() == ['1e-10', '0.001']
    assert kept['start'].tolist() == ['1', '70']

def test_input_not_modified_and_index_labels_kept():
    """
    The input is unchanged, and kept rows keep their labels, duplicates included.
    """
    df = _regions(('A', 1, 60, 1e-10), ('B', 50, 69, 1e-5), ('C', 70, 120, 1e-3))
    df.index = ['x', 'x', 'y']
    before = df.copy()
    kept = select_regions(df)
    assert df.equals(before)
    assert kept.index.tolist() == ['x', 'y']
    assert kept.dtypes.equals(df.dtypes)

def test_empty_input_keeps_columns():
    """
    No rows in, no rows out, same columns.
    """
    df = _regions(('A', 1, 2, 1e-1)).iloc[:0]
    kept = select_regions(df)
    assert kept.empty and kept.columns.tolist() == df.columns.tolist()

def test_missing_columns_raise_key_error():
    """
    Missing reference, start, end, criteria and domain columns raise KeyError.
    """
    df = _regions(('A', 1, 20, 1e-5))
    for options in ({'reference': 'nope'}, {'start': 'nope'}, {'end': 'nope'},
                    {'criteria': {'nope': True}}, {'reference': ['sequence', 'nope']}):
        with pytest.raises(KeyError):
            select_regions(df, **options)
    with pytest.raises(KeyError):
        select_regions(df.drop(columns='domain'), ignore_overlap=['A'])
    with pytest.raises(KeyError):
        select_regions(df, ignore_overlap=['A'], domain='nope')
    # No evalue column: the default ranking skips it
    assert _names(select_regions(df.drop(columns='evalue'))) == ['A']

def _reference(df, overlap, nested, ignore, criteria):
    """
    One-region-at-a-time greedy selection, written for clarity.

    Parameters
    ----------
    df : pandas.DataFrame
        Regions with numeric start, end and evalue.
    overlap : int or float
        Residues (int) or fraction (float).
    nested : bool
        Whether regions inside kept regions may be kept.
    ignore : set of str
        Domains that are always kept and never remove others.
    criteria : list of (str, bool)
        Ranking columns and whether ascending.

    Returns
    -------
    list of int
        Positions of the kept rows, in output order.
    """
    out = []
    for _, group in df.groupby('sequence', sort=True):
        rows = [(i, r) for i, r in enumerate(group.itertuples())]
        ranked = list(range(len(rows)))
        for name, ascending in reversed(criteria):
            if name == 'region_length':
                key = lambda i: rows[i][1].end - rows[i][1].start + 1
            else:
                key = lambda i, name=name: getattr(rows[i][1], name)
            ranked.sort(key=key, reverse=not ascending)
        ranked.sort(key=lambda i: rows[i][1].start)  # not the final key: redone below
        ranked = list(range(len(rows)))
        keys = [(n, a) for n, a in criteria] + [('start', True)]
        for name, ascending in reversed(keys):
            if name == 'region_length':
                key = lambda i: rows[i][1].end - rows[i][1].start + 1
            else:
                key = lambda i, name=name: getattr(rows[i][1], name)
            ranked.sort(key=key, reverse=not ascending)
        kept = []
        for i in ranked:
            a = rows[i][1]
            keep = True
            if a.domain not in ignore:
                for j in kept:
                    b = rows[j][1]
                    if b.domain in ignore:
                        continue
                    shared = min(a.end, b.end) - max(a.start, b.start) + 1
                    limit = overlap * min(a.end - a.start + 1, b.end - b.start + 1) \
                        if isinstance(overlap, float) else overlap
                    inside = b.start <= a.start and b.end >= a.end
                    if shared > limit or (not nested and inside):
                        keep = False
                        break
            if keep:
                kept.append(i)
        kept.sort(key=lambda i: (rows[i][1].start, rows[i][1].end, i))
        out.extend(group.index[kept].tolist())
    return out

@pytest.mark.parametrize("seed", range(20))
@pytest.mark.parametrize("options", [
    {'overlap': 0.4},
    {'overlap': 0.0, 'nested': False},
    {'overlap': 7},
    {'overlap': 1000, 'nested': False},
    {'overlap': 0.4, 'by_length': True},
    {'overlap': 5, 'ignore_overlap': ['A', 'C']},
])
def test_random_tables_match_one_at_a_time_reference(seed, options):
    """
    Vectorised selection equals the greedy reference on random tables.

    Parameters
    ----------
    seed : int
        Seed of the random table (40 regions in 6 sequences, many ties).
    options : dict
        Arguments of select_regions.
    """
    rng = np.random.default_rng(seed)
    n = 40
    start = rng.integers(1, 80, n)
    df = pd.DataFrame({'sequence': rng.choice(list('pqrstu'), n), 'domain': rng.choice(list('ABCD'), n),
                       'start': start, 'end': start + rng.integers(0, 40, n),
                       'evalue': rng.choice([1e-9, 1e-5, 1e-2], n)})
    df.index = rng.permutation(n)
    by_length = options.get('by_length', False)
    if by_length:
        criteria = [('region_length', False), ('start', True), ('evalue', True)]
    else:
        criteria = [('evalue', True), ('region_length', False)]
    expected = _reference(df, options['overlap'], options.get('nested', True),
                          set(options.get('ignore_overlap', ())), criteria)
    assert select_regions(df, **options).index.tolist() == expected

def _write(path, text, compress=False):
    """
    Write a table to a file, optionally gzip-compressed.

    Parameters
    ----------
    path : pathlib.Path
        Destination.
    text : str
        Contents.
    compress : bool, default False
        Compress with gzip.

    Returns
    -------
    str
        The path as text.
    """
    if compress:
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return str(path)

TABLE = ("p1\tA\t1\t60\t1e-10\n"
         "p1\tB\t50\t69\t0.5\n"
         "p2\tC\t5\t9\t0.005\n"
         "p2\tD\t10\t20\tabc\n")

def test_read_regions_columns_and_text():
    """
    Default columns, strings as written, filtered by the E-value cutoff.
    """
    df = read_regions(io.StringIO(TABLE))
    assert df.columns.tolist() == INPUT_COLUMNS
    assert df['domain'].tolist() == ['A', 'C', 'D']      # 0.5 > 0.01; 'abc' counts as 0
    assert df['evalue'].tolist() == ['1e-10', '0.005', 'abc']
    assert df.index.tolist() == [0, 1, 2]
    assert read_regions(io.StringIO(TABLE), evalue=0.001)['domain'].tolist() == ['A', 'D']

def test_read_regions_rename_exclude_and_target_cutoff():
    """
    Renaming happens first: exclude and target_cutoff use the new names.
    """
    df = read_regions(io.StringIO(TABLE), rename={'C': 'Z'}, exclude=['Z'])
    assert df['domain'].tolist() == ['A', 'D']
    df = read_regions(io.StringIO(TABLE), rename={'B': 'Y'}, target_cutoff={'Y': 1, 'A': 1e-20})
    assert df['domain'].tolist() == ['Y', 'C', 'D']

def test_read_regions_header_and_columns():
    """
    Column names from a header (per table) or from ``columns``; unknown names are dropped.
    """
    text = "evalue\tx\tend\tstart\tdomain\tsequence\n1e-5\t-\t9\t1\tA\tp1\n"
    df = read_regions(io.StringIO(text), header=True)
    assert df.iloc[0].tolist() == ['p1', 'A', '1', '9', '1e-5']
    df = read_regions(io.StringIO("p1,A,1,9,1e-5\n"), columns=['sequence', 'domain', 'start', 'end', 'evalue'],
                      delimiter=',')
    assert df.iloc[0].tolist() == ['p1', 'A', '1', '9', '1e-5']

def test_read_regions_short_lines_and_missing_evalue_column():
    """
    Missing cells are '' for sequence and domain and None for the rest.
    """
    df = read_regions(io.StringIO("p1\tA\t1\t9\n"))
    assert df['evalue'].tolist() == [None]
    df = read_regions(io.StringIO("p1\n"))
    assert df.iloc[0].tolist() == ['p1', '', None, None, None]
    df = read_regions(io.StringIO("p1\tA\t1\t9\t1\n"), columns=['sequence', 'domain', 'start', 'end'])
    assert df['evalue'].tolist() == [None]

def test_read_regions_several_files_and_compressed(tmp_path):
    """
    Several tables concatenate in file order; gzip files are read.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory.
    """
    one = _write(tmp_path / 'one.tsv', "p1\tA\t1\t9\t1e-5\n")
    two = _write(tmp_path / 'two.tsv.gz', "p2\tB\t2\t8\t1e-5\n", compress=True)
    df = read_regions([one, two])
    assert df['sequence'].tolist() == ['p1', 'p2']
    assert read_regions(two)['domain'].tolist() == ['B']
