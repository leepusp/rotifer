#!/usr/bin/env python3
"""Compare NeighborhoodDF.compact with a frozen copy of its #145 implementation.

compact_reference.py holds the old per-block implementation. Both run on
seeded random neighborhood tables and must return identical frames (dtypes
and return type included) or raise the same exception type.

The tables only use inputs on which the intended behaviour did not change:
query values 0 and 1, and blocks whose rows are contiguous. Where the
reference raises TypeError because a missing value reached the join, compact
must instead return a frame without missing values.

Run under pytest, or standalone: python test/genome/test_compact_equivalence.py
"""

import numpy as np
import pandas as pd

from rotifer.genome.data import NeighborhoodDF
from compact_reference import compact_reference

ANNOTATION = 'cluster'
EXTRA = ['lineage', 'organism', 'seq_type', 'classification']
SEEDS = range(40)
OPTIONS = [dict(strand=s, fill_unk=f)
           for s in (False, True)
           for f in (False, 'pid', ANNOTATION, 'gene')]


def random_table(rng, n_blocks, missing_query=False, nan_values=False):
    """Contiguous blocks in random block order with a non-default index.

    missing_query: some blocks have no query row.
    nan_values: some rows have no type, and some 'gene' values are missing.
    """
    rows = []
    for block_id in rng.permutation(np.arange(1, n_blocks + 1)) * 7:
        size = int(rng.integers(1, 13))
        n_queries = 0 if missing_query and rng.random() < 0.3 else int(rng.choice([1, 1, 1, 2, 3]))
        queries = set(rng.choice(size, min(n_queries, size), replace=False).tolist())
        for i in range(size):
            u = rng.random()
            types = ['tRNA', 'rRNA', 'pseudogene'] + ([None] if nan_values else [])
            rows.append(dict(
                block_id=int(block_id),
                pid=None if rng.random() < 0.05 else f'P{len(rows)}',
                cluster=None if rng.random() < 0.15 else f'C{rng.integers(0, 8)}',
                strand=np.nan if u < 0.05 else (0 if u < 0.08 else int(rng.choice([1, -1]))),
                query=1 if i in queries else 0,
                type='CDS' if rng.random() < 0.85 else rng.choice(types),
                gene=None if nan_values and rng.random() < 0.3 else f'g{rng.integers(0, 50)}',
                **{c: rng.choice(['x', 'y', None]) for c in EXTRA},
            ))
    df = pd.DataFrame(rows)
    df.index = rng.permutation(np.arange(len(df))) + 1000
    return NeighborhoodDF(df)


def table_for_seed(seed):
    rng = np.random.default_rng(seed)
    return random_table(rng, int(rng.integers(1, 16)),
                        missing_query=(seed % 5 == 0),
                        nan_values=(seed % 5 == 3))


def run(fn):
    try:
        return fn(), None
    except Exception as e:
        return None, e


def difference(table, **kwargs):
    """None if both implementations agree, otherwise a description."""
    expected, expected_error = run(lambda: compact_reference(table, annotation=ANNOTATION, **kwargs))
    actual, actual_error = run(lambda: table.compact(annotation=ANNOTATION, **kwargs))
    if isinstance(expected_error, TypeError) and actual_error is None:
        # The reference fails with TypeError when a missing type or fill_unk
        # value reaches the join; compact now shows '?' there instead.
        if actual[['compact', 'pid_compact']].isna().any().any():
            return 'compact output has missing values'
        return None
    if expected_error or actual_error:
        if type(expected_error) is type(actual_error):
            return None
        return f'reference raised {expected_error!r}, compact raised {actual_error!r}'
    if type(actual) is not type(expected):
        return f'reference returned {type(expected).__name__}, compact returned {type(actual).__name__}'
    try:
        pd.testing.assert_frame_equal(actual, expected, check_dtype=True)
    except AssertionError as e:
        return str(e)
    return None


def test_matches_reference():
    failures = []
    for seed in SEEDS:
        table = table_for_seed(seed)
        for kwargs in OPTIONS:
            diff = difference(table, **kwargs)
            if diff:
                failures.append(f'seed={seed} {kwargs}: {diff}')
    assert not failures, f'{len(failures)} mismatches:\n' + '\n'.join(failures[:5])


if __name__ == '__main__':
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith('test_') and callable(f)]
    failed = 0
    for name, test in tests:
        try:
            test()
            print(f'ok    {name}')
        except Exception as e:
            failed += 1
            print(f'FAIL  {name}: {type(e).__name__}: {e}')
    print(f'{len(tests) - failed}/{len(tests)} passed')
    raise SystemExit(1 if failed else 0)
