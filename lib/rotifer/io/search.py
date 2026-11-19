#!/usr/bin/env python3
"""
A standard table for sequence search results.

PSI-BLAST, PyHMMER and HMMER's text reports describe the same thing (a
region of a sequence matched by a query) with different column names,
coordinate fields and scores. The search functions keep returning their own
layouts; this module converts between those layouts and one standard table,
so downstream code can handle all of them in the same way.

The standard table
------------------
One row per matched region. The columns follow ``bin/hmmer2table``:

  sequence   identifier of the matched sequence (the BLAST hit, the HMMER target)
  model      identifier of the other side (the BLAST query, the HMMER profile)
  start      first aligned residue in ``sequence``
  end        last aligned residue in ``sequence``
  evalue     E-value of the region (HMMER: independent E-value)
  cov        aligned fraction of the search query, from 0 to 1
  qstart     first aligned position in ``model`` (residue or HMM state)
  qend       last aligned position in ``model``
  iteration  search round (PSI-BLAST, jackhmmer); 1 otherwise
  score      bit score of the region
  talilen    length of the region aligned in ``sequence``: end - start + 1
  qlen       full length of the search query
  estart     first residue of the envelope in ``sequence`` (BLAST: start)
  eend       last residue of the envelope in ``sequence`` (BLAST: end)
  program    program that produced the row (psiblast, hmmscan, hmmsearch, ...)
  tlen       full length of ``sequence``
  source     database or model file searched, when known

Coordinates are 1-based and inclusive: the region is ``seq[start - 1:end]``.

The search query is what the search was run with: the protein for PSI-BLAST
and hmmscan, the profile for hmmsearch. ``cov`` and ``qlen`` always refer to
it, so for hmmsearch ``cov`` is the fraction of the model that was aligned.

Differences from hmmer2table:

  * ``talilen`` is ``end - start + 1``; hmmer2table writes ``end - start``,
    one short. :func:`from_hmmer2table` corrects it and
    :func:`to_hmmer2table` writes hmmer2table's value back.
  * ``tlen`` and ``source`` are new.

Converting
----------
Every ``from_*`` function returns the standard columns first, followed by the
input's columns that have no standard equivalent, so ``to_*`` can restore
them: a native table converted to the standard table and back is unchanged.
"""

import numpy as np
import pandas as pd

# Column name -> kind. 'id' and 'text' columns are left as they are;
# 'int' columns must be complete; 'Int' columns may be missing.
SCHEMA = {
    'sequence':  'id',
    'model':     'id',
    'start':     'int',
    'end':       'int',
    'evalue':    'float',
    'cov':       'float',
    'qstart':    'Int',
    'qend':      'Int',
    'iteration': 'int',
    'score':     'float',
    'talilen':   'int',
    'qlen':      'Int',
    'estart':    'int',
    'eend':      'int',
    'program':   'text',
    'tlen':      'Int',
    'source':    'text',
}
COLUMNS = list(SCHEMA)

_DTYPES = { 'int': 'int64', 'Int': 'Int64', 'float': 'float64' }

# psiblast()'s table, as written by blast2table
BLAST_COLUMNS = ['hit', 'query', 'hitstart', 'hitend', 'evalue', 'querycoverage',
                 'querystart', 'queryend', 'iteration', 'bitscore', 'length']

# hmmer2table's columns
HMMER2TABLE_COLUMNS = ['sequence', 'model', 'start', 'end', 'evalue', 'cov', 'qstart',
                       'qend', 'iteration', 'score', 'talilen', 'qlen', 'estart', 'eend',
                       'program']

# pyhmmer_to_df's default columns plus the model coordinates. Pass it as
# columns= to hmmscan/hmmsearch to get every standard column filled.
PYHMMER_COLUMNS = ['aln_target_name', 'aln_hmm_name', 'i_evalue', 'c_evalue', 'score',
                   'env_score', 'aln_target_from', 'aln_target_to', 'aln_target_length',
                   'aln_hmm_length', 'env_from', 'env_to', 'aln_hmm_from', 'aln_hmm_to']

# What rename=True in hmmer.py does to pyhmmer_to_df's columns
PYHMMER_RENAME = {'aln_target_name': 'sequence', 'aln_hmm_name': 'model',
                  'i_evalue': 'evalue', 'env_from': 'estart', 'env_to': 'eend'}

PYHMMER_PROGRAMS = ('hmmscan', 'hmmsearch')


def empty():
    '''An empty standard table with the standard column types.'''
    return validate(pd.DataFrame({ c: [] for c in COLUMNS }))


def validate(df):
    '''
    Check that a table follows the standard and return it with the standard
    column types.

    Columns outside the standard are kept, after the standard ones.

    Raises
    ------
    ValueError
        if a standard column is missing, a column can't take its type, or a
        coordinate is below 1, a start is after its end, or ``cov`` is
        outside 0-1.
    '''
    missing = [ c for c in COLUMNS if c not in df.columns ]
    if missing:
        raise ValueError(f'missing standard columns: {", ".join(missing)}')
    df = df[COLUMNS + [ c for c in df.columns if c not in SCHEMA ]].copy()

    for column, kind in SCHEMA.items():
        if kind not in _DTYPES:
            continue
        try:
            df[column] = _numeric(df[column], kind)
        except (ValueError, TypeError) as error:
            raise ValueError(f'column {column} must be {_DTYPES[kind]}: {error}') from None

    for first, last in [('start', 'end'), ('estart', 'eend'), ('qstart', 'qend')]:
        if (df[first] < 1).any() or (df[last] < 1).any():
            raise ValueError(f'{first}/{last} must be 1-based (>= 1)')
        if (df[first] > df[last]).any():
            raise ValueError(f'{first} must not be after {last}')
    if ((df['cov'] < 0) | (df['cov'] > 1)).any():
        raise ValueError('cov must be a fraction between 0 and 1')
    return df


def from_blast(df, program='psiblast', source=None):
    '''
    Convert psiblast()'s table (blast2table's columns) to the standard table.

    blast2table doesn't keep the query's full length, so ``qlen`` is
    missing. The envelope is the aligned region.
    '''
    _require(df, BLAST_COLUMNS, 'from_blast')
    std = pd.DataFrame({
        'sequence':  df['hit'],
        'model':     df['query'],
        'start':     df['hitstart'],
        'end':       df['hitend'],
        'evalue':    df['evalue'],
        'cov':       df['querycoverage'],
        'qstart':    df['querystart'],
        'qend':      df['queryend'],
        'iteration': df['iteration'],
        'score':     df['bitscore'],
        'talilen':   df['hitend'] - df['hitstart'] + 1,
        'qlen':      pd.NA,
        'estart':    df['hitstart'],
        'eend':      df['hitend'],
        'program':   program,
        'tlen':      df['length'],
        'source':    source,
    }, index=df.index)
    return validate(_with_extra(std, df, BLAST_COLUMNS))


def to_blast(df):
    '''Convert a standard table to psiblast()'s columns (blast2table's).'''
    df = validate(df)
    out = pd.DataFrame({
        'hit':           df['sequence'],
        'query':         df['model'],
        'hitstart':      df['start'],
        'hitend':        df['end'],
        'evalue':        df['evalue'],
        'querycoverage': df['cov'],
        'querystart':    df['qstart'],
        'queryend':      df['qend'],
        'iteration':     df['iteration'],
        'bitscore':      df['score'],
        'length':        df['tlen'],
    }, index=df.index)
    return _restore_extra(out, df)


def from_pyhmmer(df, program, source=None):
    '''
    Convert hmmscan/hmmscan_linear/hmmsearch output to the standard table.

    Parameters
    ----------
    df : pandas.DataFrame
        pyhmmer_to_df's columns, renamed (rename=True) or not.
    program : {'hmmscan', 'hmmsearch'}
        Which search made the table. It decides which side is the search
        query, so what ``cov`` and ``qlen`` measure: the sequence for
        hmmscan, the model for hmmsearch.
    source : str, optional
        Model file searched. Defaults to the table's ``source`` column,
        which hmmscan and hmmsearch add.

    Columns that pyhmmer_to_df only writes on request are left missing when
    absent: ``qstart``/``qend`` need aln_hmm_from/aln_hmm_to (see
    PYHMMER_COLUMNS), and so does ``cov`` for hmmsearch.
    '''
    if program not in PYHMMER_PROGRAMS:
        raise ValueError(f'program must be one of {", ".join(PYHMMER_PROGRAMS)}, not {program!r}')
    renamed = 'sequence' in df.columns
    raw = df.rename(columns={ v: k for k, v in PYHMMER_RENAME.items() }) if renamed else df
    _require(raw, ['aln_target_name', 'aln_hmm_name', 'i_evalue', 'score',
                   'aln_target_from', 'aln_target_to', 'env_from', 'env_to'], 'from_pyhmmer')

    def column(name):
        return raw[name] if name in raw.columns else pd.Series(pd.NA, index=raw.index)

    start, end = raw['aln_target_from'], raw['aln_target_to']
    qstart, qend = column('aln_hmm_from'), column('aln_hmm_to')
    tlen, hmm_length = column('aln_target_length'), column('aln_hmm_length')
    if program == 'hmmscan':
        qlen = tlen
        cov = (end - start + 1) / _numeric(tlen, 'float')
    else:
        qlen = hmm_length
        cov = (_numeric(qend, 'float') - _numeric(qstart, 'float') + 1) / _numeric(hmm_length, 'float')

    std = pd.DataFrame({
        'sequence':  raw['aln_target_name'],
        'model':     raw['aln_hmm_name'],
        'start':     start,
        'end':       end,
        'evalue':    raw['i_evalue'],
        'cov':       cov,
        'qstart':    qstart,
        'qend':      qend,
        'iteration': 1,
        'score':     raw['score'],
        'talilen':   end - start + 1,
        'qlen':      qlen,
        'estart':    raw['env_from'],
        'eend':      raw['env_to'],
        'program':   program,
        'tlen':      tlen,
        'source':    source if source is not None else column('source'),
    }, index=raw.index)
    known = ['aln_target_name', 'aln_hmm_name', 'i_evalue', 'score', 'aln_target_from',
             'aln_target_to', 'aln_hmm_from', 'aln_hmm_to', 'aln_target_length',
             'env_from', 'env_to', 'source']
    std = _with_extra(std, raw, known)
    std.attrs['pyhmmer_renamed'] = renamed
    return validate(std)


def to_pyhmmer(df, rename=None):
    '''
    Convert a standard table to pyhmmer_to_df's columns.

    Parameters
    ----------
    rename : bool, optional
        Apply hmm.py's rename (sequence, model, evalue, estart, eend). By
        default, whatever the table had when :func:`from_pyhmmer` read it,
        and False for other tables.
    '''
    if rename is None:
        rename = df.attrs.get('pyhmmer_renamed', False)
    df = validate(df)
    hmmsearch = df['program'] == 'hmmsearch'
    out = pd.DataFrame({
        'aln_target_name':   df['sequence'],
        'aln_hmm_name':      df['model'],
        'i_evalue':          df['evalue'],
        'score':             df['score'],
        'aln_target_from':   df['start'],
        'aln_target_to':     df['end'],
        'aln_target_length': df['tlen'],
        'env_from':          df['estart'],
        'env_to':            df['eend'],
        'aln_hmm_from':      df['qstart'],
        'aln_hmm_to':        df['qend'],
    }, index=df.index)
    extra = _extra(df)
    if 'aln_hmm_length' not in extra.columns:
        out['aln_hmm_length'] = df['qlen'].where(hmmsearch, pd.NA)
    if df['source'].notna().any():
        out['source'] = df['source']
    out = _restore_extra(out, df)
    if df['qstart'].isna().all() and df['qend'].isna().all():
        out = out.drop(columns=['aln_hmm_from', 'aln_hmm_to'])
    if rename:
        out = out.rename(columns=PYHMMER_RENAME)
    return out


def from_hmmer2table(df, source=None):
    '''
    Convert hmmer2table's table to the standard table.

    ``talilen`` is corrected to ``end - start + 1``. hmmer2table doesn't keep
    the target's full length, so ``tlen`` is missing.
    '''
    _require(df, HMMER2TABLE_COLUMNS, 'from_hmmer2table')
    std = df[HMMER2TABLE_COLUMNS].copy()
    std['talilen'] = std['end'] - std['start'] + 1
    std['iteration'] = std['iteration'].fillna(1)
    std['tlen'] = pd.NA
    std['source'] = source
    return validate(_with_extra(std, df, HMMER2TABLE_COLUMNS))


def to_hmmer2table(df):
    '''
    Convert a standard table to hmmer2table's columns, with hmmer2table's
    ``talilen`` (``end - start``).
    '''
    df = validate(df)
    out = df[HMMER2TABLE_COLUMNS].copy()
    out['talilen'] = out['end'] - out['start']
    return _restore_extra(out, df)


# Helpers

def _numeric(series, kind):
    if kind == 'Int':
        series = pd.to_numeric(series.astype(object).where(series.notna(), None))
        return series.astype('Int64')
    if kind == 'int':
        if series.isna().any():
            raise ValueError('missing values')
        return pd.to_numeric(series).astype('int64')
    return pd.to_numeric(series.astype(object).where(series.notna(), np.nan)).astype('float64')


def _require(df, columns, caller):
    missing = [ c for c in columns if c not in df.columns ]
    if missing:
        raise ValueError(f'{caller}: missing columns {", ".join(missing)}')


def _with_extra(std, df, known):
    '''Append df's columns that the conversion didn't use.'''
    extra = [ c for c in df.columns if c not in known and c not in std.columns ]
    return pd.concat([std, df[extra]], axis=1) if extra else std


def _extra(df):
    return df[[ c for c in df.columns if c not in SCHEMA ]]


def _restore_extra(out, df):
    extra = _extra(df)
    extra = extra[[ c for c in extra.columns if c not in out.columns ]]
    return pd.concat([out, extra], axis=1) if len(extra.columns) else out
