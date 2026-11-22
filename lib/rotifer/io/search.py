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


# BLAST tabular output (issue #141)
#
# psiblast() used to read blast2table's report of the text alignments. The
# same information is available, in a format meant for programs, as BLAST's
# tabular output (-outfmt 7, or -outfmt 11 reformatted with blast_formatter).
# read_blast_tabular() reads it and collapse_hsps() merges its HSPs the way
# blast2table's default (non-redundant) mode does.

BLAST_TABULAR_FIELDS = 'qaccver saccver sstart send evalue qstart qend bitscore slen qlen'
BLAST_OUTFMT = '7 ' + BLAST_TABULAR_FIELDS

_TABULAR_NAMES = ['qaccver', 'saccver', 'sstart', 'send', 'evalue', 'qstart', 'qend',
                  'bitscore', 'slen', 'qlen']


def read_blast_tabular(path_or_text, source=None):
    '''
    Read BLAST tabular output into the standard table, one row per HSP.

    Parameters
    ----------
    path_or_text : str or path-like or file-like
        A file written with ``-outfmt`` :data:`BLAST_OUTFMT` (psiblast, or
        blast_formatter on an ASN.1 archive), or its text. Text is told from a
        path by being empty or containing a newline or a tab. ``.gz`` files are decompressed.
    source : str, optional
        Database searched, stored in the ``source`` column.

    Rows keep the order of the report. ``iteration`` comes from the
    ``# Iteration: N`` comment that precedes the rows (1 when there is none,
    and restarted at each ``# PSIBLAST``/``# BLAST*`` header, so queries of
    the same report don't inherit each other's rounds). A report without hits
    gives an empty standard table.

    ``cov`` is the fraction of the query the HSP spans, ``program`` is
    'psiblast' and ``talilen`` is ``end - start + 1``.
    '''
    import os
    import re

    if hasattr(path_or_text, 'read'):
        text = path_or_text.read()
    elif isinstance(path_or_text, (str, bytes)) and (not path_or_text or '\n' in path_or_text or '\t' in path_or_text):
        text = path_or_text
    else:
        path = os.fspath(path_or_text)
        if path.endswith('.gz'):
            import gzip
            with gzip.open(path, 'rt') as handle:
                text = handle.read()
        else:
            with open(path) as handle:
                text = handle.read()
    if isinstance(text, bytes):
        text = text.decode()

    rows, iteration = [], 1
    for line in text.splitlines():
        if line.startswith('#'):
            if re.match(r'#\s*(PSI)?(BLAST|TBLAST)', line):
                iteration = 1
            else:
                found = re.match(r'#\s*Iteration:\s*(\d+)', line)
                if found:
                    iteration = int(found.group(1))
            continue
        if not line.strip() or line.strip() == 'Search has CONVERGED!':
            # psiblast writes this line among the rows when it converges
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) != len(_TABULAR_NAMES):
            raise ValueError(f'read_blast_tabular: expected {len(_TABULAR_NAMES)} tab-separated '
                             f'fields ({BLAST_TABULAR_FIELDS}), found {len(fields)}: {line!r}')
        rows.append(fields + [iteration])
    if not rows:
        return empty()

    raw = pd.DataFrame(rows, columns=_TABULAR_NAMES + ['iteration'])
    for column in _TABULAR_NAMES[2:] + ['iteration']:
        raw[column] = pd.to_numeric(raw[column])
    start, end = raw['sstart'], raw['send']
    std = pd.DataFrame({
        'sequence':  raw['saccver'],
        'model':     raw['qaccver'],
        'start':     start,
        'end':       end,
        'evalue':    raw['evalue'],
        'cov':       ((raw['qend'] - raw['qstart'] + 1) / raw['qlen']).clip(upper=1.0),
        'qstart':    raw['qstart'],
        'qend':      raw['qend'],
        'iteration': raw['iteration'],
        'score':     raw['bitscore'],
        'talilen':   end - start + 1,
        'qlen':      raw['qlen'],
        'estart':    start,
        'eend':      end,
        'program':   'psiblast',
        'tlen':      raw['slen'],
        'source':    source,
    })
    return validate(std)


def blast2table_evalue(evalue):
    '''
    Round E-values the way blast2table prints them: ``%.1f`` above 1, 0 for 0
    and ``%.3e`` otherwise.

    blast2table rounds *before* it compares E-values with the inclusion
    threshold and with each other, so collapse_hsps does the same. BLAST's
    tabular output already has three significant digits, so this only changes
    values above 1.
    '''
    values = np.asarray(evalue, dtype='float64')
    out = np.empty(values.shape, dtype='float64')
    for i, value in enumerate(values.ravel()):
        if value == 0 or np.isnan(value):
            out.flat[i] = value
        elif value > 1:
            out.flat[i] = float(f'{value:.1f}')
        else:
            out.flat[i] = float(f'{value:.3e}')
    return out


def collapse_hsps(df, threshold=0.01, overlap=0.4):
    '''
    Merge redundant HSPs the way blast2table's default mode does.

    A standard table with one row per HSP, in report order (what
    :func:`read_blast_tabular` returns), becomes one with one row per region
    of each (``model``, ``sequence``) pair. Each HSP is compared with the
    regions already kept for its pair, in order; the first one that it
    overlaps *on both axes* (``sequence`` and ``model`` coordinates) is the
    one it is merged into, and later ones are not looked at:

      * if the kept region's E-value is below ``threshold`` (it was already
        significant), its coordinates are widened on each axis where the new
        HSP is longer (``start`` and ``end`` together, ``qstart`` and
        ``qend`` together). E-value, score and iteration stay the ones of
        the first significant HSP;
      * otherwise, if the new HSP has a lower E-value, it replaces the kept
        region entirely (E-value, score, iteration and coordinates);
      * otherwise the new HSP is dropped.

    A new HSP that overlaps nothing becomes a new region. Regions come out
    in the order in which they were first seen, even if they were replaced
    later, which is the order of blast2table's rows.

    Parameters
    ----------
    threshold : float
        Inclusion threshold. blast2table's default is 0.01.
    overlap : float or int
        A float is a fraction: two HSPs overlap on an axis when the shared
        length is above ``overlap`` times the length of each of them, that
        is, above that fraction of the *longer* one. An int is a number of
        residues: they overlap when they share more than that. blast2table's
        default is 0.4.

    Notes
    -----
    Which rows come out, their ``start``/``end``/``qstart``/``qend`` and their
    order follow blast2table exactly. To get the same decisions, E-values are
    rounded with :func:`blast2table_evalue` before they are compared with the
    threshold and with each other, so ties fall in the same place.

    The values reported are the correct ones, and differ from blast2table's:

      * ``evalue`` and ``score`` are the unrounded values of the HSP the
        region's E-value comes from (blast2table prints E-values with 1-2
        significant digits and rounds bit scores);
      * ``iteration`` is that HSP's true round (blast2table reports 1 with
        BioPerl's text parser);
      * ``cov`` is recomputed from the final ``qstart``/``qend`` and ``qlen``
        when ``qlen`` is known (blast2table keeps the coverage of the
        region's first HSP after widening it);
      * ``talilen``, ``estart`` and ``eend`` follow ``start`` and ``end``.

    Other columns come from the HSP that the region's E-value, score and
    iteration come from.

    Everything is computed per pair on small arrays, in pandas-version
    independent code (no ``groupby().apply``).
    '''
    df = validate(df)
    if df.empty:
        return df.reset_index(drop=True)
    if df['qstart'].isna().any() or df['qend'].isna().any():
        raise ValueError('collapse_hsps needs qstart and qend in every row')

    fraction = isinstance(overlap, float)
    ev = blast2table_evalue(df['evalue'].to_numpy())
    hs, he = df['start'].to_numpy(), df['end'].to_numpy()
    qs, qe = df['qstart'].to_numpy(dtype='int64'), df['qend'].to_numpy(dtype='int64')

    def same(a, b, c, d):
        '''Do spans (a, b) and (c, d) overlap more than allowed?'''
        shared = min(b, d) - max(a, c) + 1
        if fraction:
            return shared > overlap * (b - a + 1) and shared > overlap * (d - c + 1)
        return shared > overlap

    regions = []   # [first row, source row, hit start, hit end, query start, query end, evalue]
    for rows in df.groupby(['model', 'sequence'], sort=False).indices.values():
        kept = []
        for i in rows:
            i = int(i)
            for region in kept:
                if same(region[2], region[3], hs[i], he[i]) and same(region[4], region[5], qs[i], qe[i]):
                    if region[6] < threshold:
                        if he[i] - hs[i] > region[3] - region[2]:
                            region[2], region[3] = int(hs[i]), int(he[i])
                        if qe[i] - qs[i] > region[5] - region[4]:
                            region[4], region[5] = int(qs[i]), int(qe[i])
                    elif ev[i] < region[6]:
                        region[1:] = [i, int(hs[i]), int(he[i]), int(qs[i]), int(qe[i]), ev[i]]
                    break
            else:
                region = [i, i, int(hs[i]), int(he[i]), int(qs[i]), int(qe[i]), ev[i]]
                kept.append(region)
                regions.append(region)

    regions.sort(key=lambda region: region[0])
    table = np.array([region[1:] for region in regions], dtype='float64')
    out = df.iloc[table[:, 0].astype('int64')].reset_index(drop=True)
    out['start'], out['end'] = table[:, 1].astype('int64'), table[:, 2].astype('int64')
    out['qstart'] = pd.array(table[:, 3].astype('int64'), dtype='Int64')
    out['qend'] = pd.array(table[:, 4].astype('int64'), dtype='Int64')
    out['estart'], out['eend'] = out['start'], out['end']
    out['talilen'] = out['end'] - out['start'] + 1
    known = out['qlen'].notna()
    cov = (out['qend'] - out['qstart'] + 1).astype('float64') / out['qlen'].astype('float64')
    out['cov'] = cov.clip(upper=1.0).where(known, out['cov'])
    return validate(out)
