__doc__ = """
Domain architectures from tables of sequence regions.

A domain table lists regions found in sequences (one row per region:
sequence, domain, start, end, E-value), as written by ``hmmer2table``,
``blast2table -s`` or ``phobius2table``. For each sequence, the best set
of non-overlapping regions is its domain architecture, written compactly
as domain names joined by ``+``, with nested regions in brackets::

    ID           architecture            regions                          evalue
    NP_414543.1  AA_kinase+ACT_9         2..284&AA_kinase,319..378&ACT_9  1.2e-45&AA_kinase,...

Functions:

- ``select_regions`` picks the best non-overlapping regions of each
  sequence in a DataFrame. Every other function here uses it.
- ``read_regions`` loads domain table files into a DataFrame and applies
  the E-value, exclusion and renaming filters of ``domain2architecture``.
- ``domain2architecture`` builds the architecture table from files or a
  DataFrame; ``architecture2table`` expands it back into one row per
  region.
- ``iter_architecture_table`` and ``region_table`` are the line-based
  pieces behind the ``architecture2table`` command.
- ``PRESETS`` holds the column names and ranking rules of common tables,
  to pass as ``select_regions(df, **PRESETS['hmmer'])``.

Selection: regions are ranked by E-value (or by length, see
``by_length``), then by length and start; each region is kept unless it
overlaps a better region already kept. With a fractional ``overlap``
(default 0.4) two regions conflict when their overlap exceeds that
fraction of the shorter region; with an integer ``overlap``, when it
exceeds that number of residues. A kept region that lies inside another
is written inside brackets after it, e.g. ``A[B]``.
"""

import re
from rotifer.io import split_fields, written_as_fraction

_LEADING_NUMBER = r'^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)'

COLUMNS = ['ID', 'architecture', 'regions', 'evalue']
"""Columns of an architecture table."""

TABLE_COLUMNS = ['ID', 'domain', 'start', 'end', 'evalue']
"""Columns of an expanded architecture table (one row per region)."""

INPUT_COLUMNS = ['sequence', 'domain', 'start', 'end', 'evalue']
"""Default columns of domain tables (the first five columns)."""

DOMAIN_DELIMITER = r'[~+\[\]]+'

PRESETS = {
    'hmmer': {
        'reference': 'sequence', 'domain': 'model', 'start': 'estart', 'end': 'eend',
        'criteria': {'evalue': True, 'region_length': False}, 'overlap': 0.4,
    },
    'hhsuite': {
        'reference': 'query', 'domain': 'hit', 'start': 'qstart', 'end': 'qend',
        'criteria': {'sum_probs': False, 'probability': False, 'score': False, 'region_length': False},
        'overlap': 0.4,
    },
    'archtsv': {
        'reference': 'ID', 'domain': 'domain', 'start': 'start', 'end': 'end',
        'criteria': {'evalue': True, 'region_length': False}, 'overlap': 0.4,
    },
}
"""
Arguments of ``select_regions`` for common tables.

- ``hmmer``: ``hmmer2table`` or ``rotifer.io.search`` output, by the
  envelope coordinates (estart, eend), best E-value first.
- ``hhsuite``: HH-suite tables (query, hit, qstart, qend), largest
  sum_probs, probability and score first.
- ``archtsv``: expanded architecture tables (``TABLE_COLUMNS``).

Each value is a dict to unpack: ``select_regions(df, **PRESETS['hmmer'])``.
Keys override one another as usual, e.g.
``select_regions(df, **{**PRESETS['hmmer'], 'overlap': 10})``.
"""

def _numbers(values):
    """
    Convert values to floats the way Perl reads numbers from text.

    Numbers are kept. A string is read from its leading number (optional
    blanks, sign, digits with an optional decimal part and an optional
    exponent), so ``'1e-5'`` is 1e-05, ``'12abc'`` is 12.0 and ``'abc'``
    is 0.0. Missing values (None, NaN, ``''``) are 0.0. The former Perl
    tools compared E-values and coordinates this way, and the golden
    tests depend on it.

    Parameters
    ----------
    values : pandas.Series or array-like
        Values to convert. It is not modified.

    Returns
    -------
    numpy.ndarray
        float64 array of the same length, in the same order.

    Examples
    --------
    >>> _numbers(['1e-5', '12abc', 'abc', '', None, 3]).tolist()
    [1e-05, 12.0, 0.0, 0.0, 0.0, 3.0]
    """
    import numpy as np
    import pandas as pd
    series = values if isinstance(values, pd.Series) else pd.Series(list(values), dtype=object)
    if pd.api.types.is_numeric_dtype(series.dtype) or pd.api.types.is_bool_dtype(series.dtype):
        return np.nan_to_num(series.to_numpy(dtype=float, na_value=np.nan), nan=0.0)
    series = series.astype(object)
    result = pd.to_numeric(series, errors='coerce').to_numpy(dtype=float, copy=True)
    unread = np.isnan(result) & series.notna().to_numpy() & (series.to_numpy() != '')
    if unread.any():
        leading = series[unread].astype(str).str.extract(_LEADING_NUMBER)[0]
        result[unread] = pd.to_numeric(leading, errors='coerce').to_numpy(dtype=float)
    return np.nan_to_num(result, nan=0.0)

def select_regions(df, reference='sequence', domain=None, start='start', end='end',
                   criteria=None, by_length=False, overlap=0.4, nested=True, ignore_overlap=()):
    """
    Select the best set of non-overlapping regions of each sequence.

    Each row of ``df`` is a region (a domain hit, a transmembrane segment,
    ...) of the sequence named in the ``reference`` column(s), from
    ``start`` to ``end`` (inclusive). Regions are ranked within each
    sequence, best first, and taken in that order: a region is kept
    unless it conflicts with a region already kept (see ``overlap`` and
    ``nested``). The result is the same as checking each region against
    all better kept regions one at a time, but the work is done in
    whole-table passes, one per kept region of the sequence with the most
    kept regions.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per region. It is not modified. Coordinates and ranking
        columns may be numbers or text; text is read as numbers the way
        the former Perl tools did (``'1e-5'`` is 1e-05, missing values
        are 0), so tables read with ``dtype=str`` keep their text.
    reference : str or list of str, default 'sequence'
        Column(s) naming the sequence of each region. A single string is
        one column name. Missing values form their own group.
    domain : str, optional
        Column with the region names, used only by ``ignore_overlap``.
        None means ``'domain'`` when ``df`` has that column, else
        ``'model'`` (the name in ``hmmer2table`` and ``rotifer.io.search``
        tables).
    start, end : str, default 'start', 'end'
        Columns with the first and last position of each region.
    criteria : dict, optional
        Ranking columns, in order of precedence, mapped to True to prefer
        the smallest values or False to prefer the largest. The name
        ``'region_length'`` means ``end - start + 1`` unless ``df`` has a
        column with that name, which is then used for ranking (overlaps
        always use the coordinates). After these columns, ties are broken
        by the smallest start, then by input order. None means the order
        of ``domain2architecture``: lowest ``evalue``, then longest, then
        smallest start (``evalue`` is skipped when ``df`` has no such
        column); with ``by_length=True``, longest, then smallest start,
        then lowest ``evalue``.
    by_length : bool, default False
        Rank longer regions first (see ``criteria``). Only used when
        ``criteria`` is None.
    overlap : int, float or str, default 0.4
        How much two regions may overlap before the worse one is dropped.
        A float is a fraction: two regions conflict when they share more
        than ``overlap`` times the length of the shorter one. An int is a
        number of residues: they conflict when they share more than
        ``overlap`` positions. A string is read as written: a fraction
        when it has a decimal point or an exponent (``'0.4'``,
        ``'1e-2'``), else a number of residues (``'10'``), as in the
        ``domain2architecture`` command. Both limits are exclusive.
    nested : bool, default True
        Whether a region that lies entirely inside a kept region may be
        kept too (only possible when the overlap limit allows it, e.g.
        with a number of residues larger than the inner region, or with
        ``ignore_overlap``). False drops every region that a kept region
        contains, whatever ``overlap`` says.
    ignore_overlap : collection of str, default ()
        Region names (values of ``domain``) that are always kept and never
        cause another region to be dropped. ``'all'`` keeps every region.

    Returns
    -------
    pandas.DataFrame
        The kept rows of ``df``, with all its columns, dtypes and index
        labels, sorted by the ``reference`` column(s), then by numeric
        start and end; regions with the same coordinates keep their rank
        order, best first, as in the former Perl tool. Empty, with the
        columns of ``df``, when ``df`` is empty.

    Raises
    ------
    KeyError
        If ``reference``, ``start``, ``end`` or a ``criteria`` column is
        not in ``df``, or ``domain`` is needed (``ignore_overlap`` names a
        domain other than ``'all'``) and not found.
    ValueError
        If both ``criteria`` and ``by_length=True`` are given.

    See Also
    --------
    domain2architecture : writes the kept regions of each sequence as one
        architecture string.
    PRESETS : arguments for hmmer, HH-suite and architecture tables.

    Examples
    --------
    A is kept first (lowest E-value); B overlaps A by 11 residues, more
    than 40% of B's 20, and is dropped; C does not overlap A:

    >>> import pandas as pd
    >>> df = pd.DataFrame({'sequence': ['p1', 'p1', 'p1'], 'domain': ['A', 'B', 'C'],
    ...                    'start': [1, 50, 70], 'end': [60, 69, 120],
    ...                    'evalue': [1e-10, 1e-5, 1e-3]})
    >>> select_regions(df)
      sequence domain  start  end        evalue
    0       p1      A      1   60  1.000000e-10
    2       p1      C     70  120  1.000000e-03

    Allowing 15 shared residues keeps B as well:

    >>> select_regions(df, overlap=15)['domain'].tolist()
    ['A', 'B', 'C']
    """
    import numpy as np
    import pandas as pd
    if criteria is not None and by_length:
        raise ValueError("give either criteria or by_length=True, not both")
    names = [reference] if isinstance(reference, str) else list(reference)
    ignored = [ignore_overlap] if isinstance(ignore_overlap, str) else list(ignore_overlap)
    for name in [*names, start, end]:
        if name not in df.columns:
            raise KeyError(name)
    if criteria is None:
        has_evalue = 'evalue' in df.columns
        if by_length:
            criteria = {'region_length': False}
        else:
            criteria = {'evalue': True, 'region_length': False} if has_evalue else {'region_length': False}
    for name in criteria:
        if name not in df.columns and name != 'region_length':
            raise KeyError(name)
    if set(ignored) - {'all'}:
        if domain is None:
            domain = 'domain' if 'domain' in df.columns else 'model'
        if domain not in df.columns:
            raise KeyError(domain)
    if df.empty:
        return df.iloc[:0].copy()

    first, last = _numbers(df[start]), _numbers(df[end])
    length = last - first + 1
    ranking = []
    for name, smallest in criteria.items():
        values = length if name == 'region_length' and name not in df.columns else _numbers(df[name])
        ranking.append(values if smallest else -values)
    if by_length and 'evalue' in df.columns:
        ranking = [ranking[0], first, _numbers(df['evalue'])]
    codes = None
    for name in names:
        column, uniques = pd.factorize(df[name], sort=True, use_na_sentinel=False)
        codes = column if codes is None else pd.factorize(codes * len(uniques) + column, sort=True)[0]
    order = np.lexsort(tuple([codes, *ranking, first][::-1]))
    group, s, e, size = codes[order], first[order], last[order], length[order]

    fraction = isinstance(overlap, (float, np.floating)) or (isinstance(overlap, str) and written_as_fraction(overlap))
    limit = float(overlap)
    dropped = np.zeros(len(order), dtype=bool)
    if 'all' in ignored:
        immune = np.ones(len(order), dtype=bool)
    elif ignored:
        immune = pd.Series(df[domain].to_numpy()[order]).isin(ignored).to_numpy()
    else:
        immune = np.zeros(len(order), dtype=bool)
    undecided = np.flatnonzero(~immune)
    while undecided.size:
        groups = group[undecided]
        head = np.ones(undecided.size, dtype=bool)
        head[1:] = groups[1:] != groups[:-1]
        kept, rest = undecided[head], undecided[~head]
        if not rest.size:
            break
        better = kept[np.searchsorted(group[kept], group[rest])]
        shared = np.minimum(e[better], e[rest]) - np.maximum(s[better], s[rest]) + 1
        limits = limit * np.minimum(size[better], size[rest]) if fraction else limit
        conflict = shared > limits
        if not nested:
            conflict |= (s[better] <= s[rest]) & (e[better] >= e[rest])
        dropped[rest[conflict]] = True
        undecided = rest[~conflict]
    keep = np.flatnonzero(~dropped)
    final = np.lexsort((keep, e[keep], s[keep], group[keep]))
    return df.iloc[order[keep][final]]

def read_regions(tables, evalue=0.01, target_cutoff=None, exclude=(), rename=None,
                 columns=INPUT_COLUMNS, header=False, delimiter='\t'):
    """
    Load domain tables into one DataFrame of regions, keeping their text.

    Each line of a table is a region. Lines are split on ``delimiter``
    (dropping trailing empty fields), named by ``columns`` (or by the
    first line when ``header`` is true), and filtered as the
    ``domain2architecture`` command does: domains are renamed first, then
    a region is dropped when its E-value is above its cutoff or its
    (renamed) domain is in ``exclude``.

    Parameters
    ----------
    tables : str, list of str or file object
        Domain table file name(s), or an open table (anything
        ``rotifer.io.each_source`` accepts). Compressed files are read
        too.
    evalue : float, default 0.01
        Maximum E-value of accepted regions (inclusive).
    target_cutoff : dict, optional
        Domain-specific maximum E-values, overriding ``evalue``. Keys are
        domain names after renaming.
    exclude : collection of str, default ()
        Domains to drop (names after renaming).
    rename : dict, optional
        Map of old to new domain names.
    columns : list of str, default INPUT_COLUMNS
        Names of the input columns, in order. Must include ``sequence``,
        ``domain``, ``start`` and ``end``; ``evalue`` is optional. Other
        names are read and dropped.
    header : bool, default False
        Read column names from the first line of each table (split on
        TABs) instead of ``columns``.
    delimiter : str, default '\\t'
        Regular expression separating columns.

    Returns
    -------
    pandas.DataFrame
        Columns ``INPUT_COLUMNS`` (sequence, domain, start, end, evalue),
        all strings exactly as written in the tables, one row per accepted
        region in file order, with a 0..n-1 index. A cell missing from a
        short line is ``''`` for sequence and domain and None for the
        others; evalue is None when the tables have no evalue column. An
        E-value that is missing or not a number counts as 0 and passes
        any cutoff.

    Examples
    --------
    >>> import io
    >>> table = io.StringIO("p1\\tA\\t1\\t60\\t1e-10\\np1\\tB\\t50\\t69\\t0.5\\n")
    >>> read_regions(table)
      sequence domain start end evalue
    0       p1      A     1  60  1e-10
    """
    import pandas as pd
    from rotifer.io import each_source
    splitter = re.compile(delimiter)
    wanted = list(INPUT_COLUMNS)
    cells = {name: [] for name in wanted}
    for name, handle in each_source(tables):
        index = None if header else {name: i for i, name in enumerate(columns)}
        for line in handle:
            line = line.rstrip('\n')
            if index is None:
                index = {name: i for i, name in enumerate(re.split(r'\t+', line))}
                continue
            row = splitter.split(line)
            while row and row[-1] == '':
                row.pop()
            for name in wanted:
                i = index.get(name)
                cells[name].append(row[i] if i is not None and i < len(row) else None)
    frame = pd.DataFrame({name: pd.Series(values, dtype=object) for name, values in cells.items()})
    return _filter_regions(frame, evalue, target_cutoff, exclude, rename)

def _filter_regions(frame, evalue, target_cutoff, exclude, rename):
    """
    Rename domains and drop the regions above their E-value cutoff or in
    the excluded domains.

    Parameters
    ----------
    frame : pandas.DataFrame
        Regions with ``INPUT_COLUMNS`` (text or None). Missing sequence
        and domain names are read as ``''``.
    evalue, target_cutoff, exclude, rename
        See ``read_regions``.

    Returns
    -------
    pandas.DataFrame
        The accepted rows, as object columns with a 0..n-1 index.
    """
    import numpy as np
    import pandas as pd
    frame = frame.astype(object)
    for name in ('sequence', 'domain'):
        frame[name] = frame[name].where(frame[name].notna(), '')
    if rename:
        frame['domain'] = frame['domain'].map(lambda name: rename.get(name, name))
    domain = frame['domain'].to_numpy(dtype=object)
    cutoff = np.full(len(frame), float(evalue))
    for name, value in (target_cutoff or {}).items():
        cutoff[domain == name] = float(value)
    drop = _numbers(frame['evalue']) > cutoff
    if len(exclude):
        drop |= frame['domain'].isin(list(exclude)).to_numpy()
    return frame[~drop].reset_index(drop=True)

def _contains(a, b):
    """
    Tell whether region ``a`` contains region ``b``.

    Parameters
    ----------
    a, b : dict
        Regions with numeric ``'start'`` and ``'end'`` values.

    Returns
    -------
    bool
        True when ``a`` starts at or before ``b`` and ends at or after it.
    """
    return a['start'] <= b['start'] and a['end'] >= b['end']

def _assemble(ordered, connector):
    """
    Write regions as an architecture string, nesting regions inside the
    closest earlier region that contains them.

    Parameters
    ----------
    ordered : list of dict
        Kept regions of one sequence sorted by start and end, each with
        ``'domain'`` (str) and numeric ``'start'`` and ``'end'``.
    connector : str
        String joining domain names.

    Returns
    -------
    str
        Domain names joined by ``connector``, with each region that lies
        inside another written in brackets after it, e.g. ``A[B+C]+D``.
        This reproduces the former Perl tool, including its empty
        brackets (``A[]``) when every contained region is also contained
        by a later region.
    """
    names, inside, contains = [], {}, {}
    for i, region in enumerate(ordered):
        contained = [j for j in range(i + 1, len(ordered)) if _contains(region, ordered[j])]
        names.append(region['domain'])
        if contained:
            contains[i] = contained
            for j in contained:
                inside[j] = i
    for i in range(len(names) - 1, -1, -1):
        if i in contains:
            nested = [names[j] for j in contains[i] if inside[j] == i]
            names[i] += '[' + connector.join(nested) + ']'
    return connector.join(names[i] for i in range(len(names)) if i not in inside)

def iter_architecture_table(handle, delimiter='\t', domain_delimiter=DOMAIN_DELIMITER):
    """
    Expand the lines of an architecture table into one row per region:
    ``[ID, domain, start, end, evalue]``.
    """
    header = True
    for line in handle:
        fields = split_fields(line.rstrip('\n'), delimiter)
        if header and fields and fields[0] == 'ID':
            yield list(TABLE_COLUMNS)
            header = False
            continue
        regions = []
        pattern = domain_delimiter
        for i in range(1, len(fields)):
            for row, value in enumerate(split_fields(fields[i], pattern)):
                if row == len(regions):
                    regions.append([])
                regions[row].append(value)
            if i <= 1:
                pattern = r',+'
        for region in regions:
            values = []
            for j in range(len(fields) - 1):
                value = region[j] if j < len(region) else ''
                value = re.sub(r'&.+', '', value.replace('..', '\t', 1), count=1)
                values.append(value)
            yield [fields[0], *values]

def domain2architecture(tables, evalue=0.01, overlap='0.4', by_length=False, ignore_overlap=(),
                        target_cutoff=None, exclude=(), rename=None, connector='+'):
    """
    Build the domain architecture of each sequence in domain tables.

    Regions are loaded and filtered by ``read_regions`` (or taken from a
    DataFrame), the best non-overlapping ones are chosen by
    ``select_regions`` with nested regions allowed, and the kept regions
    of each sequence are written as one row. The output is the same as the
    ``domain2architecture`` command's, which was checked against the
    former Perl tool.

    Parameters
    ----------
    tables : str, list of str, file object or pandas.DataFrame
        Domain table file(s) (TAB-separated, columns sequence, domain,
        start, end, evalue, no header), an open table, or a DataFrame.
        DataFrame columns are taken by name when it has ``sequence``,
        ``start``, ``end`` and ``domain`` or ``model`` (``domain`` wins,
        and ``model`` fills its missing values), so tables from different
        parsers can be concatenated; otherwise its first five columns are
        used. DataFrame values are converted to text first. It is not
        modified.
    evalue, target_cutoff, exclude, rename
        Region filters, see ``read_regions``.
    overlap, by_length, ignore_overlap
        Selection rules, see ``select_regions``. ``overlap`` is read as
        written, like the command's ``-o``: ``'0.4'`` is a fraction,
        ``'10'`` a number of residues.
    connector : str, default '+'
        String joining domain names.

    Returns
    -------
    pandas.DataFrame
        Columns ``COLUMNS``, all strings, one row per sequence with at
        least one kept region, sorted by sequence, with a 0..n-1 index:

        - ``ID``: the sequence.
        - ``architecture``: kept domains in coordinate order, joined by
          ``connector``, nested regions in brackets.
        - ``regions``: ``start..end&domain`` for each kept region, joined
          by commas.
        - ``evalue``: ``evalue&domain`` for each kept region (empty
          E-value when there is none).

        Empty, with these columns, when no region passes the filters.

    See Also
    --------
    architecture2table : expands this table into one row per region.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'sequence': ['p1', 'p1', 'p1'], 'domain': ['A', 'B', 'C'],
    ...                    'start': [1, 20, 70], 'end': [60, 40, 120],
    ...                    'evalue': [1e-10, 1e-5, 1e-3]})
    >>> domain2architecture(df, overlap='30')[['ID', 'architecture']]
       ID architecture
    0  p1       A[B]+C
    """
    import pandas as pd
    options = dict(evalue=evalue, target_cutoff=target_cutoff, exclude=exclude, rename=rename)
    if isinstance(tables, pd.DataFrame):
        regions = _filter_regions(_frame_regions(tables), **options)
    else:
        regions = read_regions(tables, **options)
    return architecture_table(regions, overlap=overlap, by_length=by_length,
                              ignore_overlap=ignore_overlap, connector=connector)

def _frame_regions(df):
    """
    Regions with ``INPUT_COLUMNS`` from a DataFrame, as text. Columns are
    taken by name (sequence, domain or model, start, end, evalue) when
    present, so tables from different parsers can be concatenated;
    otherwise the first five columns are used.

    Parameters
    ----------
    df : pandas.DataFrame
        Table of regions. It is not modified.

    Returns
    -------
    pandas.DataFrame
        Object columns ``INPUT_COLUMNS``; evalue is None when missing.
    """
    import pandas as pd
    domain = df['domain'] if 'domain' in df else None
    if 'model' in df:
        domain = df['model'] if domain is None else domain.fillna(df['model'])
    if domain is not None and all(name in df for name in ('sequence', 'start', 'end')):
        columns = {'sequence': df['sequence'], 'domain': domain, 'start': df['start'], 'end': df['end']}
        if 'evalue' in df:
            columns['evalue'] = df['evalue']
    else:
        columns = {name: df.iloc[:, i] for i, name in enumerate(INPUT_COLUMNS[:df.shape[1]])}
    frame = pd.DataFrame({name: column.astype(str).astype(object) for name, column in columns.items()}).reset_index(drop=True)
    for name in INPUT_COLUMNS:
        if name not in frame:
            frame[name] = pd.Series([None] * len(frame), dtype=object)
    return frame[INPUT_COLUMNS]

def architecture_table(regions, overlap='0.4', by_length=False, ignore_overlap=(), connector='+'):
    """
    Write the architecture of each sequence from a table of regions
    already loaded and filtered (see ``read_regions``).

    Parameters
    ----------
    regions : pandas.DataFrame
        Columns ``INPUT_COLUMNS`` as text, as returned by ``read_regions``.
    overlap, by_length, ignore_overlap, connector
        See ``domain2architecture``.

    Returns
    -------
    pandas.DataFrame
        Columns ``COLUMNS``, one row per sequence, sorted by sequence.
    """
    import pandas as pd
    best = select_regions(regions, reference='sequence', domain='domain', by_length=by_length,
                          overlap=overlap, nested=True, ignore_overlap=ignore_overlap)
    sequences = best['sequence'].tolist()
    domains = best['domain'].tolist()
    starts, ends = best['start'].tolist(), best['end'].tolist()
    evalues = best['evalue'].tolist()
    numbers = list(zip(_numbers(best['start']).tolist(), _numbers(best['end']).tolist()))
    rows, first = [], 0
    for i in range(1, len(sequences) + 1):
        if i < len(sequences) and sequences[i] == sequences[first]:
            continue
        span = range(first, i)
        ordered = [{'domain': domains[j], 'start': numbers[j][0], 'end': numbers[j][1]} for j in span]
        rows.append([sequences[first], _assemble(ordered, connector),
                     ','.join(f'{starts[j]}..{ends[j]}&{domains[j]}' for j in span),
                     ','.join(f"{'' if evalues[j] is None else evalues[j]}&{domains[j]}" for j in span)])
        first = i
    return pd.DataFrame(rows, columns=COLUMNS)

def architecture2table(architectures):
    """
    Expand an architecture table into one row per region.

    Parameters
    ----------
    architectures : pandas.DataFrame, str or file object
        Architecture table (``COLUMNS``), as returned by
        ``domain2architecture``, or the name of such a file.

    Returns
    -------
    pandas.DataFrame with ``TABLE_COLUMNS`` (start and end as integers,
    evalue as float)
    """
    import pandas as pd
    from rotifer.io import each_source
    if isinstance(architectures, pd.DataFrame):
        lines = ['\t'.join(COLUMNS)] + list(architectures[COLUMNS].astype(str).agg('\t'.join, axis=1))
        rows = list(iter_architecture_table(lines))
    else:
        rows = [row for name, handle in each_source(architectures) for row in iter_architecture_table(handle)]
    # The command writes "start..end" as "start<TAB>end" in one field: split it
    return region_table([row[0], row[1], *row[2].split('\t', 1), *row[3:]]
                        for row in rows if row != TABLE_COLUMNS and len(row) == 4 and '\t' in row[2])

def region_table(rows):
    """
    DataFrame with ``TABLE_COLUMNS`` from rows of strings (start and end
    become integers, evalue a float).
    """
    import pandas as pd
    df = pd.DataFrame(list(rows), columns=TABLE_COLUMNS)
    df[['start', 'end']] = df[['start', 'end']].astype(int)
    df['evalue'] = pd.to_numeric(df['evalue'], errors='coerce')
    return df
