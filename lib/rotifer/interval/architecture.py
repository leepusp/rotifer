__doc__ = """
Domain architectures from tables of sequence regions.

A domain table lists regions found in sequences (one row per region:
sequence, domain, start, end, E-value), as written by ``hmmer2table``,
``blast2table -s`` or ``phobius2table``. For each sequence, the best set
of non-overlapping regions is its domain architecture, written compactly
as domain names joined by ``+``, with nested regions in brackets::

    ID           architecture            regions                          evalue
    NP_414543.1  AA_kinase+ACT_9         2..284&AA_kinase,319..378&ACT_9  1.2e-45&AA_kinase,...

Selection: regions are ranked by E-value (or by length, see
``by_length``), then by length and start; each region is kept unless it
overlaps a better region already kept. With a fractional ``overlap``
(default 0.4) two regions conflict when their overlap exceeds that
fraction of either region; with an integer ``overlap``, when it exceeds
that number of residues. A kept region that lies inside another is
written inside brackets after it, e.g. ``A[B]``.

``domain2architecture`` and ``architecture2table`` build and expand
architectures as DataFrames; the ``parse_*`` and ``iter_*`` functions are
the streaming versions used by the commands of the same names.
"""

import re
from functools import cmp_to_key
from rotifer.io import split_fields, written_as_fraction

COLUMNS = ['ID', 'architecture', 'regions', 'evalue']
"""Columns of an architecture table."""

TABLE_COLUMNS = ['ID', 'domain', 'start', 'end', 'evalue']
"""Columns of an expanded architecture table (one row per region)."""

INPUT_COLUMNS = ['sequence', 'domain', 'start', 'end', 'evalue']
"""Default columns of domain tables (the first five columns)."""

DOMAIN_DELIMITER = r'[~+\[\]]+'

def _number(text):
    """
    Numeric value of a string, as Perl would see it (0 if missing).
    """
    if text is None or text == '':
        return 0.0
    try:
        return float(text)
    except ValueError:
        match = re.match(r'\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?', text)
        return float(match.group(0)) if match else 0.0

def _cmp(a, b):
    return (a > b) - (a < b)

def parse_regions(handle, data=None, evalue=0.01, target_cutoff=None, exclude=(), rename=None,
                  columns=INPUT_COLUMNS, header=False, delimiter='\t'):
    """
    Load regions from the lines of a domain table into ``data``, a dict
    mapping each sequence to its list of regions.

    Parameters
    ----------
    handle : iterable of str
        Lines of a TAB-separated domain table.
    data : dict, optional
        Dictionary to extend (regions of several tables can be merged).
    evalue : float, default 0.01
        Maximum E-value of accepted regions.
    target_cutoff : dict, optional
        Domain-specific maximum E-values, overriding ``evalue``.
    exclude : collection of str
        Domains to ignore.
    rename : dict, optional
        Map of old to new domain names.
    columns : list of str
        Names of the input columns; must include sequence, domain, start,
        end and (optionally) evalue.
    header : bool, default False
        Read column names from the first line instead of ``columns``.
    delimiter : str
        Regular expression separating columns.

    Returns
    -------
    dict
    """
    data = {} if data is None else data
    target_cutoff = target_cutoff or {}
    rename = rename or {}
    index = None if header else {name: i for i, name in enumerate(columns)}
    for line in handle:
        line = line.rstrip('\n')
        if index is None:
            index = {name: i for i, name in enumerate(re.split(r'\t+', line))}
            continue
        row = split_fields(line, delimiter)
        cell = lambda name: row[index[name]] if name in index and index[name] < len(row) else None
        sequence, domain = cell('sequence') or '', cell('domain') or ''  # Perl: undef is ''
        domain = rename.get(domain, domain)
        start, end, score = cell('start'), cell('end'), cell('evalue')
        cutoff = float(target_cutoff[domain]) if domain in target_cutoff else float(evalue)
        if _number(score) > cutoff or domain in exclude:
            continue
        region = {'domain': domain, 'start': start, 'end': end,
                  'length': _number(end) - _number(start) + 1}
        if score is not None:
            region['evalue'] = score
        data.setdefault(sequence, []).append(region)
    return data

def _rank(by_length):
    """
    Sort key ranking regions, best first (a port of the Perl comparator).
    """
    def by_evalue(a, b):
        if 'evalue' in a:
            c = _cmp(_number(a['evalue']), _number(b.get('evalue')))
            if c:
                return c
        return _cmp(b['length'], a['length']) or _cmp(_number(a['start']), _number(b['start']))
    def longest(a, b):
        c = _cmp(b['length'], a['length']) or _cmp(_number(a['start']), _number(b['start']))
        if c or 'evalue' not in a:
            return c
        return _cmp(_number(a['evalue']), _number(b.get('evalue')))
    return cmp_to_key(longest if by_length else by_evalue)

def _conflict(region, better, overlap):
    overlap_length = min(_number(region['end']), _number(better['end'])) \
        - max(_number(region['start']), _number(better['start'])) + 1
    if written_as_fraction(overlap):
        fraction = float(overlap)
        return overlap_length > fraction * better['length'] or overlap_length > fraction * region['length']
    return overlap_length > float(overlap)

def _contains(a, b):
    return _number(a['start']) <= _number(b['start']) and _number(a['end']) >= _number(b['end'])

def _assemble(ordered, connector):
    """
    Write regions sorted by coordinates as an architecture string,
    nesting regions inside the closest region that contains them.
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

def select_regions(regions, overlap='0.4', by_length=False, ignore_overlap=()):
    """
    Select the best non-overlapping regions of one sequence.

    Parameters
    ----------
    regions : list of dict
        Regions of one sequence, as loaded by ``parse_regions``.
    overlap : str or number, default '0.4'
        Maximum overlap between conflicting regions: a fraction if written
        with a decimal point (or as ``1e-2``), else a number of residues.
    by_length : bool, default False
        Prefer longer regions over regions with lower E-values.
    ignore_overlap : collection of str
        Domains never discarded for overlapping others (``all`` for every
        domain).

    Returns
    -------
    list of dict
        The selected regions, sorted by start and end.
    """
    ranked = sorted(regions, key=_rank(by_length))
    best = ranked[:1]
    for region in ranked[1:]:
        keep = True
        for better in best:
            if 'all' in ignore_overlap:
                continue
            if region['domain'] in ignore_overlap:
                break
            if better['domain'] in ignore_overlap:
                continue
            if _conflict(region, better, overlap):
                keep = False
                break
        if keep:
            best.append(region)
    best.sort(key=lambda r: (_number(r['start']), _number(r['end'])))
    return best

def build_architecture(regions, overlap='0.4', by_length=False, ignore_overlap=(), connector='+'):
    """
    Write the architecture of one sequence from its best regions (see
    ``select_regions`` for the parameters).

    Parameters
    ----------
    connector : str, default '+'
        String joining domain names.

    Returns
    -------
    tuple of str
        Architecture, regions (``start..end&domain``) and E-values
        (``evalue&domain``) of the selected regions.
    """
    best = select_regions(regions, overlap=overlap, by_length=by_length, ignore_overlap=ignore_overlap)
    return (_assemble(best, connector),
            ','.join(f"{r['start']}..{r['end']}&{r['domain']}" for r in best),
            ','.join(f"{r.get('evalue', '')}&{r['domain']}" for r in best))

def iter_architectures(data, **kwargs):
    """
    Yield ``[ID, architecture, regions, evalue]`` for each sequence in
    ``data`` (see ``parse_regions``), sorted by sequence identifier.
    Keyword arguments go to ``build_architecture``.
    """
    for sequence in sorted(data):
        yield [sequence, *build_architecture(data[sequence], **kwargs)]

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

def _frame_lines(df):
    """
    Domain table lines from a DataFrame. Columns are taken by name
    (sequence, domain or model, start, end, evalue) when present, so
    tables from different parsers can be concatenated; otherwise the first
    five columns are used.
    """
    import pandas as pd
    domain = df['domain'] if 'domain' in df else None
    if 'model' in df:
        domain = df['model'] if domain is None else domain.fillna(df['model'])
    if domain is not None and all(name in df for name in ('sequence', 'start', 'end')):
        evalue = df['evalue'] if 'evalue' in df else pd.Series('', index=df.index)
        df = pd.DataFrame({'sequence': df['sequence'], 'domain': domain, 'start': df['start'],
                           'end': df['end'], 'evalue': evalue})
    return df.iloc[:, :5].astype(str).agg('\t'.join, axis=1)

def domain2architecture(tables, evalue=0.01, overlap='0.4', by_length=False, ignore_overlap=(),
                        target_cutoff=None, exclude=(), rename=None, connector='+'):
    """
    Build the domain architecture of each sequence in domain tables.

    Parameters
    ----------
    tables : str, list of str, file object or pandas.DataFrame
        Domain table file(s) (TAB-separated, columns sequence, domain,
        start, end, evalue, no header), an open table, or a DataFrame whose
        first five columns are sequence, domain, start, end and evalue.
    evalue, target_cutoff, exclude, rename
        Region filters, see ``parse_regions``.
    overlap, by_length, ignore_overlap, connector
        Selection rules, see ``build_architecture``.

    Returns
    -------
    pandas.DataFrame with ``COLUMNS``

    Examples
    --------
    >>> from rotifer.io.hmmer import read_hmmer
    >>> from rotifer.interval.architecture import domain2architecture
    >>> arch = domain2architecture(read_hmmer("seqs.hmmscan.out"), evalue=0.0101)
    """
    import pandas as pd
    from rotifer.io import each_source
    options = dict(evalue=evalue, target_cutoff=target_cutoff, exclude=exclude, rename=rename)
    data = {}
    if isinstance(tables, pd.DataFrame):
        parse_regions(_frame_lines(tables), data, **options)
    else:
        for name, handle in each_source(tables):
            parse_regions(handle, data, **options)
    rows = iter_architectures(data, overlap=overlap, by_length=by_length,
                              ignore_overlap=ignore_overlap, connector=connector)
    return pd.DataFrame(list(rows), columns=COLUMNS)

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
