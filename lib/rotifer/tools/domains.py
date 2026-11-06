__doc__ = """
Domain architecture of proteins, and counts of domains per protein.

These functions replace the shell functions that the ROTIFER Lmod module
used to define (``acc2pfam``, ``acc2profiledb``, ``dom_count``,
``dom_per_prot`` and ``distinct_dom_per_prot``); each one is also a
command of the same name.

``acc2pfam`` and ``acc2profiledb`` fetch protein sequences with
``rotifer.db.proteins`` (local database first, then NCBI) or read them
from a FASTA file, search them against Pfam with ``hmmscan``
or against a profile database with ``rpsblast``, predict signal peptides
and transmembrane segments with ``phobius``, and return the domain
architecture as one row per region (ID, domain, start, end, evalue).
"""

import os
import subprocess
import tempfile

PHOBIUS_EVALUE = '0.0101'
"""E-value given to Phobius regions, and cutoff for all regions."""

def _database(*path):
    import rotifer
    return os.path.join(rotifer.config['data'], *path)

def fetch_fasta(accessions, methods=('esl_sfetch', 'entrez')):
    """
    Fetch protein sequences as FASTA text with ``rotifer.db.proteins``:
    from the local sequence database first, then from NCBI.
    """
    from io import StringIO
    from Bio import SeqIO
    import rotifer.db
    accessions = [accessions] if isinstance(accessions, str) else list(accessions)
    records = rotifer.db.proteins(accessions, methods=list(methods), progress=False)
    fasta = StringIO()
    SeqIO.write(records, fasta, 'fasta')
    return fasta.getvalue()

def _run(command, stdout):
    with open(stdout, 'w') as out:
        subprocess.run(command, stdout=out, stderr=subprocess.DEVNULL, check=True)

def architecture_rows(region_rows, evalue=PHOBIUS_EVALUE):
    """
    Select the best regions of each sequence from domain table rows
    (lists of strings) and yield one row per selected region:
    ``[ID, domain, start, end, evalue]`` as strings, preceded by a header
    row. This is ``domain2architecture | architecture2table``.
    """
    from rotifer.interval.architecture import TABLE_COLUMNS, parse_regions, select_regions
    data = parse_regions(('\t'.join(row) for row in region_rows), evalue=float(evalue))
    yield list(TABLE_COLUMNS)
    for sequence in sorted(data):
        for region in select_regions(data[sequence]):
            yield [sequence, region['domain'], region['start'], region['end'], region.get('evalue', '')]

def pfam_rows(hmmscan_report, phobius_report):
    """
    Region rows (strings, with header) from hmmscan and Phobius reports.
    """
    from rotifer.io import open_input
    from rotifer.io.hmmer import parse_hmmer
    from rotifer.io.phobius import parse_phobius
    regions = []
    with open_input(hmmscan_report) as handle:
        regions.extend(row[:5] for row in parse_hmmer(handle, clean={'model': 'version'}))
    with open_input(phobius_report) as handle:
        regions.extend(parse_phobius(handle, evalue=PHOBIUS_EVALUE))
    return list(architecture_rows(regions))

def profiledb_rows(rpsblast_report, phobius_report):
    """
    Region rows (strings, with header) from RPS-BLAST and Phobius reports.
    """
    from rotifer.io import open_input
    from rotifer.io.blast import iter_blast
    from rotifer.io.phobius import parse_phobius
    regions = []
    with open_input(rpsblast_report) as handle:
        regions.extend(row[:5] for row in iter_blast(handle, swap=True, clean=['profiledb']))
    with open_input(phobius_report) as handle:
        regions.extend(parse_phobius(handle, evalue=PHOBIUS_EVALUE))
    return list(architecture_rows(regions))

def _search(search, fasta, accessions):
    """
    Run a search program and Phobius on sequences; return region rows.
    """
    with tempfile.TemporaryDirectory() as tmp:
        sequences = os.path.join(tmp, 'query.fa')
        with open(sequences, 'w') as out:
            out.write(fasta if fasta is not None else fetch_fasta(accessions))
        report, phobius = os.path.join(tmp, 'search.out'), os.path.join(tmp, 'phobius.out')
        rows = search(sequences, report)
        _run(['phobius', sequences], phobius)
        return rows(report, phobius)

def pfam_search_rows(accessions=None, fasta=None, database=None, cpu=4):
    """
    Region rows (strings, with header) for ``acc2pfam``.
    """
    database = database or _database('pfam', 'Pfam-A.hmm')
    def search(sequences, report):
        _run(['hmmscan', '--cpu', str(cpu), database, sequences], report)
        return pfam_rows
    return _search(search, fasta, accessions)

def profiledb_search_rows(accessions=None, fasta=None, database=None):
    """
    Region rows (strings, with header) for ``acc2profiledb``.
    """
    database = database or _database('rpsdb', 'allprofiles')
    def search(sequences, report):
        _run(['rpsblast', '-query', sequences, '-db', database], report)
        return profiledb_rows
    return _search(search, fasta, accessions)

def acc2pfam(accessions=None, fasta=None, database=None, cpu=4):
    """
    Pfam domains, signal peptides and transmembrane segments of proteins.

    Parameters
    ----------
    accessions : str or list of str, optional
        NCBI protein accessions to fetch.
    fasta : str, optional
        FASTA text to use instead of fetching sequences.
    database : str, optional
        Pfam HMM library (default: ``<rotifer data>/pfam/Pfam-A.hmm``).
    cpu : int, default 4
        Threads for hmmscan.

    Returns
    -------
    pandas.DataFrame with columns ID, domain, start, end, evalue

    Examples
    --------
    >>> from rotifer.tools.domains import acc2pfam
    >>> acc2pfam("NP_414543.1")
    """
    from rotifer.interval.architecture import region_table
    return region_table(pfam_search_rows(accessions, fasta, database, cpu)[1:])

def acc2profiledb(accessions=None, fasta=None, database=None):
    """
    Profile database domains, signal peptides and transmembrane segments
    of proteins (see ``acc2pfam``). The default database is
    ``<rotifer data>/rpsdb/allprofiles``.
    """
    from rotifer.interval.architecture import region_table
    return region_table(profiledb_search_rows(accessions, fasta, database)[1:])

def _regions(architectures):
    from rotifer.interval.architecture import architecture2table
    return architecture2table(architectures)

def dom_count(architectures):
    """
    Number of regions of each domain in each protein.

    Parameters
    ----------
    architectures : str or pandas.DataFrame
        Architecture table (file name or DataFrame), as written by
        ``domain2architecture``.

    Returns
    -------
    pandas.DataFrame with columns ID, domain, count
    """
    regions = _regions(architectures)
    return regions.groupby(['ID', 'domain'], sort=False).size().reset_index(name='count')

def dom_per_prot(architectures):
    """
    Number of regions in each protein (see ``dom_count``).

    Returns
    -------
    pandas.DataFrame with columns ID, count
    """
    regions = _regions(architectures)
    return regions.groupby('ID', sort=False).size().reset_index(name='count')

def distinct_dom_per_prot(architectures):
    """
    Number of distinct domains in each protein (see ``dom_count``).

    Returns
    -------
    pandas.DataFrame with columns ID, count
    """
    regions = _regions(architectures)
    return regions.groupby('ID', sort=False)['domain'].nunique().reset_index(name='count')
