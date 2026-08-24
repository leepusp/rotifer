"""
Access biological sequence and annotation databases.

This package provides a uniform, cursor based interface to remote
repositories (NCBI Entrez, the NCBI FTP site, the UniProt REST API)
and to local resources (indexed FASTA files, NCBI genome mirrors,
SQLite3 stores and the ETE toolkit taxonomy database).

Every data source is wrapped in a cursor class derived from
:class:`rotifer.db.core.BaseCursor`. Cursors share three access styles:

* dictionary-like access to a single entry, ``cursor[accession]``
* lazy iteration over many entries, ``fetchone``
* bulk retrieval of many entries, ``fetchall``

See :mod:`rotifer.db.ncbi` for the most complete set of cursors.

Configuration
-------------
The default path for local sequence databases is read from the user
configuration file ``~/.rotifer/etc/db.yml`` (key
``local_database_path``) and falls back to ``fadb/nr/nr`` under the
directory named by the ``ROTIFER_DATA`` environment variable.
"""

import os
import rotifer
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

# Configuration
config = loadConfig(__name__.replace('rotifer.',':'), defaults = {
    'local_database_path': os.path.join(GlobalConfig['data'],"fadb","nr","nr"),
})

# FUNCTIONS

def proteins(query, methods=['esl_sfetch','entrez'], local_database_path=config["local_database_path"], entrez_database="protein", batch_size=200, threads=None, tries=3, progress=True):
    """
    Fetch protein sequences from local or remote databases.

    Each retrieval method listed in `methods` is tried in order and
    receives only the identifiers that the previous methods could not
    resolve. Identifiers that remain unresolved after the last method
    are reported through the logging system.

    Parameters
    ----------
    query : list of str
        Sequence identifiers.
    methods : list of str, default ``['esl_sfetch', 'entrez']``
        Retrieval methods, tried in order. Supported keywords:

        ``esl_sfetch``
            Local indexed FASTA files, using
            :class:`rotifer.db.local.easel.FastaCursor`.
        ``entrez``
            NCBI Entrez, using
            :class:`rotifer.db.ncbi.entrez.FastaCursor`.
    local_database_path : str, optional
        Path to a sequence file, usually a FASTA file indexed by
        ``esl-sfetch``. Defaults to the ``local_database_path`` entry
        of the ``rotifer.db`` configuration.
    entrez_database : str, default 'protein'
        Name of the NCBI Entrez sequence database.
    batch_size : int, default 200
        Number of sequences to retrieve per batch.
    threads : int, optional
        Number of simultaneous threads used while fetching data.
        Each cursor applies its own default when not set.
    tries : int, default 3
        Maximum number of attempts to download from remote databases.
    progress : bool, default True
        Whether to print progress messages.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
        One record per sequence found.

    See Also
    --------
    rotifer.db.local.easel.FastaCursor : the local backend
    rotifer.db.ncbi.entrez.FastaCursor : the remote backend

    Notes
    -----
    This function is optimized for speed and does not provide access
    to sequence annotation, because all data is fetched as FASTA
    formatted streams. Use :class:`rotifer.db.ncbi.SequenceCursor`
    when annotations are needed.

    Examples
    --------
    Fetch two proteins, trying the local database first and falling
    back to NCBI Entrez:

    >>> import rotifer.db as rdb
    >>> seqs = rdb.proteins(["WP_063732599.1", "YP_009724395.1"])  # doctest: +SKIP
    >>> [s.id for s in seqs]  # doctest: +SKIP
    ['WP_063732599.1', 'YP_009724395.1']
    """
    from rotifer.db.local import easel
    from rotifer.db.ncbi import entrez

    result = []
    targets = set(query)
    for method in methods:
        if method == 'esl_sfetch':
            cursor = easel.FastaCursor(database_path=local_database_path, batch_size=batch_size, threads=threads, progress=progress)
        elif method == 'entrez':
            cursor = entrez.FastaCursor(database=entrez_database, batch_size=batch_size, threads=threads, tries=tries, progress=progress)
        seqs = cursor.fetchall(targets)
        if isinstance(seqs,list):
            result.extend(seqs)
        else:
            result.append(seqs)
        targets = cursor.missing
    if targets:
        logger.warn(f'A total of {len(targets)} sequences could not be found: {targets}')

    return result
