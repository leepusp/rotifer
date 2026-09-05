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
    Fetch proteins sequences from local or remote databases.
    
    Note
    ----
    This method is optimized for speed and does not provide
    access to sequence annotation data, because all data is
    fetched as FASTA formatted data streams.

    Parameters
    ----------
    query: list of string
      Sequence identifiers
    methods: list of strings
      List of keywords listing methods of data retrieval.
      The following methods are available:
      - esl_sfetch
        rotifer.db.local.EaselFastaCursor
      - entrez
        rotifer.db.ncbi.entrez.FastaCursor
    local_database_path: strin
      Path to a (possibly indexed) sequence file.
      In most cases the local file will be a FASTA file.
    entrez_database: string
      Name the Entrez sequence database
    batch_size: int, default 200
      Number of sequences to retrieve per iteration
    threads: integer
      Number of simultaneous threads to run while fetching data
    tries: integer
      Maximum number of attempts to download from remote databases
    progress: boolean, default True
      Whether to print progress messages

    Returns
    -------
      One or a list of Bio.SeqRecord objects
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
        # missing is a dataframe: iterating it yields its column names,
        # so the next method would be asked for those instead of the
        # accessions. Entries some cursor declared final are dropped
        # here too, since no later method will find them either.
        targets = cursor.missing_ids(final=False)
    if targets:
        logger.warning(f'A total of {len(targets)} sequences could not be found: {sorted(targets)}')

    return result
