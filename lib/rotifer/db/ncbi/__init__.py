# Copyright 2020 by Robson F. de Souza.  All rights reserved.
# This file is part of the Rotifer distribution and governed by 
# the "BSD 3-Clause License".
#
# Please see the LICENSE file that should have been included as part of this
# package.

r"""
=============================
Rotifer's NCBI database tools
=============================

Class attributes
----------------

config : NCBI configuration
  Automatically loaded from ~/.rotifer/etc/db/ncbi.yml
"""

# Import external modules
import os
import sys
import types
import socket
import typing
import numpy as np
import pandas as pd
from copy import deepcopy

# Import rotifer modules
import rotifer
logger = rotifer.logging.getLogger(__name__)

# Load NCBI configuration
import rotifer.db.core
import rotifer.db.methods
import rotifer.db.delegator
from rotifer.core.functions import loadConfig
config = loadConfig(__name__, defaults = {
        'local_database_path': [ os.path.join(rotifer.config['data'],"fadb","nr","nr") ],
        "entrez_database": "protein",
        "mirror": os.path.join(os.environ["ROTIFER_DATA"] if 'ROTIFER_DATA' in os.environ else "/databases","genomes"),
        'email': os.environ['USER'] + '@' + socket.gethostname() if 'USER' in os.environ  else 'Unk_user' + '@' + socket.gethostname(),
        'ftpserver': 'ftp.ncbi.nlm.nih.gov',
        'api_key': os.environ['NCBI_API_KEY'] if 'NCBI_API_KEY' in os.environ else None,
        'readers': {
            'entrez': 'rotifer.db.ncbi.entrez',
            'easel': 'rotifer.db.local.easel',
            'ete3': 'rotifer.db.local.ete3',
            'ftp': 'rotifer.db.ncbi.ftp',
            'mirror': 'rotifer.db.ncbi.mirror',
            'sqlite3': 'rotifer.db.sql.sqlite3',
        },
        'writers': {
            'sqlite3': 'rotifer.db.sql.sqlite3',
        }
    })
NcbiConfig = config # for compatibility but deprecated: to be removed!

# Classes

class SequenceCursor(rotifer.db.methods.SequenceCursor, rotifer.db.delegator.SequentialDelegatorCursor):
    """
    Fetch sequences from NCBI.

    This class loads and parses sequences in Genbank format, i.e.
    the most richly annotated format from NCBI.

    Usage
    -----
    Fetch a protein sequence
    >>> from rotifer.db import ncbi
    >>> sc = ncbi.SequenceCursor(database="protein")
    >>> seqrec = sc.fetchall("YP_009724395.1")

    Fetch several nucleotide entries
    >>> import sys
    >>> from Bio import SeqIO
    >>> from rotifer.db import ncbi
    >>> sc = ncbi.SequenceCursor(database="nucleotide")
    >>> query = ['CP084314.1', 'NC_019757.1', 'AAHROG010000026.1']
    >>> for seqrec in sc.fetchone(query):
    >>>     print(SeqIO.write(seqrec, sys.stdout, "genbank")

    Parameters
    ----------
    readers: list of strings, default ['entrez']
      List of backend reader modules
    writers: list of strings, default []
      List of backend writer modules
    database: string, default 'protein'
        Valid NCBI sequence database
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    sleep_between_tries: int, default 1
      Number of seconds to wait between download attempts
    batch_size: int, default 1
      Number of accessions per batch
    threads: integer, default 3
      Number of simultaneous threads to run

    """
    def __init__(
            self,
            readers=['entrez'],
            writers=[],
            database=config["entrez_database"],
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=None,
            threads=None,
            *args, **kwargs):
        self._shared_attributes = ['progress','tries','sleep_between_tries','batch_size','threads','database']
        self.sleep_between_tries = sleep_between_tries
        self.database = database
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

class FastaCursor(rotifer.db.methods.SequenceCursor, rotifer.db.delegator.SequentialDelegatorCursor):
    """
    Fetch sequences from NCBI.

    This class downloads FASTA files, i.e. doesn't include
    sequence annotations but is made available for fast access
    to sequence data.

    Usage
    -----
    >>> from rotifer.db import ncbi
    >>> sc = ncbi.FastaCursor(database="protein")
    >>> seqrec = sc.fetchall("YP_009724395.1")

    Parameters
    ----------
    readers: list of strings, default ['entrez']
      List of backend reader modules
    writers: list of strings, default []
      List of backend writer modules
    local_database_path: list of strings
        Path to local FASTA files indexed by esl-sfetch
    entrez_database: string, default 'protein'
        Valid NCBI sequence database
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    sleep_between_tries: int, default 1
      Number of seconds to wait between download attempts
    batch_size: int, default 1
      Number of accessions per batch
    threads: integer, default 3
      Number of simultaneous threads to run

    """
    def __init__(
            self,
            readers=['easel','entrez'],
            writers=[],
            local_database_path=config["local_database_path"],
            entrez_database=config["entrez_database"],
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=None,
            threads=None,
            *args, **kwargs):
        self._shared_attributes = ['progress','tries','sleep_between_tries','batch_size','threads','database','database_path']
        self.sleep_between_tries = sleep_between_tries
        self.database_path = local_database_path
        self.database = entrez_database
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

class IPGCursor(rotifer.db.methods.IPGCursor, rotifer.db.delegator.SequentialDelegatorCursor):
    """
    Fetch identical proteins (IPG) reports.

    Usage
    -----
    >>> from rotifer.db import ncbi
    >>> ic = ncbi.IPGCursor(database="protein")
    >>> df = ic.fetchall("YP_009724395.1")

    Parameters
    ----------
    readers: list of strings, default ['entrez']
      List of backend reader modules
    writers: list of strings, default []
      List of backend writer modules
    local_database_path: list of strings
        Path to local SQLite3 database
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    sleep_between_tries: int, default 1
      Number of seconds to wait between download attempts
    batch_size: int, default 1
      Number of accessions per batch
    threads: integer, default 3
      Number of simultaneous threads to run

    """
    def __init__(
            self,
            readers=['entrez'],
            writers=[],
            local_database_path=None,
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=None,
            threads=None,
            *args, **kwargs):
        self._shared_attributes = ['progress','tries','sleep_between_tries','batch_size','threads','path']
        self.sleep_between_tries = sleep_between_tries
        self.path = local_database_path
        if self.path != None and os.path.exists(self.path):
            readers.append("sqlite3")
            kwargs['path'] = self.path
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

    def fetchall(self, accessions):
        df = super().fetchall(accessions)
        if isinstance(df, list):
            df = pd.concat(df, ignore_index=True)
        return df

class GenomeCursor(rotifer.db.methods.GenomeCursor, rotifer.db.delegator.SequentialDelegatorCursor):
    """
    Fetch annotated genome sequences.

    Usage
    -----
    Load a sample of genomes

    >>> q = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> from rotifer.db as ncbi
    >>> gfc = ncbi.GenomeCursor(progress=True)
    >>> g = gfc.fetchall(q)

    Parameters
    ----------
    readers: list of strings, default ['entrez']
      List of backend reader modules
    writers: list of strings, default []
      List of backend writer modules
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    threads: integer, default 15
      Number of processes to run parallel downloads
    batch_size: int, default 1
      Number of accessions per batch
    mirror: string
      Path to a local mirror of the NCBI's FTP genome repository
    cache: path-like string
      Where to place temporary files

    """
    def __init__(
            self,
            readers=['mirror','ftp'],
            writers=[],
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=None,
            threads=None,
            timeout=10,
            mirror = config["mirror"],
            cache=rotifer.config['cache'],
            *args, **kwargs):
        self._shared_attributes = ['progress','tries','sleep_between_tries','batch_size','threads','cache','path']
        self.sleep_between_tries = sleep_between_tries
        self.timeout = timeout
        if 'path' in kwargs and mirror == None:
            self.path = kwargs['path']
        else:
            self.path = mirror
        self.cache = cache
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

class GenomeFeaturesCursor(rotifer.db.methods.GenomeFeaturesCursor, rotifer.db.delegator.SequentialDelegatorCursor):
    """
    Fetch genome annotation as dataframes.

    Usage
    -----
    Load a random sample of genomes

    >>> g = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> from rotifer.db as ncbi
    >>> gfc = ncbi.GenomeFeaturesCursor(progress=True)
    >>> df = gfc.fetchall(g)

    Parameters
    ----------
    readers: list of strings, default ['entrez']
      List of backend reader modules
    writers: list of strings, default []
      List of backend writer modules
    exclude_type: list of strings
      List of names for the features that must be ignored
    autopid: boolean
      Automatically set protein identifiers
    codontable: string por int, default 'Bacterial'
      Default codon table, if not set within the data
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    threads: integer, default 15
      Number of processes to run parallel downloads
    batch_size: int, default 1
      Number of accessions per batch
    path: string
      Path to a local mirror of the NCBI's FTP genome repository
    cache: path-like string
      Where to place temporary files

    """
    def __init__(
            self,
            readers=['mirror','ftp'],
            writers=[],
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=None,
            threads=None,
            timeout=10,
            path = config["mirror"],
            cache=rotifer.config['cache'],
            *args, **kwargs):
        self._shared_attributes = [
            'progress','tries','sleep_between_tries','batch_size','threads','cache','path',
            'exclude_type','autopid','codontable',
        ]
        self.sleep_between_tries = sleep_between_tries
        self.timeout = timeout
        self.path = path
        self.cache = cache
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

class GeneNeighborhoodCursor(rotifer.db.methods.GeneNeighborhoodCursor, rotifer.db.core.BaseCursor):
    """
    Fetch gene neighborhoods as dataframes
    ======================================

    This class implements the logic to search for genomic
    patches centered around a target coding gene.

    The target genes must be identified based on the accession
    numbers of their proteins.

    For the multi-query methods (fetchone and fetchall)
    results always return in a random order.

    Usage
    -----
    Using the dictionary-like interface, fetch the gene
    neighborhood around the gene encoding a target protein:

    >>> import rotifer.db a ncbi
    >>> gnc = ncbi.GeneNeighborhoodCursor()
    >>> df = gnc["EEE9598493.1"]

    Fetch all gene neighborhoods for a sample of proteins:

    >>> q = ['WP_012291365.1','WP_013208129.1','WP_122330970.1']
    >>> df = gnc.fetchall(q)

    Process gene neighborhoods while downloading:

    >>> for n in gnc.fetchone(q):
    >>>     do_something(n)

    Parameters
    ----------
    *Important note:*

    All parameters for initialization of this class are acessible
    as mutable attributes and can be modified to tune the cursor's
    behaviour.

    column : string
      Name of the column to scan for matches to the accessions
      See rotifer.genome.data.NeighborhoodDF
    before : int
      Keep at most this number of features, of the same type as the
      target, before each target
    after  : nt
      Keep at most this number of features, of the same type as the
      target, after each target
    min_block_distance : int
      Minimum distance between two consecutive blocks
    strand : string
      How to evaluate rows concerning the value of the strand column
      Possible values for this option are:

      - None : ignore strand
      - same : same strand as the targets
      -    + : positive strand features and targets only
      -    - : negative strand features and targets only

    fttype : string
      How to process feature types of neighbors
      Supported values:

      - same : consider only features of the same type as the target
      - any  : ignore feature type and count all features when
               setting neighborhood boundaries

    eukaryotes : boolean, default False
      If set to True, neighborhood data for eukaryotic genomes
    save : string, default None
      If set, the full genome annotation will be saved to a
      local SQLite3 database for later use.

      Note: currently, save won't write data but it can be
            used to access previously loaded databases

    replace : boolean, default True
      When save is set, whether to replace that file
    mirror : path
      Path to a mirror of the NCBI's FTP genomes directory.
    exclude_type: list of strings
      List of names for the features that must be ignored
    autopid: boolean
      Automatically set protein identifiers
    codontable: string por int, default 'Bacterial'
      Default codon table, if not set within the data
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    threads: integer, default 15
      Number of processes to run parallel downloads
    batch_size: int, default 1
      Number of accessions per batch
    cache: path-like string
      Where to place temporary files

    Internal state attributes
    -------------------------
    Objects of this class modify the folllowing read/write
    attributes when fetch methods are called:
    * missing
      A Pandas DataFrame describing errors and messages for
      failed attempts to download gene neighborhoods.

    """
    def __init__(
            self,
            readers = ['ftp','entrez'],
            writers = [],
            column = 'pid',
            before = 7,
            after = 7,
            min_block_distance = 0,
            strand = None,
            fttype = 'same',
            eukaryotes=False,
            save=None,
            replace=False,
            mirror=None,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=True,
            tries=3,
            batch_size=None,
            threads=None,
            cache=rotifer.config['cache'],
            *args, **kwargs
        ):
        super().__init__(
            progress = progress,
            *args, **kwargs
        )
        self.readers = readers.copy()
        self.writers = writers.copy()
        self.batch_size = batch_size
        self.threads = threads
        self.save = save

        # Setup special attributes
        self._shared_attributes = [
            'column','before','after','min_block_distance','strand','fttype','eukaryotes',
            'exclude_type','autopid','codontable',
            'progress','tries','batch_size','threads','cache',
            'giveup',
        ]

        # Loading cursors
        from rotifer.db.ncbi import ftp
        from rotifer.db.ncbi import entrez
        self.cursors = {
            'ftp': ftp.GeneNeighborhoodCursor(),
            'entrez': entrez.GeneNeighborhoodCursor(),
        }
        if mirror:
            from rotifer.db.ncbi import mirror as rdnm
            if isinstance(mirror, list):
                count=len(mirror) +1
                for mirror_path in mirror[::-1]:
                    count -= 1
                    cursor = rdnm.GeneNeighborhoodCursor(path=mirror_path, tries=1)
                    self.readers.insert(0,f'mirror_{count}')
                    self.cursors[f'mirror_{count}'] = cursor
            else:
                cursor = rdnm.GeneNeighborhoodCursor(path=mirror)
                if 'mirror' not in self.readers:
                    self.readers.insert(0,'mirror')
                self.cursors['mirror'] = cursor
        if save:
            from rotifer.db.sql import sqlite3 as rdss
            cursor = rdss.GeneNeighborhoodCursor(save, replace=replace)
            if 'sqlite3' not in self.readers:
                self.readers.insert(0,'sqlite3')
            #if 'sqlite3' not in self.writers:
            #   self.writers.insert(0,'sqlite3')
            self.cursors['sqlite3'] = cursor

        # Setup simple attributes
        self.column = column
        self.before = before
        self.after = after
        self.min_block_distance = min_block_distance
        self.strand = strand
        self.fttype = fttype
        self.eukaryotes = eukaryotes
        self.exclude_type = exclude_type.copy()
        self.autopid = autopid
        self.codontable = codontable
        self.progress = progress
        self.tries = tries
        self.cache = cache
        self.giveup.update(["HTTP Error 400"])
        # A protein with no IPG cannot be placed on a genome, and
        # every backend of this cursor needs one, so the verdict
        # binds them all
        self.final_errors.update(["no IPG","No IPG"])
        if not eukaryotes:
            self.giveup.update(["Eukaryot","eukaryot"])

    def __setattr__(self, name, value):
        super().__setattr__(name, value)
        if hasattr(self,'cursors') and hasattr(self,'_shared_attributes') and name in self._shared_attributes:
            for cursor in self.cursors.values():
                if hasattr(cursor,name) and not isinstance(value,types.NoneType):
                    cursor.__setattr__(name,value)

    def __getitem__(self, protein, ipgs=None):
        """
        Dictionary-like access to gene neighbors.

        Usage
        -----
        >>> import rotifer.db as ncbi
        >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)
        >>> n = gnc["WP_063732599.1"]

        Parameters
        ----------
        proteins: list of strings
          Database identifiers.
        ipgs : Pandas dataframe
          This parameter may be used to avoid downloading IPGs
          from NCBI several times. Example:

          >>> from rotifer.db as ncbi
          >>> ic = ncbi.IPGCursor(batch_size=1)
          >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)
          >>> i = ic.fetchall(['WP_063732599.1'])
          >>> n = gnc.__getitem__(['WP_063732599.1'], ipgs=i)

        Returns
        -------
        Generator of rotifer.genome.data.NeighborhoodDF
         """
        from rotifer.genome.utils import seqrecords_to_dataframe
        result = seqrecords_to_dataframe([])
        targets = self.parse_ids(protein)
        tried = []
        for reader in self.readers:
            if reader not in self.cursors:
                continue
            tried.append(reader)
            reader = self.cursors[reader]
            result = reader.__getitem__(targets, ipgs=ipgs)
            if not isinstance(result,types.NoneType) and len(result) > 0:
                for otherReader in tried:
                    self.cursors[otherReader].remove_missing(targets)
                #if self.save:
                #   for writer in self.writers:
                #       if writer not in self.cursors:
                #           continue
                #       self.cursors[writer].insert(result)
                break
        return result

    def fetchone(self, accessions, ipgs=None):
        """
        Fetch each gene neighborhood iteratively.

        Parameters
        ----------
        accessions: list of strings
          Database identifiers.
        ipgs : Pandas dataframe
          This parameter may be used to avoid downloading IPGs
          from NCBI several times. Example:

          >>> from rotifer.db as ncbi
          >>> ic = ncbi.IPGCursor(batch_size=1)
          >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)
          >>> i = ic.fetchall(['WP_063732599.1'])
          >>> for x in gnc.fetchone(['WP_063732599.1'], ipgs=i):
          >>>    do_something(x)

        Returns
        -------
        Generator of rotifer.genome.data.NeighborhoodDF
        """
        from rotifer.genome.utils import seqrecords_to_dataframe

        # Copy identifiers and remove redundancy
        targets = self.parse_ids(accessions)

        # Make sure we have IPGs
        if isinstance(ipgs,types.NoneType):
            from rotifer.db.ncbi import entrez
            if self.progress:
                logger.warn(f'Downloading IPGs for {len(targets)} proteins....')
            ic = entrez.IPGCursor(progress=self.progress, tries=self.tries)
            ipgs = ic.fetchall(targets)
            self.update_missing(data=ic.remove_missing())
            targets = targets - self.missing_ids()

        # Select IPGs corresponding to our queries
        ipgs = ipgs[ipgs.id.isin(ipgs[ipgs.pid.isin(targets) | ipgs.representative.isin(targets)].id)]
        missing = targets - set(ipgs.pid).union(ipgs.representative)
        if missing:
            # No IPG means no way to place the protein on a genome,
            # and every backend of this cursor needs one
            self.update_missing(missing,"Not found in IPGs",False,final=True)
            targets = targets - missing
        if len(ipgs) == 0:
            return [seqrecords_to_dataframe([])]

        # Call cursors
        tried = []
        for reader in self.readers:
            if len(targets) == 0:
                break
            if reader not in self.cursors:
                continue
            tried.append(reader)
            reader = self.cursors[reader]
            for result in reader.fetchone(targets, ipgs=ipgs):
                if self.save:
                    for writer in self.writers:
                        if writer not in self.cursors:
                            continue
                        self.cursors[writer].insert(result)
                found = targets.intersection(self.getids(result, ipgs=ipgs))
                for readerName in tried:
                    self.cursors[readerName].remove_missing(found)
                self.remove_missing(found)
                self.update_missing(data=reader._missing)
                targets = targets - found
                yield result

    def fetchall(self, proteins, ipgs=None):
        """
        Fetch all gene neighborhoods in a single dataframe.

        Parameters
        ----------
        proteins: list of strings
          Database identifiers.
        ipgs : Pandas dataframe
          This parameter may be used to avoid downloading IPGs
          from NCBI several times. Example:

          >>> from rotifer.db as ncbi
          >>> ic = ncbi.IPGCursor(batch_size=1)
          >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)
          >>> i = ic.fetchall(['WP_063732599.1'])
          >>> n = gnc.fetchall(['WP_063732599.1'], ipgs=i)

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
        """
        from rotifer.genome.utils import seqrecords_to_dataframe
        stack = []
        for df in self.fetchone(proteins, ipgs=ipgs):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])

class TaxonomyCursor(rotifer.db.delegator.SequentialDelegatorCursor):
    def __init__(self, readers=['ete3','entrez'], writers=[], progress=True, tries=3, sleep_between_tries=1, batch_size=None, threads=None, *args, **kwargs):
        self._shared_attributes = ['progress','tries','sleep_between_tries','batch_size','threads']
        self.sleep_between_tries = sleep_between_tries
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        self.taxcols = ['taxid','organism','superkingdom','lineage','classification','alternative_taxids']

    def getids(self, obj, *args, **kwargs):
        if not (isinstance(obj,list) or isinstance(obj,tuple)):
            obj = [ obj ]
        ids = set()
        for item in obj:
            ids.update(set(item.taxid.astype(str)))
            if 'alternative_taxids' in item:
                aids = item.alternative_taxids.dropna().astype(str)
                aids = aids.str.split(",").explode().dropna()
                ids.update(aids)
        return ids

    def __getitem__(self, accessions, *args, **kwargs):
        """
        Dictionary-like access to data.

        Usage
        -----
        >>> from rotifer.db import ncbi
        >>> tc = ncbi.TaxonomyCursor(progress=True)
        >>> t = tc[2599]

        Parameters
        ----------
        accessions: list of strings
          Database identifiers.

        Returns
        -------
        Pandas dataframe
        """
        result = super().__getitem__(accessions, *args, **kwargs)
        if len(result) == 0:
            return pd.DataFrame(columns=self.taxcols)
        elif isinstance(result,list):
            return pd.concat(result, ignore_index=True)
        else:
            return result

    def fetchall(self, accessions, *args, **kwargs):
        """
        Fetch data for all accessions.

        Parameters
        ----------
        accessions: list of database identifiers
          Database identifiers.

        Returns
        -------
        Pandas dataframe
        """
        df = super().fetchall(accessions, *args, **kwargs)
        if len(df) == 0:
            return pd.DataFrame(columns=self.taxcols)
        else:
            return pd.concat(df, ignore_index=True)

# FUNCTIONS

# Load NCBI assembly reports
def assemblies(baseurl=f'ftp://{config["ftpserver"]}/genomes/ASSEMBLY_REPORTS', targets=['refseq', 'genbank', 'refseq_historical', 'genbank_historical'], taxonomy=True, progress=True):
    '''
    Load a table documenting all NCBI genome assemblies.

    By default, these concatenated tables are downloaded
    from the genomes/ASSEMBLY_REPORTS directory at NCBI's
    FTP site.

    Usage
    -----
    Download from NCBI's FTP site

    >>> from rotifer.db import ncbi
    >>> a = ncbi.assembly_reports()

    Load local files at /db/ncbi

    >>> b = ncbi.assembly_reports(baseurl="/db/ncbi")
      
    If working at NIH servers use

    >>> a = ncbi.assembly_reports(baseurl="/am/ftp-genomes/ASSEMBLY_REPORTS")

    Parameters
    ----------
    baseurl: string
      URL or directory with assembly_summary_*.txt files
    targets: list
      List of genome database sections to load.
      Options are: refseq, genbank, genbank_historical, refseq_historical
    taxonomy: boolean, default True
      If set to true, taxonomy data is added to the table
    progress: boolean
      Display progress messages.

    Returns
    -------
    Pandas DataFrame

    Notes
    -----
    Some columns are added to the original table:
    source : NCBI's source database
    loaded_from : data source (same as baseurl)
    '''

    # Method dependencies
    import pandas as pd
    from glob import glob
    origLevel = logger.getEffectiveLevel()
    if progress:
        rotifer.logger.setLevel(rotifer.logging.INFO)
    logger.info(f'main: loading assembly reports...')

    # Load assembly reports
    df = list()
    for x in targets:
        if os.path.exists(baseurl): # Local file
            url = os.path.join(baseurl, f'assembly_summary_{x}.txt')
            if not os.path.exists(url):
                logger.warning(f'{__name__}: {url} not found. Ignoring...')
                continue
        else: # FTP
            url = f'{baseurl}/assembly_summary_{x}.txt'
        _ = pd.read_csv(url, sep ="\t", skiprows=[0], low_memory=False)
        _.rename({'# assembly_accession':'assembly'}, axis=1, inplace=True)
        _['source'] = x
        _['loaded_from'] = url
        df.append(_)
        logger.info(f'{url}, {len(_)} rows, {len(df)} loaded')
    df = pd.concat(df, ignore_index=True)
    df.taxid = df.taxid.astype(str)
    logger.info(f'loaded {len(df)} assembly summaries.')

    # Make sure the ftp_path columns refers to the ftp site as we expect
    if 'ftp_path' in df.columns:
        df.ftp_path = df.ftp_path.str.replace('https','ftp')

    # Add taxonomy
    if taxonomy:
        cursor = TaxonomyCursor(progress=progress)
        taxonomy = cursor.fetchall(df.taxid.unique().tolist())
        taxonomy['_same'] = (taxonomy.taxid == taxonomy.alternative_taxids).astype(int)
        taxonomy.sort_values(['taxid','_same'], ascending=True, inplace=True)
        taxonomy.drop('_same', axis=1, inplace=True)
        taxonomy.drop_duplicates('taxid', keep='first', inplace=True)
        df = df.merge(taxonomy, left_on='taxid', right_on='taxid', how='left')
        logger.info(f'{len(df)} df left-merged with taxonomy dataframe.')

    # Reset ncbi object, update missing list and return
    logger.info(f'main: {len(df)} assembly reports loaded!')
    if progress:
        rotifer.logger.setLevel(origLevel)
    return df

# END
if __name__ == '__main__':
    pass
