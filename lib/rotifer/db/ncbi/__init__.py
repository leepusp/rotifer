# Copyright 2020 by Robson F. de Souza.  All rights reserved.
# This file is part of the Rotifer distribution and governed by 
# the "BSD 3-Clause License".
#
# Please see the LICENSE file that should have been included as part of this
# package.

r"""
Access NCBI sequence, genome and taxonomy data.

This package is the main entry point for retrieving data from the
National Center for Biotechnology Information (NCBI). The cursors
defined here are delegators: each one combines several backends
(Entrez E-utilities, the NCBI FTP site, local genome mirrors, local
indexed FASTA files, the ETE toolkit taxonomy database and SQLite3
stores) and tries them in order until every requested identifier is
resolved.

Network, authentication and rate limits
---------------------------------------
Entrez based backends send the user email registered in the
configuration and honor the ``NCBI_API_KEY`` environment variable.
Without an API key, NCBI limits clients to 3 requests per second and
the Entrez backends restrict themselves to 3 simultaneous threads;
with a key, up to 10 threads are used. FTP backends download files
to the rotifer cache directory (``rotifer.config['cache']``).

Configuration
-------------
The module level ``config`` dictionary is loaded from
``~/.rotifer/etc/db/ncbi.yml`` when that file exists. Its keys
include the Entrez email and API key, the FTP server address, the
path of a local genome mirror and the mapping of backend names to
reader and writer modules.
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
    Fetch annotated sequences from NCBI.

    Sequences are downloaded and parsed in GenBank format, the most
    richly annotated format NCBI provides.

    Parameters
    ----------
    readers : list of str, default ``['entrez']``
        Backend reader modules, tried in order.
    writers : list of str, default []
        Backend writer modules.
    database : str, default 'protein'
        A valid NCBI sequence database name, such as ``protein`` or
        ``nucleotide``.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of simultaneous threads.

    See Also
    --------
    FastaCursor : faster download without annotations
    rotifer.db.ncbi.entrez.SequenceCursor : the default backend

    Examples
    --------
    Fetch a protein sequence:

    >>> from rotifer.db import ncbi
    >>> sc = ncbi.SequenceCursor(database="protein")  # doctest: +SKIP
    >>> seqrec = sc.fetchall("YP_009724395.1")  # doctest: +SKIP

    Fetch several nucleotide entries:

    >>> import sys
    >>> from Bio import SeqIO
    >>> query = ['CP084314.1', 'NC_019757.1', 'AAHROG010000026.1']
    >>> sc = ncbi.SequenceCursor(database="nucleotide")  # doctest: +SKIP
    >>> for seqrec in sc.fetchone(query):  # doctest: +SKIP
    ...     SeqIO.write(seqrec, sys.stdout, "genbank")
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
    Fetch sequences from local FASTA files or NCBI, without
    annotations.

    Sequences are retrieved as FASTA data, which is much faster than
    GenBank format but carries no annotation. By default a local
    ``esl-sfetch`` indexed database is tried before NCBI Entrez.

    Parameters
    ----------
    readers : list of str, default ``['easel', 'entrez']``
        Backend reader modules, tried in order.
    writers : list of str, default []
        Backend writer modules.
    local_database_path : list of str, optional
        Paths to local FASTA files indexed by ``esl-sfetch``.
        Defaults to the ``local_database_path`` configuration entry.
    entrez_database : str, default 'protein'
        A valid NCBI sequence database name.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of simultaneous threads.

    See Also
    --------
    SequenceCursor : slower download with full annotations
    rotifer.db.local.easel.FastaCursor : the local backend
    rotifer.db.ncbi.entrez.FastaCursor : the remote backend

    Examples
    --------
    >>> from rotifer.db import ncbi
    >>> fc = ncbi.FastaCursor()  # doctest: +SKIP
    >>> seqrec = fc.fetchall("YP_009724395.1")  # doctest: +SKIP
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
    Fetch identical protein group (IPG) reports.

    An IPG report lists every protein sequence identical to the
    query known to NCBI, together with the nucleotide sequences and
    genome assemblies encoding them. When a local SQLite3 database
    path is given and exists, it is queried before NCBI Entrez.

    Parameters
    ----------
    readers : list of str, default ``['entrez']``
        Backend reader modules, tried in order.
    writers : list of str, default []
        Backend writer modules.
    local_database_path : str, optional
        Path to a local SQLite3 database. Appended to the readers
        when the file exists.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of simultaneous threads.

    See Also
    --------
    rotifer.db.ncbi.entrez.IPGCursor : the remote backend
    rotifer.db.sql.sqlite3.IPGCursor : the local backend

    Examples
    --------
    >>> from rotifer.db import ncbi
    >>> ic = ncbi.IPGCursor()  # doctest: +SKIP
    >>> df = ic.fetchall("YP_009724395.1")  # doctest: +SKIP
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
        """
        Fetch the IPG reports of all accessions as one dataframe.

        Parameters
        ----------
        accessions : list of str
            NCBI protein accessions.

        Returns
        -------
        pandas.DataFrame
            The concatenated IPG reports.
        """
        df = super().fetchall(accessions)
        if isinstance(df, list):
            df = pd.concat(df, ignore_index=True)
        return df

class GenomeCursor(rotifer.db.methods.GenomeCursor, rotifer.db.delegator.SequentialDelegatorCursor):
    """
    Fetch annotated genome sequences.

    Genomes are located by their assembly accession and parsed from
    GenBank flat files into Bio.SeqRecord objects. A local mirror of
    the NCBI genomes repository, when available, is tried before the
    NCBI FTP site.

    Parameters
    ----------
    readers : list of str, default ``['mirror', 'ftp']``
        Backend reader modules, tried in order.
    writers : list of str, default []
        Backend writer modules.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of processes used for parallel downloads.
    timeout : int, default 10
        Maximum time, in seconds, to wait for a server connection.
    mirror : str, optional
        Path to a local mirror of the NCBI FTP genomes repository.
        Defaults to the ``mirror`` configuration entry.
    cache : str, optional
        Directory for temporary files. Defaults to the rotifer cache
        directory.

    See Also
    --------
    GenomeFeaturesCursor : genome annotation as dataframes
    rotifer.db.ncbi.ftp.GenomeCursor : the FTP backend
    rotifer.db.ncbi.mirror.GenomeCursor : the local mirror backend

    Examples
    --------
    Load a sample of genomes:

    >>> from rotifer.db import ncbi
    >>> q = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> gc = ncbi.GenomeCursor(progress=True)  # doctest: +SKIP
    >>> g = gc.fetchall(q)  # doctest: +SKIP
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
    Fetch genome annotation as feature tables.

    Genomes are located by their assembly accession and converted to
    dataframes with one row per annotated feature. A local mirror of
    the NCBI genomes repository, when available, is tried before the
    NCBI FTP site.

    Parameters
    ----------
    readers : list of str, default ``['mirror', 'ftp']``
        Backend reader modules, tried in order.
    writers : list of str, default []
        Backend writer modules.
    exclude_type : list of str, default ``['source', 'gene', 'mRNA']``
        Feature types to ignore.
    autopid : bool, default False
        Automatically set protein identifiers.
    codontable : str or int, default 'Bacterial'
        Codon table used when the data does not define one.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of processes used for parallel downloads.
    timeout : int, default 10
        Maximum time, in seconds, to wait for a server connection.
    path : str, optional
        Path to a local mirror of the NCBI FTP genomes repository.
        Defaults to the ``mirror`` configuration entry.
    cache : str, optional
        Directory for temporary files. Defaults to the rotifer cache
        directory.

    See Also
    --------
    GenomeCursor : genomes as annotated sequence records
    GeneNeighborhoodCursor : only the regions around target genes

    Examples
    --------
    Load the feature tables of two genomes:

    >>> from rotifer.db import ncbi
    >>> g = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> gfc = ncbi.GenomeFeaturesCursor(progress=True)  # doctest: +SKIP
    >>> df = gfc.fetchall(g)  # doctest: +SKIP
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
    Fetch gene neighborhoods as dataframes.

    This class searches for genomic patches centered around target
    coding genes, identified by the accession numbers of their
    protein products. Backends are tried in order: a local SQLite3
    store (when `save` is set), local genome mirrors (when `mirror`
    is set), the NCBI FTP site and NCBI Entrez.

    For the multi query methods (:meth:`fetchone` and
    :meth:`fetchall`) results are returned in random order.

    Parameters
    ----------
    readers : list of str, default ``['ftp', 'entrez']``
        Backend reader modules, tried in order. The `mirror` and
        `save` parameters prepend their backends to this list.
    writers : list of str, default []
        Backend writer modules.
    column : str, default 'pid'
        Name of the column to scan for matches to the accessions.
        See :class:`rotifer.genome.data.NeighborhoodDF`.
    before : int, default 7
        Keep at most this number of features, of the same type as
        the target, before each target.
    after : int, default 7
        Keep at most this number of features, of the same type as
        the target, after each target.
    min_block_distance : int, default 0
        Minimum distance between two consecutive blocks.
    strand : str, optional
        How to evaluate rows concerning the value of the strand
        column. Supported values:

        * ``None`` : ignore strand
        * ``same`` : same strand as the targets
        * ``+`` : positive strand features and targets only
        * ``-`` : negative strand features and targets only
    fttype : {'same', 'any'}, default 'same'
        How to process feature types of neighbors. With ``same``,
        only features of the same type as the target are considered.
        With ``any``, all features count when setting neighborhood
        boundaries.
    eukaryotes : bool, default False
        Whether to process eukaryotic genomes.
    save : str, optional
        Path to a local SQLite3 database. When set, the database is
        queried before any remote source. Note: currently the save
        backend does not write new data; it only reads previously
        loaded databases.
    replace : bool, default False
        When `save` is set, whether to replace that file.
    mirror : str or list of str, optional
        Path(s) to local mirrors of the NCBI FTP genomes directory,
        queried before any remote source.
    exclude_type : list of str, default ``['source', 'gene', 'mRNA']``
        Feature types to ignore.
    autopid : bool, default False
        Automatically set protein identifiers.
    codontable : str or int, default 'Bacterial'
        Codon table used when the data does not define one.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of processes used for parallel downloads.
    cache : str, optional
        Directory for temporary files. Defaults to the rotifer cache
        directory.

    Notes
    -----
    All initialization parameters are accessible as mutable
    attributes and may be modified to tune the cursor's behavior.
    Changes are propagated to the backend cursors. Every fetch
    method also updates the ``missing`` dataframe, which describes
    the errors of failed download attempts.

    Examples
    --------
    Using the dictionary-like interface, fetch the gene neighborhood
    around the gene encoding a target protein:

    >>> from rotifer.db import ncbi
    >>> gnc = ncbi.GeneNeighborhoodCursor()  # doctest: +SKIP
    >>> df = gnc["EEE9598493.1"]  # doctest: +SKIP

    Fetch all gene neighborhoods for a sample of proteins:

    >>> q = ['WP_012291365.1', 'WP_013208129.1', 'WP_122330970.1']
    >>> df = gnc.fetchall(q)  # doctest: +SKIP

    Process gene neighborhoods while downloading:

    >>> for n in gnc.fetchone(q):  # doctest: +SKIP
    ...     print(n.block_id.nunique())
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
        self.giveup.update(["no IPG","No IPG"])
        if not eukaryotes:
            self.giveup.update(["Eukaryot","eukaryot"])

    def __setattr__(self, name, value):
        """
        Set an attribute, propagating shared ones to the backends.

        Attributes named in ``_shared_attributes`` are also assigned
        on every backend cursor that already defines them, so that
        retuning the delegator keeps its backends in sync. ``None``
        values are never propagated.

        Parameters
        ----------
        name : str
            Attribute name.
        value : object
            Value to assign. Forwarded to the backends only when it
            is not None.
        """
        super().__setattr__(name, value)
        if hasattr(self,'cursors') and hasattr(self,'_shared_attributes') and name in self._shared_attributes:
            for cursor in self.cursors.values():
                if hasattr(cursor,name) and not isinstance(value,types.NoneType):
                    cursor.__setattr__(name,value)

    def __getitem__(self, protein, ipgs=None):
        """
        Fetch gene neighborhoods for one or more proteins,
        dictionary style.

        Backends are tried in order until one returns data.

        Parameters
        ----------
        protein : str or iterable of str
            NCBI protein accessions.
        ipgs : pandas.DataFrame, optional
            Precomputed identical protein group reports, used to
            avoid downloading IPGs several times.

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
            The neighborhoods found. Empty when nothing could be
            retrieved.

        Examples
        --------
        >>> from rotifer.db import ncbi
        >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)  # doctest: +SKIP
        >>> n = gnc["WP_063732599.1"]  # doctest: +SKIP

        Reuse previously downloaded IPGs:

        >>> ic = ncbi.IPGCursor(batch_size=1)  # doctest: +SKIP
        >>> i = ic.fetchall(['WP_063732599.1'])  # doctest: +SKIP
        >>> n = gnc.__getitem__(['WP_063732599.1'], ipgs=i)  # doctest: +SKIP
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
        Iterate over gene neighborhoods as they are retrieved.

        Backends are tried in order; each one only receives the
        identifiers that previous backends could not resolve.

        Parameters
        ----------
        accessions : list of str
            NCBI protein accessions.
        ipgs : pandas.DataFrame, optional
            Precomputed identical protein group reports, used to
            avoid downloading IPGs several times.

        Yields
        ------
        rotifer.genome.data.NeighborhoodDF
            One dataframe per retrieved gene neighborhood.

        Examples
        --------
        >>> from rotifer.db import ncbi
        >>> ic = ncbi.IPGCursor(batch_size=1)  # doctest: +SKIP
        >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)  # doctest: +SKIP
        >>> i = ic.fetchall(['WP_063732599.1'])  # doctest: +SKIP
        >>> for x in gnc.fetchone(['WP_063732599.1'], ipgs=i):  # doctest: +SKIP
        ...     print(len(x))
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
            self.update_missing(missing,"Not found in IPGs",False)
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
        Fetch all gene neighborhoods as a single dataframe.

        Parameters
        ----------
        proteins : list of str
            NCBI protein accessions.
        ipgs : pandas.DataFrame, optional
            Precomputed identical protein group reports, used to
            avoid downloading IPGs several times.

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
            The concatenated neighborhoods. Empty when nothing could
            be retrieved.

        Examples
        --------
        >>> from rotifer.db import ncbi
        >>> ic = ncbi.IPGCursor(batch_size=1)  # doctest: +SKIP
        >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)  # doctest: +SKIP
        >>> i = ic.fetchall(['WP_063732599.1'])  # doctest: +SKIP
        >>> n = gnc.fetchall(['WP_063732599.1'], ipgs=i)  # doctest: +SKIP
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
    """
    Fetch NCBI Taxonomy data as dataframes.

    Taxonomy identifiers are first searched in the local ETE toolkit
    copy of the NCBI Taxonomy database and only the identifiers not
    found there are sent to NCBI Entrez.

    Parameters
    ----------
    readers : list of str, default ``['ete3', 'entrez']``
        Backend reader modules, tried in order.
    writers : list of str, default []
        Backend writer modules.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of simultaneous threads.

    Attributes
    ----------
    taxcols : list of str
        Columns of the returned dataframes: ``taxid``, ``organism``,
        ``superkingdom``, ``lineage``, ``classification`` and
        ``alternative_taxids``.

    See Also
    --------
    rotifer.db.local.ete3.TaxonomyCursor : the local backend
    rotifer.db.ncbi.entrez.TaxonomyCursor : the remote backend

    Examples
    --------
    >>> from rotifer.db import ncbi
    >>> tc = ncbi.TaxonomyCursor(progress=False)  # doctest: +SKIP
    >>> t = tc.fetchall([2599])  # doctest: +SKIP
    """

    def __init__(self, readers=['ete3','entrez'], writers=[], progress=True, tries=3, sleep_between_tries=1, batch_size=None, threads=None, *args, **kwargs):
        self._shared_attributes = ['progress','tries','sleep_between_tries','batch_size','threads']
        self.sleep_between_tries = sleep_between_tries
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        self.taxcols = ['taxid','organism','superkingdom','lineage','classification','alternative_taxids']

    def getids(self, obj, *args, **kwargs):
        """
        Extract taxonomy identifiers from taxonomy dataframes.

        Parameters
        ----------
        obj : pandas.DataFrame or list of pandas.DataFrame
            Taxonomy dataframes produced by the cursor.

        Returns
        -------
        set of str
            All identifiers in the ``taxid`` and
            ``alternative_taxids`` columns.
        """
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
        Fetch taxonomy data for one or more taxids, dictionary
        style.

        Parameters
        ----------
        accessions : int, str or iterable
            NCBI Taxonomy identifiers.

        Returns
        -------
        pandas.DataFrame
            One row per taxon, with the columns listed in
            ``taxcols``.

        Examples
        --------
        >>> from rotifer.db import ncbi
        >>> tc = ncbi.TaxonomyCursor(progress=False)  # doctest: +SKIP
        >>> t = tc[2599]  # doctest: +SKIP
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
        Fetch taxonomy data for all taxids as one dataframe.

        Parameters
        ----------
        accessions : list
            NCBI Taxonomy identifiers.

        Returns
        -------
        pandas.DataFrame
            One row per taxon, with the columns listed in
            ``taxcols``.
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

    By default, the concatenated assembly summary tables are
    downloaded from the ``genomes/ASSEMBLY_REPORTS`` directory of the
    NCBI FTP site. A local directory containing
    ``assembly_summary_*.txt`` files may be used instead.

    Parameters
    ----------
    baseurl : str, optional
        URL or directory with ``assembly_summary_*.txt`` files.
        Defaults to the NCBI FTP site named in the configuration.
    targets : list of str, optional
        Genome database sections to load. Options are ``refseq``,
        ``genbank``, ``refseq_historical`` and
        ``genbank_historical``. All four are loaded by default.
    taxonomy : bool, default True
        Whether to add taxonomy data to the table, using
        :class:`TaxonomyCursor`.
    progress : bool, default True
        Display progress messages.

    Returns
    -------
    pandas.DataFrame
        One row per assembly.

    Notes
    -----
    Two columns are added to the original table: ``source`` (the
    genome database section) and ``loaded_from`` (the URL or path
    each row was read from). The ``ftp_path`` column is rewritten to
    use the ``ftp`` scheme.

    Examples
    --------
    Download from the NCBI FTP site:

    >>> from rotifer.db import ncbi
    >>> a = ncbi.assemblies()  # doctest: +SKIP

    Load local copies stored in ``/db/ncbi``:

    >>> b = ncbi.assemblies(baseurl="/db/ncbi")  # doctest: +SKIP
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
