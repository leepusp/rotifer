"""
Access genome data from the NCBI FTP site.

Cursors in this module locate genome assemblies under the ``genomes``
directory of the NCBI FTP server (``ftp.ncbi.nlm.nih.gov`` by
default), download their GenBank flat files and parse them into
sequence records, feature tables or gene neighborhoods.

Network and caching
-------------------
Connections are anonymous; no account or API key is needed. Every
downloaded file is verified against the ``md5checksums.txt`` published
next to it and stored in the rotifer cache directory
(``rotifer.config['cache']``), from which it is removed once parsed.
The FTP server address can be changed through the ``ftpserver`` entry
of the :mod:`rotifer.db.ncbi` configuration.
"""

import os
import sys
import types
import socket
import typing
import numpy as np
import pandas as pd
from tqdm import tqdm
from ftplib import FTP
from copy import deepcopy

import rotifer
import rotifer.db.core
import rotifer.db.parallel
from rotifer.db.ncbi import NcbiConfig
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

# Configuration
_defaults = {
    'batch_size': 4,
    "maxgetitem": 1,
    "threads": 15,
}
config = loadConfig(__name__, defaults = _defaults)

# Classes

class connection():
    """
    A live connection to the NCBI FTP server.

    Parameters
    ----------
    url : str, optional
        NCBI FTP server address. Defaults to the ``ftpserver``
        configuration entry.
    tries : int, default 3
        Maximum number of attempts to download a file.
    timeout : int, default 50
        Maximum time, in seconds, to wait for a server connection.
    cache : str, optional
        Directory for temporary files. Defaults to the rotifer
        cache directory.

    Examples
    --------
    >>> from rotifer.db.ncbi import ftp as ncbiftp
    >>> conn = ncbiftp.connection(tries=3)  # doctest: +SKIP
    >>> localpath = conn.ftp_get('genomes/README.txt')  # doctest: +SKIP
    """

    def __init__(self, url=NcbiConfig['ftpserver'], tries=3, timeout=50, cache=rotifer.config['cache']):
        self.url = url
        self.tries = tries
        self.timeout = timeout
        self.cache = cache
        self.connect()

    def connect(self):
        """
        Connect or reconnect to the server.

        An existing connection is probed with a ``NOOP`` command and
        replaced when stale. Connection attempts are repeated up to
        ``tries`` times, sleeping one second after each timeout.
        """
        import time
        attempt = 0
        while attempt < self.tries:
            try:
                if hasattr(self,"connection"):
                    self.connection.sendcmd("NOOP")
                else:
                    self.connection = FTP(self.url, timeout=self.timeout)
                    self.connection.login()
                break
            except:
                try:
                    self.connection = FTP(self.url, timeout=self.timeout)
                    self.connection.login()
                    break
                except socket.timeout:
                    time.sleep(1)
                except TimeoutError:
                    time.sleep(1)
            attempt += 1

    # Download files
    def ftp_get(self, target, avoid_collision=False, outdir=None):
        '''
        Download a file from the NCBI FTP site.

        Parameters
        ----------
        target : str
            URL or path of the file, relative to the server root.
        avoid_collision : bool, default False
            Avoid name collisions by adding a random suffix to the
            local file name.
        outdir : str, optional
            Output directory. Defaults to the connection's cache
            directory, which is created when missing.

        Returns
        -------
        str
            Path to the downloaded file.

        Raises
        ------
        IOError
            If the output directory cannot be created or the
            download fails.

        Examples
        --------
        >>> from rotifer.db.ncbi import ftp as ncbiftp
        >>> conn = ncbiftp.connection(tries=3)  # doctest: +SKIP
        >>> localpath = conn.ftp_get('genomes/README.txt')  # doctest: +SKIP
        '''
        from tempfile import NamedTemporaryFile

        # Create output directory, if necessary
        if not outdir:
            outdir = self.cache
        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except:
                raise IOError(f'failed to create download directory {outdir}')

        # Retrieve contents for each folder
        # Prepare local file handle
        if avoid_collision:
            parts = os.path.splitext(os.path.basename(target))
            prefix = parts[0] + '.' 
            suffix = None if parts[1] == '' else parts[1]
            outfh = NamedTemporaryFile(mode='w+b', suffix=suffix, prefix=prefix, dir=outdir, delete=False)
            outfile = outfh.name
        else:
            outfile = os.path.join(outdir, os.path.basename(target))
            outfh   = open(outfile,'wb')

        # To avoid problems with very long names, I change to the target
        # directory and, later, back to /
        p = target.replace('ftp://',"")
        p = p.replace(self.url,"")
        p = p[1:] if p[0] == "/" else p
        p = p.split("/")
        self.connect()
        for i in list(range(len(p)-1)):
            self.connection.cwd(p[i])
        try:
            self.connection.retrbinary("RETR " + p[-1], outfh.write)
        except:
            raise IOError(f'Unable to write stream to {target}')
        self.connection.cwd("/")
        outfh.close()

        # Return pandas object
        logger.debug(f'Download complete {self.url}/{target}')
        return outfile

    # List files in ftp directory
    def ftp_ls(self, targets):
        '''
        List the contents of one or more FTP directories.

        Parameters
        ----------
        targets : str or list of str
            Path of one or more directories.

        Returns
        -------
        pandas.DataFrame
            One row per directory entry, with the server's file
            facts plus ``target`` (the listed directory) and
            ``name`` (the entry name).

        Raises
        ------
        FileNotFoundError
            If a directory listing cannot be retrieved.

        Examples
        --------
        >>> from rotifer.db.ncbi import ftp as ncbiftp
        >>> conn = ncbiftp.connection()  # doctest: +SKIP
        >>> contents = conn.ftp_ls('genomes')  # doctest: +SKIP
        '''
        import pandas as pd

        # Process targets
        d = []
        if not (isinstance(targets,list) or isinstance(targets,tuple)):
            targets = [targets]
        self.connect()
        for target in targets:
            try:
                for x in self.connection.mlsd(target):
                    if x[0] == "." or x[0] == "..":
                        continue
                    x[1]["target"] = target
                    x[1]["name"] = x[0]
                    d.append(x[1])
            except:
                raise FileNotFoundError(f'''Could not retrieve list for directory {target} at the NCBI's FTP site.''')
        d = pd.DataFrame(d)
        return d

    # Mimick opening of local files
    def ftp_open(self, target,  mode='rt', avoid_collision=True, delete=True):
        '''
        Open a file stored at the NCBI FTP site.

        The file is first downloaded to the cache directory and then
        opened. Compressed files are uncompressed on the fly.

        Parameters
        ----------
        target : str
            Path of the file, relative to the server root.
        mode : str, default 'rt'
            Read mode used to open the file, such as ``r``, ``rt``
            or ``rb``.
        avoid_collision : bool, default True
            Avoid name collisions by adding a random suffix to the
            local file name.
        delete : bool, default True
            Whether the local copy is removed when the returned
            stream is closed. When False, the file remains in the
            cache directory.

        Returns
        -------
        file-like
            The open data stream.

        Examples
        --------
        Open a genome GBFF file:

        >>> from rotifer.db.ncbi import ftp as ncbiftp
        >>> conn = ncbiftp.connection()  # doctest: +SKIP
        >>> path = "genomes/all/GCA/900/547/725/GCA_900547725.1_UMGS1014/"
        >>> path += "GCA_900547725.1_UMGS1014_genomic.gbff.gz"
        >>> fh = conn.ftp_open(path, mode="rt")  # doctest: +SKIP
        '''
        from tempfile import _TemporaryFileWrapper
        import rotifer.core.functions as rcf
        self.connect()
        outfile = self.ftp_get(target, avoid_collision=avoid_collision)
        return _TemporaryFileWrapper(rcf.open_compressed(outfile, mode), outfile, delete)

class GenomeCursor(rotifer.db.methods.GenomeCursor, rotifer.db.parallel.SimpleParallelProcessCursor):
    """
    Fetch genome sequences from the NCBI FTP site.

    Genomes are located by their assembly accession, downloaded as
    GenBank flat files, verified against the published MD5 checksums
    and parsed into Bio.SeqRecord objects.

    Parameters
    ----------
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    batch_size : int, default 4
        Number of accessions per batch. The default may be changed
        by the module configuration.
    threads : int, default 15
        Number of processes used for parallel downloads.
    timeout : int, default 10
        Maximum time, in seconds, to wait for a server connection.
    cache : str, optional
        Directory for temporary files. Defaults to the rotifer
        cache directory.

    See Also
    --------
    rotifer.db.ncbi.GenomeCursor : delegator that also checks local
        mirrors

    Examples
    --------
    Load a sample of genomes:

    >>> from rotifer.db.ncbi import ftp
    >>> g = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> gc = ftp.GenomeCursor()  # doctest: +SKIP
    >>> genomes = gc.fetchall(g)  # doctest: +SKIP
    """
    def __init__(
            self,
            progress=True,
            tries=3,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
            timeout=10,
            cache=rotifer.config['cache'],
            *args, **kwargs
        ):
        threads = threads or _defaults['threads']
        super().__init__(progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        self.timeout = timeout
        self.cache = cache

    def open_genome(self, accession, assembly_reports=None):
        """
        Open the GBFF file of a genome from the NCBI FTP site.

        The file's MD5 checksum is downloaded first and the genome
        download is repeated until the checksum matches, at most
        ``tries`` times.

        Parameters
        ----------
        accession : str
            Genome assembly accession.
        assembly_reports : pandas.DataFrame, optional
            NCBI assembly summary table, as loaded by
            :func:`rotifer.db.ncbi.assemblies`. When given, the
            genome path is read from its ``ftp_path`` column instead
            of being searched on the server.

        Returns
        -------
        file-like or None
            The open data stream, with the assembly accession in an
            ``assembly`` attribute, or None when the genome is not
            found.

        Raises
        ------
        IOError
            If the checksum file cannot be parsed or the GBFF file
            cannot be opened.

        Examples
        --------
        Open and parse a genome:

        >>> from rotifer.db.ncbi import ftp
        >>> from Bio import SeqIO
        >>> gc = ftp.GenomeCursor()  # doctest: +SKIP
        >>> fh = gc.open_genome("GCA_900547725.1")  # doctest: +SKIP
        >>> s = [x for x in SeqIO.parse(fh, "genbank")]  # doctest: +SKIP
        >>> fh.close()  # doctest: +SKIP
        """
        import rotifer.core.functions as rcf
        from rotifer.db.ncbi import ftp as ncbiftp
        ftp = ncbiftp.connection(tries=self.tries, timeout=self.timeout, cache=self.cache)

        # find genome and download
        path = self.genome_path(accession, assembly_reports=assembly_reports)
        if len(path) == 0:
           return None

        # Download checksum
        md5url = "/".join([path[0],"md5checksums.txt"])
        ftp.connect()
        for attempt in range(0,self.tries):
            md5 = ftp.ftp_open(md5url, mode='rt', avoid_collision=True, delete=True)
            try:
                md5 = pd.read_csv(md5, sep=' +', names=['md5','filename'], engine="python")
                md5 = md5[md5.filename.fillna("_").str.contains('_genomic.gbff.gz')]
                md5 = md5.md5.iloc[0]
                if md5:
                    break
            except:
                raise IOError(f'''Parsing of checksum for {accession} failed.''')
        if not md5:
            return None

        # Download genome
        gz = None
        ftp.connect()
        for attempt in range(0,self.tries):
            gz = ftp.ftp_open("/".join(path),  mode='rt', avoid_collision=True, delete=True)
            try:
                md5gz = rcf.md5(gz.name)
            except:
                raise IOError(f'''Could not open GBFF for {accession}.''')
            if md5 == md5gz:
                gz.assembly = accession
                break

        # Return file object
        return gz

    def genome_path(self, accession, assembly_reports=None):
        """
        Find the path of a genome at the NCBI FTP site.

        When the accession matches several versions, the newest one
        is chosen.

        Parameters
        ----------
        accession : str
            Genome assembly accession.
        assembly_reports : pandas.DataFrame, optional
            NCBI assembly summary table, as loaded by
            :func:`rotifer.db.ncbi.assemblies`. When given, the path
            is read from its ``ftp_path`` column instead of being
            searched on the server.

        Returns
        -------
        tuple of str
            The directory and file name of the genome's GBFF file.
            Empty when the genome is not found.

        Examples
        --------
        >>> from rotifer.db.ncbi import ftp
        >>> gc = ftp.GenomeCursor()  # doctest: +SKIP
        >>> path = gc.genome_path("GCA_900547725.1")  # doctest: +SKIP
        >>> print("/".join(path))  # doctest: +SKIP
        """
        from rotifer.db.ncbi import NcbiConfig
        from rotifer.db.ncbi import ftp as ncbiftp
        ftp = ncbiftp.connection(tries=self.tries, timeout=self.timeout, cache=self.cache)
        path = ()

        # Extract genome path from assembly reports
        if isinstance(assembly_reports, pd.DataFrame) and not assembly_reports.empty:
            path = assembly_reports.query(f'assembly == "{accession}"')
            if not path.empty:
                path = path.ftp_path.iloc[0]
                path = path.replace(f'ftp://{NcbiConfig["ftpserver"]}','')
                path = (path,os.path.basename(path) + "_genomic.gbff.gz")

        # Retrieve genome path for newest version
        if len(path) == 0:
            path = accession[0:accession.find(".")].replace("_","")
            path = "/".join([ path[i : i + 3] for i in range(0, len(path), 3) ])
            path = f'/genomes/all/{path}'
            path = ftp.ftp_ls(path)
            if path.empty:
                return ()
            path = path.query(f'name.str.contains("{accession}")')
            path = path.sort_values(['name'], ascending=False).iloc[0]
            path = path.target + "/" + path['name']

            # Retrieve GBFF path
            if not path:
                return ()
            path = ftp.ftp_ls(path)
            if path.empty:
                return ()
            path = path[path['name'].str.contains(".gbff.gz")]
            if path.empty:
                return ()
            path = (path.target.iloc[0], path['name'].iloc[0])

        return path

    def genome_report(self, accession):
        """
        Fetch and parse a genome's assembly report.

        Parameters
        ----------
        accession : str
            Genome assembly accession.

        Returns
        -------
        tuple
            A pair ``(contigs, properties)``. ``contigs`` is a
            pandas.DataFrame listing the assembly's sequences.
            ``properties`` is a pandas.DataFrame, indexed by
            property name, holding the assembly metadata, similar
            to a row of :func:`rotifer.db.ncbi.assemblies`. When
            the genome is not found, ``contigs`` is an empty list.

        Examples
        --------
        >>> from rotifer.db.ncbi import ftp
        >>> gc = ftp.GenomeCursor()  # doctest: +SKIP
        >>> contigs, assembly = gc.genome_report("GCA_900547725.1")  # doctest: +SKIP
        """
        from rotifer.db.ncbi import ftp as ncbiftp
        ftp = ncbiftp.connection(tries=self.tries, timeout=self.timeout, cache=self.cache)

        # Column names
        arcolumn = f"""                Assembly name : assembly_name
                                       Organism name : organism_name
                                             Isolate : isolate
                                               Taxid : taxid
                                           BioSample : biosample
                                          BioProject : bioproject
                                           Submitter : submitter
                                                Date : submission_date
                                       Assembly type : assembly_type
                                        Release type : release_type
                                      Assembly level : assembly_level
                               Genome representation : representative
                                         WGS project : wgs
                                     Assembly method : assembly_method
                                     Genome coverage : genome_coverage
                               Sequencing technology : sequencing
                                     RefSeq category : refseq_category
                          GenBank assembly accession : genbank
                           RefSeq assembly accession : refseq
                                Excluded from RefSeq : excluded_from_refseq
    RefSeq assembly and GenBank assemblies identical : identical""".split("\n")
        arcolumn = [ x.strip().split(" : ") for x in arcolumn ]
        arcolumn = { x[0]:x[1] for x in arcolumn }

        # Opening data file
        ar = [['column','value'], ['assembly', accession]]
        sc = []
        path = self.genome_path(accession)
        if len(path):
            report = ftp.ftp_ls(path[0])
        else:
            return ([],pd.DataFrame(columns=ar[0]))
        report = report[report.name.str.contains("_assembly_report.txt")]
        report = path[0] + "/" + report.name.iloc[0]
        report = ftp.ftp_open(report)

        inar = True
        for row in report:
            row = row.strip()
            if row == "#" or row[0:2] == "##":
                inar = False
            elif row[0:15] == "# Sequence-Name":
                sc.append(['assembly'] + row[2:].split("\t"))
            elif inar and row[0:2] == "# ":
                ar.append(row[2:].split(":", maxsplit=1))
            elif row[0] != '#':
                sc.append([accession] + row.split("\t"))

        ar = pd.DataFrame(ar[1:], columns=ar[0])
        ar.value = ar.value.str.strip()
        ar.column = ar.column.replace(arcolumn)
        sc = pd.DataFrame(sc[1:], columns=sc[0])
        ar = ar.set_index("column")

        return sc, ar

class GenomeFeaturesCursor(rotifer.db.methods.GenomeFeaturesCursor, GenomeCursor):
    """
    Fetch genome annotation from the NCBI FTP site as dataframes.

    Genomes are downloaded like :class:`GenomeCursor` does, then
    converted to feature tables with one row per annotated feature.

    Parameters
    ----------
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
    batch_size : int, default 4
        Number of accessions per batch.
    threads : int, default 15
        Number of processes used for parallel downloads.
    timeout : int, default 10
        Maximum time, in seconds, to wait for a server connection.
    cache : str, optional
        Directory for temporary files. Defaults to the rotifer
        cache directory.

    Examples
    --------
    Load the feature tables of two genomes:

    >>> from rotifer.db.ncbi import ftp
    >>> g = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> gfc = ftp.GenomeFeaturesCursor()  # doctest: +SKIP
    >>> df = gfc.fetchall(g)  # doctest: +SKIP
    """
    def __init__(
            self,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=True,
            tries=3,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
            timeout=10,
            cache=rotifer.config['cache'],
            *args, **kwargs
        ):
        threads = threads or _defaults['threads']
        super().__init__(progress=progress, tries=tries, batch_size=batch_size, threads=threads, timeout=timeout, cache=cache, *args, **kwargs)
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable

class GeneNeighborhoodCursor(rotifer.db.methods.GeneNeighborhoodCursor, rotifer.db.parallel.GeneNeighborhoodCursor, GenomeFeaturesCursor):
    """
    Fetch gene neighborhoods from genomes at the NCBI FTP site.

    Target proteins are resolved to genome assemblies through
    identical protein group (IPG) reports, the genomes are
    downloaded from the FTP site and the annotated regions around
    each target gene are returned as dataframes.

    Parameters
    ----------
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
    batch_size : int, default 4
        Number of accessions per batch.
    threads : int, default 15
        Number of processes used for parallel downloads.
    timeout : int, default 10
        Maximum time, in seconds, to wait for a server connection.
    cache : str, optional
        Directory for temporary files. Defaults to the rotifer
        cache directory.

    See Also
    --------
    rotifer.db.ncbi.GeneNeighborhoodCursor : delegator that combines
        this backend with mirrors and Entrez

    Examples
    --------
    >>> from rotifer.db.ncbi import ftp
    >>> gnc = ftp.GeneNeighborhoodCursor(progress=True)  # doctest: +SKIP
    >>> df = gnc.fetchall(["EEE9598493.1"])  # doctest: +SKIP
    """
    def __init__(
            self,
            column = 'pid',
            before = 7,
            after = 7,
            min_block_distance = 0,
            strand = None,
            fttype = 'same',
            eukaryotes=False,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=True,
            tries=3,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
            timeout=10,
            cache=rotifer.config['cache'],
            *args, **kwargs
        ):

        threads = threads or _defaults['threads']

        super().__init__(
            column = column,
            before = before,
            after = after,
            min_block_distance = min_block_distance,
            strand = strand,
            fttype = fttype,
            eukaryotes = eukaryotes,
            exclude_type = exclude_type,
            autopid = autopid,
            codontable = codontable,
            progress = progress,
            tries = tries,
            batch_size = batch_size,
            threads = threads,
            *args, **kwargs
        )
        self.timeout = timeout
        self.cache = cache

# Is this library being used as a script?
if __name__ == '__main__':
    pass
