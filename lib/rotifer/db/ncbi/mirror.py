"""
Access genome data from a local mirror of the NCBI FTP site.

Cursors in this module behave like their
:mod:`rotifer.db.ncbi.ftp` counterparts but read GenBank flat files
from a local directory tree that mirrors the ``genomes`` section of
the NCBI FTP site. No network connection is used.

Configuration
-------------
The mirror location defaults to the ``mirror`` entry of the
:mod:`rotifer.db.ncbi` configuration, which falls back to the
``genomes`` directory under ``ROTIFER_DATA``.
"""

import os
import sys
import types
import typing
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
from copy import deepcopy

import rotifer
import rotifer.db.parallel
import rotifer.db.methods
from rotifer import GlobalConfig
from rotifer.db.ncbi import config as NcbiConfig
from rotifer.db.ncbi import utils as rdnu
from rotifer.core.functions import loadConfig
from rotifer.genome.utils import seqrecords_to_dataframe
logger = rotifer.logging.getLogger(__name__)

# Configuration
from rotifer.core.functions import loadConfig
_defaults = {
    "path": NcbiConfig['mirror'] or os.path.join(rotifer.config['data'],"genomes"),
    "batch_size": None,
    "threads": int(np.floor(os.cpu_count()/2)),
}
config = loadConfig(__name__, defaults = _defaults)

# Classes

class GenomeCursor(rotifer.db.methods.GenomeCursor, rotifer.db.parallel.SimpleParallelProcessCursor):
    """
    Fetch genome sequences from a local mirror of the NCBI genomes
    repository.

    Genomes are located by their assembly accession and parsed from
    the mirror's GenBank flat files into Bio.SeqRecord objects.

    Parameters
    ----------
    progress : bool, default False
        Whether to print a progress bar.
    tries : int, default 1
        Accepted for interface compatibility; local reads are
        attempted only once.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of processes used for parallel reads. Defaults to
        half of the CPU count.
    path : str, optional
        Path to a mirror of the genomes section of the NCBI FTP
        site. Contents are expected to be the same as, or a subset
        of, the ``genomes`` directory.

    See Also
    --------
    rotifer.db.ncbi.ftp.GenomeCursor : the same interface over FTP

    Examples
    --------
    >>> from rotifer.db.ncbi import mirror
    >>> gc = mirror.GenomeCursor()  # doctest: +SKIP
    >>> genomes = gc.fetchall(['GCA_900547725.1'])  # doctest: +SKIP
    """
    def __init__(
            self,
            progress=False,
            tries=1,
            batch_size = config["batch_size"],
            threads = config["threads"] or _defaults['threads'],
            path = config["path"],
            *args, **kwargs):
        threads = threads or _defaults['threads']
        super().__init__(progress=progress, tries=1, batch_size=batch_size, threads=threads)
        self.path = path

    def open_genome(self, accession, assembly_reports=None):
        """
        Open the GBFF file of a genome from the local mirror.

        Compressed files are uncompressed on the fly.

        Parameters
        ----------
        accession : str
            Genome assembly accession.
        assembly_reports : pandas.DataFrame, optional
            NCBI assembly summary table, as loaded by
            :func:`rotifer.db.ncbi.assemblies`. When given, the
            genome path is derived from its ``ftp_path`` column
            instead of being searched in the mirror tree.

        Returns
        -------
        file-like or None
            The open data stream, with the assembly accession in an
            ``assembly`` attribute, or None when the genome is not
            found.

        Examples
        --------
        Open and parse a genome:

        >>> from rotifer.db.ncbi import mirror
        >>> from Bio import SeqIO
        >>> gc = mirror.GenomeCursor()  # doctest: +SKIP
        >>> fh = gc.open_genome("GCA_900547725.1")  # doctest: +SKIP
        >>> s = [x for x in SeqIO.parse(fh, "genbank")]  # doctest: +SKIP
        >>> fh.close()  # doctest: +SKIP
        """
        import rotifer.core.functions as rcf

        # find genome and download
        path = self.genome_path(accession, assembly_reports=assembly_reports)
        if len(path) == 0:
           return None

        # Download genome
        path = os.path.join(*path)
        gz = rcf.open_compressed(path,  mode='rt')
        gz.assembly = accession

        # Return file object
        return gz

    def genome_path(self, accession, assembly_reports=None):
        """
        Find the path of a genome in the local NCBI mirror.

        When the accession matches several versions, the newest one
        is chosen.

        Parameters
        ----------
        accession : str
            Genome assembly accession.
        assembly_reports : pandas.DataFrame, optional
            NCBI assembly summary table, as loaded by
            :func:`rotifer.db.ncbi.assemblies`. When given, the path
            is derived from its ``ftp_path`` column instead of being
            searched in the mirror tree.

        Returns
        -------
        tuple
            The directory and file name of the genome's GBFF file.
            The file name is None when no GBFF file exists in the
            genome's directory.

        Raises
        ------
        FileNotFoundError
            If the mirror has no directory for the accession.
        IOError
            If the genome's directory cannot be read.

        Examples
        --------
        >>> from rotifer.db.ncbi import mirror
        >>> gc = mirror.GenomeCursor()  # doctest: +SKIP
        >>> path = gc.genome_path("GCA_900547725.1")  # doctest: +SKIP
        >>> print("/".join(path))  # doctest: +SKIP
        """
        from rotifer.db.ncbi import NcbiConfig
        path = ()

        # Extract genome path from assembly reports
        if isinstance(assembly_reports, pd.DataFrame) and not assembly_reports.empty:
            path = assembly_reports.query(f'assembly == "{accession}"')
            if not path.empty:
                path = path.ftp_path.iloc[0]
                path = path.replace(f'ftp://{NcbiConfig["ftpserver"]}/genomes/','')
                path = (os.path.join(self.path,path),os.path.basename(path) + "_genomic.gbff.gz")

        # Retrieve genome path for newest version
        if len(path) == 0:
            path = accession[0:accession.find(".")].replace("_","")
            path = [ path[i : i + 3] for i in range(0, len(path), 3) ]
            path = os.path.join(self.path,'all',*path)
            if os.path.exists(path):
                ls = os.listdir(path)
            else:
                raise FileNotFoundError(f'No directory {path} for {accession}')
            ls = [ x for x in sorted(ls) if accession in x and os.path.isdir(os.path.join(path,x)) ]
            if len(ls):
                ls = ls[-1] # Expected to be the latest version of the target genome
            else:
                raise FileNotFoundError(f'Empty directory for {accession} in {path}')
            path = os.path.join(path,ls)

            # Retrieve GBFF path
            try:
                ls = os.listdir(path)
            except:
                raise IOError(f'Unable to read directory for {accession} in {path}')
            ls = [ x for x in sorted(ls) if '.gbff.gz' in x ]
            if len(ls):
                ls = ls[0] # Only one GBFF is expected
            else:
                ls = None
            path = (path, ls)

        return path

    def genome_report(self, accession):
        """
        Fetch and parse a genome's assembly report from the mirror.

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
            the genome or its report is not found, ``contigs`` is
            an empty list.

        Examples
        --------
        >>> from rotifer.db.ncbi import mirror
        >>> gc = mirror.GenomeCursor()  # doctest: +SKIP
        >>> contigs, assembly = gc.genome_report("GCA_900547725.1")  # doctest: +SKIP
        """

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
        ar = [['column','value'], ['assembly', accession]]
        sc = []

        # Find report 
        path = self.genome_path(accession)
        if len(path):
            report = os.listdir(path[0])
        else:
            return ([],pd.DataFrame(columns=ar[0]))
        report = [ x for x in report if "_assembly_report.txt" in x ]
        if len(report):
            report = os.path.join(path[0],report[0])
        else:
            return ([],pd.DataFrame(columns=ar[0]))

        # Parse report
        report = open(report)
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
    Fetch genome annotation from a local mirror as dataframes.

    Genomes are read like :class:`GenomeCursor` does, then converted
    to feature tables with one row per annotated feature.

    Parameters
    ----------
    path : str, optional
        Path to a mirror of the genomes section of the NCBI FTP
        site.
    exclude_type : list of str, default ``['source', 'gene', 'mRNA']``
        Feature types to ignore.
    autopid : bool, default False
        Automatically set protein identifiers.
    codontable : str or int, default 'Bacterial'
        Codon table used when the data does not define one.
    progress : bool, default False
        Whether to print a progress bar.
    tries : int, default 1
        Accepted for interface compatibility; local reads are
        attempted only once.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of processes used for parallel reads.

    Examples
    --------
    Load the feature tables of two genomes:

    >>> from rotifer.db.ncbi import mirror
    >>> g = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> gfc = mirror.GenomeFeaturesCursor()  # doctest: +SKIP
    >>> df = gfc.fetchall(g)  # doctest: +SKIP
    """
    def __init__(
            self,
            path = config["path"],
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=1,
            batch_size = config["batch_size"],
            threads = config["threads"] or _defaults['threads'],
            *args, **kwargs
        ):
        threads = threads or _defaults['threads']
        super().__init__(progress=progress, tries=1, batch_size=batch_size, threads=threads, path=path, *args, **kwargs)
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable

class GeneNeighborhoodCursor(rotifer.db.methods.GeneNeighborhoodCursor, rotifer.db.parallel.GeneNeighborhoodCursor, GenomeFeaturesCursor):
    """
    Fetch gene neighborhoods from genomes in a local mirror.

    Target proteins are resolved to genome assemblies through
    identical protein group (IPG) reports, the genomes are read from
    the mirror and the annotated regions around each target gene are
    returned as dataframes.

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
    path : str, optional
        Path to a mirror of the genomes section of the NCBI FTP
        site.
    exclude_type : list of str, default ``['source', 'gene', 'mRNA']``
        Feature types to ignore.
    autopid : bool, default False
        Automatically set protein identifiers.
    codontable : str or int, default 'Bacterial'
        Codon table used when the data does not define one.
    progress : bool, default False
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to process each batch.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of processes used for parallel reads.

    See Also
    --------
    rotifer.db.ncbi.GeneNeighborhoodCursor : delegator that combines
        this backend with the FTP and Entrez backends

    Examples
    --------
    >>> from rotifer.db.ncbi import mirror
    >>> gnc = mirror.GeneNeighborhoodCursor(progress=True)  # doctest: +SKIP
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
            path = config["path"],
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=3,
            batch_size = config["batch_size"],
            threads = config["threads"] or _defaults['threads'],
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
        self.path = path

# Is this library being used as a script?
if __name__ == '__main__':
    pass
