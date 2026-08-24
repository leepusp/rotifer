"""
Access NCBI databases through the Entrez E-utilities.

Cursors in this module download sequences, identical protein group
(IPG) reports, taxonomy records, nucleotide feature tables and gene
neighborhoods with Biopython's :mod:`Bio.Entrez` interface to the
NCBI E-utilities web service.

Network, authentication and rate limits
---------------------------------------
Every request sends the user email registered in the
:mod:`rotifer.db.ncbi` configuration and the API key from the
``NCBI_API_KEY`` environment variable, when set. Without an API key
NCBI allows 3 requests per second and the cursors cap themselves at
3 simultaneous threads; with a key the cap is 10 threads. No data is
cached locally.
"""

# Import external modules
import os
import sys
import types
import socket
import typing
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from copy import deepcopy

# Import submodules
import rotifer
from rotifer import GlobalConfig
import rotifer.db.core
import rotifer.db.parallel
import rotifer.db.methods
from rotifer.db.ncbi import config as NcbiConfig
from rotifer.db.ncbi import utils as rdnu
from rotifer.core.functions import loadConfig
from rotifer.genome.utils import seqrecords_to_dataframe
logger = rotifer.logging.getLogger(__name__)

# Configuration
_defaults = {
    'batch_size': 20,
    "maxgetitem": 200,
    "threads": 10,
}
config = loadConfig(__name__, defaults = _defaults)

class SequenceCursor(rotifer.db.methods.SequenceCursor, rotifer.db.parallel.SimpleParallelProcessCursor):
    """
    Fetch annotated sequences from NCBI with the E-utilities.

    Sequences are downloaded with EFetch and parsed from GenBank
    format, the most richly annotated format NCBI provides.

    Parameters
    ----------
    database : str, default 'nucleotide'
        A valid NCBI sequence database name, such as ``protein`` or
        ``nucleotide``.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, default 20
        Number of accessions per batch. The default may be changed
        by the module configuration.
    threads : int, default 10
        Number of simultaneous threads. Capped at 3 without an NCBI
        API key and at 10 with one.

    See Also
    --------
    FastaCursor : faster download without annotations
    rotifer.db.ncbi.SequenceCursor : delegator that combines backends

    Examples
    --------
    Fetch a protein sequence:

    >>> from rotifer.db.ncbi import entrez
    >>> sc = entrez.SequenceCursor(database="protein")  # doctest: +SKIP
    >>> seqrec = sc.fetchall("YP_009724395.1")  # doctest: +SKIP

    Fetch several nucleotide entries:

    >>> import sys
    >>> from Bio import SeqIO
    >>> query = ['CP084314.1', 'NC_019757.1', 'AAHROG010000026.1']
    >>> sc = entrez.SequenceCursor(database="nucleotide")  # doctest: +SKIP
    >>> for seqrec in sc.fetchone(query):  # doctest: +SKIP
    ...     SeqIO.write(seqrec, sys.stdout, "genbank")
    """
    def __init__(
            self,
            database="nucleotide",
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
            *args, **kwargs):
        super().__init__(progress=progress, *args, **kwargs)
        self.maxgetitem = config['maxgetitem']
        self._tries = tries
        self.tries = 1
        self.database = database
        self.sleep_between_tries = sleep_between_tries
        self.batch_size = batch_size
        self.threads = threads or _defaults['threads']
        if self.threads > 3:
            if NcbiConfig['api_key']:
                if self.threads > 10:
                    self.threads = 10
            else:
                self.threads = 3

        # Register rules for giving up
        self.giveup.update(["HTTP Error 400"])

        # Private attributes (may be overloaded by children)
        self._rettype = 'gbwithparts'
        self._format = 'genbank'
        self._retmode = 'text'

    def parser(self, stream, accession):
        """
        Parse an EFetch stream into sequence records.

        Parameters
        ----------
        stream : file-like
            Open stream returned by :meth:`fetcher`.
        accession : str or iterable of str
            Database entry identifiers, kept for interface
            compatibility.

        Returns
        -------
        list of Bio.SeqRecord.SeqRecord
            The parsed records.
        """
        stack = []
        for s in SeqIO.parse(stream, self._format):
            stack.append(s)
        return stack

    def fetcher(self, accession):
        """
        Run EFetch and return the response stream.

        Parameters
        ----------
        accession : str or iterable of str
            Database entry identifiers.

        Returns
        -------
        file-like
            The open EFetch response.
        """
        from Bio import Entrez
        Entrez.email = NcbiConfig["email"]
        Entrez.api_key = NcbiConfig["api_key"]
        targets = self.parse_ids(accession, as_string=True)
        return Entrez.efetch(
                db=self.database,
                rettype=self._rettype,
                retmode=self._retmode,
                id=",".join(targets),
                max_tries=self._tries,
                sleep_between_tries=self.sleep_between_tries,
        )

    def __getitem__(self, accession):
        """
        Download data for one batch of entries, dictionary style.

        Entries that fail with a recoverable error are retried once
        in a second pass.

        Parameters
        ----------
        accession : str or iterable of str
            NCBI database entry identifiers. A string may contain
            several accessions separated by commas.

        Returns
        -------
        Bio.SeqRecord.SeqRecord or list of Bio.SeqRecord.SeqRecord
            A single record when a single accession yields a single
            record, otherwise a list.
        """
        targets = sorted(list(self.parse_ids(accession)))
        objlist = []
        batch = [",".join(targets)]
        for attempt in range(0,2):
            #logger.debug(f'Process {os.getpid()}, Attempt: {attempt}, processing {len(targets)} accessions by sending {len(batch)} strings') 
            for accs in batch:
                result = super().__getitem__(accs)
                if isinstance(result, types.NoneType):
                    continue
                elif isinstance(result,list):
                    objlist.extend(result)
                else:
                    objlist.append(result)
            batch = [ k for k, v in self._missing.items() if v[-1] == True ]
            if len(batch) == 0:
                break
        if len(targets) == 1 and len(objlist) == 1:
            objlist = objlist[0]
        return objlist

class FastaCursor(SequenceCursor):
    """
    Fetch sequences from NCBI as FASTA, without annotations.

    Identical to :class:`SequenceCursor` except that data is
    downloaded in FASTA format, which is faster but carries no
    annotation.

    Parameters
    ----------
    database : str, default 'protein'
        A valid NCBI sequence database name.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, default 20
        Number of accessions per batch.
    threads : int, default 10
        Number of simultaneous threads. Capped at 3 without an NCBI
        API key and at 10 with one.

    Examples
    --------
    >>> from rotifer.db.ncbi import entrez
    >>> fc = entrez.FastaCursor()  # doctest: +SKIP
    >>> seqs = fc.fetchall(["YP_009724395.1"])  # doctest: +SKIP
    """

    def __init__(self,
            database="protein",
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
            *args, **kwargs):
        threads = threads or _defaults['threads']
        super().__init__(database=database, progress=progress, tries=tries, sleep_between_tries=sleep_between_tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        self._rettype = "fasta"
        self._format = 'fasta'

class IPGCursor(rotifer.db.methods.IPGCursor, SequenceCursor):
    """
    Fetch identical protein group (IPG) reports from NCBI.

    IPG reports are downloaded with EFetch from the ``ipg`` database
    and returned as dataframes, one row per identical sequence, with
    the genomic coordinates of the encoding nucleotide sequences.

    Parameters
    ----------
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, default 20
        Number of accessions per batch.
    threads : int, default 10
        Number of simultaneous threads. Capped at 3 without an NCBI
        API key and at 10 with one.

    Notes
    -----
    Three columns are added to the original NCBI report:
    ``order`` (position of each row within its IPG), ``is_query``
    (whether the row's protein was part of the query) and
    ``representative`` (the query protein that represents the IPG).

    Examples
    --------
    >>> from rotifer.db.ncbi import entrez
    >>> ic = entrez.IPGCursor()  # doctest: +SKIP
    >>> ipgs = ic.fetchall(["WP_063732599.1"])  # doctest: +SKIP
    """

    def __init__(self,
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
        ):
        threads = threads or _defaults['threads']
        super().__init__(database="ipg", progress=progress, tries=tries, sleep_between_tries=sleep_between_tries, batch_size=batch_size, threads=threads)
        self._rettype = "ipg"
        self._columns = ['id','ipg_source','nucleotide','start','stop','strand','pid','description','ipg_organism','strain','assembly']
        self._added_columns = ['order','is_query','representative']
        self.giveup.update(["no IPG","No IPG"])

    def _seqrecords_to_ipg(self, seqrecords):
        """
        Build an IPG-like dataframe from GenPept records.

        Records with identical sequences are grouped under negative
        IPG identifiers and coordinates are read from their
        ``coded_by`` qualifiers.

        Parameters
        ----------
        seqrecords : iterable of Bio.SeqRecord.SeqRecord
            GenPept records.

        Returns
        -------
        pandas.DataFrame
            A dataframe with the same columns as a parsed IPG
            report.
        """
        ipg = dict()
        representative = dict()
        order = dict()
        source = "genpept"
        ipgFromGenPept = []
        for x in seqrecords:
            seq = str(x.seq)
            if seq not in ipg:
                order[seq] = 0
                representative[seq] = x.id
                if len(ipg) == 0:
                    ipg[seq] = -1
                else:
                    ipg[seq] = min(ipg.values()) - 1
            else:
                order[seq] += 1
            desc = x.description.split("[")
            if len(desc) > 1:
                org = desc[1].replace("]","")
                desc = desc[0]
            else:
                org = np.nan
            strain = np.nan
            for f in x.features:
                for k in f.qualifiers:
                    if "strain" in f.qualifiers:
                        strain = f.qualifiers['strain'][0]
                    if not (f.type == "CDS" and k == "coded_by"):
                        continue
                    for v in f.qualifiers[k]:
                        strand =  "-" if "complement" in v else "+"
                        coord = v.replace('complement(',"").replace('join(',"").replace(")","").replace("..",":").replace(",",":").split(":")
                        acc = ",".join(pd.Series([ y.strip() for y in coord if "." in y ]).unique())
                        coord = [ int(y.replace(">","").replace("<","")) for y in coord if "." not in y ]
                        ipgFromGenPept.append([ipg[seq],source,acc,min(coord),max(coord),strand,x.id, desc,org,strain,np.nan,1,order[seq],representative[seq]])
        ipgFromGenPept = pd.DataFrame(ipgFromGenPept, columns=self._columns + self._added_columns)
        return ipgFromGenPept

    def parser(self, stream, accession):
        """
        Parse an EFetch stream into an IPG dataframe.

        Rows with malformed identifiers are dropped, query proteins
        are flagged in the ``is_query`` column, the original row
        order is recorded and each IPG is annotated with the query
        protein that represents it.

        Parameters
        ----------
        stream : file-like
            Open stream returned by the fetcher.
        accession : str or iterable of str
            The queried protein accessions.

        Returns
        -------
        pandas.DataFrame
            The parsed report.

        Raises
        ------
        ValueError
            If no valid IPG rows remain after error removal.
        """
        targets = self.parse_ids(accession)
        ipg = pd.read_csv(stream, sep='\t', names=self._columns, header=0).drop_duplicates()

        # Make sure all IPG ids are numbers
        numeric = pd.to_numeric(ipg.id, errors="coerce")
        errors = numeric.isna()
        if errors.any():
            logger.debug(f'Errors in IPG for accessions {", ".join(accession)}:\n{ipg[errors].to_string()}')
        ipg.id = numeric
        ipg = ipg[~errors]
        if ipg.empty:
            error = f'After removing errors, no IPG reports were found!'
            raise ValueError(error)

        # Register query proteins
        ipg['is_query'] = ipg.pid.isin(targets).astype(int).to_list()

        # Register original order of the table's rows
        o = pd.Series(range(1, len(ipg) + 1))
        c = ipg.id.map(o.groupby(ipg.id).min().to_dict())
        #c = pd.Series(np.where(ipg.id != ipg.id.shift(1), o.values, pd.NA)).ffill()
        ipg['order'] = (o - c).values

        # Annotate representatives
        if len(targets) > 1: # Many queries
            #  Register first query protein as representative
            rep = ipg.query('is_query == 1').drop_duplicates('id', keep='first')
            rep = rep.set_index('id').pid.to_dict()
            ipg['representative'] = ipg['id'].map(rep)
            # Remove IPGs with no known query
            ipg = ipg[ipg.representative.notna()]
        else: # One query
            ipg['representative'] = list(targets)[0]

        # Register all accessions found and return sliced DataFrame
        #return [ x[1].copy() for x in ipg.groupby('id') ]
        return ipg

    def fetchone(self,accessions):
        """
        Iterate over IPG reports as they are retrieved.

        Reports whose accessions were all seen in previously yielded
        reports are skipped.

        Parameters
        ----------
        accessions : list of str
            NCBI protein accessions.

        Yields
        ------
        pandas.DataFrame
            One IPG report per batch.
        """
        seen = set()
        for ipg in super().fetchone(accessions):
            if len(ipg) == 0:
                continue
            ids = self.getids(ipg).intersection(accessions)
            if seen.issuperset(ids):
                continue
            seen.update(self.getids(ipg))
            yield ipg

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
            The concatenated reports. Empty, but with the expected
            columns, when nothing could be retrieved.
        """
        targets = self.parse_ids(accessions)
        df = list(self.fetchone(targets))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = pd.DataFrame(columns = self._columns + self._added_columns)
        return df

class TaxonomyCursor(SequenceCursor):
    """
    Fetch NCBI Taxonomy records with the E-utilities.

    Taxonomy entries are downloaded as XML with EFetch and returned
    as dataframes.

    Parameters
    ----------
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, default 20
        Number of accessions per batch.
    threads : int, default 10
        Number of simultaneous threads. Capped at 3 without an NCBI
        API key and at 10 with one.

    Attributes
    ----------
    columns : list of str
        Columns of the returned dataframes: ``taxid``, ``organism``,
        ``superkingdom``, ``lineage``, ``classification`` and
        ``alternative_taxids``.

    See Also
    --------
    rotifer.db.ncbi.TaxonomyCursor : delegator that prefers the
        local ETE toolkit database

    Examples
    --------
    >>> from rotifer.db.ncbi import entrez
    >>> tc = entrez.TaxonomyCursor()  # doctest: +SKIP
    >>> t = tc.fetchall([2599])  # doctest: +SKIP
    """

    def __init__(self,
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
        ):
        threads = threads or _defaults['threads']
        super().__init__(database="taxonomy",progress=progress,tries=tries,sleep_between_tries=sleep_between_tries,batch_size=batch_size,threads=threads)
        self._rettype = "full"
        self._retmode = 'xml'
        self.columns = ['taxid','organism','superkingdom','lineage','classification','alternative_taxids']
        self.giveup.update(["no taxonomy"])

    def getids(self,obj):
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
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            ids.update(set(o.taxid))
            ids.update(set(o.alternative_taxids.str.split(",").explode().dropna()))
        return ids

    def parser(self, stream, accession):
        """
        Parse an EFetch XML stream into a taxonomy dataframe.

        Parameters
        ----------
        stream : file-like
            Open stream returned by the fetcher.
        accession : str or iterable of str
            The queried taxonomy identifiers, kept for interface
            compatibility.

        Returns
        -------
        list of pandas.DataFrame
            A single element list holding the parsed dataframe.

        Raises
        ------
        ValueError
            If the stream contains no taxonomy records.
        """
        from Bio import Entrez
        from rotifer.taxonomy.utils import lineage
        taxdf = [ x for x in Entrez.parse(stream) ]
        if len(taxdf) == 0:
            raise ValueError(f'Empty data stream: no taxonomy')
        taxdf = pd.DataFrame(taxdf)
        if "AkaTaxIds" in taxdf.columns:
            taxdf["alternative_taxids"] = taxdf["AkaTaxIds"].fillna("").map(lambda x: ",".join(x))
            taxdf.alternative_taxids = np.where(taxdf.alternative_taxids == "", taxdf.TaxId, taxdf.alternative_taxids)
        elif "TaxId" in taxdf:
            taxdf["alternative_taxids"] = taxdf.TaxId
        taxdf['superkingdom'] = taxdf.Lineage.str.replace("cellular organisms; ","").str.split("; ", expand=True)[0]
        taxdf.rename({'Lineage':'classification', 'TaxId':'taxid', 'ScientificName':'organism'}, axis=1, inplace=1)
        taxdf['lineage'] = lineage(taxdf.classification)
        taxdf = taxdf.loc[:,self.columns].applymap(lambda x: str(x))
        stream.close()
        return [taxdf]

    def fetchall(self, accessions):
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
            ``columns``.
        """
        df = list(self.fetchone(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = pd.DataFrame(self.columns)
        return df

class NucleotideFeaturesCursor(SequenceCursor):
    """
    Fetch nucleotide annotation as feature tables.

    Nucleotide entries are downloaded in GenBank format and
    converted to dataframes with one row per annotated feature.

    Parameters
    ----------
    exclude_type : list of str, default ``['source', 'gene', 'mRNA']``
        Feature types to ignore.
    autopid : bool, default False
        Automatically set protein identifiers.
    assembly : str, optional
        Assembly accession assigned to every parsed feature.
    codontable : str or int, default 'Bacterial'
        Codon table used when the data does not define one.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 3
        Number of attempts to download data.
    sleep_between_tries : int, default 1
        Number of seconds to wait between download attempts.
    batch_size : int, default 20
        Number of accessions per batch.
    threads : int, default 10
        Number of simultaneous threads. Capped at 3 without an NCBI
        API key and at 10 with one.

    Examples
    --------
    >>> from rotifer.db.ncbi import entrez
    >>> nfc = entrez.NucleotideFeaturesCursor()  # doctest: +SKIP
    >>> df = nfc.fetchall(['CP084314.1'])  # doctest: +SKIP
    """

    def __init__(
            self,
            exclude_type = ['source','gene','mRNA'],
            autopid = False,
            assembly = None,
            codontable= 'Bacterial',
            progress = True,
            tries = 3,
            sleep_between_tries=1,
            batch_size = config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
        ):
        threads = threads or _defaults['threads']
        super().__init__(
                database='nucleotide',
                progress=progress,
                tries=tries,
                sleep_between_tries=sleep_between_tries,
                batch_size=batch_size,
                threads=threads
        )
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.assembly = assembly
        self.codontable = codontable

    def getids(self,obj):
        """
        Extract nucleotide accessions from feature tables.

        Parameters
        ----------
        obj : pandas.DataFrame or list of pandas.DataFrame
            Feature tables produced by the cursor.

        Returns
        -------
        set of str
            The values of the ``nucleotide`` column.
        """
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            ids.update(o.nucleotide.dropna().to_list())
        return ids

    def parser(self, stream, accession):
        """
        Parse an EFetch GenBank stream into feature tables.

        Parameters
        ----------
        stream : file-like
            Open stream returned by the fetcher.
        accession : str or iterable of str
            The queried nucleotide accessions, kept for interface
            compatibility.

        Returns
        -------
        list of pandas.DataFrame
            One feature table per nucleotide sequence.
        """
        stream = SeqIO.parse(stream, self._format)
        stream = seqrecords_to_dataframe(stream, exclude_type=self.exclude_type, autopid=self.autopid, assembly=self.assembly, codontable=self.codontable)
        stream = [ x[1].copy() for x in stream.groupby('nucleotide') ]
        return stream

    def fetchall(self, accessions):
        """
        Fetch the feature tables of all entries as one dataframe.

        Parameters
        ----------
        accessions : list of str
            NCBI nucleotide accessions.

        Returns
        -------
        pandas.DataFrame
            The concatenated feature tables. Empty when nothing
            could be retrieved.
        """
        df = list(self.fetchone(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = seqrecords_to_dataframe([])
        return df

class GeneNeighborhoodCursor(rotifer.db.methods.GeneNeighborhoodCursor, NucleotideFeaturesCursor):
    """
    Fetch gene neighborhoods from nucleotide annotation.

    Target proteins are first resolved to nucleotide sequences
    through identical protein group (IPG) reports, and the annotated
    regions around each target are then extracted from the
    nucleotide entries downloaded with EFetch. This backend does not
    require a genome assembly, so it also covers proteins encoded on
    unassembled sequences.

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
        Whether to process eukaryotic nucleotide sequences.
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
    batch_size : int, default 20
        Number of accessions per batch.
    threads : int, default 10
        Number of simultaneous threads. Capped at 3 without an NCBI
        API key and at 10 with one.

    See Also
    --------
    rotifer.db.ncbi.GeneNeighborhoodCursor : delegator that combines
        this backend with the FTP and mirror backends

    Examples
    --------
    >>> from rotifer.db.ncbi import entrez
    >>> gnc = entrez.GeneNeighborhoodCursor(progress=True)  # doctest: +SKIP
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
            eukaryotes = False,
            exclude_type = ['source','gene','mRNA'],
            autopid = False,
            codontable = 'Bacterial',
            progress = True,
            tries = 3,
            sleep_between_tries = 1,
            batch_size = config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
        ):

        threads = threads or _defaults['threads']

        super().__init__(
            exclude_type = exclude_type,
            autopid = autopid,
            codontable = codontable,
            progress = progress,
            tries = tries,
            sleep_between_tries=sleep_between_tries,
            batch_size = batch_size,
            threads = threads,
        )
        self.column = column
        self.before = before
        self.after = after
        self.min_block_distance = min_block_distance
        self.strand = strand
        self.fttype = fttype
        self.eukaryotes = eukaryotes
        self.giveup.update(["HTTP Error 400"])
        self.giveup.update(["no IPG","No IPG"])
        if not eukaryotes:
            self.giveup.update(["Eukaryot","eukaryot"])

    def getids2(self, obj, *args, **kwargs):
        """
        List protein identifiers per nucleotide sequence.

        Parameters
        ----------
        obj : pandas.DataFrame
            A neighborhood dataframe.

        Returns
        -------
        list of str
            The identifiers found in the ``pid`` and, when present,
            ``replaced`` columns.
        """
        columns = ['pid']
        if 'replaced' in obj.columns:
            columns.append('replaced')
        ids = obj.melt(id_vars=['nucleotide'], value_vars=columns, value_name='id', var_name='type')
        ids.drop('type', axis=1, inplace=True)
        ids.set_index('nucleotide', inplace=True)
        ids.drop_duplicates(inplace=True)
        return ids.id.tolist()

    def __getitem__(self, accessions, ipgs=None):
        """
        Find gene neighborhoods around one or more target proteins.

        Parameters
        ----------
        accessions : str or iterable of str
            NCBI protein accessions.
        ipgs : pandas.DataFrame, optional
            Precomputed identical protein group reports. When not
            given, they are downloaded with :class:`IPGCursor`.

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
            The neighborhoods found. Empty when nothing could be
            retrieved.
        """
        objlist = seqrecords_to_dataframe([])

        # Make sure no identifiers are used twice
        targets = self.parse_ids(accessions)

        if isinstance(ipgs,types.NoneType):
            from rotifer.db.ncbi import entrez
            ic = entrez.IPGCursor(progress=False, tries=self.tries, batch_size=self.batch_size, threads=self.threads)
            ipgs = ic.fetchall(targets)
            targets = targets - ic.missing_ids()
            self.update_missing(data=ic.remove_missing())
        ipgids = set(ipgs[ipgs.pid.isin(targets) | ipgs.representative.isin(targets)].id)
        ipgs = ipgs[ipgs.id.isin(ipgids) & (ipgs.assembly.notna() | ipgs.nucleotide.notna())]
        best = rdnu.best_ipgs(ipgs)
        best = best[best.nucleotide.notna()]
        ipgs = ipgs[ipgs.nucleotide.isin(best.nucleotide)]
        missing = targets - self.getids(ipgs)
        if missing:
            self.update_missing(missing, error="No IPGs", retry=False)
            targets = targets - missing
            if len(targets) == 0:
                return objlist

        # Identify DNA data
        nucleotides = ipgs.filter(['nucleotide','pid','representative'])
        nucleotides = nucleotides.drop_duplicates(ignore_index=True)
        nucleotides = nucleotides.groupby('nucleotide').apply(lambda x: x.set_index('pid').representative.to_dict())
        nucleotides = nucleotides.to_dict()
        #assemblies, nucleotides = rdnu.ipgs_to_dicts(ipgs)

        # Download and parse
        objlist = []
        for accession in nucleotides.keys():
            expected = set([ y for x in nucleotides[accession].items() for y in x ])
            expected = targets.intersection(expected)

            obj = None
            for attempt in range(0,self.tries):
                # Download and open data file
                error = None
                stream = None
                try:
                    stream = self.fetcher(accession)
                except RuntimeError:
                    error = f'Runtime error for nucleotide {accession}: {sys.exc_info()[1]}'
                    logger.debug(error)
                    continue
                except ValueError:
                    error = f'Value error for nucleotide {accession}: {sys.exc_info()[1]}'
                    logger.debug(error)
                    break
                except:
                    error = f'Failed to download nucleotide {accession}: {sys.exc_info()[1]}'
                    logger.debug(error)
                    continue

                if error or isinstance(stream, types.NoneType):
                    if self.update_missing(expected, error):
                        continue
                    else:
                        break

                # Use parser to process results
                try:
                    obj = self.parser(stream, accession, nucleotides[accession])
                    break
                except ValueError:
                    error = f'Value error for nucleotide {accession}: {sys.exc_info()[1]}'
                    logger.debug(error)
                    break
                except:
                    error = f"Failed to parse nucleotide {accession}: {sys.exc_info()[1]}"
                    logger.debug(error)

                # See if the error indicates we should give up
                if error:
                    if self.update_missing(expected, error):
                        continue
                    else:
                        break

            if isinstance(obj, types.NoneType):
                self.update_missing(expected, error)
            elif len(obj) == 0:
                error = f'No anchors in nucleotide sequence {accession}'
                self.update_missing(expected, error)
            else:
                objlist.extend(obj)

        # No data?
        if len(objlist) == 0:
            return seqrecords_to_dataframe([])

        # Concatenate and evaluate
        objlist = pd.concat(objlist, ignore_index=True)

        # Return data
        if len(objlist) > 0:
            self.remove_missing(self.getids(objlist))
        return objlist

    def parser(self, stream, accession, proteins):
        """
        Extract gene neighborhoods from a nucleotide stream.

        Parameters
        ----------
        stream : file-like
            Open stream returned by the fetcher.
        accession : str
            Nucleotide accession.
        proteins : dict
            Mapping of each target protein accession to its
            identical protein group representative.

        Returns
        -------
        list of rotifer.genome.data.NeighborhoodDF
            The neighborhoods found, one dataframe per nucleotide
            sequence, each with a ``replaced`` column mapping
            proteins to the queries they represent.

        Raises
        ------
        ValueError
            If the sequence is eukaryotic and ``eukaryotes`` is
            False.
        """
        stack = []
        for df in super().parser(stream, accession):
            taxonomy = df.classification.fillna("").iloc[0].split(";")
            if (not self.eukaryotes) and "Eukaryota" in taxonomy:
                raise ValueError(f"Eukaryotic nucleotide sequence {accession} ignored.")
            df = df.neighbors(
                df[self.column].isin(proteins.keys()),
                before = self.before,
                after = self.after,
                min_block_distance = self.min_block_distance,
                strand = self.strand,
                fttype = self.fttype,
            )
            df['replaced'] = df.pid.replace(proteins)
            stack.append(df)
        return stack

    def worker(self, chunk):
        """
        Process one batch of nucleotide sequences in a worker
        process.

        Parameters
        ----------
        chunk : list of tuple
            Batch entries produced by :meth:`splitter`. Each entry
            holds a set of target protein accessions and the IPG
            rows of one nucleotide sequence.

        Returns
        -------
        dict
            A dictionary with two keys: ``result``, a list of
            neighborhood dataframes (one per gene block), and
            ``missing``, the registry of entries that could not be
            retrieved.
        """
        result = []
        for args in chunk:
            df = self.__getitem__(*args)
            if len(df) == 0:
                continue
            for x in df.groupby('block_id'):
                result.append(x[1])
        return {"result":result,"missing":self.remove_missing()}

    def splitter(self, accessions, ipgs):
        """
        Group targets and their IPG rows into per-nucleotide
        batches.

        Parameters
        ----------
        accessions : set of str
            Target protein accessions.
        ipgs : pandas.DataFrame
            Identical protein group reports restricted to the
            selected nucleotide sequences.

        Returns
        -------
        list of list of tuple
            Batches of ``(proteins, ipg_rows)`` pairs, one pair per
            nucleotide sequence, at most ``batch_size`` pairs per
            batch.
        """
        size = self.batch_size
        if size == None or size == 0:
            size = max(int(ipgs.nucleotide.nunique()/self.threads),1)
        batch = []
        for x, y in ipgs.groupby('nucleotide'):
            proteins = accessions.intersection(self.getids(ipgs))
            batch.append((proteins, y.copy()))
        batch = [ batch[x:x+size] for x in range(0,len(batch),size) ]
        return batch

    def nucleotide_ids(self, obj):
        """
        Extract nucleotide accessions from neighborhood dataframes.

        Parameters
        ----------
        obj : pandas.DataFrame or list of pandas.DataFrame
            Neighborhood dataframes.

        Returns
        -------
        set of str
            The unique values of the ``nucleotide`` column.
        """
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            ids.update(o.nucleotide.unique().tolist())
        return ids

    def fetchone(self, accessions, ipgs=None):
        """
        Iterate over gene neighborhoods as downloads complete.

        Results are yielded in completion order, not input order.

        Parameters
        ----------
        accessions : list of str
            NCBI protein identifiers.
        ipgs : pandas.DataFrame, optional
            Precomputed identical protein group reports, used to
            avoid downloading IPGs again.

        Yields
        ------
        rotifer.genome.data.NeighborhoodDF
            One dataframe per retrieved gene neighborhood.
        """
        from concurrent.futures import ProcessPoolExecutor, as_completed

        # Make sure no identifiers are used twice
        targets = self.parse_ids(accessions)

        # Make sure we have IPGs
        if isinstance(ipgs,types.NoneType):
            from rotifer.db.ncbi import entrez
            if self.progress:
                logger.warning(f'Downloading IPGs for {len(targets)} proteins...')
            ic = entrez.IPGCursor(progress=self.progress, tries=self.tries, threads=self.threads)
            ipgs = ic.fetchall(targets)
            targets = targets - ic.missing_ids()
            self.update_missing(data=ic.remove_missing())
        ipgs = ipgs[ipgs.pid.isin(targets) | ipgs.representative.isin(targets)]
        if len(ipgs) == 0:
            self.update_missing(targets,"No IPGs to match a nucleotide sequence")
            return [seqrecords_to_dataframe([])]
        missing = targets - self.getids(ipgs)
        if missing:
            self.update_missing(missing,"No IPGs")
            targets = targets - missing
        nucleotides = rdnu.best_ipgs(ipgs)
        nucleotides = nucleotides[nucleotides.nucleotide.notna()]
        nucleotides = ipgs[ipgs.nucleotide.isin(nucleotides.nucleotide)]
        if len(nucleotides) == 0:
            return [seqrecords_to_dataframe([])]

        # Split jobs and execute
        todo = set(nucleotides.nucleotide.unique())
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            if self.progress:
                pids = set(nucleotides.pid).union(nucleotides.representative)
                pids = len(pids.intersection(targets))
                logger.warn(f'Downloading {len(todo)} nucleotides for {pids} proteins...')
                p = tqdm(total=len(todo), initial=0)
            tasks = []
            missing = self.remove_missing()
            for chunk in self.splitter(targets, nucleotides):
                tasks.append(executor.submit(self.worker, chunk))
            self.update_missing(data=missing)
            completed = set()
            for x in as_completed(tasks):
                data = x.result()
                for acc in completed.intersection(data['missing'].keys()):
                    data['missing'].pop(acc, None)
                self.update_missing(data=data['missing'])
                for obj in data['result']:
                    found = targets.intersection(self.getids(obj))
                    self.remove_missing(found)
                    done = todo.intersection(self.nucleotide_ids(obj)) - completed
                    if  len(done) > 0:
                        completed.update(done)
                        if self.progress:
                            p.update(len(done))
                    todo = todo - done
                    yield obj

def elink(accessions, dbfrom="protein", dbto="taxonomy", linkname=None):
    """
    Find related database entries with the ELink E-utility.

    Parameters
    ----------
    accessions : str or list of str
        NCBI accessions to search links for.
    dbfrom : str, default 'protein'
        Name of the source database.
    dbto : str, default 'taxonomy'
        Name of the target database.
    linkname : str, optional
        Type of link between `dbfrom` and `dbto`. When not set,
        ``{dbfrom}_{dbto}`` is used.

    Returns
    -------
    pandas.DataFrame
        One row per link, with columns ``qacc``, ``quid``,
        ``dbfrom``, ``linkname``, ``dbto`` and ``tuid``.

    Examples
    --------
    Find the taxonomy entry of a protein:

    >>> from rotifer.db.ncbi import entrez
    >>> links = entrez.elink("YP_009724395.1")  # doctest: +SKIP
    """
    from Bio import Entrez
    Entrez.email = NcbiConfig["email"]
    Entrez.api_key = NcbiConfig["api_key"]

    # Fix input
    if not isinstance(accessions,list):
        accessions = [accessions]
    if not linkname:
        linkname = dbfrom + "_" + dbto

    data = []
    for acc in accessions:
        try:
            raw = list(Entrez.read(Entrez.elink(dbfrom=dbfrom, linkname=linkname, id=acc)))
        except:
            logger.info(f'Entrez.elink failed for accession {acc}, dbfrom: {dbfrom}, dbto: {dbto}. Error: '+str(sys.exc_info()[0]))
            continue
        for d in raw:
            for x in d["LinkSetDb"]:
                for y in x["Link"]:
                    data.append([acc, d["IdList"][0], d["DbFrom"], x["LinkName"], dbto, y["Id"]])
    data = pd.DataFrame(data, columns=["qacc", "quid", "dbfrom", "linkname", "dbto", "tuid"])

    return data

def nucleotide2assembly(nucids):
    """
    Map nucleotide accessions to their genome assemblies.

    Links are resolved with ELink from the ``nuccore`` database to
    the ``assembly`` database, and the assembly accessions are then
    read from the assembly document summaries.

    Parameters
    ----------
    nucids : str or list of str
        NCBI nucleotide accessions.

    Returns
    -------
    pandas.DataFrame
        One row per link, with columns ``nucleotide``, ``nuid``,
        ``auid`` and ``assembly``. The ``assembly`` column is NaN
        for links whose summary could not be resolved.

    Examples
    --------
    >>> from rotifer.db.ncbi import entrez
    >>> t = entrez.nucleotide2assembly(["CP084314.1"])  # doctest: +SKIP
    """
    from Bio import Entrez
    Entrez.email = NcbiConfig["email"]
    Entrez.api_key = NcbiConfig["api_key"]
    t = elink(nucids, dbfrom="nuccore", dbto="assembly")
    t.rename({'qacc':'nucleotide','quid':'nuid','tuid':'auid'}, axis=1, inplace=True)
    t.drop(['dbfrom','linkname','dbto'], axis=1, inplace=True)
    t['assembly'] = np.nan
    if len(t) == 0:
        return t
    auids = t.auid.to_list()
    fh = Entrez.efetch(db="assembly", rettype="docsum", retmode="xml", id=",".join(auids))
    data = Entrez.read(fh)
    if "DocumentSummarySet" not in data:
        return t
    if "DocumentSummary" not in data["DocumentSummarySet"]:
        return t
    for x in data["DocumentSummarySet"]["DocumentSummary"]:
        auid = x.attributes['uid']
        if auid in auids:
            t.loc[t.auid == auid,"assembly"] = x["AssemblyAccession"]
    return t
