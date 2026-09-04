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
    Fetch sequences from NCBI using accession numbers and EUtilities.

    This class downloads and parses sequences in Genbank format, i.e.
    the most richly annotated format from NCBI.

    Usage
    -----
    Fetch a protein sequence
    >>> from rotifer.db.ncbi import entrez
    >>> sc = entrez.SequenceCursor(database="protein")
    >>> seqrec = sc.fetchall("YP_009724395.1")

    Fetch several nucleotide entries
    >>> import sys
    >>> from Bio import SeqIO
    >>> from rotifer.db.ncbi import entrez
    >>> sc = entrez.SequenceCursor(database="nucleotide")
    >>> query = ['CP084314.1', 'NC_019757.1', 'AAHROG010000026.1']
    >>> for seqrec in sc.fetchone(query):
    >>>     print(SeqIO.write(seqrec, sys.stdout, "genbank")

    Parameters
    ----------
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
        stack = []
        for s in SeqIO.parse(stream, self._format):
            stack.append(s)
        return stack

    def fetcher(self, accession):
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
        Run EFetch to download NCBI data.

        Parameters
        ----------
        accession: string
          Comma separated list of NCBI database entries.

        Returns
        -------
        List of Bio.SeqRecord objects

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
        targets = self.parse_ids(accessions)
        df = list(self.fetchone(targets))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = pd.DataFrame(columns = self._columns + self._added_columns)
        return df

class TaxonomyCursor(SequenceCursor):
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
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            ids.update(set(o.taxid))
            ids.update(set(o.alternative_taxids.str.split(",").explode().dropna()))
        return ids

    def parser(self, stream, accession):
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
        df = list(self.fetchone(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = pd.DataFrame(self.columns)
        return df

class NucleotideFeaturesCursor(SequenceCursor):
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
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            ids.update(o.nucleotide.dropna().to_list())
        return ids

    def parser(self, stream, accession):
        stream = SeqIO.parse(stream, self._format)
        stream = seqrecords_to_dataframe(stream, exclude_type=self.exclude_type, autopid=self.autopid, assembly=self.assembly, codontable=self.codontable)
        stream = [ x[1].copy() for x in stream.groupby('nucleotide') ]
        return stream

    def fetchall(self, accessions):
        df = list(self.fetchone(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = seqrecords_to_dataframe([])
        return df

class GeneNeighborhoodCursor(rotifer.db.methods.GeneNeighborhoodCursor, NucleotideFeaturesCursor):
    """
    Fetch gene neighbors from nucleotide annotation.

    Usage
    -----
    >>> from rotifer.db.ncbi import ftp
    >>> gfc = ftp.GeneNeighborhoodCursor(progress=True)
    >>> df = gfc.fetchall(["EEE9598493.1"])

    Parameters
    ----------
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
      If set to True, neighborhood data for eukaryotic nucleotide sequences
    exclude_type: list of strings
      List of names for the features that must be ignored
    autopid: boolean
      Automatically set protein identifiers
    codontable: string por int, default 'Bacterial'
      Default codon table, if not set within the data
    progress: boolean, deafult True
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    threads: integer, default 15
      Number of processes to run parallel downloads
    batch_size: int, default 1
      Number of accessions per batch

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
        # A protein with no IPG cannot be placed on a genome, and
        # every backend of this cursor needs one, so the verdict
        # binds them all
        self.final_errors.update(["no IPG","No IPG"])
        if not eukaryotes:
            self.giveup.update(["Eukaryot","eukaryot"])

    def getids2(self, obj, *args, **kwargs):
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
        Find gene neighborhoods in a nucleotide sequence.

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
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
            # No IPG means no way to place the protein on a genome,
            # and every backend of this cursor needs one
            self.update_missing(missing, error="No IPGs", retry=False, final=True)
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
        result = []
        for args in chunk:
            df = self.__getitem__(*args)
            if len(df) == 0:
                continue
            for x in df.groupby('block_id'):
                result.append(x[1])
        return {"result":result,"missing":self.remove_missing()}

    def splitter(self, accessions, ipgs):
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
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            ids.update(o.nucleotide.unique().tolist())
        return ids

    def fetchone(self, accessions, ipgs=None):
        """
        Asynchronously fetch gene neighborhoods from NCBI.

        Parameters
        ----------
        accessions: list of strings
          NCBI protein identifiers

        Returns
        -------
        A generator for rotifer.genome.data.NeighborhoodDF objects
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
    Find related database entries via NCBI's EUtilities.

    Usage:
      # download from NCBI's FTP site
      from rotifer.db.ncbi import entrez
      eutils = entrez.cursor()
      a = eutils.elink("YP_009724395.1")

    Returns:
      Pandas DataFrame

    Parameters:
      accessions : string or list of strings
        NCBI accessions to search links for
      to: string
        Name of the target database
      linkname: string
        Type of link between dbfrom and dbto
        If not set, {dbfrom}_{dbto} is used
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
