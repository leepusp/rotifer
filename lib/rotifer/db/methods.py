"""
Mixin classes shared by cursors of the same data type.

Each mixin in this module implements behavior that depends on the kind
of data a cursor returns (sequences, identical protein reports, genome
records, feature tables or gene neighborhoods) but not on where the
data comes from. Concrete cursors combine one of these mixins with a
transport class such as
:class:`rotifer.db.parallel.SimpleParallelProcessCursor` or
:class:`rotifer.db.delegator.SequentialDelegatorCursor`.
"""

import types
import typing
import pandas as pd
from copy import deepcopy
import rotifer
logger = rotifer.logging.getLogger(__name__)

class SequenceCursor:
    """
    Mixin for cursors that return Bio.SeqRecord objects.
    """

    def getids(self, obj):
        """
        Extract accessions from sequence records.

        Parameters
        ----------
        obj : Bio.SeqRecord.SeqRecord or list of Bio.SeqRecord.SeqRecord
            Sequence records produced by the cursor.

        Returns
        -------
        set of str
            The record identifiers.
        """
        import typing
        if isinstance(obj,list) or isinstance(obj,tuple):
            return set([ x.id for x in obj ])
        else:
            return {obj.id}

class IPGCursor:
    """
    Mixin for cursors that return identical protein group reports.

    Identical protein group (IPG) reports are tables that list, for
    each queried protein, all identical sequences known to NCBI
    together with their genomic coordinates.
    """

    @property
    def columns(self):
        """
        Column names of the IPG dataframe.

        Returns
        -------
        list of str
            The original NCBI report columns followed by the columns
            added by rotifer (``order``, ``is_query`` and
            ``representative``).
        """
        return self._columns + self._added_columns

    def getids(self,obj):
        """
        Extract protein accessions from IPG dataframes.

        Parameters
        ----------
        obj : pandas.DataFrame or list of pandas.DataFrame
            IPG reports produced by the cursor.

        Returns
        -------
        set of str
            All accessions in the ``pid`` and ``representative``
            columns.
        """
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            ids.update(set(o.pid))
            ids.update(set(o.representative))
        return ids

class GenomeCursor:
    """
    Mixin for cursors that return annotated genome sequences.

    Genomes are retrieved as GenBank flat files and parsed into
    Bio.SeqRecord objects, one per contig, each carrying the source
    assembly accession in an ``assembly`` attribute.
    """

    def getids(self,obj):
        """
        Extract assembly accessions from genome records.

        Parameters
        ----------
        obj : Bio.SeqRecord.SeqRecord, list, set or None
            Sequence records produced by the cursor, or a set of
            assembly accessions, which is returned unchanged.

        Returns
        -------
        set of str
            Assembly accessions, read from each record's ``assembly``
            attribute or from ``Assembly:`` entries in ``dbxrefs``.
        """
        if isinstance(obj,types.NoneType):
            return set()
        elif isinstance(obj,set):
            return deepcopy(obj)
        elif not isinstance(obj,list):
            obj = [obj]
        assemblies = set()
        for s in obj:
            if hasattr(s,"assembly"):
                if isinstance(s.assembly,str):
                    assemblies.add(s.assembly)
                else:
                    logger.warn(f'Unknown assembly type {type(assembly)}: {assembly}')
            elif hasattr(s,"dbxrefs") and isinstance(s.dbxrefs,list):
                for x in s.dbxrefs:
                    if 'Assembly:' in x:
                        assemblyID = x.split(':')[-1]
                        assemblies.add(assemblyID)
        return assemblies

    def fetcher(self, accession):
        """
        Open the data streams for one or more genomes.

        Parameters
        ----------
        accession : str or iterable of str
            Assembly accessions.

        Returns
        -------
        list
            Open file-like objects, one per genome found. Each
            stream carries the assembly accession in an ``assembly``
            attribute.
        """
        tries = self.tries
        targets = self.parse_ids(accession)
        stream = []
        for acc in targets:
            fh = self.open_genome(acc)
            if fh == None:
                continue
            fh.assembly = acc
            stream.append(fh)
        self.tries = tries
        return stream

    def parser(self, stream, accession):
        """
        Parse GenBank streams into sequence records.

        Parameters
        ----------
        stream : list
            Open file-like objects returned by :meth:`fetcher`.
        accession : str or iterable of str
            Assembly accessions, kept for interface compatibility.

        Returns
        -------
        list of Bio.SeqRecord.SeqRecord
            One record per contig, each annotated with the source
            assembly accession in the ``assembly`` attribute.
        """
        from Bio import SeqIO
        stack = []
        for fh in stream:
            if isinstance(fh,types.NoneType):
                continue
            for s in SeqIO.parse(fh,"genbank"):
                setattr(s,"assembly",fh.assembly)
                stack.append(s)
            fh.close()
        return stack

class GenomeFeaturesCursor(GenomeCursor):
    """
    Mixin for cursors that return genome annotation as dataframes.

    Genomes are parsed into feature tables, one row per annotated
    genomic feature, using
    :func:`rotifer.genome.utils.seqrecords_to_dataframe`.
    """

    def getids(self,obj):
        """
        Extract assembly accessions from feature tables.

        Parameters
        ----------
        obj : pandas.DataFrame, list, set or None
            Feature tables produced by the cursor, a list of objects
            with an ``assembly`` attribute, or a set of assembly
            accessions, which is returned unchanged.

        Returns
        -------
        set of str
            Assembly accessions.

        Raises
        ------
        TypeError
            If `obj` is of an unsupported type.
        """
        if isinstance(obj,types.NoneType):
            return set()
        elif isinstance(obj,set):
            return deepcopy(obj)
        elif isinstance(obj,list):
            return set([ x.assembly for x in obj ])
        elif isinstance(obj,pd.DataFrame) and "assembly" in obj.columns:
            return set(obj.assembly)
        else:
            raise TypeError(f'Unknown object type {type(obj)}: {obj}')

    def parser(self, stream, accession):
        """
        Parse GenBank streams into a feature table.

        Parameters
        ----------
        stream : file-like or list of file-like
            Open data streams returned by :meth:`GenomeCursor.fetcher`.
        accession : str or iterable of str
            Assembly accessions, kept for interface compatibility.

        Returns
        -------
        pandas.DataFrame
            Concatenated feature table for all input genomes. The
            ``exclude_type``, ``autopid`` and ``codontable``
            attributes of the cursor control the conversion.
        """
        from Bio import SeqIO
        from rotifer.genome.utils import seqrecords_to_dataframe
        if not isinstance(stream, list):
            stream = [stream]
        data = []
        for fh in stream:
            datum = SeqIO.parse(fh,"genbank")
            datum = seqrecords_to_dataframe(
                datum,
                exclude_type = self.exclude_type,
                autopid = self.autopid,
                assembly = fh.assembly,
                codontable = self.codontable,
            )
            data.append(datum)
            fh.close()
        if len(data) > 0:
            data = pd.concat(data)
        else:
            data = seqrecords_to_dataframe([])
        return data

    def fetchall(self, accessions, *args, **kwargs):
        """
        Fetch all accessions as a single dataframe.

        Parameters
        ----------
        accessions : list of str
            Assembly accessions.

        Returns
        -------
        pandas.DataFrame
            The concatenated feature table for every genome found.
            Empty when nothing could be retrieved.
        """
        from rotifer.genome.utils import seqrecords_to_dataframe
        stack = []
        for df in self.fetchone(accessions):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])

class GeneNeighborhoodCursor:
    """
    Mixin for cursors that return gene neighborhood dataframes.

    A gene neighborhood is the set of annotated features located
    around a target gene, identified by the accession of its protein
    product. Neighborhood tables follow the layout of
    :class:`rotifer.genome.data.NeighborhoodDF`.
    """

    def getids(self, obj, ipgs=None):
        """
        Extract protein accessions from neighborhood data.

        Parameters
        ----------
        obj : pandas.DataFrame, str or iterable
            A neighborhood dataframe, a single accession or an
            iterable of accessions. For dataframes, the columns
            named by the cursor's ``column`` attribute are scanned,
            plus the ``pid``, ``replaced`` and ``representative``
            columns when proteins are being searched.
        ipgs : pandas.DataFrame, optional
            Identical protein group reports. When given, accessions
            that share an IPG with the identifiers found in `obj`
            are also returned.

        Returns
        -------
        set of str
            Protein accessions.
        """
        import types

        # extract ids from dataframe
        if isinstance(obj,pd.DataFrame):
            # Load columns
            columns = self.column
            if not isinstance(columns,typing.Iterable) or isinstance(columns,str):
                columns = [columns]
            else:
                columns = list(columns)

            # when searching for proteins, ensure all columns with protein IDs are used
            pids = ['pid','replaced','representative']
            if set(columns).intersection(pids):
                columns += [ x for x in pids if x not in columns ]

            # Load identifiers from object
            ids = set()
            for col in columns:
                if col in obj.columns:
                    ids.update(set(obj[col].dropna().drop_duplicates()))

        # If obj is not a Pandas Dataframe
        elif not isinstance(obj,typing.Iterable) or isinstance(obj,str):
            ids = {obj}
        elif isinstance(obj,typing.Iterable):
            ids = set(obj)
        else:
            logger.error(f'Unknown type {type(obj)}')

        # Add synonyms from IPGs
        if not isinstance(ipgs,types.NoneType):
            ipgids = ipgs[ipgs.pid.isin(ids) | ipgs.representative.isin(ids)].id
            ipgids = ipgs[ipgs.id.isin(ipgids)]
            ids.update(ipgids.pid.dropna())
            ids.update(ipgids.representative.dropna())

        return ids

    def ipgs_to_dict(self, ipgs, column='assembly'):
        """
        Group IPG rows and map each protein to its representative.

        Parameters
        ----------
        ipgs : pandas.DataFrame
            Identical protein group reports.
        column : str, default 'assembly'
            Column to group by, usually ``assembly`` or
            ``nucleotide``.

        Returns
        -------
        dict
            Keys are the values of `column`; values are dictionaries
            mapping each protein accession (``pid``) to its IPG
            representative.
        """
        d = { k: v.set_index('pid').representative.to_dict() for k,v in ipgs.groupby(column) }
        return d

    def ipg_proteins(self, ipgs):
        """
        List every protein accession mentioned in IPG reports.

        Parameters
        ----------
        ipgs : pandas.DataFrame
            Identical protein group reports.

        Returns
        -------
        set of str
            The union of the ``pid`` and ``representative`` columns.
        """
        allipgids = set(ipgs.pid).union(ipgs.representative.drop_duplicates())
        return allipgids

    def genome_ids(self, obj):
        """
        Extract genome identifiers from neighborhood dataframes.

        Parameters
        ----------
        obj : pandas.DataFrame or list of pandas.DataFrame
            Neighborhood dataframes.

        Returns
        -------
        set of str
            The unique values of the cursor's target column, usually
            assembly accessions.
        """
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            ids.update(o[self._target_column].drop_duplicates().dropna())
        return ids

    def fetchall(self, accessions, ipgs=None):
        """
        Fetch all gene neighborhoods as a single dataframe.

        Parameters
        ----------
        accessions : list of str
            NCBI protein identifiers.
        ipgs : pandas.DataFrame, optional
            Precomputed identical protein group reports, passed to
            ``fetchone`` to avoid repeated downloads.

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
            The concatenated neighborhoods.
        """
        stack = []
        for df in self.fetchone(accessions, ipgs=ipgs):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])

