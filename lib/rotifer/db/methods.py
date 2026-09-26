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
                    logger.warning(f'Unknown assembly type {type(assembly)}: {assembly}')
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


class MappingCursor:
    """
    Mixin for cursors that translate identifiers through UniProt.

    Every such query has the same shape, whichever direction it runs
    in: identifiers of one or more databases are matched, their
    UniProtKB accessions are found, and the identifiers those
    accessions carry in one or more other databases are returned.
    Asking for the cross-references of an accession, asking which
    accession an external identifier belongs to, and translating
    between two external databases are that one query with different
    ends pinned, so cursors mixing this class in implement it once.

    The two ends are named by :meth:`fetchall` and :meth:`fetchone`
    rather than by the constructor, because they describe a question
    rather than a data source: one cursor can answer many of them.

    Attributes
    ----------
    columns : list of str
        ``['source', 'source_type', 'accession', 'target',
        'target_type']``: the queried identifier, the database it
        belongs to, the UniProtKB accession linking the two ends, the
        identifier found, and the database that one belongs to.
    """

    #: Database of sequence checksums in ``idmapping.dat``. A sequence
    #: is looked up by its checksum, which is what makes an identical
    #: sequence findable whatever it has been called.
    CHECKSUM = 'CRC64'

    #: Name standing for a UniProtKB accession where a database name
    #: is expected. It is not a row of the mapping table but its join
    #: key, so both ends accept it and mean the accession itself.
    UNIPROTKB = 'UniProtKB-AC'

    _columns = ['source','source_type','accession','target','target_type']

    @property
    def columns(self):
        """
        Column names of the mapping dataframe.

        Returns
        -------
        list of str
        """
        return list(self._columns)

    @staticmethod
    def is_sequence(obj):
        """
        Find whether an object carries sequences rather than names one.

        Parameters
        ----------
        obj : object

        Returns
        -------
        bool
            True for a Biopython record and for rotifer's own sequence
            object, which holds many at once.
        """
        if hasattr(obj, 'seq'):
            return True
        frame = getattr(obj, 'df', None)
        return isinstance(frame, pd.DataFrame) and 'sequence' in frame.columns

    @classmethod
    def checksum(cls, sequence):
        """
        The checksum UniProt knows a sequence by.

        UniProt identifies a sequence by a CRC64 of its residues, and
        publishes those in ``idmapping.dat`` like any other identifier.
        Biopython computes the same value, prefixed, so the prefix is
        dropped.

        Parameters
        ----------
        sequence : str
            The residues.

        Returns
        -------
        str
            Sixteen uppercase hexadecimal digits.
        """
        from Bio.SeqUtils.CheckSum import crc64

        return crc64(str(sequence).strip().upper()).replace('CRC-', '')

    @classmethod
    def sequence_checksums(cls, obj):
        """
        Read the sequences out of an object and check each one.

        Parameters
        ----------
        obj : Bio.SeqRecord.SeqRecord, rotifer sequence, or iterable
            One or more sequences, in any of the shapes rotifer passes
            them around in.

        Returns
        -------
        dict
            Checksum to the names the sequences carrying it were given.
            Several names can share one checksum, which is the point:
            the same sequence under two names is one sequence.
        """
        found = {}

        def record(name, residues):
            if not residues:
                return
            found.setdefault(cls.checksum(residues), set()).add(str(name))

        if hasattr(obj, 'seq'):
            record(getattr(obj, 'id', ''), obj.seq)
            return found
        frame = getattr(obj, 'df', None)
        if isinstance(frame, pd.DataFrame) and 'sequence' in frame.columns:
            names = frame['id'] if 'id' in frame.columns else frame.index
            for name, residues in zip(names, frame['sequence']):
                record(name, residues)
            return found
        if isinstance(obj, typing.Iterable) and not isinstance(obj, str):
            for item in obj:
                for key, names in cls.sequence_checksums(item).items():
                    found.setdefault(key, set()).update(names)
        return found

    def parse_ids(self, accessions, as_string=True):
        """
        Accept sequences wherever identifiers are accepted.

        A sequence is not an identifier, but it has one: its checksum,
        which UniProt publishes like any other. Turning them into
        checksums here means every cursor and every access style takes
        sequences without knowing it, and the names they came under
        are remembered so a caller can find their way back from a
        result.

        Parameters
        ----------
        accessions : str, iterable, sequence object or mixture
            Identifiers, sequences, or both.
        as_string : bool, default True
            Passed through.

        Returns
        -------
        set of str
        """
        if not isinstance(accessions, (list, tuple, set)):
            accessions = [accessions]
        plain, checksums = [], {}
        for item in accessions:
            if self.is_sequence(item):
                for key, names in self.sequence_checksums(item).items():
                    checksums.setdefault(key, set()).update(names)
            else:
                plain.append(item)
        if checksums:
            if not hasattr(self, '_checksums'):
                self._checksums = {}
            for key, names in checksums.items():
                self._checksums.setdefault(key, set()).update(names)
            plain.extend(checksums)
        return super().parse_ids(plain, as_string=as_string)

    @property
    def checksums(self):
        """
        What the sequences given so far were called.

        A result names a sequence by its checksum, since that is what
        the data holds. This says which of the sequences handed in
        carried it.

        Returns
        -------
        pandas.DataFrame
            Columns ``source`` and ``sequence``, ready to be merged
            onto a result by its ``source`` column.
        """
        rows = [ {'source': key, 'sequence': name}
                 for key, names in getattr(self, '_checksums', {}).items()
                 for name in sorted(names) ]
        return pd.DataFrame(rows, columns=['source','sequence'])

    @staticmethod
    def parse_databases(databases):
        """
        Normalize a ``source`` or ``target`` argument to a list.

        Parameters
        ----------
        databases : str, iterable of str or None
            One database name, several, or None for every database.

        Returns
        -------
        list of str or None
            None when no filter was asked for, which every backend
            reads as "any database".
        """
        if isinstance(databases, types.NoneType):
            return None
        if isinstance(databases, str) or not isinstance(databases, typing.Iterable):
            databases = [databases]
        names = [ str(x) for x in databases if not isinstance(x, types.NoneType) ]
        return names or None

    def databases(self):
        """
        Name the databases this cursor can map to and from.

        A backend that knows its own vocabulary lets a delegator ask
        it only for the databases it could answer for, and record the
        rest as unanswerable here rather than as absent everywhere.
        The distinction matters: a database missing from a backend is
        a gap in that copy of the data, while a database no backend
        supports is a gap in the answer.

        Returns
        -------
        set of str or None
            None when the cursor cannot enumerate them, which is read
            as "any", so nothing is narrowed on its account.
        """
        return None

    def unsupported(self, databases):
        """
        Pick the databases this cursor cannot answer for.

        Parameters
        ----------
        databases : list of str or None
            Databases asked for, or None for every database.

        Returns
        -------
        list of str
            Empty when the cursor supports them all, or cannot say.
        """
        known = self.databases()
        if isinstance(known, types.NoneType) or isinstance(databases, types.NoneType):
            return []
        return [ x for x in databases if x != self.UNIPROTKB and x not in known ]

    def supported(self, databases):
        """
        Narrow a list of databases to the ones this cursor can serve.

        Parameters
        ----------
        databases : list of str or None
            Databases asked for, or None for every database.

        Returns
        -------
        list of str or None
            None is passed through, meaning every database this
            cursor has.
        """
        known = self.databases()
        if isinstance(known, types.NoneType) or isinstance(databases, types.NoneType):
            return databases
        return [ x for x in databases if x == self.UNIPROTKB or x in known ]

    def empty(self):
        """
        Build an empty mapping dataframe.

        Returns
        -------
        pandas.DataFrame
            No rows, and the columns listed in :attr:`columns`.
        """
        return pd.DataFrame([], columns=self.columns)

    def getids(self, obj, *args, **kwargs):
        """
        Report which queried identifiers were translated.

        Only the ``source`` column is read. It holds what the caller
        asked about, so what a delegator still has to look for
        elsewhere is exactly what is missing from it; the identifiers
        found are answers, not queries, and counting them would mark
        the wrong things as done.

        Parameters
        ----------
        obj : pandas.DataFrame, set, str, iterable or None
            Mappings produced by the cursor, or a collection of
            identifiers, returned as a set.

        Returns
        -------
        set of str
        """
        if isinstance(obj, types.NoneType):
            return set()
        elif isinstance(obj, set):
            return deepcopy(obj)
        elif isinstance(obj, pd.DataFrame):
            if obj.empty or 'source' not in obj.columns:
                return set()
            return set(obj['source'].dropna().astype(str))
        elif isinstance(obj, str):
            return {obj}
        elif isinstance(obj, typing.Iterable):
            return set([ str(x) for x in obj ])
        else:
            raise TypeError(f'Unknown object type {type(obj)}: {obj}')

    def _with_checksums(self, accessions, source):
        """
        Make sure a query built from sequences looks where they live.

        Naming source databases and then passing a sequence would look
        past it: a sequence is only ever found under
        :attr:`CHECKSUM`, so that is added rather than the query
        quietly returning nothing.

        Parameters
        ----------
        accessions : object
            Whatever was handed to the query.
        source : list of str or None
            The source databases asked for.

        Returns
        -------
        list of str or None
            None is left alone, since it already looks everywhere.
        """
        source = self.parse_databases(source)
        if isinstance(source, types.NoneType):
            return source
        items = accessions if isinstance(accessions, (list, tuple, set)) else [accessions]
        if any(self.is_sequence(x) for x in items) and self.CHECKSUM not in source:
            return list(source) + [self.CHECKSUM]
        return source

    def fetchall(self, accessions, source=None, target=None, *args, **kwargs):
        """
        Translate every identifier, as a single dataframe.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to translate.
        source : str or list of str, optional
            Databases the queried identifiers belong to. None accepts
            any, and :attr:`UNIPROTKB` means they are UniProtKB
            accessions.
        target : str or list of str, optional
            Databases to translate into. None returns every database,
            and :attr:`UNIPROTKB` returns the accession itself.

        Returns
        -------
        pandas.DataFrame
            The columns listed in :attr:`columns`. Empty, with those
            columns, when nothing is found.
        """
        source = self._with_checksums(accessions, source)
        stack = []
        for df in self.fetchone(accessions, source=source, target=target, *args, **kwargs):
            stack.append(df)
        if not stack:
            return self.empty()
        # Sources overlap: a database several of them carry yields the
        # same row from each, and a mapping stated twice is still one
        # mapping.
        return pd.concat(stack, ignore_index=True).drop_duplicates().reset_index(drop=True)
