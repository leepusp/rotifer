"""
Convert Biopython sequence records into gene neighbourhood tables.

:func:`seqrecords_to_dataframe` is the bridge between genome files and
:class:`rotifer.genome.data.NeighborhoodDF`: it takes the records that
:func:`rotifer.genome.io.parse` or :func:`Bio.SeqIO.parse` yield and
writes one row per feature location. Every NCBI, SQLite and parallel
cursor in :mod:`rotifer.db` builds its tables with it.

The private ``_parse_*`` functions are ready-made column extractors for
its ``parsers`` argument.
"""

import types

import numpy as np
import pandas as pd

import rotifer
logger = rotifer.logging.getLogger(__name__)

# Columns of an empty result, in the order seqrecords_to_dataframe writes them
_columns = ['nucleotide', 'start', 'end', 'strand', 'nlen', 'block_id', 'query', 'pid', 'type', 'plen', 'locus', 'seq_type', 'assembly', 'gene', 'origin', 'topology', 'product', 'taxid', 'organism', 'lineage', 'classification', 'feature_order', 'internal_id', 'is_fragment']

def _parse_assembly_from_dbxrefs(column, feature, seqrecord):
    """
    Return the assembly accession from a record's ``Assembly:`` dbxref.

    A column extractor for the ``parsers`` argument of
    :func:`seqrecords_to_dataframe`.

    Parameters
    ----------
    column : str
        Name of the column being filled (unused).
    feature : Bio.SeqFeature.SeqFeature
        The feature of the current row (unused).
    seqrecord : Bio.SeqRecord.SeqRecord
        The record the feature belongs to.

    Returns
    -------
    str or float
        The text after the last ``:`` of the first dbxref containing
        ``Assembly:``, or NaN when there is none.
    """
    for x in seqrecord.dbxrefs:
        if 'Assembly:' in x:
            return x.split(':')[-1]
    return np.nan

def _parse_annotations(column, feature, seqrecord):
    """
    Return the record annotation named like the column.

    A column extractor for the ``parsers`` argument of
    :func:`seqrecords_to_dataframe`, for record-level values such as
    ``topology`` or ``data_file_division``.

    Parameters
    ----------
    column : str
        Annotation key to read, also the column name.
    feature : Bio.SeqFeature.SeqFeature or None
        The feature of the current row (unused).
    seqrecord : Bio.SeqRecord.SeqRecord
        The record whose ``annotations`` are read.

    Returns
    -------
    object
        ``seqrecord.annotations[column]``, or NaN when the key is absent.
    """
    if column in seqrecord.annotations:
        return seqrecord.annotations[column]
    return np.nan

def _parse_qualifiers(column, feature, seqrecord):
    """
    Return the feature qualifier named like the column.

    A column extractor for the ``parsers`` argument of
    :func:`seqrecords_to_dataframe`; ``bin/rgaparser`` uses it for the
    ``gembase`` qualifier.

    Parameters
    ----------
    column : str
        Qualifier key to read, also the column name.
    feature : Bio.SeqFeature.SeqFeature
        The feature whose ``qualifiers`` are read.
    seqrecord : Bio.SeqRecord.SeqRecord
        The record the feature belongs to (unused).

    Returns
    -------
    str or float
        The single value of the qualifier; several values joined with
        ``\\x01`` (the separator NCBI uses for the descriptions of
        identical sequences, which :mod:`rotifer.db.local.easel` splits
        on); NaN when the feature lacks the qualifier.
    """
    if column in feature.qualifiers:
        if len(feature.qualifiers[column]) == 1:
            return feature.qualifiers[column][0]
        else:
            return "\x01".join(feature.qualifiers[column])
    return np.nan

def seqrecords_to_dataframe(seqrecs=None, exclude_type=[], autopid=False, assembly=None, codontable='Bacterial', parsers={}):
    '''
    Write the features of sequence records into a NeighborhoodDF.

    Each feature location becomes one row; a feature that runs through
    the origin of a circular replicon becomes one row per part, with
    ``origin`` set to 1. Features are numbered in the order of their
    start and end coordinates. Records that are not SeqRecords are
    logged as errors and skipped; records without features are skipped.
    The input records are not modified, except that each record's
    feature list is sorted by location.

    Parameters
    ----------
    seqrecs : Bio.SeqRecord.SeqRecord, list or iterator of them, optional
        The records, for example from :func:`rotifer.genome.io.parse` or
        :func:`Bio.SeqIO.parse`. A single record is treated as a list of
        one. None gives an empty table.
    exclude_type : list of str, default []
        Feature types to leave out (for example ``['gene', 'source']``).
        Pseudogenes are typed ``PSE``.
    autopid : {False, 'auto', 'copy'}, default False
        How to fill ``pid`` for CDSs without a ``protein_id`` or ``ID``
        qualifier: ``'auto'`` writes ``<locus>.p<cds number>``,
        ``'copy'`` copies the locus tag, False leaves it missing.
    assembly : str, optional
        Assembly accession for records that have no ``Assembly:``
        dbxref and whose source file name gives none. When None, the
        record id is used.
    codontable : int or str, default 'Bacterial'
        Genetic code used to compute the protein length of CDSs without
        a ``translation`` qualifier; a ``transl_table`` qualifier takes
        precedence. A CDS that fails to translate is typed ``PSE``.
    parsers : dict, default {}
        Extra columns: each key is a column name and each value a
        function called as ``f(column, feature, record)`` for every
        row. See the ``_parse_*`` functions of this module.

    Returns
    -------
    rotifer.genome.data.NeighborhoodDF
        One row per feature location, with the columns listed in the
        Genomes section of the documentation (``nucleotide``,
        ``start`` and ``end`` 1-based and inclusive, ``strand`` 1 or
        -1, ``block_id`` set to the record id, ``query`` 0, ...)
        followed by one column per key of ``parsers``. ``lineage`` is
        computed from ``classification``. The index is a RangeIndex.
        Empty, with the same columns, when no record has features.

    Examples
    --------
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
    >>> from Bio.SeqRecord import SeqRecord
    >>> rec = SeqRecord(Seq("ATGGCTTAA" * 4), id="NC_1",
    ...                 annotations={"molecule_type": "DNA"})
    >>> rec.features = [
    ...     SeqFeature(FeatureLocation(0, 9, 1), type="CDS",
    ...                qualifiers={"protein_id": ["P1.1"], "locus_tag": ["L1"]}),
    ...     SeqFeature(FeatureLocation(18, 27, -1), type="CDS",
    ...                qualifiers={"protein_id": ["P2.1"], "locus_tag": ["L2"]}),
    ... ]
    >>> df = seqrecords_to_dataframe([rec])
    >>> df[["nucleotide", "start", "end", "strand", "pid", "plen", "locus"]]
      nucleotide  start  end  strand   pid  plen locus
    0       NC_1      1    9       1  P1.1     2    L1
    1       NC_1     19   27      -1  P2.1     3    L2
    '''
    import re
    import os
    import Bio
    from math import ceil, log10
    from rotifer.genome.data import NeighborhoodDF
    ncbiaccre = re.compile(r'^GC[AF]_[0-9]+\.[0-9]+') # Regular expression to extract NCBI assembly IDs

    # Make sure first argument is a list
    if not (isinstance(seqrecs,types.GeneratorType) or isinstance(seqrecs,list) or isinstance(seqrecs,Bio.SeqIO.Interfaces.SequenceIterator)):
        seqrecs = [ seqrecs ]

    # Process each seqrecord
    data = []
    optionalFound = []
    for seqrecord in seqrecs:
        # Check input
        if not isinstance(seqrecord, Bio.SeqRecord.SeqRecord):
            logger.error(f'Element is not a seqrecord: {seqrecord}')
            continue
        if not seqrecord.features:
            #logger.error(f'Bio.Reqrecord {seqrecord.id} contains has features to extract')
            continue
        digits = max(6,ceil(log10(len(seqrecord.features))))
        assemblyID = None

        # Extract SeqRecord data
        nlen = len(seqrecord)
        if seqrecord.annotations:
            annotations = seqrecord.annotations
            topology = annotations['topology'] if 'topology' in annotations else 'linear'
            organism = annotations['organism'] if 'organism' in annotations else 'Unknown'
            taxonomy = '; '.join(annotations['taxonomy']) if 'taxonomy' in annotations else ''
            for x in seqrecord.dbxrefs:
                if 'Assembly:' in x:
                    assemblyID = x.split(':')[-1]
        else:
            topology = 'linear'
            organism = 'Unknown'
            taxonomy = pd.NA

        # SeqRecord doesn't contain a reference to its assembly ID
        # Let's assume all sequences in the current file belong to
        # the same assembly
        if assemblyID == None:
            if hasattr(seqrecs,"stream"): # Bio.SeqIO.Interfaces.SequenceIterator!
                if hasattr(seqrecs.stream,"filename") and callable(seqrecs.stream.filename):
                    assemblyID = seqrecs.stream.filename() # fileinput.FileInput
                elif hasattr(seqrecs.stream,"name"):
                    assemblyID = seqrecs.stream.name # Another _io class
                if assemblyID != None:
                    try:
                        assemblyID = os.path.basename(assemblyID)
                        if ncbiaccre.match(assemblyID):
                            assemblyID = ncbiaccre.search(assemblyID).group(0)
                    except:
                        logger.debug(f"Assembly {str(assemblyID)} doesn't parse as a string")
        if assemblyID == None:
            if assembly:
                assemblyID = assembly
            else:
                assemblyID = seqrecord.id

        # Process source feature
        taxid = np.nan
        seq_type = 'chromosome'
        ftsource = [ x for x in seqrecord.features if x.type == "source" ]
        if len(ftsource):
            ftsource = ftsource[0]
            if "db_xref" in ftsource.qualifiers:
                for dbx in ftsource.qualifiers['db_xref']:
                    if "taxon:" == dbx[0:6]:
                        taxid = int(dbx[6:])
            if 'plasmid' in ftsource.qualifiers:
                seq_type = 'plasmid'

        # Extract data for each feature
        feature_order = {}
        internal_id = -1
        seqrecord.features.sort(key=lambda x: (x.location.start, x.location.end))
        for ft in seqrecord.features:
            internal_id += 1 # Internal ID depends solely on the features annotated in the data source
            qualifiers = ft.qualifiers

            # Feature type
            feature_type = ft.type

            # Pseudogenes
            if feature_type not in feature_order:
                feature_order[feature_type] = 0
            if feature_type == 'pseudogene':
                feature_type = 'PSE'
                if 'PSE' not in feature_order:
                    feature_order['PSE'] = 0
            if exclude_type and feature_type in exclude_type:
                continue

            # PID
            locus = qualifiers['locus_tag'][0] if 'locus_tag' in qualifiers else seqrecord.id + f'.{internal_id:0{digits}}'
            plen = sum([ len(x) for x in ft.location.parts ])
            pid = np.nan
            if feature_type == 'CDS':
                if 'pseudo' in qualifiers:
                    feature_type = 'PSE'
                    if 'PSE' not in feature_order:
                        feature_order['PSE'] = 0
                else:
                    if 'translation' in qualifiers:
                        plen = len(qualifiers['translation'][0])
                    else:
                        selectedTable = codontable
                        if 'transl_table' in qualifiers:
                            selectedTable = int(qualifiers['transl_table'][0])
                        try:
                            plen = len(ft.translate(seqrecord, table=selectedTable, cds=False, to_stop=True))
                        except:
                            feature_type = 'PSE'
                            if 'PSE' not in feature_order:
                                feature_order['PSE'] = 0
                if feature_type == 'CDS':
                    if 'protein_id' in qualifiers:
                        pid = qualifiers['protein_id'][0]
                    elif 'ID' in qualifiers:
                        pid = qualifiers['ID'][0]
                    elif autopid == 'auto':
                        pid = f'{locus}.p{feature_order["CDS"]:0{digits}}'
                    elif autopid == 'copy':
                        pid = locus

            # Other feature attributes
            gene = qualifiers['gene'][0] if 'gene' in qualifiers else np.nan
            product = np.nan
            for tag in ['product','inference']:
                if tag in qualifiers:
                    product = qualifiers[tag][0]
                    break

            # Strand
            try:
                strand = int(ft.location.strand)
            except:
                strand = 1

            # Feature location: check if multiple locations imply running through the origin
            origin = 0
            location = []
            for l in ft.location.parts:
                start = l.start+1
                end = int(l.end)
                if location:
                    if ((end - location[-1][1]) * strand) <= 0:
                        if end == nlen and location[-1][0] == 1:
                            location.insert(-1,[start,end])
                        else:
                            location.append([start,end])
                        origin = 1
                    else:
                        location[-1][0] = min(location[-1][0],start)
                        location[-1][1] = max(location[-1][1],end)
                else:
                    location.append([start,end])

            # Optional qualifiers
            optionalDict = {}
            for column in parsers.keys():
                optionalDict[column] = parsers[column](column, ft, seqrecord)

            # Store each location
            for l in location:
                datum = {
                    'nucleotide': seqrecord.id,
                    'start': l[0],
                    'end': l[1],
                    'strand': strand,
                    'nlen': nlen,
                    'block_id': seqrecord.id,
                    'query':0,
                    'pid': pid,
                    'type': feature_type,
                    'plen': plen,
                    'locus': locus,
                    'seq_type': seq_type,
                    'assembly': assemblyID,
                    'gene':gene,
                    'origin':origin,
                    'topology': topology,
                    'product': product,
                    'taxid': taxid,
                    'organism': organism,
                    'lineage':pd.NA,
                    'classification': taxonomy,
                    'feature_order':feature_order[feature_type],
                    'internal_id':internal_id,
                    'is_fragment':0
                }
                datum.update(optionalDict)
                data.append(datum)

            # Increment feature counters
            feature_order[feature_type] += 1

        #END: for ft in seqrecord.features
    #END: for seqrecord in seqrecs

    # Build dataframe and return
    if len(data) > 0:
        df = pd.DataFrame(data)
        df = NeighborhoodDF(df, update_lineage='classification')
        return(df)
    else:
        return(NeighborhoodDF(pd.DataFrame(columns = _columns + list(parsers.keys()))))

