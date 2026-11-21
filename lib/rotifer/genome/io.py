"""
Read genome files into Biopython sequence records.

The public functions are:

:func:`parse`
    Iterate over the :class:`~Bio.SeqRecord.SeqRecord` objects of one or
    more genome files (GenBank, EMBL, FASTA or anything else
    :func:`Bio.SeqIO.parse` reads, plus GFF). Records get a usable
    sequence id, a description and an ``Assembly:`` entry in ``dbxrefs``.
:func:`gff`
    Read a GFF file, with or without a ``##FASTA`` section, and add the
    annotations that GenBank records carry (date, organism, taxonomy,
    topology, keywords). :func:`parse` calls it when ``informat`` is
    ``'gff'``.
:func:`GenbankAssemblyFromFilename`, :func:`SequenceNameToAssembly`
    Ready-made rules for the ``assembly`` argument of :func:`parse`.

A typical session turns files into a table of gene neighbourhoods::

    from rotifer.genome.io import parse
    from rotifer.genome.utils import seqrecords_to_dataframe
    from rotifer.genome.data import NeighborhoodDF

    records = parse(['GCF_000005845.2_ASM584v2_genomic.gbff'], 'genbank')
    ndf = NeighborhoodDF(seqrecords_to_dataframe(records))

Nothing here contacts a remote service, except :func:`gff` when it is asked
to look up an organism by NCBI taxid.
"""
import datetime
import os
import re
import sys
from contextlib import contextmanager
from io import StringIO

from Bio import SeqIO


# Bio.SeqIO formats that must be read from binary files
_BINARY_FORMATS = {"abi", "gck", "nib", "seqxml", "sff", "sff-trim", "snapgene",
                   "twobit", "uniprot-xml", "xdna"}


@contextmanager
def _passthrough(source):
    """
    Yield ``source`` unchanged, as a context manager that closes nothing.

    Parameters
    ----------
    source : object
        A file name or handle that a reader opens and closes itself.

    Yields
    ------
    object
        ``source``.
    """
    yield source


@contextmanager
def _open_for_seqio(source, informat):
    """
    Open a file name for Bio.SeqIO and close it afterwards.

    Parameters
    ----------
    source : str, os.PathLike or file-like
        A file name or a handle.
    informat : str
        A Bio.SeqIO format name, used to choose text or binary mode.

    Yields
    ------
    file-like
        A file opened here when ``source`` is a file name, and closed on
        exit. A handle is yielded as is and left open. For ``'gff'`` the
        name is passed on untouched, because :func:`gff` does its own
        opening.
    """
    if hasattr(source, "read") or informat == "gff":
        yield source
    else:
        with open(source, "rb" if informat in _BINARY_FORMATS else "r") as handle:
            yield handle


def _get_basename(path):
    """
    Name a file without its directory and last extension.

    Parameters
    ----------
    path : str, os.PathLike or file-like
        A file name, or an open file handle. The object is not modified.

    Returns
    -------
    str
        The base name of the file without its last extension. For handles
        the name comes from their ``name`` attribute. Two placeholders
        stand for handles that have no file name: ``'__stdin__'`` for
        ``sys.stdin`` and ``'__stream__'`` for any other handle without a
        usable name (``io.StringIO``, a pipe opened by file descriptor).

    Raises
    ------
    TypeError
        If ``path`` is neither a file name nor a file-like object.

    Examples
    --------
    >>> _get_basename('/data/GCF_000005845.2_ASM584v2_genomic.gbff')
    'GCF_000005845.2_ASM584v2_genomic'
    >>> import io
    >>> _get_basename(io.StringIO('LOCUS'))
    '__stream__'
    """
    if path is sys.stdin:
        return "__stdin__"
    if hasattr(path, "read"):
        name = getattr(path, "name", None)
        if isinstance(name, (str, bytes, os.PathLike)) and name:
            path = name
        else:
            return "__stream__"
    try:
        path = os.fsdecode(path)
    except TypeError:
        raise TypeError("Unknown file-like object") from None
    return os.path.splitext(os.path.basename(path))[0]


def SequenceNameToAssembly(path, s):
    """
    Use the sequence name as the assembly identifier.

    A rule for the ``assembly`` argument of :func:`parse`.

    Parameters
    ----------
    path : str or file-like
        The file being read. Not used.
    s : Bio.SeqRecord.SeqRecord
        A sequence record. Not modified.

    Returns
    -------
    str
        ``s.name``.

    Examples
    --------
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Seq import Seq
    >>> SequenceNameToAssembly('x.gbff', SeqRecord(Seq('A'), id='NC_1', name='NC_1'))
    'NC_1'
    """
    return s.name


def GenbankAssemblyFromFilename(iostream, s):
    r"""
    Read the assembly identifier from the file name.

    A rule for the ``assembly`` argument of :func:`parse`, and its
    default. File names that start like ``GCF_000005845.2_...`` or
    ``GCA_...`` (the NCBI genomes FTP convention, such as
    ``GCF_000005845.2_ASM584v2_genomic.gbff.gz``) give the leading
    ``GC[AF]_<digits>.<version>`` part. Any other name is returned
    without its last extension.

    Parameters
    ----------
    iostream : str, os.PathLike or file-like
        The file being read, as a name or an open handle. Not modified.
    s : Bio.SeqRecord.SeqRecord
        A sequence record. Not used, and not modified.

    Returns
    -------
    str
        The assembly identifier.

    Examples
    --------
    >>> GenbankAssemblyFromFilename('GCF_000005845.2_ASM584v2_genomic.gbff', None)
    'GCF_000005845.2'
    >>> GenbankAssemblyFromFilename('/tmp/my_genome.gbk', None)
    'my_genome'
    """
    aid = _get_basename(iostream)
    if re.match(r'^(GC[AF]_\d+\.\d+)', aid):
        aid = re.sub(r'^(GC[AF]_\d+\.\d+).*', r'\1', aid)
    return aid


def parse(handle, informat, *args, assembly=GenbankAssemblyFromFilename, **kwargs):
    """
    Iterate over the sequence records of one or more genome files.

    Each record is changed before it is returned: an id made only of
    digits becomes ``<file base name>:<id>``, an empty description becomes
    ``"<name>, whole genome shotgun sequence"``, and the ``Assembly:``
    entry of ``dbxrefs`` is set as asked by ``assembly``.

    Assembly identifiers are not part of the INSDC standard but are found
    in some NCBI entries under the DBLINK tag::

        DBLINK    BioProject: PRJNA224116
                  BioSample: SAMEA3138382
                  Assembly: GCF_000210855.2 <-- here

    Parameters
    ----------
    handle : str, os.PathLike, file-like or list of these
        One input or a list of inputs. Files that are opened here are
        closed when the iteration ends; handles given by the caller are
        read but left open.
    informat : str
        ``'gff'`` or a format name accepted by :func:`Bio.SeqIO.parse`.
    *args
        Passed to :func:`gff` when ``informat`` is ``'gff'``, ignored
        otherwise.
    assembly : callable, str or None, default GenbankAssemblyFromFilename
        Rule for the assembly identifier. A callable is called with the
        input (as given in ``handle``) and the record, and its result is
        converted to ``str``. A string is used as is for every record.
        A false value (``None``, ``''``) keeps what the record has: the
        ``Assembly:`` entry already in ``dbxrefs`` or, if none, the
        sequence name. See :func:`GenbankAssemblyFromFilename` and
        :func:`SequenceNameToAssembly`.
    **kwargs
        Passed to :func:`gff` when ``informat`` is ``'gff'``, ignored
        otherwise.

    Yields
    ------
    Bio.SeqRecord.SeqRecord
        The records, in file order. Their ``dbxrefs`` list always has an
        ``Assembly:`` entry.

    Raises
    ------
    TypeError
        When the first record is requested, if ``assembly`` is not a
        callable, a string or a false value.

    Examples
    --------
    >>> import tempfile, os
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio import SeqIO
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = os.path.join(d, 'GCF_000005845.2_ASM584v2_genomic.gbk')
    ...     rec = SeqRecord(Seq('ATGC'), id='NC_1', name='NC_1', description='toy',
    ...                     annotations={'molecule_type': 'DNA'})
    ...     n = SeqIO.write([rec], p, 'genbank')
    ...     [r.dbxrefs for r in parse(p, 'genbank')]
    [['Assembly:GCF_000005845.2']]
    """
    if assembly and not (callable(assembly) or isinstance(assembly, str)):
        raise TypeError(
            "assembly must be a callable, a string or None, "
            f"not {type(assembly).__name__}")

    if not isinstance(handle, list):
        handle = [handle]

    # Process each file, return and wait for next request
    for eachhandle in handle:
        # Parse file basename
        bn = _get_basename(eachhandle)

        # Files opened here are closed when the loop ends or is abandoned.
        # gff() opens and closes its own input.
        if informat == "gff":
            opener = _passthrough(eachhandle)
        else:
            opener = _open_for_seqio(eachhandle, informat)
        with opener as source:
            # Choose parser
            if informat == "gff":
                generator = gff(source, *args, **kwargs)
            else:
                generator = SeqIO.parse(source, informat)

            # Process each Bio.SeqRecord
            for s in generator:
                # Change number-only sequence names
                if re.match(r'^\d+$', s.id):
                    s.id = bn + ":" + s.id
                    s.name = s.id

                # Fix DESCRIPTION
                if not s.description:
                    s.description = s.name + ", whole genome shotgun sequence"

                # Find assembly index
                assemblyStr = s.name # If s has no assembly accession, use nucleotide accession
                assemblyIndex = None
                for i in range(0, len(s.dbxrefs)):
                    if "Assembly:" == s.dbxrefs[i][0:9]:
                        assemblyIndex = i
                        assemblyStr = s.dbxrefs[i]
                        break

                # Change assembly if requested
                if assembly:
                    if callable(assembly):
                        assemblyStr = str(assembly(eachhandle, s))
                    else:
                        assemblyStr = assembly
                    assemblyStr = f'Assembly:{assemblyStr}'

                # Fix assembly
                if assemblyIndex is None:
                    s.dbxrefs.append(assemblyStr)
                else:
                    s.dbxrefs[assemblyIndex] = assemblyStr

                # Return and wait next call
                yield s


def _taxonomy_from_taxid(organism):
    """
    Look up a taxid in the NCBI Taxonomy.

    Parameters
    ----------
    organism : int or str
        An NCBI taxid, as a number or a string of digits.

    Returns
    -------
    tuple of (str, list of str)
        The scientific name and the classification, from the highest
        rank to the lowest.

    Raises
    ------
    ValueError
        If ``organism`` is not a taxid (organism names cannot be looked
        up) or the taxid is not found.
    """
    if isinstance(organism, bool) or not (
            isinstance(organism, int) or
            (isinstance(organism, str) and re.match(r'^\d+$', organism.strip()))):
        raise ValueError(
            f"organism must be an NCBI taxid, not {organism!r}: "
            "organism names cannot be looked up")
    from rotifer.db.ncbi import TaxonomyCursor
    result = TaxonomyCursor(progress=False).fetchall([int(organism)])
    if len(result) == 0:
        raise ValueError(f"taxid {organism} not found in the NCBI Taxonomy")
    name = result.loc[result.index[0], 'organism']
    classification = result.loc[result.index[0], 'classification']
    taxonomy = [x for x in str(classification).split("; ") if x] if isinstance(classification, str) else []
    return name, taxonomy


def gff(path, seqid=None, organism=None, strain=None, description=None, topology='linear', keywords=None, ete3=None):
    """
    Read a GFF file as sequence records.

    The file may end with a ``##FASTA`` section holding the sequences. When
    it does not, the records have no sequence, only features. GenBank-style
    annotations are added to each record: ``date`` (the modification time
    of the file, or today when the input is not a file on disk),
    ``topology``, and optionally ``organism``, ``source``, ``taxonomy`` and
    ``keywords``. Sequences in the ``##FASTA`` section that have no
    features are returned too.

    Parameters
    ----------
    path : str, os.PathLike or file-like
        GFF file name or text handle. A file opened here is closed before
        returning; a handle given by the caller is read to the end and
        left open.
    seqid : callable, optional
        Renames sequences. Called as ``seqid(record, path)`` and the result
        becomes the record id.
    organism : int or str, optional
        NCBI taxid (a number or a string of digits). Sets the ``organism``,
        ``source`` and ``taxonomy`` annotations. Names of organisms are
        not supported.
    strain : str or callable, optional
        Text appended to the organism name, but not to its taxonomy. A
        callable is called as ``strain(path, record)``. Only used with
        ``organism``.
    description : str or callable, optional
        Sets the record description. A callable is called as
        ``description(path, record)``.
    topology : str, default 'linear'
        ``'linear'`` or ``'circular'``, used when the file does not say.
    keywords : list of str, optional
        Keywords added to the ``keywords`` annotation. The list is not
        modified.
    ete3 : optional
        Not used. Kept so that old calls do not break.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
        One record per sequence, in file order. Empty when the input is
        empty.

    Raises
    ------
    ValueError
        If ``organism`` is not a taxid or is not found.

    Examples
    --------
    >>> import tempfile, os
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = os.path.join(d, 'toy.gff')
    ...     with open(p, 'w') as f:
    ...         n = f.write("##gff-version 3\\n")
    ...         n = f.write("chr1\\tsrc\\tgene\\t1\\t9\\t.\\t+\\t.\\tID=g1\\n")
    ...     recs = gff(p)
    >>> [(r.id, len(r.features), r.annotations['topology']) for r in recs]
    [('chr1', 1, 'linear')]
    """
    from BCBio import GFF
    from collections import OrderedDict

    if hasattr(path, "read"):
        lines = path.read().splitlines(keepends=True)
        name = getattr(path, "name", None)
    else:
        with open(path) as h:
            lines = h.readlines()
        name = path

    if path is not sys.stdin and isinstance(name, (str, os.PathLike)) and os.path.exists(name):
        date = datetime.datetime.fromtimestamp(os.path.getmtime(name)).strftime("%d-%b-%Y").upper()
    else:
        date = datetime.date.today().strftime("%d-%b-%Y").upper()

    if len(lines) == 0: # Input is empty
        return []

    # Loading taxonomy
    organismName = None
    if organism:
        organismName, taxonomy = _taxonomy_from_taxid(organism)

    # Split the tabular section from the sequences, if any
    fasta = [x for x, line in enumerate(lines) if line.rstrip("\r\n") == "##FASTA"]
    if fasta:
        n = fasta[0]
        sio = SeqIO.parse(StringIO("".join(lines[n+1:])), "fasta")
        seqrecords = OrderedDict({x.id: x for x in sio})
        base_dict = seqrecords
    else:
        n = len(lines)
        seqrecords = OrderedDict()
        base_dict = None

    # Process each Bio.Seqrecord
    for s in GFF.parse(StringIO("".join(lines[0:n])), base_dict=base_dict):
        # Change sequence ID
        seqrecordID = s.id
        if callable(seqid):
            s.id = seqid(s, path)

        # Annotate date
        s.annotations["date"] = date

        # Annotate taxonomy
        if organismName is not None:
            s.annotations["organism"] = organismName
            if strain:
                if callable(strain):
                    strainStr = strain(path, s)
                else:
                    strainStr = strain
                s.annotations["organism"] = s.annotations["organism"] + strainStr
            s.annotations["source"] = s.annotations["organism"]
            s.annotations["taxonomy"] = list(taxonomy)

        # Add or set keywords
        if keywords:
            if 'keywords' in s.annotations:
                s.annotations["keywords"].extend(keywords)
            else:
                s.annotations["keywords"] = list(keywords)

        # Set description
        if description:
            if callable(description):
                descriptionStr = description(path, s)
            else:
                descriptionStr = description
            s.description = descriptionStr

        # Set topology
        if 'topology' not in s.annotations:
            s.annotations['topology'] = topology

        # Update seqrecord
        if 'gff-version' in s.annotations:
            del(s.annotations['gff-version'])
        if 'sequence-region' in s.annotations:
            del(s.annotations['sequence-region'])
        seqrecords[seqrecordID] = s

    return list(seqrecords.values())
