__doc__ = """
Fetch sequences from local FASTA files with Easel.

This module wraps the ``esl-sfetch`` program from the Easel toolkit
(distributed with HMMER) to retrieve sequences from indexed FASTA
files. The target files must be indexed; missing ``.ssi`` indices
are built automatically on first use, which requires write access to
the database directory. No network connection is used.
"""

import re
import os
import types
import typing
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from io import StringIO

import rotifer
import rotifer.db.parallel
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

# Defaults
_defaults = {
    'local_database_path': [ os.path.join(rotifer.config['data'],"fadb","nr","nr") ],
    "batch_size": 200,
    "threads": int(np.floor(os.cpu_count()/2)),
}
config = loadConfig(__name__.replace('rotifer.',':'), defaults = _defaults)

class FastaCursor(rotifer.db.parallel.SimpleParallelProcessCursor):
    """
    Fetch biomolecular sequences with Easel's ``esl-sfetch``.

    Each database file is checked at construction time: files that
    do not exist are dropped and missing ``.ssi`` indices are built
    by running ``esl-sfetch --index``.

    Parameters
    ----------
    database_path : str or list of str, optional
        Paths to FASTA files. Defaults to the
        ``local_database_path`` configuration entry.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, default 1
        Accepted for interface compatibility; local reads are
        attempted only once.
    batch_size : int, default 200
        Number of accessions per batch.
    threads : int, optional
        Number of simultaneous processes. Defaults to half of the
        CPU count.

    See Also
    --------
    rotifer.db.ncbi.FastaCursor : delegator that falls back to NCBI

    Examples
    --------
    >>> from rotifer.db.local import easel
    >>> fc = easel.FastaCursor("/databases/fadb/nr/nr")  # doctest: +SKIP
    >>> seqs = fc.fetchall(["WP_063732599.1"])  # doctest: +SKIP
    """
    def __init__(
            self,
            database_path=config["local_database_path"],
            progress=True,
            tries=1,
            batch_size = config['batch_size'],
            threads = config['threads'] or _defaults['threads'],
            *args, **kwargs):
        threads = threads or _defaults['threads']
        super().__init__(progress=progress, tries=1, batch_size=batch_size, threads=threads, *args, **kwargs)
        self.executable = "esl-sfetch"
        self.maxgetitem = 1
        if isinstance(database_path,str) or not isinstance(database_path,typing.Iterable):
            database_path = [ database_path ]
        self.path = []
        for p in database_path:
            if not os.path.exists(p):
                logger.error(f'{p}: no such file!')
                continue
            if not os.path.exists(p + ".ssi"):
                logger.warn(f'Building {self.executable} index for {p}...')
                try:
                    subprocess.run([self.executable,"--index",p])
                except:
                    logger.error("Unable to create index for file {p} ({self.executable})")
                    continue
            self.path.append(p)
        if len(self.path) == 0:
            logger.critical("No index or database for executable {self.executable} in {self.path}")

    def _clean_description(self, seqrec):
        """
        Strip the identifier and record separators from a description.

        Parameters
        ----------
        seqrec : Bio.SeqRecord.SeqRecord
            Record to clean. Modified in place.

        Returns
        -------
        Bio.SeqRecord.SeqRecord
            The same record, with its description reduced to the
            text following the identifier and truncated at the
            first ``\\x01`` byte, which NCBI uses to join the
            descriptions of identical sequences.
        """
        seqrec.description = re.sub("\x01.+", "", seqrec.description.replace(seqrec.id, "").lstrip())
        return seqrec

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

    def fetcher(self, accession, *args, **kwargs):
        """
        Run ``esl-sfetch`` and collect its FASTA output.

        Each configured database file is searched in turn until
        every accession is found or every file was tried.
        Accessions not found in any database are registered as
        missing.

        Parameters
        ----------
        accession : str or iterable of str
            Sequence accessions.

        Returns
        -------
        io.StringIO
            The concatenated FASTA output.
        """
        import tempfile
        from subprocess import Popen, PIPE
        targets = self.parse_ids(accession)
        data = ""
        for db in self.path:
            missing = set()
            while len(targets) > 0:
                if len(targets) == 1:
                    p = Popen([self.executable,db,next(iter(targets))], stderr=PIPE, stdout=PIPE, text=True)
                    o, e = p.communicate()
                else:
                    f = tempfile.NamedTemporaryFile(mode="w+t", delete=True)
                    print("\n".join(list(targets)), file=f)
                    f.flush()
                    p = Popen([self.executable,"-f",db,f.name], stderr=PIPE, stdout=PIPE, text=True)
                    o, e = p.communicate()
                if len(e) > 0:
                    missing.update({e.split(" ")[1]})
                found = set()
                if len(o) > 0:
                    data += o
                    found = set([ x.replace(">","").split(" ")[0] for x in o.split("\n") if len(x) > 0 and x[0] == ">" ])
                    missing.discard(found)
                targets = targets - found - missing
            targets = missing
            if not targets:
                break
        if targets:
            self.update_missing(targets, "Not found")
        return StringIO(data)

    def parser(self, stream, accession, *args, **kwargs):
        """
        Parse a FASTA stream into sequence records.

        NCBI style concatenated descriptions are trimmed to the
        first entry.

        Parameters
        ----------
        stream : file-like
            Stream returned by :meth:`fetcher`.
        accession : str or iterable of str
            Sequence accessions, kept for interface compatibility.

        Returns
        -------
        list of Bio.SeqRecord.SeqRecord
            The parsed records.
        """
        sequence = []
        for seq in SeqIO.parse(stream,"fasta"):
            seq = self._clean_description(seq)
            sequence.append(seq)
        stream.close()
        return sequence
